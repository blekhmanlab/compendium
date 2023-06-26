"""
This module provides helper functions for interacting with a SQLite database
and loading external data into it.
"""
import sqlite3
import time
import xml.etree.ElementTree as ET

import requests

import config

class Connection(object):
    """Data type holding the data required to maintain a database
    connection and perform queries.

    """
    def __init__(self):
        """Stores db connection info in memory and initiates a
        connection to the specified db."""
        try:
            self.db = sqlite3.connect(config.db_path)
        except sqlite3.Error as ex:
            print(f'FATAL: {ex.sqlite_errorname}')
            exit(1)
        print('Connected!')
        self.setup_tables()

    def write(self, query, params=None):
        """
        Executes a query against the database that
        """
        cursor = self.db.cursor()
        if params is not None:
            if isinstance(params, tuple):
                cursor.execute(query, params)
            elif isinstance(params, list):
                cursor.executemany(query, params)
            else:
                raise Exception('Parameters must either be in a tuple (execute) or list (executemany)')

        else:
            cursor.execute(query)
        self.db.commit()
        results = []
        for result in cursor:
            results.append(result)
        cursor.close()
        return results

    def read(self, query, params=None):
        """Helper function that converts results returned stored in a
        sqlite3 cursor into a less temperamental list format.

        Arguments:
            - query: The SQL query to be executed.
            - params: Any parameters to be substituted into the query.
                sqlite3 handles this better than Python does.
        Returns:
            - A list of tuples, one for each row of results.

        """

        results = []
        try:
            cursor = self.db.cursor()
            if params is not None:
                cursor.execute(query, params)
            else:
                cursor.execute(query)

            for result in cursor:
                results.append(result)

            cursor.close()

            return results

        except sqlite3.Error as ex:
            print(f'ERROR with db query execution: {ex.sqlite_errorname}')

    def setup_tables(self):
        """
        Makes sure all required tables are created in the specified database.
        """
        self.write("""
            CREATE TABLE IF NOT EXISTS samples(
                srs TEXT PRIMARY KEY,
                project TEXT,
                taxon TEXT,
                srr TEXT,
                library_strategy TEXT,
                library_source TEXT,
                pubdate TEXT,
                total_bases INTEGER,
                exported INTEGER DEFAULT 0
            )
        """)

        self.write("""
            CREATE TABLE IF NOT EXISTS tags (
                tagid INTEGER PRIMARY KEY,
                srs TEXT,
                tag TEXT,
                value TEXT
            )
        """)

        self.write("""
            CREATE TABLE IF NOT EXISTS tags (
                tagid INTEGER PRIMARY KEY,
                srs TEXT,
                tag TEXT,
                value TEXT
            )
        """)

        self.write("""
            CREATE TABLE IF NOT EXISTS status (
                project TEXT PRIMARY KEY,
                status TEXT NOT NULL,
                rerun_as_single_end INTEGER DEFAULT 0,
                paired INTEGER,
                note1 TEXT,
                note2 TEXT
            )
        """)

        # Results:
        self.write("""
            CREATE TABLE IF NOT EXISTS asv_counts (
                entryid INTEGER PRIMARY KEY,
                project TEXT NOT NULL,
                sample TEXT NOT NULL,
                asv TEXT NOT NULL,
                count INTEGER NOT NULL
            )
        """)

        self.write("""
            CREATE TABLE IF NOT EXISTS asv_sequences (
                asv_id INTEGER PRIMARY KEY,
                project TEXT NOT NULL,
                asv TEXT NOT NULL,
                seq TEXT
            )
        """)

        self.write("""
            CREATE TABLE IF NOT EXISTS asv_assignments (
                asv_id INTEGER PRIMARY KEY,
                tdatabase TEXT,
                kingdom TEXT,
                phylum TEXT,
                tclass TEXT,
                torder TEXT,
                family TEXT,
                genus TEXT
            )
        """)

    def __del__(self):
        """Closes the database connection when the Connection object
        is destroyed."""
        if self.db is not None:
            self.db.close()

def load_xml(taxon, filename, save_samples=True, save_tags=False):
    """
    Loads the "full text XML" exported from a search of BioSamples and adds
    them to the database.

    Inputs:
        - taxon: The taxon ID from the NCBI taxonomy browser associated with the samples.
    """
    connection = Connection()
    print(f'\n\n\n===================\nProcessing XML for taxon {taxon}\n==========\n\n')

    # load the XML file
    print('loading xml...')
    tree = ET.parse(filename)
    biosamples = tree.getroot()
    print('processing samples!')
    # iterate through each entry in the file
    done = -1
    skipped = 0
    for sample in biosamples:
        done += 1
        if done % 10000 == 0:
            print(f'   {done} out of {len(biosamples)} complete.')
        # find SRA ID of sample
        # example: <BioSample> <Ids> <Id db="SRA">SRS5588834</Id> </Ids> </BioSample>
        sra = None
        for entry in sample.iter('Id'):
            if 'db' in entry.attrib.keys() and entry.attrib['db'] == 'SRA':
                sra = entry.text
        if sra is None:
            skipped += 1
            if skipped % 1000 == 0:
                print(f'Skipped {skipped} samples so far.')
            continue # skip samples without an SRA sample
        #  NOTE: we used to check for BioProject ID here,
        #  but for some reason half the samples don't list a bioproject
        #  even if they have one.

        if save_samples:
            # write sample into table
            connection.write('INSERT INTO samples (srs, taxon) VALUES (?, ?);', (sra, taxon))

        if save_tags:
            # go through all the attributes and tally them
            all_tags = {}
            for tag in sample.iter('Attribute'):
                text = tag.text.lower()
                if 'harmonized_name' in tag.attrib.keys():
                    all_tags[tag.attrib['harmonized_name']] = text
                elif 'attribute_name' in tag.attrib.keys():
                    all_tags[tag.attrib['attribute_name']] = text
            # add all the tags to the tag table
            sql = 'INSERT INTO tags (srs, tag, value) VALUES (?,?,?);'
            params = [(sra, tag, value) for (tag, value) in all_tags.items()]
            connection.write(sql, params)

    print(f'{len(biosamples)} total samples evaluated, {skipped} skipped')

def find_runs(count, per_query=80):
    """
    Queries the NCBI eUtils API to use sample IDs ("SRS" codes)
    to get information about runs ("SRR" codes) that can then
    be downloaded as FASTQ files.

    Inputs:
        - count: int. The upper limit for how many entries to search in total.
        - per_query: int. The number of entries to request in each web request
    """
    connection = Connection()

    todo = connection.read("""
        SELECT srs FROM samples
        WHERE srr IS NULL
        LIMIT ?""", (count,)
    )

    todo = [x[0] for x in todo] # each ID is nested inside a tuple of length 1
    print(f'Found {len(todo)} samples to process')
    cursor = 0
    multiple_runs = 0
    while cursor < len(todo):
        if cursor > 0 and cursor % (per_query * 100) == 0:
            time.sleep(10) # pause for a bit every 100 requests
        if cursor % 1000 == 0:
            print(f'COMPLETE: {cursor}')

        url = config.esearch_url
        # Build up the URL project by project:
        for _ in range(0, per_query):
            url += f'{todo[cursor]}[accn] or '
            cursor += 1
            if cursor == len(todo):
                break # in case the total isn't a multiple of "per_query"
        url = url[:-4] # trim off trailing " or "
        if len(url) >1950:
            print(url)
            print('\n\n\nURL IS TOO LONG! Bailing to avoid cutting off request.')
            exit(1)
        try:
            req = requests.get(url, timeout=config.timeout)
        except requests.exceptions.HTTPError:
            print('ERROR: Error sending request for webenv data. Skipping.')
            time.sleep(1)
            continue

        try:
            tree = ET.fromstring(req.text)
        except ET.ParseError:
            print(f'ERROR: Couldnt parse response retrieving webenv data: {req.text}')
            print('Skipping.')
            time.sleep(1)
            continue

        webenv = tree.find('WebEnv')
        if webenv is None:
            print('\n---------\n')
            print(req.text)
            print("WARNING: Got response without a 'webenv' field. Moving on.")
            print('\n---\n')
            time.sleep(1)
            continue
        time.sleep(0.5)
        url = f'{config.efetch_url}&WebEnv={webenv.text}'
        if len(url) >1950:
            print(url)
            print('\n\n\nURL IS TOO LONG! Bailing to avoid cutting off request.')
            exit(1)

        req = requests.get(url, timeout=config.timeout)
        try:
            tree = ET.fromstring(req.text)
        except ET.ParseError:
            print("WARNING: Misformed response from call to eFetch. Skipping.")
            time.sleep(10)
            continue
        multiple_runs += _record_data(tree)

    print(f"\n\nTOTAL SAMPLES WITH MULTIPLE RUNS: {multiple_runs}.\n\n")

def _record_data(data):
    """Parses a response from the efetch endpoint that has info about
    all the samples in the query."""
    connection = Connection()

    multiple_runs = 0

    for package in data.findall('EXPERIMENT_PACKAGE'):
        sample = None
        tosave = {'run': []}

        for entry in package.iter('SAMPLE'):
            if 'accession' in entry.attrib.keys():
                sample = entry.attrib['accession']
        for entry in package.iter('RUN'):
            if 'accession' in entry.attrib.keys():
                tosave['run'].append(entry.attrib['accession'])
            if 'published' in entry.attrib.keys():
                tosave['pubdate'] = entry.attrib['published']
            if 'total_bases' in entry.attrib.keys():
                tosave['total_bases'] = entry.attrib['total_bases']
        for entry in package.iter('EXTERNAL_ID'):
            if 'namespace' in entry.attrib.keys():
                if entry.attrib['namespace'] == 'BioProject':
                    tosave['project'] = entry.text
                    break
        for entry in package.iter('LIBRARY_STRATEGY'):
            tosave['library_strategy'] = entry.text
        for entry in package.iter('LIBRARY_SOURCE'):
            tosave['library_source'] = entry.text
        for entry in package.iter('INSTRUMENT_MODEL'):
            tosave['instrument'] = entry.text

        # If we found multiple runs, combine them into a single string
        if len(tosave['run']) == 0:
            tosave['run'] = None
        elif len(tosave['run']) == 1:
            tosave['run'] = tosave['run'][0]
        else:
            print(f"MULTIPLE RUNS! {len(tosave['run'])}")
            multiple_runs += 1
            delim = ';'
            tosave['run'] = delim.join(tosave['run'])

        # If there is no SRA run identified, SKIP this entry.
        # Sometimes a sample will have multiple entries, one with
        # a run (and lots of metadata) and another without any info
        # but DIFFERENT metadata. We only want ones that have a run.
        if tosave.get('run') is None:
            continue

        connection.write("""
            UPDATE samples
            SET srr=?
            WHERE srs=?
        """, (tosave.get('run'), sample))


        if tosave.get('project') is not None:
            connection.write("""
                UPDATE samples
                SET project=?
                WHERE srs=?
            """, (tosave.get('project'), sample))
        if tosave.get('library_strategy') is not None:
            connection.write("""
                UPDATE samples
                SET library_strategy=?
                WHERE srs=?
            """, (tosave.get('library_strategy'), sample))
        if tosave.get('library_source') is not None:
            connection.write("""
                UPDATE samples
                SET library_source=?
                WHERE srs=?
            """, (tosave.get('library_source'), sample))
        if tosave.get('pubdate') is not None:
            connection.write("""
                UPDATE samples
                SET pubdate=?
                WHERE srs=?
            """, (tosave.get('pubdate'), sample))
        if tosave.get('total_bases') is not None:
            connection.write("""
                UPDATE samples
                SET total_bases=?
                WHERE srs=?
            """, (tosave.get('total_bases'), sample))
        if tosave.get('instrument') is not None:
            connection.write("""
                UPDATE samples
                SET instrument=?
                WHERE srs=?
            """, (tosave.get('instrument'), sample))
    return multiple_runs
