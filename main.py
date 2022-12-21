from collections import defaultdict
import os
import sys # for the command-line params
import time
import xml.etree.ElementTree as ET

import requests

import config
import db
import results

connection = db.Connection()

def load_xml(taxon, filename, save_samples=True, save_tags=False):
    """
    Loads the "full text XML" exported from a search of BioSamples and adds
    them to the database.

    Inputs:
        - taxon: The taxon ID from the NCBI taxonomy browser associated with the samples.
    """
    print(f'\n\n\n===================\nProcessing XML for taxon {taxon}\n==========\n\n')

    # load the XML file
    print('loading xml...')
    tree = ET.parse(filename)
    biosamples = tree.getroot()
    print('processing samples!')
    # iterate through each entry in the file
    done = -1
    for sample in biosamples:
        done += 1
        if done % 10000 == 0:
            print(f'   {done} out of {len(biosamples)} complete.')
        # find SRA ID of sample
        # example: <BioSample> <Ids> <Id db="SRA">SRS5588834</Id> </Ids> </BioSample>
        sra = None
        for x in sample.iter('Id'):
            if 'db' in x.attrib.keys() and x.attrib['db'] == 'SRA':
                sra = x.text
        if sra is None:
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
            with connection.db.cursor() as cursor:
                sql = 'INSERT INTO tags (srs, tag, value) VALUES (?,?,?);'
                params = [(sra, tag, value) for (tag, value) in all_tags.items()]
                cursor.executemany(sql, params)
                connection.db.commit()

    print(f'{len(biosamples)} total samples')

def find_runs(count, per_query=80):
    """
    Queries the NCBI eUtils API to use sample IDs ("SRS" codes)
    to get information about runs ("SRR" codes) that can then
    be downloaded as FASTQ files.

    Inputs:
        - count: int. The upper limit for how many entries to search in total.
        - per_query: int. The number of entries to request in each web request
    """

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

        url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?tool={config.Tool}&email={config.Email}&api_key={config.Key}&db=sra&usehistory=y&term='
        for x in range(0, per_query):
            url += f'{todo[cursor]}[accn] or '
            cursor += 1
            if cursor == len(todo): break # in case the total isn't a multiple of "per_query"
        url = url[:-4] # trim off trailing " or "
        if len(url) >1950:
            print(url)
            print('\n\n\nURL IS TOO LONG! Bailing to avoid cutting off request.')
            exit(1)
        try:
            r = requests.get(url)
        except:
            print('ERROR: Error sending request for webenv data. Skipping.')
            time.sleep(1)
            continue

        try:
            tree = ET.fromstring(r.text)
        except:
            print(f'ERROR: Couldnt parse response retrieving webenv data: {r.text}')
            print('Skipping.')
            time.sleep(1)
            continue
        
        webenv = tree.find('WebEnv')
        if webenv is None:
            print('\n---------\n')
            print(r.text)
            print("WARNING: Got response without a 'webenv' field. Moving on.")
            print('\n---\n')
            time.sleep(1)
            continue
        time.sleep(0.5)
        url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?tool={config.Tool}&email={config.Email}&db=sra&query_key=1&WebEnv={webenv.text}'
        if len(url) >1950:
            print(url)
            print('\n\n\nURL IS TOO LONG! Bailing to avoid cutting off request.')
            exit(1)

        r = requests.get(url)
        try:
            tree = ET.fromstring(r.text)
        except ET.ParseError:
            print("WARNING: Misformed response from call to eFetch. Skipping.")
            time.sleep(10)
            continue
        multiple_runs += _record_data(tree)

    print(f"\n\nTOTAL SAMPLES WITH MULTIPLE RUNS: {multiple_runs}.\n\n")

def _record_data(data):
    """Parses a response from the efetch endpoint that has info about
    all the samples in the query."""
    multiple_runs = 0

    for package in data.findall('EXPERIMENT_PACKAGE'):
        sample = None
        tosave = {'run': []}

        for x in package.iter('SAMPLE'):
            if 'accession' in x.attrib.keys():
                sample = x.attrib['accession']
        for x in package.iter('RUN'):
            if 'accession' in x.attrib.keys():
                tosave['run'].append(x.attrib['accession'])
            if 'published' in x.attrib.keys():
                tosave['pubdate'] = x.attrib['published']
            if 'total_bases' in x.attrib.keys():
                tosave['total_bases'] = x.attrib['total_bases']
        for x in package.iter('EXTERNAL_ID'):
            if 'namespace' in x.attrib.keys():
                if x.attrib['namespace'] == 'BioProject':
                    tosave['project'] = x.text
                    break
        for x in package.iter('LIBRARY_STRATEGY'):
            tosave['library_strategy'] = x.text
        for x in package.iter('LIBRARY_SOURCE'):
            tosave['library_source'] = x.text
        for x in package.iter('INSTRUMENT_MODEL'):
            tosave['instrument'] = x.text

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

def write_lists(min_samples=10):
    """
    Fetches a list of SRA projects from the local database and generates a list
    of samples for each project.

    Inputs:
        - min_samples: int. The minimum number of samples that a project needs to
            have to get a list generated.
    """

    todo = connection.read("""
        SELECT samplecount.project
        FROM (
            SELECT project, COUNT(srs) AS tally
            FROM samples
            WHERE srr IS NOT NULL AND project IS NOT NULL
                AND library_source IN ('GENOMIC','METAGENOMIC')
                AND library_strategy='AMPLICON'
            GROUP BY 1
            ORDER BY 2 ASC
        ) AS samplecount
        WHERE samplecount.tally >= ?""", (min_samples,))
    todo = [x[0] for x in todo] # each ID is nested inside a tuple of length 1

    project_samples = []

    for project in todo:
        samples = connection.read("""
            SELECT s.srr
            FROM SAMPLES s
            WHERE srr IS NOT NULL
                AND library_source IN ('GENOMIC','METAGENOMIC')
                AND library_strategy='AMPLICON'
            AND project=?""", (project,))
        samples = [x[0] for x in samples]
        project_samples.append((project, len(samples)))

        if os.path.exists(f'accession_lists/{project}/SraAccList.txt'):
            print(f'Project already recorded: {project}')
            continue

        os.mkdir(f'accession_lists/{project}')

        with open(f'accession_lists/{project}/SraAccList.txt','w') as f:
            for sample in samples:
                f.write(f'{sample}\n')

        sql = """
        UPDATE samples SET exported=true
        WHERE srr IN (
            SELECT s.srr
            FROM SAMPLES s
            WHERE srr IS NOT NULL
                AND library_source IN ('GENOMIC','METAGENOMIC')
                AND library_strategy='AMPLICON'
            AND project=?
        )
        """
        connection.write(sql, (project,))

    with open(f'samples_per_project.csv','a') as f:
        for x in project_samples:
            f.write(f'{x[0]}, {x[1]}\n')

if __name__ == "__main__":
    #load_xml('txid408170', 'txid408170.221220.xml', save_samples=True, save_tags=False)
    #load_xml('txid408170', save_samples=False, save_tags=True)

    # only command-line param is how many to do in this session
    if sys.argv[1] == 'runs':
        todo = 2000 if len(sys.argv) < 3 else sys.argv[2]
        find_runs(todo, per_query=80)
    elif sys.argv[1] == 'results':
        x = results.Load_asv_data('PRJNA842201')
        print(x[0:2])
    #write_lists(min_samples=50)

