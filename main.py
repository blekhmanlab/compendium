# This script opens an XML file of exported BioSample search results
# (see notes.txt), parses the samples that include what we're looking for,
# and saves them to a database.

from collections import defaultdict
import sys # for the command-line params
import time
import xml.etree.ElementTree as ET

import requests

import db

connection = db.Connection()

def load_xml(xmlfile):
    """
    Loads the "full text XML" exported from a search of BioSamples and adds
    them to the database.

    Inputs:
        - xmlfile: String. Path to the file to be read.
    """
    # for counting entries in each field:
    source = defaultdict(int)
    host = defaultdict(int)
    # for tracking which tags were used for the "host" and "source" fields:
    used_host = defaultdict(int)
    used_source = defaultdict(int)

    source_tags = [
            'host_body_product','body_product','host_body_habitat',
            'sample_type','isolation_source','env_medium'
        ]

    host_tags = ['host_taxid','host']
    # load the XML file
    print('loading xml...')
    tree = ET.parse(xmlfile)
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
        # even if they have one.
        
        # go through all the attributes and tally them
        filterdata = {}
        all_tags = {}
        for tag in sample.iter('Attribute'):
            text = tag.text.lower()
            if 'harmonized_name' in tag.attrib.keys():
                if tag.attrib['harmonized_name'] in (source_tags + host_tags):
                    filterdata[tag.attrib['harmonized_name']] = text
                else:
                    all_tags[tag.attrib['harmonized_name']] = text
            elif 'attribute_name' in tag.attrib.keys():
                all_tags[tag.attrib['attribute_name']] = text
        # Then wade through the tags we found:
        samplehost = None
        samplesource = None
        for key in host_tags:
            if filterdata.get(key) is not None:
                used_host[key] += 1
                samplehost = filterdata[key]
                break
        if samplehost is None: continue

        host[samplehost] += 1
        # iterate through the potential values for "source" and get the first
        # one that has a value
        for key in source_tags:
            if filterdata.get(key) is not None:
                used_source[key] += 1
                samplesource = filterdata[key]
                break
        
        if samplesource is None: continue
        source[samplesource] += 1
        
        # write sample into table
        with connection.db.cursor() as cursor:
            cursor.execute('INSERT INTO samples (srs, host, source) VALUES (%s, %s, %s);', (sra, samplehost, samplesource))


        # add all the random tags to the tag table
        with connection.db.cursor() as cursor:
            sql = 'INSERT INTO tags (srs, tag, value) VALUES (%s, %s, %s);'
            params = [(sra, tag, value) for (tag, value) in filterdata.items()]
            cursor.executemany(sql, params)
        with connection.db.cursor() as cursor:
            sql = 'INSERT INTO tags (srs, tag, value) VALUES (%s, %s, %s);'
            params = [(sra, tag, value) for (tag, value) in all_tags.items()]
            cursor.executemany(sql, params)

    print(f'{len(biosamples)} total samples')

def get_entrez(count):
    """
    Queries the NCBI eUtils API to turn sample IDs ("SRS" codes)
    into Entrez ID numbers, which can be used for other searches later.
    
    Inputs:
        - count: int. The upper limit for how many entries to search.
    """
    todo = connection.read("""SELECT s.srs FROM samples s
    INNER JOIN acceptable_hosts ah ON ah.host=s.host
    INNER JOIN acceptable_sources asa ON asa.source=s.source
    WHERE ah.keep AND asa.keep
        AND entrez_id IS NULL
    LIMIT %s""", (count,))
    todo = [x[0] for x in todo] # each ID is nested inside a tuple of length 1
    per_query = 5 # how many to send at one time
    cursor = 0
    while cursor < len(todo):
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term='
        for x in range(0, per_query):
            url += f'{todo[cursor]}[accn] or '
            cursor += 1
            if cursor == len(todo): break # in case the total isn't a multiple of "per_query"
        url = url[:-4] # trim off trailing " or "
        print(url)
        r = requests.get(url)
        tree = ET.fromstring(r.text)
        found = []
        for x in tree.iter('Id'):
            print(f'   ID: {x.text}')
            found.append(x.text)
        if len(found) != per_query:
            print(f'   WAIT! Asked for {per_query} IDs but got {len(found)}\n\n')
        time.sleep(3)
        
if __name__ == "__main__":
    #load_xml(sys.argv[1])
    get_entrez(50)