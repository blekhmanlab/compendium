# This script opens an XML file of exported BioSample search results
# (see notes.txt), parses the samples that include what we're looking for,
# and saves them to a database.

from collections import defaultdict
import sys # for the command-line params
import xml.etree.ElementTree as ET

#import db

#connection = db.Connection()

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
tree = ET.parse(sys.argv[1])
biosamples = tree.getroot()

# iterate through each entry in the file
for sample in biosamples:
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
            samplesource = filterdata[key]
            break
    
    for key in source_tags:
        if filterdata.get(key) is not None:
            used_source[key] += 1
            samplesource = filterdata[key]
            break
    
    # TODO: write sample into table

    # TODO: add all the random tags to the tag table

print('\n\n\n\n\n\n\n\n\n\n\n\n\n---!!------------HOST:\n')
print(host)

print('\n\n\n\n\n\n\n\n\n\n\n\n\n----!!-----------HOST KEY:\n')
print(used_host)

print('\n----!!-----------SOURCE KEY:\n')
print(used_source)

print(f'{len(biosamples)} total samples')
