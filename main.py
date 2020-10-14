from collections import defaultdict
import os
import sys
import xml.etree.ElementTree as ET

isolation_source = defaultdict(int)
env_material = defaultdict(int)
source = defaultdict(int)
host = defaultdict(int)
tags = defaultdict(int)

tree = ET.parse(sys.argv[1])
root = tree.getroot()
biosamples = [x for x in root]

noproject = 0
nosra = 0
for sample in biosamples:
    # find SRA ID of sample
    # example: <BioSample> <Ids> <Id db="SRA">SRS5588834</Id> </Ids> </BioSample>
    sra = None
    for x in sample.iter('Id'):
        if 'db' in x.attrib.keys() and x.attrib['db'] == 'SRA':
            sra = x.text
    if sra is None:
        nosra += 1
        continue # skip samples without an SRA sample
    #  NOTE: we used to check for BioProject ID here,
    #  but for some reason half the samples don't list a bioproject
    # even if they have one.
    
    # go through all the attributes and tally them
    sampleenv = None
    samplehost = None
    sampleiso = None
    for tag in sample.iter('Attribute'):
        text = tag.text.lower()
        if 'attribute_name' in tag.attrib.keys():
            tags[tag.attrib['attribute_name']] += 1
            if tag.attrib['attribute_name'] == 'isolation_source':
                #isolation_source[text] += 1
                sampleiso = text
            if tag.attrib['attribute_name'] == 'env_material':
                #env_material[text] += 1
                sampleenv = text
            if tag.attrib['attribute_name'] == 'host':
                host[text] += 1
                samplehost = text
    # load up the manually curated list of acceptable values for "host":
    host_tokeep = []
    with open('tokeep_host.txt','r') as f:
        for line in f:
            host_tokeep.append(line[:-1]) # chop off newline character

    if samplehost in host_tokeep:
        if sampleiso is not None:
            source[sampleiso] += 1
        elif sampleenv is not None:
            source[sampleenv] += 1
#print('\n\n\n\n\n\n\n\n\n\n\n\n\n---!!------------HOST:\n')
#print(host)

#print('\n\n\n\n\n\n\n\n\n\n\n\n\n----!!-----------SOURCE:\n')
#print(source)
print(f'{len(biosamples)} total samples ({noproject} without a project, {nosra} without SRA)')
