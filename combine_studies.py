from collections import defaultdict
import csv

unique_taxa = []

# TODO: This is where we can strip out blanks that got included with the other samples.
# If we know the approximate read depth for the real samples, anything with way fewer
# reads can get chucked. There are studies where most samples have 50,000 reads and
# 3 samples have literally 35 reads, for example.

sampledata = defaultdict(lambda: defaultdict(int))
samples = []
print('Loading count data')
for inputfile in ['PRJDB5310', 'PRJDB5316', 'PRJDB5420', 'PRJDB7128', 'PRJDB8820', 'PRJDB9608', 'PRJEB1799', 'PRJEB34610', 'PRJEB39038', 'PRJEB39707', 'PRJNA279828', 'PRJNA282013', 'PRJNA297268', 'PRJNA306884', 'PRJNA311147', 'PRJNA314806', 'PRJNA321051', 'PRJNA322554', 'PRJNA324147', 'PRJNA324825', 'PRJNA327713', 'PRJNA342805', 'PRJNA353658', 'PRJNA362530', 'PRJNA362944', 'PRJNA379120', 'PRJNA383300', 'PRJNA386260', 'PRJNA388732', 'PRJNA389481', 'PRJNA403824', 'PRJNA413706', 'PRJNA418903', 'PRJNA419373', 'PRJNA421288', 'PRJNA422125', 'PRJNA427170', 'PRJNA428226', 'PRJNA445346', 'PRJNA450340', 'PRJNA450540', 'PRJNA450690', 'PRJNA471615', 'PRJNA472768', 'PRJNA480312', 'PRJNA487119', 'PRJNA489693', 'PRJNA496539', 'PRJNA509882', 'PRJNA510713', 'PRJNA511459', 'PRJNA513489', 'PRJNA515810', 'PRJNA521738', 'PRJNA522306', 'PRJNA522626', 'PRJNA528344', 'PRJNA528960', 'PRJNA529487', 'PRJNA535518', 'PRJNA545546', 'PRJNA547558', 'PRJNA547806', 'PRJNA564636', 'PRJNA574920', 'PRJNA579560', 'PRJNA589036', 'PRJNA591924', 'PRJNA596333', 'PRJNA608722', 'PRJNA608934', 'PRJNA622267', 'PRJNA637202', 'PRJNA642975', 'PRJNA661679', 'PRJNA666641']:
    print(f'Processing file {inputfile}')
    with open(f'results/taxa_files/{inputfile}_consolidated.tsv', 'r') as f:
        reader = csv.reader(f, dialect='excel-tab')
        study_samples = next(reader)[1:] # get (ordered) list of samples
        study_samples = [f'{inputfile}_{x}' for x in study_samples] # add study ID to beginning to prevent collisions
        samples += study_samples # add samples to master list
        for row in reader:
            for subj, count in enumerate(row[1:]): # skip the first item, which is the name
                sampledata[study_samples[subj]][row[0]] = int(count)
                # make sure we have this taxon in the list
                if row[0] not in unique_taxa:
                    unique_taxa.append(row[0])
# write the data
print('Writing count data')
with open('results/taxa_files/studies_consolidated.tsv','w') as out:
    writer = csv.writer(out, dialect='excel-tab')
    # columns are SAMPLES:
    # writer.writerow([''] + samples)
    # for taxon in unique_taxa:
    #     row = [taxon]
    #     for sample in samples:
    #         row.append(sampledata[sample][taxon])
    #     writer.writerow(row)
    # columns are TAXA:
    writer.writerow([''] + unique_taxa)
    for sample in samples:
        row = [sample]
        for taxon in unique_taxa:
            row.append(sampledata[sample][taxon])
        writer.writerow(row)

print(f'Done. {len(unique_taxa)} taxa found in {len(samples)} samples.')
