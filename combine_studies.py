from collections import defaultdict
import os
import csv

unique_taxa = []

# TODO: This is where we can strip out blanks that got included with the other samples.
# If we know the approximate read depth for the real samples, anything with way fewer
# reads can get chucked. There are studies where most samples have 50,000 reads and
# 3 samples have literally 35 reads, for example.

sampledata = defaultdict(lambda: defaultdict(int))
samples = []
print('Loading count data')

files = os.listdir('../results/taxa_files')

stats = []

with open('../results/taxa_files/studies_consolidated_LOG.csv','w') as out:
    writer = csv.writer(out)
    writer.writerow(['study','samples','taxa'])
    for inputfile in files:
        print(f'Processing file {inputfile}')
        with open(f'../results/taxa_files/{inputfile}', 'r') as f:
            reader = csv.reader(f, dialect='excel-tab')
            try:
                study_samples = next(reader)[1:] # get (ordered) list of samples
            except StopIteration:
                # this happens if the file is empty
                print(f'Empty file?? {inputfile}')
                continue
            study_samples = [f'{inputfile}_{x}' for x in study_samples] # add study ID to beginning to prevent collisions
            samples += study_samples # add samples to master list

            study_taxa = []
            for row in reader:
                for subj, count in enumerate(row[1:]): # skip the first item, which is the name
                    sampledata[study_samples[subj]][row[0]] = int(count)
                    study_taxa.append(row[0])
                    # make sure we have this taxon in the list
                    if row[0] not in unique_taxa:
                        unique_taxa.append(row[0])

        writer.writerow((inputfile, len(study_samples), len(study_taxa)))

# write the data
print('Writing count data')
with open('../results/taxa_files/studies_consolidated.tsv','w') as out:
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
