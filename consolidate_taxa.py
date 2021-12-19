from collections import defaultdict
import os
import csv


all_projects = os.listdir('../results')
done = os.listdir('../results/taxa_files')
done = [x.split('_')[0] for x in done] # trim off filenames

todo = [x for x in all_projects if x not in done]
print(f'PROJECTS TO PROCESS: {len(todo)}')
for project in todo:
    if project[0:3] != 'PRJ':
        print(f'Entry {project} doesnt match pattern. Skipping.')
        continue
    print(f'Starting project {project}')
    asv_taxa = {}
    unique_taxa = []
    print('  Loading taxonomy data')
    # load the taxonomic assignments for each ASV
    with open(f'../results/{project}/ASVs_taxonomy.tsv', 'r') as f:
        reader = csv.reader(f, dialect='excel-tab')
        next(reader) # skip the header
        for row in reader:
            taxon = row[1:]
            asv_taxa[row[0]] = taxon
            if taxon not in unique_taxa:
                unique_taxa.append(taxon)

    # load the data
    subjectdata = defaultdict(dict)
    print('  Loading count data')
    with open(f'../results/{project}/ASVs_counts.tsv', 'r') as f:
        reader = csv.reader(f, dialect='excel-tab')
        subjects = next(reader)[1:] # get (ordered) list of subjects
        for row in reader:
            for subj, count in enumerate(row[1:]): # skip the first item, which is the name
                subjectdata[subjects[subj]][row[0]] = int(count)
    # write the data
    print('  Writing count data')
    with open(f'../results/taxa_files/{project}_consolidated.tsv','w') as out:
        writer = csv.writer(out, dialect='excel-tab')
        writer.writerow([''] + subjects)
        for taxon in unique_taxa:
            row = [' '.join(taxon)]
            for subject in subjects:
                result = 0
                for asv, classification in asv_taxa.items():
                    if classification == taxon:
                        result += subjectdata[subject][asv]
                row.append(result)
            writer.writerow(row)
    print('  Done')
