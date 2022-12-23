import config

def Load_counts(project):
    """Loads a tab-delimited file in which each column
    is a sample and each row is an ASV, with the cells
    indicating read counts."""

    to_write = []

    with open(f'results/{project}/ASVs_counts.tsv', 'r') as file:
        samples = file.readline().split('\t')[1:] # first entry is blank
        samples[-1] = samples[-1][:-1] # strip newline
        for line in file:
            line = line.split('\t')
            line[-1] = line[-1][:-1]
            asv = [line[0]] * len(samples)
            entries = list(zip(samples, asv, line[1:]))
            # convert counts to ints
            entries = [
                (x[0], x[1], int(x[2]))
                for x in entries
            ]
            to_write += [x for x in entries if x[2] > 0]
    return(to_write)

def Load_asv_data(project):
    """Loads a tab-delimited file in which each row is
    a numbered ASV, associated with its inferred taxonomic
    source."""

    asvs = {}
    # Get exact sequences
    with open(f'results/{project}/ASVs.fa', 'r') as file:
        while file:
            asv = file.readline()
            if asv == '': # empty line with no newline
                break
            if len(asv) > 2:
                asv = asv[1:-1] # strip off leading '>' and trailing newline
            seq = file.readline()
            if len(seq) > 1:
                seq = seq[0:-1] # strip trailing newline
            asvs[asv]=[seq]

    # The get taxonomic assignments
    with open(f'results/{project}/ASVs_taxonomy.tsv', 'r') as file:
        file.readline() # skip header
        for line in file:
            line = line.split('\t')
            line[-1] = line[-1][:-1]
            asvs[line[0]] += line[1:]

    # example entry: ('ASV_1', 'CCTACGGG', 'Bacteria', 'Bacteroidota', 'Bacteroidia', 'Bacteroidales',
    #    'Bacteroidaceae', 'Bacteroides')
    return([tuple([asv]+values) for asv, values in asvs.items()])
