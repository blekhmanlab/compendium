
def Load_counts():
    """Loads a tab-delimited file in which each column
    is a sample and each row is an ASV, with the cells
    indicating read counts."""

    to_write = []

    with open('results/ASVs_counts.tsv', 'r') as file:
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
