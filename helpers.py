def filter_sources():
    """
    Given a list of "sources" for samples (pulled from various tags),
    records a new text file indicating whether each source fits our string
    criteria for inclusion in the dataset.
    """
    potentialsources = {}
    tokeep = ['fec', 'faec', 'stool', '2003','meconium'] 
    totoss = ['cultur','swab','incubat','rectum','rectal', 'virus','tissue','soil']
    kept = 0
    with open('sources.txt','r') as f:
        for line in f:
            test = line[:-1] # chop off newline character
            if any(x in test for x in tokeep) and not any(x in test for x in totoss):
                potentialsources[test] = 'yes'
                kept += 1
            else:
                potentialsources[test] = 'no'
    with open('processed_sources.txt','w') as f:
        for key, decision in potentialsources.items():
            f.write(f'{key}\t{decision}\n')
    print(f'Kept {kept} out of {len(potentialsources.keys())}')
