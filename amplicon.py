from collections import defaultdict
import statistics # for mean

import skbio

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2562909/
BOUNDARIES = {
    'v1': (69,99),
    'v2': (137,242),
    'v3': (433,497),
    'v4': (576, 682),
    'v5': (822, 879),
    'v6': (986, 1043),
    'v7': (1117, 1173),
    'v8': (1243, 1294),
    'v9': (1435, 1465)
}

# https://www.ncbi.nlm.nih.gov/nuccore/J01859
whole16s = 'aaattgaagagtttgatcatggctcagattgaacgctggcggcaggcctaacacatgcaagtcgaacggtaacaggaagaagcttgctctttgctgacgagtggcggacgggtgagtaatgtctgggaaactgcctgatggagggggataactactggaaacggtagctaataccgcataacgtcgcaagaccaaagagggggaccttcgggcctcttgccatcggatgtgcccagatgggattagctagtaggtggggtaacggctcacctaggcgacgatccctagctggtctgagaggatgaccagccacactggaactgagacacggtccagactcctacgggaggcagcagtggggaatattgcacaatgggcgcaagcctgatgcagccatgccgcgtgtatgaagaaggccttcgggttgtaaagtactttcagcggggaggaagggagtaaagttaatacctttgctcattgacgttacccgcagaagaagcaccggctaactccgtgccagcagccgcggtaatacggagggtgcaagcgttaatcggaattactgggcgtaaagcgcacgcaggcggtttgttaagtcagatgtgaaatccccgggctcaacctgggaactgcatctgatactggcaagcttgagtctcgtagaggggggtagaattccaggtgtagcggtgaaatgcgtagagatctggaggaataccggtggcgaaggcggccccctggacgaagactgacgctcaggtgcgaaagcgtggggagcaaacaggattagataccctggtagtccacgccgtaaacgatgtcgacttggaggttgtgcccttgaggcgtggcttccggagctaacgcgttaagtcgaccgcctggggagtacggccgcaaggttaaaactcaaatgaattgacgggggcccgcacaagcggtggagcatgtggtttaattcgatgcaacgcgaagaaccttacctggtcttgacatccacggaagttttcagagatgagaatgtgccttcgggaaccgtgagacaggtgctgcatggctgtcgtcagctcgtgttgtgaaatgttgggttaagtcccgcaacgagcgcaacccttatcctttgttgccagcggtccggccgggaactcaaaggagactgccagtgataaactggaggaaggtggggatgacgtcaagtcatcatggcccttacgaccagggctacacacgtgctacaatggcgcatacaaagagaagcgacctcgcgagagcaagcggacctcataaagtgcgtcgtagtccggattggagtctgcaactcgactccatgaagtcggaatcgctagtaatcgtggatcagaatgccacggtgaatacgttcccgggccttgtacacaccgcccgtcacaccatgggagtgggttgcaaaagaagtaggtagcttaaccttcgggagggcgcttaccactttgtgattcatgactggggtgaagtcgtaacaaggtaaccgtaggggaacctgcggttggatcacctcctta'

def find_region(location, direction='f'):
    if direction == 'f': # forward
        for v, coords in BOUNDARIES.items():
            if location < coords[0]:
                return(v)
            # If the read starts WITHIN a hypervariable region:
            if location > coords[0] and location < coords[1]:
                hyperlen = coords[1] - coords[0]
                covered = coords[1] - location
                # It still counts if at least half of
                # the region is covered
                if covered / hyperlen >= 0.5:
                    return(v)

    elif direction == 'r': # reverse
        for v, coords in reversed(BOUNDARIES.items()):
            if location > coords[1]:
                return(v)
            # If the read ends WITHIN a hypervariable region:
            if location > coords[0] and location < coords[1]:
                hyperlen = coords[1] - coords[0]
                covered = location - coords[0]
                # It still counts if more than half
                # the region is covered
                if covered / hyperlen >= 0.5:
                    return(v)
    else:
        raise ValueError('direction parameter must be "f" or "r".')

def process_project(asvs, verbose=False):
    results = {}
    lengths = []
    for asv in asvs:
        lengths.append(len(asv))
    avglength = statistics.mean(lengths)

    asv_align = skbio.alignment.StripedSmithWaterman(whole16s)
    forwards = defaultdict(int)
    reverses = defaultdict(int)
    start = None
    end = None
    evaluated = 0
    for asv in asvs:
        evaluated += 1
        results[asv] = asv_align(asv)
        if verbose:
            print(f'{results[asv].query_begin} / {results[asv].query_end} / {results[asv].optimal_alignment_score}')
        align_length = results[asv].query_end - results[asv].query_begin

        # Only keep matches where more than 70% of bases count toward the score:
        if align_length <= len(asv) * 0.7:
            continue
        if start is None:
            region = find_region(results[asv].query_begin, direction='f')
            forwards[region] += 1
            if forwards[region] > len(asvs) / 2:
                start = region
                if verbose:
                    print(f'\n\n !!!!\nDetermined start region! {start}')
        if end is None:
            region = find_region(results[asv].query_end, direction='r')
            reverses[region] += 1
            if reverses[region] > len(asvs) / 2:
                end = region
                if verbose:
                    print(f'\n\n !!!!\nDetermined end region! {end}')
        if start is not None and end is not None:
            break

    # If we have the start but not the end:
    if start is not None and end is None:
        if verbose:
            print(f'Could not determine end region. Using start region and ASV length.')
        startpoint = BOUNDARIES[start][0]
        endpoint = startpoint + avglength
        end = find_region(endpoint, direction='r')

    # If we have the END but not the start:
    if start is None and end is not None:
        if verbose:
            print(f'Could not determine start region. Using end region and ASV length.')
        endpoint = BOUNDARIES[end][1]

        startpoint = endpoint - avglength
        print(f'endpoint is {endpoint}, start is {startpoint}')
        start = find_region(startpoint, direction='f')
    # In some wonky studies with multiple amplicons, we end
    # up with a "start" region that's after the "end" region.
    # We don't want these.
    if start is not None and end is not None and start > end:
        print(f'Start region {start} is after end region {end}. Throwing out result.')
        start = None
        end = None
    # Print!
    assignment = f'{start}{f"-{end}" if end != start else ""}'
    if verbose:
        print('FORWARD:')
        for v, count in forwards.items():
            print(f'{v}: {count}')
        print('REVERSE:')
        for v, count in reverses.items():
            print(f'{v}: {count}')

        print(f'Evaluated {evaluated} of {len(asvs)} ASVs in {proj}')
        print(f'Our guess: {assignment}')
        print(f'Average length is {avglength} +/- {statistics.stdev(lengths)}')
    return((assignment, avglength))
