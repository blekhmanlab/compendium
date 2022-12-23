import os
from pathlib import Path
import shutil
import sys

import pytest

@pytest.fixture
def Dada_dir():
    """
    Generates a directory that looks like the results of a successful
    analysis pipeline.
    """

    proj = 'PRRJA12345'
    real_files = [
        'ASVs.fa','ASVs_counts.tsv',
        'ASVs_taxonomy.tsv', 'SraAccList.txt',
        'summary.tsv'
    ]
    fake_files = [ 'filtered_out.rds',
        'forward_error_model.pdf',
        'reverse_error_model.pdf',
        'err_forward_reads.rds',
        'err_reverse_reads.rds',
        'ASV.tsv', 'asv.rds',
    ]
    os.mkdir(proj)
    os.mkdir(f'{proj}/fastq')
    os.mkdir(f'{proj}/intermediate')

    for x in real_files:
        shutil.copyfile(f'tests/mocks/{x}', f'{proj}/{x}')
    for x in fake_files:
        Path(f'{proj}/{x}').touch()
    with open(f'{proj}/SraAccList.txt', 'r') as f:
        for sample in f.readlines():
            for ext in [1,2]:
                Path(f'{proj}/fastq/{sample[:-1]}_{ext}.fastq').touch()
                Path(f'{proj}/intermediate/{sample[:-1]}.R{ext}.fastq.gz').touch()
    yield True
    # cleanup
    try:
       shutil.rmtree(proj)
    except OSError as e:
        print(f'Error deleting dir tests/fastq: {e.strerror}')
