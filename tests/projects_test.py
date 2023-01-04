import os
import sys

import pytest

from fixtures import Dada_dir

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import projects.parsing

def test_load_summary(Dada_dir):
    proj = projects.parsing.Project('PRRJA12345')
    proj.Load_results_summary()
    assert len(proj.samples) == 6

def test_Project_remove_reverse_reads(Dada_dir):
    proj = projects.parsing.Project('PRRJA12345')
    proj.Load_results_summary()
    # make sure we're starting with paired-end files
    assert len(os.listdir('PRRJA12345/fastq')) == 12

    # Get ready to re-run as single end by deleting the reverse reads
    proj._remove_reverse_reads()

    # Make sure half the files are gone
    fastqs = os.listdir('PRRJA12345/fastq')
    assert len(fastqs) == 6
    # Make sure all the files are forward reads
    for f in fastqs:
        assert f[-8:] == '_1.fastq'

def test_Project_remove_previous_dada(Dada_dir):
    proj = projects.parsing.Project('PRRJA12345')
    proj.Load_results_summary()

    files = os.listdir('PRRJA12345')
    # make sure the summary file is there
    assert 'summary.tsv' in files
    assert 'previous_summary.tsv' not in files
    proj._remove_previous_dada()

    files = os.listdir('PRRJA12345')
    # make sure the summary file has been renamed
    assert 'summary.tsv' not in files
    assert 'previous_summary.tsv' in files

    # do it again, and make sure the summary file gets re-renamed
    proj._remove_previous_dada()

    files = os.listdir('PRRJA12345')
    # make sure the summary file has been renamed
    assert 'previous_summary.tsv' not in files
    assert 'previous_previous_summary.tsv' in files
