import os
import sys

from fixtures import Dada_dir

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import results

def test_answer(Dada_dir):
    x = results.Load_asv_data('PRRJA12345')
    assert len(x) > 5000
