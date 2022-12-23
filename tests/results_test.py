import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import results.parsing

def inc(x):
    return x + 1

def test_answer():
    x = results.parsing.Load_asv_data('tests')
    assert len(x) > 5000