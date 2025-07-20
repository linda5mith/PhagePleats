import pandas as pd
import numpy as np
import pytest
from PhagePleats.phagepleats import read_search, map_query_to_genome, classify_novelty

def test_read_search(tmp_path):
    # Create a dummy Foldseek TSV
    f = tmp_path / "dummy.tsv"
    f.write_text("q1.pdb\tt1.pdb\t100\t10\t0\t0\t1\t10\t1\t10\t1e-5\t50\n")
    df = read_search(f)
    assert df.shape[0] == 1
    assert df['query'].iloc[0] == "q1"
    assert df['target'].iloc[0] == "t1"

def test_map_query_to_genome():
    search = pd.DataFrame({'query': ['p1', 'p2']})
    metadata = pd.DataFrame({'protein': ['p1'], 'genome': ['g1']})
    mapped = map_query_to_genome(search, metadata)
    assert mapped['query_genome'].isna().sum() == 1
    assert mapped['query_genome'].iloc[0] == 'g1'

@pytest.mark.parametrize("z_s, z_d, expected", [
    (None, 1.2, "Unknown"),
    (-3.2, 0.5, "Potential new order"),
    (-2.1, 1.5, "Potential new family"),
    (-0.5, 3.1, "Potential new order"),
    (0.1, 0.3, "Likely member"),
])
def test_classify_novelty(z_s, z_d, expected):
    assert classify_novelty(z_s, z_d, rank="Order") == expected
