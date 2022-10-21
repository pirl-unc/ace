from .data import get_data_path
from acelib.main import *


def test_generate_vanilla():
    n_peptides = 25
    n_peptides_per_pool = 10
    n_coverage = 3
    df_configuration = run_ace_generate(
        n_peptides=n_peptides,
        n_peptides_per_pool=n_peptides_per_pool,
        n_coverage=n_coverage,
        peptide_ids=[],
        disallowed_peptide_pairs=[],
        num_processes=1
    )
