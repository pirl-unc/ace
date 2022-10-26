from .data import get_data_path
from acelib.main import *


def test_generate_disallowed_pairs():
    n_peptides = 25
    n_peptides_per_pool = 5
    n_coverage = 3
    df_configuration = run_ace_generate(
        n_peptides=n_peptides,
        n_peptides_per_pool=n_peptides_per_pool,
        n_coverage=n_coverage,
        peptide_ids=[],
        disallowed_peptide_pairs=[('peptide_6', 'peptide_7'), ('peptide_3', 'peptide_4')],
        num_processes=1
    )
