import pandas as pd
import pytest
from .data import get_data_path
from acelib.main import run_ace_generate
from golfy import init, is_valid, optimize
from acelib.utilities import convert_golfy_results


@pytest.fixture
def small_elispot_configuration():
    data = {
        'peptide_id': [],
        'peptide_sequence': []
    }
    for i in range(1, 26):
        data['peptide_id'].append('peptide_%i' % i)
        data['peptide_sequence'].append('')
    df_peptides = pd.DataFrame(data)
    status, df_configuration = run_ace_generate(
        df_peptides=df_peptides,
        num_peptides_per_pool=5,
        num_coverage=3,
        num_processes=1,
        random_seed=1
    )
    return df_configuration


@pytest.fixture
def large_elispot_configuration():
    golfy_solution = init(
        num_peptides=120,
        peptides_per_pool=12,
        num_replicates=3
    )
    optimize(golfy_solution, max_iters=1000)
    df_configuration = convert_golfy_results(golfy_assignment=golfy_solution.assignments)
    return df_configuration

