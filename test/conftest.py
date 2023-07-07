import pandas as pd
import pytest
from .data import get_data_path
from acelib.main import run_ace_golfy, run_ace_sat_solver
from golfy import init, is_valid, optimize
from acelib.utilities import convert_golfy_results


@pytest.fixture
def small_golfy_elispot_configuration():
    data = {
        'peptide_id': [],
        'peptide_sequence': []
    }
    for i in range(1, 26):
        data['peptide_id'].append('peptide_%i' % i)
        data['peptide_sequence'].append('')
    df_peptides = pd.DataFrame(data)
    is_valid, df_configuration = run_ace_golfy(
        df_peptides=df_peptides,
        num_peptides_per_pool=5,
        num_coverage=3,
        random_seed=1,
        max_iters=2000,
        init_mode='greedy'
    )
    assert is_valid, 'With 2000 max iterations, we should have had a valid solution.'
    return df_configuration


@pytest.fixture
def small_sat_solver_elispot_configuration():
    data = {
        'peptide_id': [],
        'peptide_sequence': []
    }
    for i in range(1, 26):
        data['peptide_id'].append('peptide_%i' % i)
        data['peptide_sequence'].append('')
    df_peptides = pd.DataFrame(data)
    df_configuration = run_ace_sat_solver(
        df_peptides=df_peptides,
        num_peptides_per_pool=5,
        num_coverage=3,
        num_peptides_per_batch=100,
        random_seed=1,
        num_processes=1,
        is_first_coverage=True
    )
    return df_configuration


@pytest.fixture
def large_golfy_elispot_configuration():
    data = {
        'peptide_id': [],
        'peptide_sequence': []
    }
    for i in range(1, 121):
        data['peptide_id'].append('peptide_%i' % i)
        data['peptide_sequence'].append('')
    df_peptides = pd.DataFrame(data)
    is_valid, df_configuration = run_ace_golfy(
        df_peptides=df_peptides,
        num_peptides_per_pool=12,
        num_coverage=3,
        random_seed=1,
        max_iters=2000,
        init_mode='greedy'
    )
    assert is_valid, 'With 2000 max iterations, we should have had a valid solution.'
    return df_configuration
