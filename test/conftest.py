import pandas as pd
import pytest
from .data import get_data_path
from src.acelib.main import run_ace_generate


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
        num_processes=1
    )
    return df_configuration

