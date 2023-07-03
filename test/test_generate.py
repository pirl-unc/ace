import pandas as pd
from .data import get_data_path
from acelib.main import run_ace_generate


def test_generate_small_configuration(small_elispot_configuration):
    print(len(small_elispot_configuration))

def test_generate_large_configuration(large_elispot_configuration):
    print(len(large_elispot_configuration))

def test_generate_disallowed_peptide_pairs():
    data = {
        'peptide_id': [],
        'peptide_sequence': []
    }
    for i in range(1, 26):
        data['peptide_id'].append('peptide_%i' % i)
        data['peptide_sequence'].append('')
    df_peptides = pd.DataFrame(data)
    disallowed_peptide_pairs = [('peptide_6', 'peptide_7'), ('peptide_3', 'peptide_4')]
    enforced_peptide_pairs = [('peptide_1', 'peptide_2')]
    status, df_configuration = run_ace_generate(
        df_peptides=df_peptides,
        num_peptides_per_pool=5,
        num_coverage=3,
        num_processes=1,
        random_seed=1,
        disallowed_peptide_pairs=disallowed_peptide_pairs,
        enforced_peptide_pairs=enforced_peptide_pairs
    )
