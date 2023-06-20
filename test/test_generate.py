import pandas as pd
from .data import get_data_path
from acelib.main import run_ace_generate


def test_generate_small_configuration(small_elispot_configuration):
    print(len(small_elispot_configuration))

def test_generate_medium_configuration(medium_elispot_configuration):
    print(len(medium_elispot_configuration))

def test_generate_disallowed_peptide_pairs():
    data = {
        'peptide_id': [],
        'peptide_sequence': []
    }
    for i in range(1, 26):
        data['peptide_id'].append('peptide_%i' % i)
        data['peptide_sequence'].append('')
    df_peptides = pd.DataFrame(data)

    def func(df):
        return [('peptide_6', 'peptide_7'), ('peptide_3', 'peptide_4')]

    status, df_configuration = run_ace_generate(
        df_peptides=df_peptides,
        num_peptides_per_pool=5,
        num_coverage=3,
        num_processes=1,
        dissimilarity_inference_func=func
    )

def test_generate_recycle_configuration():
    df_template_configuration = get_data_path(name='100peptides_5perpool_3x.csv')
