from .data import get_data_path
from acelib.main import *


def test_verify():
    tsv_file = get_data_path(name='elispot_configuration_disallowed_pairs_sample.tsv')
    df_configuration = pd.read_csv(tsv_file, sep='\t')
    is_valid = run_ace_verify(
        df_configuration=df_configuration,
        n_peptides=25,
        n_peptides_per_pool=5,
        n_coverage=3
    )

    assert is_valid, \
        "'elispot_configuration_sample.tsv' is a valid ELIspot configuration."
