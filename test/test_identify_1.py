from .data import get_data_path
from acelib.main import *
from acelib.default_parameters import *


def test_visualize():
    readout_tsv_file = get_data_path(name='elispot_readout_sample.tsv')
    configuration_tsv_file = get_data_path(name='elispot_configuration_sample.tsv')
    df_readout = pd.read_csv(readout_tsv_file, sep='\t')
    df_configuration = pd.read_csv(configuration_tsv_file, sep='\t')
    df_hits = run_ace_identify(
        df_readout=df_readout,
        df_configuration=df_configuration,
        min_positive_spot_count=MIN_POSITIVE_SPOT_COUNT
    )

    assert 'peptide_1' in df_hits['hit_peptide_id'].values.tolist(), \
        "'peptide_1' is a positive hit peptide but is not included in the DataFrame."
    assert 'peptide_20' in df_hits['hit_peptide_id'].values.tolist(), \
        "'peptide_20' is a positive hit peptide but is not included in the DataFrame."
    assert 'peptide_25' in df_hits['hit_peptide_id'].values.tolist(), \
        "'peptide_25' is a positive hit peptide but is not included in the DataFrame."