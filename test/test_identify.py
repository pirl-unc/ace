import pandas as pd
from .data import get_data_path
from acelib.main import run_ace_identify
from acelib.aid_plate_reader import AidPlateReader


def test_identify_small_configuration(small_elispot_configuration):
    hit_pool_ids = ['pool_4', 'pool_5', 'pool_6', 'pool_7', 'pool_11', 'pool_14']
    df_hits = run_ace_identify(
        hit_pool_ids=hit_pool_ids,
        df_configuration=small_elispot_configuration
    )
    df_hits = df_hits[df_hits['num_coverage'] == 3]
    assert 'peptide_1' in df_hits['peptide_id'].unique(), 'peptide_1 is a hit.'
    assert 'peptide_10' in df_hits['peptide_id'].unique(), 'peptide_10 is a hit.'

def test_identify_small_configuration_aid_version():
    configuration_file = get_data_path(name='25peptides_5perpool_3x.csv')
    df_configuration = pd.read_csv(configuration_file)
    excel_file = get_data_path(name='25peptides_5perpool_3x_aid_plate_reader_simulated_results.xlsx')
    df_hits = AidPlateReader.load_readout_file(excel_file=excel_file, plate_id=1)
    df_readout = pd.merge(df_configuration, df_hits, on=['plate_id', 'well_id'])
    hit_pool_ids = list(df_readout.loc[df_readout['spot_count'] >= 300, 'pool_id'].unique())
    df_hits = run_ace_identify(
        hit_pool_ids=hit_pool_ids,
        df_configuration=df_configuration
    )
    df_hits = df_hits[df_hits['num_coverage'] == 3]
    assert 'peptide_1' in df_hits['peptide_id'].unique(), 'peptide_1 is a hit.'
    assert 'peptide_10' in df_hits['peptide_id'].unique(), 'peptide_10 is a hit.'

def test_identify_medium_configuration(medium_elispot_configuration):
    hit_pool_ids = ['pool_1', 'pool_10', 'pool_20', 'pool_29', 'pool_35', 'pool_40', 'pool_60']
    df_hits = run_ace_identify(
        hit_pool_ids=hit_pool_ids,
        df_configuration=medium_elispot_configuration
    )
    df_hits = df_hits[df_hits['num_coverage'] == 3]
    assert 'peptide_1' in df_hits['peptide_id'].unique(), 'peptide_1 is a hit.'
    assert 'peptide_10' in df_hits['peptide_id'].unique(), 'peptide_10 is a hit.'
    assert 'peptide_100' in df_hits['peptide_id'].unique(), 'peptide_100 is a hit.'
