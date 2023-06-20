import pandas as pd
from .data import get_data_path
from acelib.constants import DeconvolutionResults
from acelib.main import run_ace_identify
from acelib.aid_plate_reader import AIDPlateReader


def test_identify_small_configuration(small_elispot_configuration):
    hit_pool_ids = ['pool_4', 'pool_5', 'pool_6', 'pool_7', 'pool_11', 'pool_14']
    df_hits_max = run_ace_identify(
        hit_pool_ids=hit_pool_ids,
        df_configuration=small_elispot_configuration
    )
    expected_peptide_ids = ['peptide_1', 'peptide_10'].sort()
    hit_peptide_ids = df_hits_max.loc[df_hits_max['deconvolution_result'] == DeconvolutionResults.HIT, 'peptide_id'].values.tolist().sort()
    assert hit_peptide_ids == expected_peptide_ids, 'peptide_1 and peptide_10 are hits.'

def test_identify_small_configuration_aid_version():
    configuration_file = get_data_path(name='25peptides_5perpool_3x.csv')
    df_configuration = pd.read_csv(configuration_file)
    excel_file = get_data_path(name='25peptides_5perpool_3x_aid_plate_reader_simulated_results.xlsx')
    df_hits = AIDPlateReader.load_readout_file(excel_file=excel_file, plate_id=1)
    df_readout = pd.merge(df_configuration, df_hits, on=['plate_id', 'well_id'])
    hit_pool_ids = list(df_readout.loc[df_readout['spot_count'] >= 300, 'pool_id'].unique())
    df_hits_max = run_ace_identify(
        hit_pool_ids=hit_pool_ids,
        df_configuration=df_configuration
    )
    expected_peptide_ids = ['peptide_1', 'peptide_10'].sort()
    hit_peptide_ids = df_hits_max.loc[df_hits_max['deconvolution_result'] == DeconvolutionResults.HIT, 'peptide_id'].values.tolist().sort()
    assert hit_peptide_ids == expected_peptide_ids, 'peptide_1 and peptide_10 are hits.'

def test_identify_small_configuration_pool_id_version():
    configuration_file = get_data_path(name='25peptides_5perpool_3x.csv')
    df_configuration = pd.read_csv(configuration_file)
    csv_file = get_data_path(name='25peptides_5perpool_3x_simulated_results.csv')
    df_readout = pd.read_csv(csv_file)
    hit_pool_ids = list(df_readout.loc[df_readout['spot_count'] >= 300, 'pool_id'].unique())
    df_hits_max = run_ace_identify(
        hit_pool_ids=hit_pool_ids,
        df_configuration=df_configuration
    )
    expected_peptide_ids = ['peptide_1', 'peptide_7', 'peptide_14', 'peptide_25'].sort()
    expected_second_round_peptide_ids = ['peptide_9']
    hit_peptide_ids = df_hits_max.loc[df_hits_max['deconvolution_result'] == DeconvolutionResults.HIT, 'peptide_id'].values.tolist().sort()
    second_round_assay_peptide_ids = df_hits_max.loc[df_hits_max['deconvolution_result'] == DeconvolutionResults.CANDIDATE_HIT, 'peptide_id'].values.tolist()
    assert hit_peptide_ids == expected_peptide_ids, 'peptide_1, peptide_7, peptide_14, peptide_25 are hits.'
    assert second_round_assay_peptide_ids == expected_second_round_peptide_ids, 'peptide_9 is a candidate hit.'

def test_identify_medium_configuration(medium_elispot_configuration):
    hit_pool_ids = ['pool_1', 'pool_10', 'pool_20', 'pool_29', 'pool_35', 'pool_40', 'pool_60']
    df_hits_max = run_ace_identify(
        hit_pool_ids=hit_pool_ids,
        df_configuration=medium_elispot_configuration
    )
    expected_peptide_ids = ['peptide_1', 'peptide_10', 'peptide_100'].sort()
    hit_peptide_ids = df_hits_max.loc[df_hits_max['deconvolution_result'] == DeconvolutionResults.HIT, 'peptide_id'].values.tolist().sort()
    assert hit_peptide_ids == expected_peptide_ids, 'peptide_1 and peptide_10 are hits.'
