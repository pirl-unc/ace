import pandas as pd
from .data import get_data_path
from src.acelib.constants import DeconvolutionResults
from src.acelib.main import run_ace_identify
from src.acelib.aid_plate_reader import AIDPlateReader


def test_identify_small_configuration(small_elispot_configuration):
    ground_truth_hit_peptide_ids = ['peptide_1','peptide_10']
    hit_pool_ids = list(small_elispot_configuration.loc[small_elispot_configuration['peptide_id'].isin(ground_truth_hit_peptide_ids), 'pool_id'].unique())
    df_hits_max = run_ace_identify(
        hit_pool_ids=hit_pool_ids,
        df_configuration=small_elispot_configuration
    )
    hit_peptide_ids = sorted(df_hits_max.loc[df_hits_max['deconvolution_result'] == DeconvolutionResults.HIT, 'peptide_id'].values.tolist())
    assert hit_peptide_ids == ground_truth_hit_peptide_ids, 'peptide_1 and peptide_10 are hits.'

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
    expected_peptide_ids = sorted(['peptide_1', 'peptide_10'])
    hit_peptide_ids = sorted(df_hits_max.loc[df_hits_max['deconvolution_result'] == DeconvolutionResults.HIT, 'peptide_id'].values.tolist())
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
    expected_peptide_ids = sorted(['peptide_1', 'peptide_7', 'peptide_14', 'peptide_25'])
    expected_second_round_peptide_ids = ['peptide_9']
    hit_peptide_ids = sorted(df_hits_max.loc[df_hits_max['deconvolution_result'] == DeconvolutionResults.HIT, 'peptide_id'].values.tolist())
    second_round_assay_peptide_ids = sorted(df_hits_max.loc[df_hits_max['deconvolution_result'] == DeconvolutionResults.CANDIDATE_HIT, 'peptide_id'].values.tolist())
    assert hit_peptide_ids == expected_peptide_ids, 'peptide_1, peptide_7, peptide_14, peptide_25 are hits.'
    assert second_round_assay_peptide_ids == expected_second_round_peptide_ids, 'peptide_9 is a candidate hit.'