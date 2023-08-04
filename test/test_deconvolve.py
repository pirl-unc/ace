import pandas as pd
from .data import get_data_path
from src.acelib.constants import DeconvolutionLabels
from acelib.block_assignment import BlockAssignment
from acelib.block_design import BlockDesign
from acelib.deconvolution import empirical_deconvolve
from acelib.plate_readout import PlateReadout


def test_deconvolve_golfy_assignment_1(golfy_assignment_1):
    ground_truth_hit_peptide_ids = ['peptide_1','peptide_10']
    df_assignment = golfy_assignment_1.to_dataframe()
    hit_pool_ids = list(df_assignment.loc[df_assignment['peptide_id'].isin(ground_truth_hit_peptide_ids), 'pool_id'].unique())
    deconvolution_result = empirical_deconvolve(
        hit_pool_ids=hit_pool_ids,
        df_assignment=golfy_assignment_1.to_dataframe(),
        min_coverage=3
    )
    df_deconvolution = deconvolution_result.to_dataframe()
    hit_peptide_ids = sorted(df_deconvolution.loc[df_deconvolution['deconvolution_result'] == DeconvolutionLabels.CONFIDENT_HIT, 'peptide_id'].values.tolist())
    assert hit_peptide_ids == ground_truth_hit_peptide_ids, 'peptide_1 and peptide_10 are hits.'


def test_deconvolve_aid_plate_reader_readout():
    block_assignment = BlockAssignment.read_excel_file(
        excel_file=get_data_path(name='25peptides_5perpool_3x_noseqsim_sat-solver.xlsx'),
        sheet_name='block_assignment'
    )
    plate_readout = PlateReadout.read_aid_plate_reader_file(
        excel_file=get_data_path(name='25peptides_5perpool_3x_aid-plate-reader_readout.xlsx'),
        plate_id=1
    )
    df_readout = pd.merge(block_assignment.to_dataframe(), plate_readout.to_dataframe(), on=['plate_id', 'well_id'])
    hit_pool_ids = list(df_readout.loc[df_readout['spot_count'] >= 300, 'pool_id'].unique())
    deconvolution_result = empirical_deconvolve(
        hit_pool_ids=hit_pool_ids,
        df_assignment=block_assignment.to_dataframe(),
        min_coverage=3
    )
    df_deconvolution = deconvolution_result.to_dataframe()
    expected_peptide_ids = sorted(['peptide_1', 'peptide_10'])
    hit_peptide_ids = sorted(df_deconvolution.loc[df_deconvolution['deconvolution_result'] == DeconvolutionLabels.CONFIDENT_HIT, 'peptide_id'].values.tolist())
    assert hit_peptide_ids == expected_peptide_ids, 'peptide_1 and peptide_10 are hits.'


def test_deconvolve_pool_id_readout():
    block_assignment = BlockAssignment.read_excel_file(
        excel_file=get_data_path(name='25peptides_5perpool_3x_noseqsim_sat-solver.xlsx')
    )
    df_readout = pd.read_excel(
        get_data_path(name='25peptides_5perpool_3x_pool-id_readout.xlsx')
    )
    hit_pool_ids = list(df_readout.loc[df_readout['spot_count'] >= 300, 'pool_id'].unique())
    deconvolution_result = empirical_deconvolve(
        hit_pool_ids=hit_pool_ids,
        df_assignment=block_assignment.to_dataframe(),
        min_coverage=3
    )
    df_deconvolution = deconvolution_result.to_dataframe()
    expected_peptide_ids = sorted(['peptide_1', 'peptide_10'])
    hit_peptide_ids = sorted(df_deconvolution.loc[df_deconvolution['deconvolution_result'] == DeconvolutionLabels.CONFIDENT_HIT, 'peptide_id'].values.tolist())
    assert hit_peptide_ids == expected_peptide_ids, 'peptide_1 and peptide_10 are hits.'

