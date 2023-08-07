import pandas as pd
from .data import get_data_path
from acelib.block_assignment import BlockAssignment
from acelib.main import run_ace_deconvolve
from acelib.constants import DeconvolutionLabels, DeconvolveModes
from acelib.plate_readout import PlateReadout


def test_deconvolve_golfy_assignment_25pep5per3x_empirical(golfy_assignment_25pep5per3x):
    ground_truth_hit_peptide_ids = ['peptide_1','peptide_10']
    df_assignment = golfy_assignment_25pep5per3x.to_dataframe()
    hit_pool_ids = list(df_assignment.loc[df_assignment['peptide_id'].isin(ground_truth_hit_peptide_ids), 'pool_id'].unique())
    df_readout = pd.DataFrame({
        'pool_id': hit_pool_ids,
        'spot_count': [300] * len(hit_pool_ids)
    })
    deconvolution_result = run_ace_deconvolve(
        df_readout=df_readout,
        block_assignment=golfy_assignment_25pep5per3x,
        mode=DeconvolveModes.EMPIRICAL,
        empirical_min_coverage=3,
        empirical_min_spot_count=300,
        statistical_min_peptide_activity=1.0
    )
    df_deconvolution = deconvolution_result.to_dataframe()
    hit_peptide_ids = sorted(df_deconvolution.loc[df_deconvolution['deconvolution_result'] == DeconvolutionLabels.CONFIDENT_HIT, 'peptide_id'].values.tolist())
    assert hit_peptide_ids == ground_truth_hit_peptide_ids, 'peptide_1 and peptide_10 are hits.'


def test_deconvolve_golfy_assignment_25pep5per3x_em(golfy_assignment_25pep5per3x):
    ground_truth_hit_peptide_ids = ['peptide_1','peptide_10']
    df_assignment = golfy_assignment_25pep5per3x.to_dataframe()
    df_readout = pd.DataFrame({
        'pool_id': list(df_assignment['pool_id'].unique()),
        'spot_count': [0] * len(list(df_assignment['pool_id'].unique()))
    })
    hit_pool_ids = list(df_assignment.loc[df_assignment['peptide_id'].isin(ground_truth_hit_peptide_ids), 'pool_id'].unique())
    df_readout.loc[df_readout['pool_id'].isin(hit_pool_ids), 'spot_count'] = 300
    deconvolution_result = run_ace_deconvolve(
        df_readout=df_readout,
        block_assignment=golfy_assignment_25pep5per3x,
        mode=DeconvolveModes.EM,
        empirical_min_coverage=3,
        empirical_min_spot_count=300,
        statistical_min_peptide_activity=1.0
    )
    df_deconvolution = deconvolution_result.to_dataframe()
    hit_peptide_ids = sorted(df_deconvolution.loc[df_deconvolution['deconvolution_result'] == DeconvolutionLabels.CONFIDENT_HIT, 'peptide_id'].values.tolist())
    assert hit_peptide_ids == ground_truth_hit_peptide_ids, 'peptide_1 and peptide_10 are hits.'


def test_deconvolve_golfy_assignment_25pep5per3x_lasso(golfy_assignment_25pep5per3x):
    ground_truth_hit_peptide_ids = ['peptide_1','peptide_10']
    df_assignment = golfy_assignment_25pep5per3x.to_dataframe()
    df_readout = pd.DataFrame({
        'pool_id': list(df_assignment['pool_id'].unique()),
        'spot_count': [0] * len(list(df_assignment['pool_id'].unique()))
    })
    hit_pool_ids = list(df_assignment.loc[df_assignment['peptide_id'].isin(ground_truth_hit_peptide_ids), 'pool_id'].unique())
    df_readout.loc[df_readout['pool_id'].isin(hit_pool_ids), 'spot_count'] = 300
    deconvolution_result = run_ace_deconvolve(
        df_readout=df_readout,
        block_assignment=golfy_assignment_25pep5per3x,
        mode=DeconvolveModes.LASSO,
        empirical_min_coverage=3,
        empirical_min_spot_count=300,
        statistical_min_peptide_activity=1.0
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
    deconvolution_result = run_ace_deconvolve(
        df_readout=df_readout,
        block_assignment=block_assignment,
        mode=DeconvolveModes.EMPIRICAL,
        empirical_min_coverage=3,
        empirical_min_spot_count=300,
        statistical_min_peptide_activity=1.0
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
    deconvolution_result = run_ace_deconvolve(
        df_readout=df_readout,
        block_assignment=block_assignment,
        mode=DeconvolveModes.EMPIRICAL,
        empirical_min_coverage=3,
        empirical_min_spot_count=300,
        statistical_min_peptide_activity=1.0
    )
    df_deconvolution = deconvolution_result.to_dataframe()
    expected_peptide_ids = sorted(['peptide_1', 'peptide_10'])
    hit_peptide_ids = sorted(df_deconvolution.loc[df_deconvolution['deconvolution_result'] == DeconvolutionLabels.CONFIDENT_HIT, 'peptide_id'].values.tolist())
    assert hit_peptide_ids == expected_peptide_ids, 'peptide_1 and peptide_10 are hits.'

