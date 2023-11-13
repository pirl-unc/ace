import pandas as pd
from .data import get_data_path
from acelib.block_assignment import BlockAssignment
from acelib.main import run_ace_deconvolve
from acelib.constants import DeconvolutionLabels, DeconvolutionMethods
from acelib.plate_readout import PlateReadout


def test_deconvolve_golfy_assignment_25pep5per3x_em(golfy_assignment_25pep5per3x):
    ground_truth_hit_peptide_ids = ['peptide_1','peptide_10']
    df_assignment = golfy_assignment_25pep5per3x.to_dataframe()
    hit_pool_ids = list(df_assignment.loc[df_assignment['peptide_id'].isin(ground_truth_hit_peptide_ids), 'pool_id'].unique())
    data = {
        'pool_id': [],
        'spot_count': []
    }
    for pool_id in df_assignment['pool_id'].unique():
        data['pool_id'].append(pool_id)
        if pool_id in hit_pool_ids:
            data['spot_count'].append(300)
        else:
            data['spot_count'].append(0)
    df_readout = pd.DataFrame(data)
    deconvolution_result = run_ace_deconvolve(
        df_readout=df_readout,
        block_assignment=golfy_assignment_25pep5per3x,
        method=DeconvolutionMethods.EM,
        min_coverage=3,
        min_pool_spot_count=300
    )
    df_deconvolution = deconvolution_result.to_dataframe()
    hit_peptide_ids = sorted(df_deconvolution.loc[df_deconvolution['deconvolution_result'] == DeconvolutionLabels.CONFIDENT_HIT, 'peptide_id'].values.tolist())
    assert hit_peptide_ids == ground_truth_hit_peptide_ids, 'peptide_1 and peptide_10 are hits.'


def test_deconvolve_golfy_assignment_25pep5per3x_empirical(golfy_assignment_25pep5per3x):
    ground_truth_hit_peptide_ids = ['peptide_1','peptide_10']
    df_assignment = golfy_assignment_25pep5per3x.to_dataframe()
    hit_pool_ids = list(df_assignment.loc[df_assignment['peptide_id'].isin(ground_truth_hit_peptide_ids), 'pool_id'].unique())
    data = {
        'pool_id': [],
        'spot_count': []
    }
    for pool_id in df_assignment['pool_id'].unique():
        data['pool_id'].append(pool_id)
        if pool_id in hit_pool_ids:
            data['spot_count'].append(300)
        else:
            data['spot_count'].append(0)
    df_readout = pd.DataFrame(data)
    deconvolution_result = run_ace_deconvolve(
        df_readout=df_readout,
        block_assignment=golfy_assignment_25pep5per3x,
        method=DeconvolutionMethods.EMPIRICAL,
        min_coverage=3,
        min_pool_spot_count=300
    )
    df_deconvolution = deconvolution_result.to_dataframe()
    hit_peptide_ids = sorted(df_deconvolution.loc[df_deconvolution['deconvolution_result'] == DeconvolutionLabels.CONFIDENT_HIT, 'peptide_id'].values.tolist())
    assert hit_peptide_ids == ground_truth_hit_peptide_ids, 'peptide_1 and peptide_10 are hits.'


def test_deconvolve_golfy_assignment_25pep5per3x_cem(golfy_assignment_25pep5per3x):
    ground_truth_hit_peptide_ids = ['peptide_1','peptide_10']
    df_assignment = golfy_assignment_25pep5per3x.to_dataframe()
    hit_pool_ids = list(df_assignment.loc[df_assignment['peptide_id'].isin(ground_truth_hit_peptide_ids), 'pool_id'].unique())
    data = {
        'pool_id': [],
        'spot_count': []
    }
    for pool_id in df_assignment['pool_id'].unique():
        data['pool_id'].append(pool_id)
        if pool_id in hit_pool_ids:
            data['spot_count'].append(300)
        else:
            data['spot_count'].append(0)
    df_readout = pd.DataFrame(data)
    deconvolution_result = run_ace_deconvolve(
        df_readout=df_readout,
        block_assignment=golfy_assignment_25pep5per3x,
        method=DeconvolutionMethods.CONSTRAINED_EM,
        min_coverage=3,
        min_pool_spot_count=300
    )
    df_deconvolution = deconvolution_result.to_dataframe()


def test_deconvolve_golfy_assignment_25pep5per3x_lasso(golfy_assignment_25pep5per3x):
    ground_truth_hit_peptide_ids = ['peptide_1','peptide_10']
    df_assignment = golfy_assignment_25pep5per3x.to_dataframe()
    hit_pool_ids = list(df_assignment.loc[df_assignment['peptide_id'].isin(ground_truth_hit_peptide_ids), 'pool_id'].unique())
    data = {
        'pool_id': [],
        'spot_count': []
    }
    for pool_id in df_assignment['pool_id'].unique():
        data['pool_id'].append(pool_id)
        if pool_id in hit_pool_ids:
            data['spot_count'].append(300)
        else:
            data['spot_count'].append(0)
    df_readout = pd.DataFrame(data)
    deconvolution_result = run_ace_deconvolve(
        df_readout=df_readout,
        block_assignment=golfy_assignment_25pep5per3x,
        method=DeconvolutionMethods.LASSO,
        min_coverage=3,
        min_pool_spot_count=300
    )
    df_deconvolution = deconvolution_result.to_dataframe()


def test_deconvolve_aid_plate_reader_readout():
    block_assignment = BlockAssignment.read_excel_file(
        excel_file=get_data_path(name='25peptides_5perpool_3x_configuration.xlsx'),
    )
    plate_readout = PlateReadout.read_aid_plate_reader_file(
        excel_file=get_data_path(name='25peptides_5perpool_3x_aid-plate-reader_readout.xlsx'),
        plate_id=1
    )
    df_readout = pd.merge(block_assignment.to_dataframe(), plate_readout.to_dataframe(), on=['plate_id', 'well_id'])
    deconvolution_result = run_ace_deconvolve(
        df_readout=df_readout,
        block_assignment=block_assignment,
        method=DeconvolutionMethods.CONSTRAINED_EM,
        min_coverage=3,
        min_pool_spot_count=300
    )
    df_deconvolution = deconvolution_result.to_dataframe()
    expected_peptide_ids = sorted(['peptide_1', 'peptide_10'])
    hit_peptide_ids = sorted(df_deconvolution.loc[df_deconvolution['deconvolution_result'] == DeconvolutionLabels.CONFIDENT_HIT, 'peptide_id'].values.tolist())
    assert hit_peptide_ids == expected_peptide_ids, 'peptide_1 and peptide_10 are hits.'


def test_deconvolve_pool_id_readout():
    block_assignment = BlockAssignment.read_excel_file(
        excel_file=get_data_path(name='25peptides_5perpool_3x_configuration.xlsx')
    )
    df_assignment = block_assignment.to_dataframe()
    df_readout = pd.read_excel(
        get_data_path(name='25peptides_5perpool_3x_pool-id_readout.xlsx')
    )
    pool_ids = []
    for index, row in df_readout.iterrows():
        curr_plate_id = row['plate_id']
        curr_well_id = row['well_id']
        curr_pool_id = df_assignment.loc[
            (df_assignment['plate_id'] == curr_plate_id) &
            (df_assignment['well_id'] == curr_well_id),'pool_id'].values[0]
        pool_ids.append(curr_pool_id)
    df_readout['pool_id'] = pool_ids
    deconvolution_result = run_ace_deconvolve(
        df_readout=df_readout,
        block_assignment=block_assignment,
        method=DeconvolutionMethods.EM,
        min_coverage=3,
        min_pool_spot_count=300
    )
    df_deconvolution = deconvolution_result.to_dataframe()
    expected_peptide_ids = sorted(['peptide_1', 'peptide_10'])
    hit_peptide_ids = sorted(df_deconvolution.loc[df_deconvolution['deconvolution_result'] == DeconvolutionLabels.CONFIDENT_HIT, 'peptide_id'].values.tolist())
    assert hit_peptide_ids == expected_peptide_ids, 'peptide_1 and peptide_10 are hits.'

