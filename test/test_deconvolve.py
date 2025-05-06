import pandas as pd
from acelib.block_assignment import BlockAssignment
from acelib.main import run_ace_deconvolve
from acelib.constants import DeconvolutionMethod
from acelib.plate_readout import PlateReadout
from .data import get_data_path


def test_deconvolve_1():
    block_assignment = BlockAssignment.read_excel_file(excel_file=get_data_path(name='25peptides_5perpool_3x_configuration.xlsx'))
    plate_readout = PlateReadout.read_aid_plate_reader_file(
        excel_file=get_data_path(name='25peptides_5perpool_3x_readout_aid-plate-reader.xlsx'),
        plate_id=1
    )
    plate_readout.assign_pool_ids(block_assignment=block_assignment)
    df_readout = plate_readout.to_dataframe()
    deconvolved_peptide_set = run_ace_deconvolve(
        df_readout=df_readout,
        block_assignment=block_assignment,
        method=DeconvolutionMethod.CONSTRAINED_EM,
        min_coverage=3,
        min_pool_spot_count=200
    )

    df_deconvolution = deconvolved_peptide_set.to_dataframe()

    hit_peptide_ids = df_deconvolution.loc[df_deconvolution['deconvolution_result'] != 'not_a_hit','peptide_id'].values.tolist()

    assert len(hit_peptide_ids) == 2
    assert 'peptide_1' in hit_peptide_ids
    assert 'peptide_10' in hit_peptide_ids
    assert df_deconvolution.loc[df_deconvolution['peptide_id'] == 'peptide_1','deconvolution_result'].values[0] == 'confident_hit'
    assert df_deconvolution.loc[df_deconvolution['peptide_id'] == 'peptide_10','deconvolution_result'].values[0] == 'confident_hit'


def test_deconvolve_2():
    block_assignment = BlockAssignment.read_excel_file(excel_file=get_data_path(name='25peptides_5perpool_3x_configuration.xlsx'))
    plate_readout = PlateReadout.read_aid_plate_reader_file(
        excel_file=get_data_path(name='25peptides_5perpool_3x_readout_aid-plate-reader.xlsx'),
        plate_id=1
    )
    plate_readout.assign_pool_ids(block_assignment=block_assignment)
    df_readout = plate_readout.to_dataframe()
    deconvolved_peptide_set = run_ace_deconvolve(
        df_readout=df_readout,
        block_assignment=block_assignment,
        method=DeconvolutionMethod.EMPIRICAL,
        min_coverage=3,
        min_pool_spot_count=200
    )

    df_deconvolution = deconvolved_peptide_set.to_dataframe()

    hit_peptide_ids = df_deconvolution.loc[df_deconvolution['deconvolution_result'] != 'not_a_hit','peptide_id'].values.tolist()

    assert len(hit_peptide_ids) == 2
    assert 'peptide_1' in hit_peptide_ids
    assert 'peptide_10' in hit_peptide_ids
    assert df_deconvolution.loc[df_deconvolution['peptide_id'] == 'peptide_1','deconvolution_result'].values[0] == 'confident_hit'
    assert df_deconvolution.loc[df_deconvolution['peptide_id'] == 'peptide_10','deconvolution_result'].values[0] == 'confident_hit'


def test_deconvolve_3():
    block_assignment = BlockAssignment.read_excel_file(excel_file=get_data_path(name='25peptides_5perpool_3x_configuration.xlsx'))
    plate_readout = PlateReadout.read_aid_plate_reader_file(
        excel_file=get_data_path(name='25peptides_5perpool_3x_readout_aid-plate-reader.xlsx'),
        plate_id=1
    )
    plate_readout.assign_pool_ids(block_assignment=block_assignment)
    df_readout = plate_readout.to_dataframe()
    deconvolved_peptide_set = run_ace_deconvolve(
        df_readout=df_readout,
        block_assignment=block_assignment,
        method=DeconvolutionMethod.EM,
        min_coverage=3,
        min_pool_spot_count=200
    )

    df_deconvolution = deconvolved_peptide_set.to_dataframe()

    hit_peptide_ids = df_deconvolution.loc[df_deconvolution['deconvolution_result'] != 'not_a_hit','peptide_id'].values.tolist()

    assert len(hit_peptide_ids) == 2
    assert 'peptide_1' in hit_peptide_ids
    assert 'peptide_10' in hit_peptide_ids
    assert df_deconvolution.loc[df_deconvolution['peptide_id'] == 'peptide_1','deconvolution_result'].values[0] == 'confident_hit'
    assert df_deconvolution.loc[df_deconvolution['peptide_id'] == 'peptide_10','deconvolution_result'].values[0] == 'confident_hit'


def test_deconvolve_4():
    block_assignment = BlockAssignment.read_excel_file(excel_file=get_data_path(name='25peptides_5perpool_3x_configuration.xlsx'))
    plate_readout = PlateReadout.read_aid_plate_reader_file(
        excel_file=get_data_path(name='25peptides_5perpool_3x_readout_aid-plate-reader.xlsx'),
        plate_id=1
    )
    plate_readout.assign_pool_ids(block_assignment=block_assignment)
    df_readout = plate_readout.to_dataframe()
    deconvolved_peptide_set = run_ace_deconvolve(
        df_readout=df_readout,
        block_assignment=block_assignment,
        method=DeconvolutionMethod.LASSO,
        min_coverage=3,
        min_pool_spot_count=200
    )

    df_deconvolution = deconvolved_peptide_set.to_dataframe()

    hit_peptide_ids = df_deconvolution.loc[df_deconvolution['deconvolution_result'] != 'not_a_hit','peptide_id'].values.tolist()

    assert len(hit_peptide_ids) == 2
    assert 'peptide_1' in hit_peptide_ids
    assert 'peptide_10' in hit_peptide_ids
    assert df_deconvolution.loc[df_deconvolution['peptide_id'] == 'peptide_1','deconvolution_result'].values[0] == 'confident_hit'
    assert df_deconvolution.loc[df_deconvolution['peptide_id'] == 'peptide_10','deconvolution_result'].values[0] == 'confident_hit'


def test_deconvolve_5():
    block_assignment = BlockAssignment.read_excel_file(excel_file=get_data_path(name='25peptides_5perpool_3x_configuration_without_coverage_ids.xlsx'))
    plate_readout = PlateReadout.read_aid_plate_reader_file(
        excel_file=get_data_path(name='25peptides_5perpool_3x_readout_aid-plate-reader.xlsx'),
        plate_id=1
    )
    plate_readout.assign_pool_ids(block_assignment=block_assignment)
    df_readout = plate_readout.to_dataframe()
    deconvolved_peptide_set = run_ace_deconvolve(
        df_readout=df_readout,
        block_assignment=block_assignment,
        method=DeconvolutionMethod.CONSTRAINED_EM,
        min_coverage=3,
        min_pool_spot_count=200
    )

    df_deconvolution = deconvolved_peptide_set.to_dataframe()

    hit_peptide_ids = df_deconvolution.loc[df_deconvolution['deconvolution_result'] != 'not_a_hit','peptide_id'].values.tolist()

    assert len(hit_peptide_ids) == 2
    assert 'peptide_1' in hit_peptide_ids
    assert 'peptide_10' in hit_peptide_ids
    assert df_deconvolution.loc[df_deconvolution['peptide_id'] == 'peptide_1','deconvolution_result'].values[0] == 'confident_hit'
    assert df_deconvolution.loc[df_deconvolution['peptide_id'] == 'peptide_10','deconvolution_result'].values[0] == 'confident_hit'



