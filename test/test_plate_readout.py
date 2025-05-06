from acelib.block_assignment import BlockAssignment
from acelib.plate_readout import PlateReadout
from .data import get_data_path


def test_plate_readout_1():
    block_assignment = BlockAssignment.read_excel_file(excel_file=get_data_path(name='25peptides_5perpool_3x_configuration.xlsx'))
    plate_readout = PlateReadout.read_aid_plate_reader_file(
        excel_file=get_data_path(name='25peptides_5perpool_3x_readout_aid-plate-reader.xlsx'),
        plate_id=1
    )
    plate_readout.assign_pool_ids(block_assignment=block_assignment)
    df_readout = plate_readout.to_dataframe()

    assert len(df_readout) == 15
    assert len(df_readout['pool_id'].unique()) == 15
    assert df_readout.loc[df_readout['well_id'] == 'A1','spot_count'].values[0] == 300
    assert df_readout.loc[df_readout['well_id'] == 'A3','spot_count'].values[0] == 300
    assert df_readout.loc[df_readout['well_id'] == 'A6','spot_count'].values[0] == 300
    assert df_readout.loc[df_readout['well_id'] == 'A8','spot_count'].values[0] == 300
    assert df_readout.loc[df_readout['well_id'] == 'B2','spot_count'].values[0] == 300
    assert df_readout.loc[df_readout['well_id'] == 'B3','spot_count'].values[0] == 300


def test_plate_readout_2():
    block_assignment = BlockAssignment.read_excel_file(excel_file=get_data_path(name='25peptides_5perpool_3x_configuration.xlsx'))
    plate_readout = PlateReadout.read_pool_id_file(
        file=get_data_path(name='25peptides_5perpool_3x_readout_pool-ids.xlsx')
    )
    plate_readout.assign_pool_ids(block_assignment=block_assignment)
    df_readout = plate_readout.to_dataframe()

    assert len(df_readout) == 15
    assert len(df_readout['pool_id'].unique()) == 15
    assert df_readout.loc[df_readout['well_id'] == 'A1','spot_count'].values[0] == 300
    assert df_readout.loc[df_readout['well_id'] == 'A3','spot_count'].values[0] == 300
    assert df_readout.loc[df_readout['well_id'] == 'A6','spot_count'].values[0] == 300
    assert df_readout.loc[df_readout['well_id'] == 'A8','spot_count'].values[0] == 300
    assert df_readout.loc[df_readout['well_id'] == 'B2','spot_count'].values[0] == 300
    assert df_readout.loc[df_readout['well_id'] == 'B3','spot_count'].values[0] == 300