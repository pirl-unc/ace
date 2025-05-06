from acelib.block_design import BlockDesign
from acelib.constants import GenerateMode, NumPlateWells
from acelib.main import run_ace_generate
from acelib.peptide import Peptide
from .data import get_data_path


def test_block_design_1():
    peptides = []
    for i in range(1, 26):
        peptide = Peptide(id='peptide_%i' % i, sequence='')
        peptides.append(peptide)

    block_assignment, block_design = run_ace_generate(
        peptides=peptides,
        num_peptides_per_pool=5,
        num_coverage=3,
        trained_model_file='',
        cluster_peptides=False,
        mode=GenerateMode.CPSAT_SOLVER,
        num_plate_wells=NumPlateWells.WELLS_96,
        cpsat_solver_max_peptides_per_block=25,
        verbose=False
    )

    df_design = block_design.metadata_dataframe

    assert set(block_assignment.coverage_ids) == {1, 2, 3}
    assert block_assignment.num_peptides == 25
    assert block_assignment.num_pools == 15
    assert block_assignment.num_violations == 0
    assert len(block_assignment.plate_map.keys()) == 15
    assert block_assignment.is_optimal(num_coverage=3, num_peptides_per_pool=5) == True

    assert len(block_design.all_peptide_ids) == len(peptides)
    assert block_design.num_peptides == len(peptides)
    assert block_design.num_total_peptides == len(peptides)
    assert df_design['num_peptides'].values[0] == len(peptides)
    assert df_design['num_peptides_per_pool'].values[0] == 5
    assert df_design['num_coverage'].values[0] == 3
    assert df_design['num_plate_wells'].values[0] == 96


def test_block_design_2():
    excel_file = get_data_path(name='25peptides_5perpool_3x_configuration.xlsx')

    block_design = BlockDesign.read_excel_file(excel_file=excel_file)

    df_design = block_design.metadata_dataframe

    assert len(block_design.all_peptide_ids) == 100 # maximum peptides per block was 100
    assert block_design.num_peptides == 25
    assert block_design.num_total_peptides == 100
    assert df_design['num_peptides'].values[0] == 25
    assert df_design['num_peptides_per_pool'].values[0] == 5
    assert df_design['num_coverage'].values[0] == 3
    assert df_design['num_plate_wells'].values[0] == 96


