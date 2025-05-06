from acelib.constants import GenerateMode, NumPlateWells
from acelib.block_assignment import compute_transitive_neighbors, BlockAssignment
from acelib.main import run_ace_generate
from acelib.peptide import Peptide
from .data import get_data_path


def test_compute_transitive_neighbors_1():
    peptide_pairs = [
        ('peptide_1','peptide_2', 1.0),
        ('peptide_2','peptide_3', 1.0),
        ('peptide_4', 'peptide_5', 1.0),
        ('peptide_5', 'peptide_6', 1.0)
    ]
    peptide_neighbors = compute_transitive_neighbors(
        peptide_pairs=peptide_pairs
    )

    for cluster in peptide_neighbors:
        if 'peptide_1' in cluster:
            assert 'peptide_2' in cluster
            assert 'peptide_3' in cluster
        if 'peptide_4' in cluster:
            assert 'peptide_5' in cluster
            assert 'peptide_6' in cluster


def test_block_assignment_1():
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
        verbose=False
    )

    block_assignment.assign_well_ids(num_plate_wells=NumPlateWells.WELLS_96)

    assert set(block_assignment.coverage_ids) == {1, 2, 3}
    assert block_assignment.num_peptides == 25
    assert block_assignment.num_pools == 15
    assert block_assignment.num_violations == 0
    assert len(block_assignment.plate_map.keys()) == 15
    assert block_assignment.is_optimal(num_coverage=3, num_peptides_per_pool=5) == True


def test_block_assignment_2():
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
        verbose=False
    )

    df_assignment = block_assignment.to_dataframe()

    block_assignment = BlockAssignment.load_from_dataframe(df_assignment=df_assignment)

    assert set(block_assignment.coverage_ids) == {1, 2, 3}
    assert block_assignment.num_peptides == 25
    assert block_assignment.num_pools == 15
    assert block_assignment.num_violations == 0
    assert len(block_assignment.plate_map.keys()) == 15
    assert block_assignment.is_optimal(num_coverage=3, num_peptides_per_pool=5) == True


def test_block_assignment_3():
    excel_file = get_data_path(name='25peptides_5perpool_3x_configuration.xlsx')

    block_assignment = BlockAssignment.read_excel_file(excel_file=excel_file)

    df_assignment = block_assignment.to_dataframe()

    assert set(block_assignment.coverage_ids) == {1, 2, 3}
    assert block_assignment.num_peptides == 25
    assert block_assignment.num_pools == 15
    assert block_assignment.num_violations == 0
    assert len(block_assignment.plate_map.keys()) == 15
    assert block_assignment.is_optimal(num_coverage=3, num_peptides_per_pool=5) == True
    assert all(len(group['pool_id'].unique()) == 3 for _, group in df_assignment.groupby('peptide_id'))
    assert all(len(group['coverage_id'].unique()) == 3 for _, group in df_assignment.groupby('peptide_id'))


def test_block_assignment_4():
    excel_file = get_data_path(name='25peptides_5perpool_3x_configuration.xlsx')

    block_assignment = BlockAssignment.read_excel_file(excel_file=excel_file)

    golfy_design = block_assignment.to_golfy_design()

    assert golfy_design[0].num_peptides == 25
    assert golfy_design[0].max_peptides_per_pool == 5
    assert golfy_design[0].num_replicates == 3
    assert len(golfy_design[0].assignments.keys()) == 3
    assert len(golfy_design[0].assignments[0].keys()) == 5
    assert len(golfy_design[0].assignments[1].keys()) == 5
    assert len(golfy_design[0].assignments[2].keys()) == 5


def test_block_assignment_5():
    excel_file = get_data_path(name='25peptides_5perpool_3x_configuration_without_coverage_ids.xlsx')

    block_assignment = BlockAssignment.read_excel_file(excel_file=excel_file)

    golfy_design = block_assignment.to_golfy_design()

    assert golfy_design[0].num_peptides == 25
    assert golfy_design[0].max_peptides_per_pool == 5
    assert golfy_design[0].num_replicates == 3
    assert len(golfy_design[0].assignments.keys()) == 3
    assert len(golfy_design[0].assignments[0].keys()) == 5
    assert len(golfy_design[0].assignments[1].keys()) == 5
    assert len(golfy_design[0].assignments[2].keys()) == 5
