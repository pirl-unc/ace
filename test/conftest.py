import pytest
from acelib.block_assignment import BlockAssignment
from acelib.main import run_ace_generate


@pytest.fixture
def golfy_assignment_25pep5per3x() -> BlockAssignment:
    peptides = []
    for i in range(1, 26):
        peptides.append(('peptide_%i' % i, ''))
    block_assignment, block_design = run_ace_generate(
        peptides=peptides,
        num_peptides_per_pool=5,
        num_coverage=3,
        trained_model_file='',
        cluster_peptides=False,
        mode='golfy',
        golfy_random_seed=1
    )
    block_assignment.is_optimal(
        num_coverage=3,
        num_peptides_per_pool=5
    )
    return block_assignment


@pytest.fixture
def sat_solver_assignment_25pep5per3x() -> BlockAssignment:
    peptides = []
    for i in range(1, 26):
        peptides.append(('peptide_%i' % i, ''))
    block_assignment, block_design = run_ace_generate(
        peptides=peptides,
        num_peptides_per_pool=5,
        num_coverage=3,
        cluster_peptides=False,
        sequence_similarity_threshold=0.7,
        sequence_similarity_function='euclidean',
        trained_model_file='',
        mode='cpsat_solver',
        golfy_random_seed=1
    )
    is_optimal = block_assignment.is_optimal(
        num_coverage=3,
        num_peptides_per_pool=5
    )
    assert is_optimal, 'We should have had an optimal solution.'
    return block_assignment


@pytest.fixture
def golfy_assignment_120pep12per3x() -> BlockAssignment:
    peptides = []
    for i in range(1, 121):
        peptides.append(('peptide_%i' % i, ''))
    block_assignment, block_design = run_ace_generate(
        peptides=peptides,
        num_peptides_per_pool=12,
        num_coverage=3,
        trained_model_file='',
        cluster_peptides=False,
        mode='golfy',
        golfy_random_seed=1
    )
    return block_assignment

