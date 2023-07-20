import pandas as pd
import pytest
from .data import get_data_path
from acelib.block_assignment import BlockAssignment
from acelib.block_design import BlockDesign
from acelib.main import run_ace_golfy, run_ace_sat_solver
from golfy import init, is_valid, optimize


@pytest.fixture
def golfy_assignment_1() -> BlockAssignment:
    peptides = []
    for i in range(1, 26):
        peptides.append(('peptide_%i' % i, ''))
    block_design = BlockDesign(
        peptides=peptides,
        num_peptides_per_pool=5,
        num_coverage=3,
        max_peptides_per_block=25
    )
    block_assignment = run_ace_golfy(
        block_design=block_design,
        random_seed=1,
        max_iters=2000,
        init_mode='greedy'
    )
    is_optimal = block_assignment.is_optimal(
        num_coverage=3,
        num_peptides_per_pool=5
    )
    return block_assignment


@pytest.fixture
def sat_solver_assignment_1() -> BlockAssignment:
    peptides = []
    for i in range(1, 26):
        peptides.append(('peptide_%i' % i, ''))
    block_design = BlockDesign(
        peptides=peptides,
        num_peptides_per_pool=5,
        num_coverage=3,
        max_peptides_per_block=25
    )
    block_assignment = run_ace_sat_solver(
        block_design=block_design,
        max_peptides_per_pool=10,
        num_processes=1
    )
    is_optimal = block_assignment.is_optimal(
        num_coverage=3,
        num_peptides_per_pool=5
    )
    assert is_optimal, 'We should have had an optimal solution.'
    return block_assignment


@pytest.fixture
def golfy_assignment_2() -> BlockAssignment:
    peptides = []
    for i in range(1, 121):
        peptides.append(('peptide_%i' % i, ''))
    block_design = BlockDesign(
        peptides=peptides,
        num_peptides_per_pool=12,
        num_coverage=3,
        max_peptides_per_block=100
    )
    block_assignment = run_ace_golfy(
        block_design=block_design,
        random_seed=1,
        max_iters=2000,
        init_mode='greedy'
    )
    is_optimal = block_assignment.is_optimal(
        num_coverage=3,
        num_peptides_per_pool=12
    )
    return block_assignment

