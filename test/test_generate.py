import pandas as pd
import pkg_resources
import torch
from .data import get_data_path
from acelib.block_assignment import BlockAssignment
from acelib.block_design import BlockDesign
from acelib.main import run_ace_golfy, run_ace_sat_solver, run_ace_generate
from acelib.utilities import convert_dataframe_to_peptides
from transformers import AutoTokenizer, AutoModelForMaskedLM
from acelib.sequence_features import AceNeuralEngine


def test_generate_golfy_assignment_1(golfy_assignment_1):
    print(golfy_assignment_1.num_pools)


def test_generate_sat_solver_assignment_1(sat_solver_assignment_1):
    print(sat_solver_assignment_1.num_pools)


def test_generate_golfy_assignment_2(golfy_assignment_2):
    print(golfy_assignment_2.num_pools)


def test_generate_golfy_assignment_with_preferred_peptide_pairs():
    # Step 1. Load peptide information
    excel_file = get_data_path(name='25peptide_sequences.xlsx')
    df_peptides = pd.read_excel(excel_file)
    peptides = convert_dataframe_to_peptides(df_peptides=df_peptides)

    # Step 2. Generate a block design
    block_assignment, block_design = run_ace_generate(
        peptides=peptides,
        num_peptides_per_pool=5,
        num_coverage=3,
        trained_model_file=pkg_resources.resource_filename('acelib', 'resources/models/trained_model5.pt'),
        mode='golfy',
        golfy_random_seed=1
    )
    block_assignment.is_optimal(
        num_coverage=3,
        num_peptides_per_pool=5
    )

    # Step 3. Check if the block assignment is optimal
    is_optimal = block_assignment.is_optimal(
        num_coverage=3,
        num_peptides_per_pool=12
    )


def test_generate_sat_solver_assignment_with_preferred_peptide_pairs():
    # Step 1. Load peptide information
    excel_file = get_data_path(name='25peptide_sequences.xlsx')
    df_peptides = pd.read_excel(excel_file)
    peptides = convert_dataframe_to_peptides(df_peptides=df_peptides)

    # Step 2. Generate a block design
    block_assignment, block_design = run_ace_generate(
        peptides=peptides,
        num_peptides_per_pool=5,
        num_coverage=3,
        trained_model_file=pkg_resources.resource_filename('acelib', 'resources/models/trained_model5.pt'),
        mode='cpsat_solver'
    )
    is_optimal = block_assignment.is_optimal(
        num_coverage=3,
        num_peptides_per_pool=5
    )

    assert is_optimal, 'We should have had an optimal solution.'

