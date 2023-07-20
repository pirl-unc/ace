import pandas as pd
import pkg_resources
import torch
from .data import get_data_path
from acelib.block_assignment import BlockAssignment
from acelib.block_design import BlockDesign
from acelib.main import run_ace_golfy, run_ace_sat_solver
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

    # Step 2. Load sequence similarity model
    trained_model_file = pkg_resources.resource_filename('acelib', 'resources/models/seq_sim_trained_model.pt')
    ESM2_TOKENIZER = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
    ESM2_MODEL = AutoModelForMaskedLM.from_pretrained("facebook/esm2_t6_8M_UR50D", return_dict=True, output_hidden_states=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    ace_eng = AceNeuralEngine(ESM2_MODEL, ESM2_TOKENIZER, device)
    ace_eng.load_weights(trained_model_file)

    # Step 3. Predict preferred peptide pairs
    preferred_peptide_pairs = ace_eng.find_paired_peptides(
        peptide_ids=[p[0] for p in peptides],
        peptide_sequences=[p[1] for p in peptides],
        sim_fxn='euclidean',
        threshold=0.7
    )
    preferred_peptide_pairs = [(p1, p2) for p1, p2, score in preferred_peptide_pairs]

    # Step 4. Generate a block design
    block_design = BlockDesign(
        peptides=peptides,
        num_peptides_per_pool=5,
        num_coverage=3,
        max_peptides_per_block=25,
        preferred_peptide_pairs=preferred_peptide_pairs
    )

    # Step 5. Generate a block assignment
    block_assignment = run_ace_golfy(
        block_design=block_design,
        random_seed=1,
        max_iters=2000,
        init_mode='greedy'
    )

    # Step 6. Check if the block assignment is optimal
    is_optimal = block_assignment.is_optimal(
        num_coverage=3,
        num_peptides_per_pool=12
    )


def test_generate_sat_solver_assignment_with_preferred_peptide_pairs():
    # Step 1. Load peptide information
    excel_file = get_data_path(name='25peptide_sequences.xlsx')
    df_peptides = pd.read_excel(excel_file)
    peptides = convert_dataframe_to_peptides(df_peptides=df_peptides)

    # Step 2. Load sequence similarity model
    trained_model_file = pkg_resources.resource_filename('acelib', 'resources/models/seq_sim_trained_model.pt')
    ESM2_TOKENIZER = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
    ESM2_MODEL = AutoModelForMaskedLM.from_pretrained("facebook/esm2_t6_8M_UR50D", return_dict=True, output_hidden_states=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    ace_eng = AceNeuralEngine(ESM2_MODEL, ESM2_TOKENIZER, device)
    ace_eng.load_weights(trained_model_file)

    # Step 3. Predict preferred peptide pairs
    preferred_peptide_pairs = ace_eng.find_paired_peptides(
        peptide_ids=[p[0] for p in peptides],
        peptide_sequences=[p[1] for p in peptides],
        sim_fxn='euclidean',
        threshold=0.7
    )
    preferred_peptide_pairs = [(p1, p2) for p1, p2, score in preferred_peptide_pairs]

    # Step 4. Generate a block design
    block_design = BlockDesign(
        peptides=peptides,
        num_peptides_per_pool=5,
        num_coverage=3,
        max_peptides_per_block=25,
        preferred_peptide_pairs=preferred_peptide_pairs
    )

    # Step 5. Run SAT solver
    block_assignment = run_ace_sat_solver(
        block_design=block_design,
        max_peptides_per_pool=10,
        num_processes=1
    )

    # Step 6. Check if the block assignment is optimal
    is_optimal = block_assignment.is_optimal(
        num_coverage=3,
        num_peptides_per_pool=5
    )
    assert is_optimal, 'We should have had an optimal solution.'

