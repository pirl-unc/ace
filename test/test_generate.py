import pandas as pd
import pkg_resources
import torch
from .data import get_data_path
from acelib.main import run_ace_golfy, run_ace_sat_solver
from transformers import AutoTokenizer, AutoModelForMaskedLM
from acelib.sequence_features import AceNeuralEngine


def test_generate_small_golfy_configuration(small_golfy_elispot_configuration):
    print(len(small_golfy_elispot_configuration))


def test_generate_small_sat_solver_configuration(small_sat_solver_elispot_configuration):
    print(len(small_sat_solver_elispot_configuration))


def test_generate_large_golfy_configuration(large_golfy_elispot_configuration):
    print(len(large_golfy_elispot_configuration))


def test_generate_golfy_preferred_peptide_pairs():
    csv_file = get_data_path(name='25peptide_sequences.csv')
    df_peptides = pd.read_csv(csv_file)
    trained_model_file = pkg_resources.resource_filename('acelib', 'resources/models/trained_model3.pt')
    ESM2_TOKENIZER = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
    ESM2_MODEL = AutoModelForMaskedLM.from_pretrained("facebook/esm2_t6_8M_UR50D", return_dict=True, output_hidden_states=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    ace_eng = AceNeuralEngine(ESM2_MODEL, ESM2_TOKENIZER, device)
    ace_eng.load_weights(trained_model_file)
    preferred_peptide_pairs = ace_eng.find_paired_peptides(
        peptide_ids=df_peptides['peptide_id'].values.tolist(),
        peptide_sequences=df_peptides['peptide_sequence'].values.tolist(),
        sim_fxn='euclidean',
        threshold=0.6
    )
    is_valid, df_configuration = run_ace_golfy(
        df_peptides=df_peptides,
        num_peptides_per_pool=5,
        num_coverage=3,
        random_seed=1,
        max_iters=2000,
        init_mode='greedy',
        preferred_peptide_pairs=preferred_peptide_pairs
    )
    assert is_valid, 'With 2000 max iterations, we should have had a valid solution.'

