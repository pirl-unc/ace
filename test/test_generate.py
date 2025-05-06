import pandas as pd
from acelib.constants import GenerateMode
from acelib.main import run_ace_generate
from acelib.peptide import Peptide
from acelib.utilities import convert_dataframe_to_peptides
from importlib import resources
from .data import get_data_path


def test_generate_golfy_1():
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
        mode=GenerateMode.GOLFY,
        golfy_random_seed=1,
        verbose=False
    )

    df_assignment = block_assignment.to_dataframe()

    assert len(df_assignment['peptide_id'].unique()) == 25
    assert len(df_assignment['pool_id'].unique()) == 15
    assert len(df_assignment['coverage_id'].unique()) == 3
    assert all(len(group['pool_id'].unique()) == 3 for _, group in df_assignment.groupby('peptide_id'))


def test_generate_golfy_2():
    peptides = []
    for i in range(1, 26):
        peptide = Peptide(id='peptide_%i' % i, sequence='')
        peptides.append(peptide)

    block_assignment, block_design = run_ace_generate(
        peptides=peptides,
        num_peptides_per_pool=5,
        num_coverage=4,
        trained_model_file='',
        cluster_peptides=False,
        mode=GenerateMode.GOLFY,
        golfy_random_seed=1,
        verbose=False
    )

    df_assignment = block_assignment.to_dataframe()

    assert len(df_assignment['peptide_id'].unique()) == 25
    assert len(df_assignment['pool_id'].unique()) == 20
    assert len(df_assignment['coverage_id'].unique()) == 4
    assert all(len(group['pool_id'].unique()) == 4 for _, group in df_assignment.groupby('peptide_id'))


def test_generate_golfy_3():
    peptides = []
    for i in range(1, 121):
        peptide = Peptide(id='peptide_%i' % i, sequence='')
        peptides.append(peptide)

    block_assignment, block_design = run_ace_generate(
        peptides=peptides,
        num_peptides_per_pool=12,
        num_coverage=3,
        trained_model_file='',
        cluster_peptides=False,
        mode=GenerateMode.GOLFY,
        golfy_random_seed=1,
        verbose=False
    )

    df_assignment = block_assignment.to_dataframe()

    assert len(df_assignment['peptide_id'].unique()) == 120
    assert len(df_assignment['pool_id'].unique()) == 30
    assert len(df_assignment['coverage_id'].unique()) == 3
    assert all(len(group['pool_id'].unique()) == 3 for _, group in df_assignment.groupby('peptide_id'))


def test_generate_golfy_4():
    peptides = []
    for i in range(1, 121):
        peptide = Peptide(id='peptide_%i' % i, sequence='')
        peptides.append(peptide)

    block_assignment, block_design = run_ace_generate(
        peptides=peptides,
        num_peptides_per_pool=12,
        num_coverage=4,
        trained_model_file='',
        cluster_peptides=False,
        mode=GenerateMode.GOLFY,
        golfy_random_seed=1,
        verbose=False
    )

    df_assignment = block_assignment.to_dataframe()

    assert len(df_assignment['peptide_id'].unique()) == 120
    assert len(df_assignment['pool_id'].unique()) == 40
    assert len(df_assignment['coverage_id'].unique()) == 4
    assert all(len(group['pool_id'].unique()) == 4 for _, group in df_assignment.groupby('peptide_id'))


def test_generate_golfy_5():
    peptides = []
    for i in range(1, 401):
        peptide = Peptide(id='peptide_%i' % i, sequence='')
        peptides.append(peptide)

    block_assignment, block_design = run_ace_generate(
        peptides=peptides,
        num_peptides_per_pool=20,
        num_coverage=3,
        trained_model_file='',
        cluster_peptides=False,
        mode=GenerateMode.GOLFY,
        golfy_random_seed=1,
        verbose=False
    )

    df_assignment = block_assignment.to_dataframe()

    assert len(df_assignment['peptide_id'].unique()) == 400
    assert len(df_assignment['pool_id'].unique()) == 60
    assert len(df_assignment['coverage_id'].unique()) == 3
    assert all(len(group['pool_id'].unique()) == 3 for _, group in df_assignment.groupby('peptide_id'))


def test_generate_golfy_6():
    peptides = []
    for i in range(1, 401):
        peptide = Peptide(id='peptide_%i' % i, sequence='')
        peptides.append(peptide)

    block_assignment, block_design = run_ace_generate(
        peptides=peptides,
        num_peptides_per_pool=20,
        num_coverage=4,
        trained_model_file='',
        cluster_peptides=False,
        mode=GenerateMode.GOLFY,
        golfy_random_seed=1,
        verbose=False
    )

    df_assignment = block_assignment.to_dataframe()

    assert len(df_assignment['peptide_id'].unique()) == 400
    assert len(df_assignment['pool_id'].unique()) == 80
    assert len(df_assignment['coverage_id'].unique()) == 4
    assert all(len(group['pool_id'].unique()) == 4 for _, group in df_assignment.groupby('peptide_id'))


def test_generate_golfy_7():
    excel_file = get_data_path(name='25peptide_sequences.xlsx')
    df_peptides = pd.read_excel(excel_file)
    peptides = convert_dataframe_to_peptides(df_peptides=df_peptides)

    block_assignment, block_design = run_ace_generate(
        peptides=peptides,
        num_peptides_per_pool=5,
        num_coverage=3,
        trained_model_file=str(resources.path('acelib.resources.models', 'trained_model_w_data_augmentation_b3000.pt')),
        cluster_peptides=True,
        mode=GenerateMode.GOLFY,
        golfy_random_seed=1,
        verbose=False
    )

    df_assignment = block_assignment.to_dataframe()

    assert len(df_assignment['peptide_id'].unique()) == 25
    assert len(df_assignment['pool_id'].unique()) == 15
    assert len(df_assignment['coverage_id'].unique()) == 3
    assert all(len(group['pool_id'].unique()) == 3 for _, group in df_assignment.groupby('peptide_id'))


def test_generate_golfy_8():
    excel_file = get_data_path(name='100peptide_sequences.xlsx')
    df_peptides = pd.read_excel(excel_file)
    peptides = convert_dataframe_to_peptides(df_peptides=df_peptides)

    block_assignment, block_design = run_ace_generate(
        peptides=peptides,
        num_peptides_per_pool=10,
        num_coverage=3,
        trained_model_file=str(resources.path('acelib.resources.models', 'trained_model_w_data_augmentation_b3000.pt')),
        cluster_peptides=True,
        mode=GenerateMode.GOLFY,
        golfy_random_seed=1,
        verbose=False
    )

    df_assignment = block_assignment.to_dataframe()

    assert len(df_assignment['peptide_id'].unique()) == 100
    assert len(df_assignment['pool_id'].unique()) == 30
    assert len(df_assignment['coverage_id'].unique()) == 3
    assert all(len(group['pool_id'].unique()) == 3 for _, group in df_assignment.groupby('peptide_id'))


def test_generate_cpsat_1():
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

    assert block_assignment.is_optimal(num_coverage=3, num_peptides_per_pool=5) == True
    assert len(df_assignment['peptide_id'].unique()) == 25
    assert len(df_assignment['pool_id'].unique()) == 15
    assert len(df_assignment['coverage_id'].unique()) == 3
    assert all(len(group['pool_id'].unique()) == 3 for _, group in df_assignment.groupby('peptide_id'))


def test_generate_cpsat_2():
    peptides = []
    for i in range(1, 101):
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

    assert block_assignment.is_optimal(num_coverage=3, num_peptides_per_pool=5) == True
    assert len(df_assignment['peptide_id'].unique()) == 100
    assert len(df_assignment['pool_id'].unique()) == 60
    assert len(df_assignment['coverage_id'].unique()) == 3
    assert all(len(group['pool_id'].unique()) == 3 for _, group in df_assignment.groupby('peptide_id'))


def test_generate_cpsat_3():
    excel_file = get_data_path(name='25peptide_sequences.xlsx')
    df_peptides = pd.read_excel(excel_file)
    peptides = convert_dataframe_to_peptides(df_peptides=df_peptides)

    block_assignment, block_design = run_ace_generate(
        peptides=peptides,
        num_peptides_per_pool=5,
        num_coverage=3,
        trained_model_file=str(resources.path('acelib.resources.models', 'trained_model_w_data_augmentation_b3000.pt')),
        cluster_peptides=True,
        mode=GenerateMode.CPSAT_SOLVER,
        golfy_random_seed=1,
        verbose=False
    )

    df_assignment = block_assignment.to_dataframe()

    assert block_assignment.is_optimal(num_coverage=3, num_peptides_per_pool=5) == True
    assert len(df_assignment['peptide_id'].unique()) == 25
    assert len(df_assignment['pool_id'].unique()) == 15
    assert len(df_assignment['coverage_id'].unique()) == 3
    assert all(len(group['pool_id'].unique()) == 3 for _, group in df_assignment.groupby('peptide_id'))
