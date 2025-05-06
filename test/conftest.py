# import pytest
# from acelib.block_assignment import BlockAssignment
# from acelib.constants import GenerateMode
# from acelib.main import run_ace_generate
# from acelib.peptide import Peptide
#
# from playground.test import df_assignment


# @pytest.fixture
# def golfy_assignment_25pep5per3x() -> BlockAssignment:
#     peptides = []
#     for i in range(1, 26):
#         peptide = Peptide(id='peptide_%i' % i, sequence='')
#         peptides.append(peptide)
#
#     block_assignment, block_design = run_ace_generate(
#         peptides=peptides,
#         num_peptides_per_pool=5,
#         num_coverage=3,
#         trained_model_file='',
#         cluster_peptides=False,
#         mode=GenerateMode.GOLFY,
#         golfy_random_seed=1
#     )
#
#     df_assignment = block_assignment.to_dataframe()
#
#     assert len(df_assignment['peptide_id'].unique()) == 25
#     assert len(df_assignment['pool_id'].unique()) == 15
#
#     return block_assignment


# @pytest.fixture
# def sat_solver_assignment_25pep5per3x() -> BlockAssignment:
#     peptides = []
#     for i in range(1, 26):
#         peptides.append(('peptide_%i' % i, ''))
#     block_assignment, block_design = run_ace_generate(
#         peptides=peptides,
#         num_peptides_per_pool=5,
#         num_coverage=3,
#         cluster_peptides=False,
#         sequence_similarity_threshold=0.7,
#         sequence_similarity_function='euclidean',
#         trained_model_file='',
#         mode='cpsat_solver',
#         golfy_random_seed=1
#     )
#     is_optimal = block_assignment.is_optimal(
#         num_coverage=3,
#         num_peptides_per_pool=5
#     )
#     assert is_optimal, 'We should have had an optimal solution.'
#     return block_assignment
#
#
# @pytest.fixture
# def golfy_assignment_120pep12per3x() -> BlockAssignment:
#     peptides = []
#     for i in range(1, 121):
#         peptides.append(('peptide_%i' % i, ''))
#     block_assignment, block_design = run_ace_generate(
#         peptides=peptides,
#         num_peptides_per_pool=12,
#         num_coverage=3,
#         trained_model_file='',
#         cluster_peptides=False,
#         mode='golfy',
#         golfy_random_seed=1
#     )
#     return block_assignment
#
