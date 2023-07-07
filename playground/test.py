import pandas as pd
from acelib.elispot import ELISpot
from acelib.main import run_ace_sat_solver


data = {
    'peptide_id': ['peptide_%i' % i for i in list(range(1, 26))],
    'peptide_sequence': [''] * 25
}
df_peptides = pd.DataFrame(data)

peptide_clusters = [
    ['peptide_12', 'peptide_14', 'peptide_15', 'peptide_17'],
    ['peptide_7', 'peptide_10']
]

df_configuration_first_coverage = ELISpot.generate_first_coverage_configuration(
    df_peptides=df_peptides,
    num_peptides_per_pool=5,
    peptide_clusters=peptide_clusters
)

disallowed_peptide_pairs = ELISpot.compute_disallowed_peptide_pairs(
    df_configuration=df_configuration_first_coverage
)

df_configuration = run_ace_sat_solver(
    df_peptides=df_peptides,
    num_peptides_per_pool=5,
    num_coverage=2,
    num_peptides_per_batch=100,
    random_seed=1,
    num_processes=4,
    is_first_coverage=False,
    disallowed_peptide_pairs=disallowed_peptide_pairs
)
df_configuration = pd.concat([df_configuration_first_coverage, df_configuration])
df_configuration.to_csv('test_sat_solver.csv', index=False)