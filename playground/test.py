from acelib.elispot import ELIspot

elispot = ELIspot(
    num_peptides_per_pool=5,
    num_coverage=3,
    num_processes=4,
    peptide_ids=list(range(0, 25))
)

df_config = elispot.generate_configuration()