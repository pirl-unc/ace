ace deconvolve \
  --readout-file-type pool_id \
  --readout-files 25peptides_5perpool_3x_pool-id_readout.xlsx \
  --assignment-excel-file 25peptides_5perpool_3x_seqsim_sat-solver.xlsx \
  --min-spot-count 300 \
  --output-excel-file 25peptides_5perpool_3x_pool-id_readout_deconvolved_empirical.xlsx

ace deconvolve \
  --readout-file-type pool_id \
  --readout-files 25peptides_5perpool_3x_pool-id_readout.xlsx \
  --assignment-excel-file 25peptides_5perpool_3x_seqsim_sat-solver.xlsx \
  --mode em \
  --output-excel-file 25peptides_5perpool_3x_pool-id_readout_deconvolved_em.xlsx

ace deconvolve \
  --readout-file-type pool_id \
  --readout-files 25peptides_5perpool_3x_pool-id_readout.xlsx \
  --assignment-excel-file 25peptides_5perpool_3x_seqsim_sat-solver.xlsx \
  --mode lasso \
  --output-excel-file 25peptides_5perpool_3x_pool-id_readout_deconvolved_lasso.xlsx
