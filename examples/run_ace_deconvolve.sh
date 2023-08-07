echo "Deconvolution Example 1: Empirical (Pool ID)."
ace deconvolve \
  --readout-file-type pool_id \
  --readout-files ../test/data/25peptides_5perpool_3x_pool-id_readout.xlsx \
  --assignment-excel-file ../test/data/25peptides_5perpool_3x_noseqsim_sat-solver.xlsx \
  --min-spot-count 300 \
  --output-excel-file outputs/25peptides_5perpool_3x_pool-id_readout_deconvolved_empirical.xlsx

echo "Deconvolution Example 2: Expectation Maximization (Pool ID)."
ace deconvolve \
  --readout-file-type pool_id \
  --readout-files ../test/data/25peptides_5perpool_3x_pool-id_readout.xlsx \
  --assignment-excel-file ../test/data/25peptides_5perpool_3x_noseqsim_sat-solver.xlsx \
  --mode em \
  --output-excel-file outputs/25peptides_5perpool_3x_pool-id_readout_deconvolved_em.xlsx

echo "Deconvolution Example 3: LASSO (Pool ID)."
ace deconvolve \
  --readout-file-type pool_id \
  --readout-files ../test/data/25peptides_5perpool_3x_pool-id_readout.xlsx \
  --assignment-excel-file ../test/data/25peptides_5perpool_3x_noseqsim_sat-solver.xlsx \
  --mode lasso \
  --output-excel-file outputs/25peptides_5perpool_3x_pool-id_readout_deconvolved_lasso.xlsx

echo "Deconvolution Example 4: Empirical (AID Plate Reader)."
ace deconvolve \
  --readout-file-type aid_plate_reader \
  --readout-files ../test/data/25peptides_5perpool_3x_aid-plate-reader_readout.xlsx \
  --assignment-excel-file ../test/data/25peptides_5perpool_3x_noseqsim_sat-solver.xlsx \
  --min-spot-count 300 \
  --output-excel-file outputs/25peptides_5perpool_3x_aid-plate-reader_readout_deconvolved_empirical.xlsx

echo "Deconvolution Example 5: Expectation Maximization (AID Plate Reader)."
ace deconvolve \
  --readout-file-type aid_plate_reader \
  --readout-files ../test/data/25peptides_5perpool_3x_aid-plate-reader_readout.xlsx \
  --assignment-excel-file ../test/data/25peptides_5perpool_3x_noseqsim_sat-solver.xlsx \
  --mode em \
  --output-excel-file outputs/25peptides_5perpool_3x_aid-plate-reader_readout_deconvolved_em.xlsx

echo "Deconvolution Example 6: LASSO (AID Plate Reader)."
ace deconvolve \
  --readout-file-type aid_plate_reader \
  --readout-files ../test/data/25peptides_5perpool_3x_aid-plate-reader_readout.xlsx \
  --assignment-excel-file ../test/data/25peptides_5perpool_3x_noseqsim_sat-solver.xlsx \
  --mode lasso \
  --output-excel-file outputs/25peptides_5perpool_3x_aid-plate-reader_readout_deconvolved_lasso.xlsx
