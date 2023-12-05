echo "Deconvolution Example 1: Any Reader"
ace deconvolve \
  --readout-file-type pool_id \
  --readout-files ../test/data/25peptides_5perpool_3x_pool-id_readout.xlsx \
  --assignment-excel-file ../test/data/25peptides_5perpool_3x_configuration.xlsx \
  --min-positive-pool-spot-count 300 \
  --output-excel-file outputs/25peptides_5perpool_3x_pool-id_readout_deconvolved.xlsx

echo "Deconvolution Example 2: AID Plate Reader"
ace deconvolve \
  --readout-file-type aid_plate_reader \
  --readout-files ../test/data/25peptides_5perpool_3x_aid-plate-reader_readout.xlsx \
  --assignment-excel-file ../test/data/25peptides_5perpool_3x_configuration.xlsx \
  --min-positive-pool-spot-count 300 \
  --output-excel-file outputs/25peptides_5perpool_3x_aid-plate-reader_readout_deconvolved.xlsx
