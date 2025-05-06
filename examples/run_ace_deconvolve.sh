ace deconvolve \
  --readout-file-type pool_id \
  --readout-files ../test/data/25peptides_5perpool_3x_readout_pool-ids.xlsx \
  --assignment-excel-file ../test/data/25peptides_5perpool_3x_configuration.xlsx \
  --min-pool-spot-count 300 \
  --output-excel-file outputs/25peptides_5perpool_3x_readout_pool-ids_deconvolved.xlsx

ace deconvolve \
  --readout-file-type aid_plate_reader \
  --readout-files ../test/data/25peptides_5perpool_3x_readout_aid-plate-reader.xlsx \
  --assignment-excel-file ../test/data/25peptides_5perpool_3x_configuration.xlsx \
  --min-pool-spot-count 300 \
  --output-excel-file outputs/25peptides_5perpool_3x_readout_aid-plate-reader_deconvolved.xlsx
