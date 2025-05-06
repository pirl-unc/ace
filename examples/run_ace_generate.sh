ace generate \
  --peptides-file ../test/data/25peptide_sequences.xlsx \
  --num-peptides-per-pool 5 \
  --num-coverage 3 \
  --mode golfy \
  --output-excel-file outputs/25peptides_5perpool_3x_configuration.xlsx

ace generate \
  --peptides-file ../test/data/100peptide_sequences.xlsx \
  --num-peptides-per-pool 5 \
  --num-coverage 3 \
  --mode golfy \
  --output-excel-file outputs/100peptides_5perpool_3x_configuration.xlsx

ace generate \
  --peptides-file ../test/data/100peptide_sequences.xlsx \
  --num-peptides-per-pool 10 \
  --num-coverage 3 \
  --mode golfy \
  --output-excel-file outputs/100peptides_10perpool_3x_configuration.xlsx

ace generate \
  --peptides-file ../test/data/100peptide_sequences.xlsx \
  --num-peptides-per-pool 5 \
  --num-coverage 4 \
  --mode golfy \
  --output-excel-file outputs/100peptides_5perpool_4x_configuration.xlsx
