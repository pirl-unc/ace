# Golfy with peptide sequences
ace generate \
  --peptides-csv-file ../test/data/25peptide_sequences.csv \
  --num-peptides-per-pool 5 \
  --num-coverage 3 \
  --num-processes 1 \
  --mode golfy \
  --sequence-similarity-threshold 0.7 \
  --output-csv-file 25peptides_5perpool_3x_golfy.csv
