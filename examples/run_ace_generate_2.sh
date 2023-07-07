# SAT solver without peptide sequences
ace generate \
  --num-peptides 100 \
  --num-peptides-per-pool 5 \
  --num-coverage 3 \
  --num-processes 4 \
  --mode sat_solver \
  --output-csv-file 100peptides_5perpool_3x.csv \
  --assign-well-ids 1 \
  --plate-type 96-well_plate