ace generate \
  --num-peptides 120 \
  --num-peptides-per-pool 12 \
  --num-coverage 3 \
  --num-processes 1 \
  --output-csv-file 120peptides_12perpool_3x.csv \
  --assign-well-ids 1 \
  --plate-type 96-well_plate

ace generate \
  --num-peptides 100 \
  --num-peptides-per-pool 10 \
  --num-coverage 3 \
  --num-processes 4 \
  --output-csv-file 100peptides_10perpool_3x.csv \
  --assign-well-ids 1 \
  --plate-type 96-well_plate
