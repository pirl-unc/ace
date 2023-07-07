#ace generate \
#  --peptides-csv-file /Users/leework/Documents/Research/projects/project_ace/ace/test/data/25peptide_sequences.csv \
#  --num-peptides-per-pool 5 \
#  --num-coverage 3 \
#  --num-processes 1 \
#  --sequence-similarity-threshold 0.7 \
#  --output-csv-file test_golfy.csv \
#  --assign-well-ids 1 \
#  --plate-type 96-well_plate

ace generate \
  --peptides-csv-file /Users/leework/Documents/Research/projects/project_ace/ace/test/data/25peptide_sequences.csv \
  --num-peptides-per-pool 5 \
  --num-coverage 3 \
  --num-processes 1 \
  --mode sat_solver \
  --sequence-similarity-threshold 0.7 \
  --output-csv-file 25peptides_5perpool_3x_sat_solver.csv \
  --assign-well-ids 1 \
  --plate-type 96-well_plate

