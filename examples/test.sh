#ace generate \
#  --num-peptides 120 \
#  --num-peptides-per-pool 12 \
#  --num-coverage 3 \
#  --num-processes 6 \
#  --mode golfy \
#  --output-excel-file test_golfy.xlsx \
#  --assign-well-ids 1 \
#  --plate-type 96-well_plate

ace generate \
  --num-peptides 120 \
  --num-peptides-per-pool 12 \
  --num-coverage 3 \
  --num-processes 6 \
  --mode golfy \
  --golfy-allow-extra-pools False \
  --output-excel-file test_golfy.xlsx \
  --assign-well-ids 1 \
  --plate-type 96-well_plate

echo ""

ace generate \
  --num-peptides 120 \
  --num-peptides-per-pool 12 \
  --num-coverage 3 \
  --num-processes 6 \
  --mode sat_solver \
  --output-excel-file test_sat.xlsx \
  --shuffle-iters 1000 \
  --max-peptides-per-block 100 \
  --max-peptides-per-pool 10 \
  --assign-well-ids 1 \
  --plate-type 96-well_plate

#ace generate \
#  --num-peptides 100 \
#  --num-peptides-per-pool 5 \
#  --num-coverage 3 \
#  --num-processes 6 \
#  --mode sat_solver \
#  --output-excel-file test.xlsx \
#  --shuffle-iters 0 \
#  --max-peptides-per-block 100 \
#  --max-peptides-per-pool 10 \
#  --assign-well-ids 1 \
#  --plate-type 96-well_plate

