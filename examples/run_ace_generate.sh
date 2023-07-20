echo "Example 1. Golfy without peptide sequences (120/12/3x)"
ace generate \
  --num-peptides 120 \
  --num-peptides-per-pool 12 \
  --num-coverage 3 \
  --mode golfy \
  --output-excel-file 120peptides_12perpool_3x_noseqsim_golfy.xlsx \
  --assign-well-ids 1 \
  --plate-type 96-well_plate
echo ""

echo "Example 2. SAT solver without peptide sequences (120/12/3x)"
ace generate \
  --num-peptides 120 \
  --num-peptides-per-pool 12 \
  --num-coverage 3 \
  --num-processes 6 \
  --mode sat_solver \
  --output-excel-file 120peptides_12perpool_3x_noseqsim_sat-solver.xlsx \
  --assign-well-ids 1 \
  --plate-type 96-well_plate
echo ""

echo "Example 3. Golfy solver without peptide sequences (100/5/3x)"
ace generate \
  --num-peptides 100 \
  --num-peptides-per-pool 5 \
  --num-coverage 3 \
  --mode golfy \
  --output-excel-file 100peptides_5perpool_3x_noseqsim_golfy.xlsx \
  --assign-well-ids 1 \
  --plate-type 96-well_plate
echo ""

echo "Example 4. SAT solver without peptide sequences (100/5/3x)"
ace generate \
  --num-peptides 100 \
  --num-peptides-per-pool 5 \
  --num-coverage 3 \
  --num-processes 6 \
  --mode sat_solver \
  --output-excel-file 100peptides_5perpool_3x_noseqsim_sat-solver.xlsx \
  --assign-well-ids 1 \
  --plate-type 96-well_plate
echo ""

echo "Example 5. Golfy without peptide sequences (120/8/3x)"
ace generate \
  --num-peptides 120 \
  --num-peptides-per-pool 8 \
  --num-coverage 3 \
  --mode golfy \
  --output-excel-file 120peptides_8perpool_3x_noseqsim_golfy.xlsx \
  --assign-well-ids 1 \
  --plate-type 96-well_plate
echo ""

echo "Example 6. SAT solver without peptide sequences (120/8/3x)"
ace generate \
  --num-peptides 120 \
  --num-peptides-per-pool 8 \
  --num-coverage 3 \
  --num-processes 6 \
  --max-peptides-per-block 64 \
  --mode sat_solver \
  --output-excel-file 120peptides_8perpool_3x_noseqsim_sat-solver.xlsx \
  --assign-well-ids 1 \
  --plate-type 96-well_plate
echo ""

echo "Example 7. Golfy with peptide sequences (25/5/3x)"
ace generate \
  --peptides-excel-file ../test/data/25peptide_sequences.xlsx \
  --num-peptides-per-pool 5 \
  --num-coverage 3 \
  --num-processes 1 \
  --mode golfy \
  --sequence-similarity-threshold 0.7 \
  --output-excel-file 25peptides_5perpool_3x_seqsim_golfy.xlsx
echo ""

echo "Example 8. SAT solver with peptide sequences (25/5/3x)"
ace generate \
  --peptides-excel-file ../test/data/25peptide_sequences.xlsx \
  --num-peptides-per-pool 5 \
  --num-coverage 3 \
  --num-processes 6 \
  --mode sat_solver \
  --sequence-similarity-threshold 0.7 \
  --output-excel-file 25peptides_5perpool_3x_seqsim_sat-solver.xlsx
echo ""

echo "Example 9. Golfy without peptide sequences (90/9/3x)"
ace generate \
  --num-peptides 90 \
  --num-peptides-per-pool 9 \
  --num-coverage 3 \
  --mode golfy \
  --output-excel-file 90peptides_9perpool_3x_noseqsim_golfy.xlsx \
  --assign-well-ids 1 \
  --plate-type 96-well_plate
echo ""

echo "Example 10. SAT solver without peptide sequences (90/9/3x)"
ace generate \
  --num-peptides 90 \
  --num-peptides-per-pool 9 \
  --num-coverage 3 \
  --num-processes 6 \
  --mode sat_solver \
  --output-excel-file 90peptides_9perpool_3x_noseqsim_sat-solver.xlsx \
  --assign-well-ids 1 \
  --plate-type 96-well_plate
echo ""

echo "Example 11. Golfy without peptide sequences (220/11/3x)"
ace generate \
  --num-peptides 220 \
  --num-peptides-per-pool 11 \
  --num-coverage 3 \
  --mode golfy \
  --output-excel-file 220peptides_11perpool_3x_noseqsim_golfy.xlsx \
  --assign-well-ids 1 \
  --plate-type 96-well_plate
echo ""

echo "Example 12. SAT solver without peptide sequences (220/11/3x)"
ace generate \
  --num-peptides 220 \
  --num-peptides-per-pool 11 \
  --num-coverage 3 \
  --num-processes 6 \
  --mode sat_solver \
  --output-excel-file 220peptides_11perpool_3x_noseqsim_sat-solver.xlsx \
  --assign-well-ids 1 \
  --plate-type 96-well_plate
echo ""

echo "Example 13. Golfy without peptide sequences (240/12/3x)"
ace generate \
  --num-peptides 240 \
  --num-peptides-per-pool 12 \
  --num-coverage 3 \
  --mode golfy \
  --output-excel-file 240peptides_12perpool_3x_noseqsim_golfy.xlsx \
  --assign-well-ids 1 \
  --plate-type 96-well_plate
echo""

echo "Example 14. SAT solver without peptide sequences (240/12/3x)"
ace generate \
  --num-peptides 240 \
  --num-peptides-per-pool 12 \
  --num-coverage 3 \
  --num-processes 6 \
  --mode sat_solver \
  --output-excel-file 240peptides_12perpool_3x_noseqsim_sat-solver.xlsx \
  --assign-well-ids 1 \
  --plate-type 96-well_plate
echo""

