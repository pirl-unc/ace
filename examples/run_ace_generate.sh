# Golfy
echo "Golfy Example 1: 90 Peptides, 9 Peptides per Pool, 3x Coverage."
ace generate \
  --num-peptides 90 \
  --num-peptides-per-pool 9 \
  --num-coverage 3 \
  --mode golfy \
  --output-excel-file outputs/90peptides_9perpool_3x_noseqsim_golfy.xlsx \
  --assign-well-ids 1 \
  --num-plate-wells 96 \

echo "Golfy Example 2: 100 Peptides, 5 Peptides per Pool, 3x Coverage."
ace generate \
  --num-peptides 100 \
  --num-peptides-per-pool 5 \
  --num-coverage 3 \
  --mode golfy \
  --output-excel-file outputs/100peptides_5perpool_3x_noseqsim_golfy.xlsx \
  --assign-well-ids 1 \
  --num-plate-wells 96 \

echo "Golfy Example 3: 120 Peptides, 8 Peptides per Pool, 3x Coverage."
ace generate \
  --num-peptides 120 \
  --num-peptides-per-pool 8 \
  --num-coverage 3 \
  --mode golfy \
  --output-excel-file outputs/120peptides_8perpool_3x_noseqsim_golfy.xlsx \
  --assign-well-ids 1 \
  --num-plate-wells 96

echo "Golfy Example 4: 120 Peptides, 12 Peptides per Pool, 3x Coverage."
ace generate \
  --num-peptides 120 \
  --num-peptides-per-pool 12 \
  --num-coverage 3 \
  --mode golfy \
  --output-excel-file outputs/120peptides_12perpool_3x_noseqsim_golfy.xlsx \
  --assign-well-ids 1 \
  --num-plate-wells 96

echo "Golfy Example 5: 220 Peptides, 11 Peptides per Pool, 3x Coverage."
ace generate \
  --num-peptides 220 \
  --num-peptides-per-pool 11 \
  --num-coverage 3 \
  --mode golfy \
  --output-excel-file outputs/220peptides_11perpool_3x_noseqsim_golfy.xlsx \
  --assign-well-ids 1 \
  --num-plate-wells 96

echo "Golfy Example 6: 240 Peptides, 12 Peptides per Pool, 3x Coverage."
ace generate \
  --num-peptides 240 \
  --num-peptides-per-pool 12 \
  --num-coverage 3 \
  --mode golfy \
  --output-excel-file outputs/240peptides_12perpool_3x_noseqsim_golfy.xlsx \
  --assign-well-ids 1 \
  --num-plate-wells 96

echo "Golfy Example 7: 400 Peptides, 5 Peptides per Pool, 3x Coverage."
ace generate \
  --num-peptides 400 \
  --num-peptides-per-pool 5 \
  --num-coverage 3 \
  --mode golfy \
  --output-excel-file outputs/240peptides_12perpool_3x_noseqsim_golfy_384plates.xlsx \
  --assign-well-ids 1 \
  --num-plate-wells 384

echo "Golfy Example 8: 25 Peptides, 5 Peptides per Pool, 3x Coverage (with Sequences)."
ace generate \
  --peptides-excel-file ../test/data/25peptide_sequences.xlsx \
  --num-peptides-per-pool 5 \
  --num-coverage 3 \
  --mode golfy \
  --output-excel-file outputs/25peptides_5perpool_3x_seqsim_golfy.xlsx

# CP-SAT Solver
echo "CP-SAT Solver Example 1: 90 Peptides, 9 Peptides per Pool, 3x Coverage."
ace generate \
  --num-peptides 90 \
  --num-peptides-per-pool 9 \
  --num-coverage 3 \
  --cpsat-solver-num-processes 6 \
  --mode cpsat_solver \
  --output-excel-file outputs/90peptides_9perpool_3x_noseqsim_sat-solver.xlsx \
  --assign-well-ids 1 \
  --num-plate-wells 96

echo "CP-SAT Solver Example 2: 100 Peptides, 5 Peptides per Pool, 3x Coverage."
ace generate \
  --num-peptides 100 \
  --num-peptides-per-pool 5 \
  --num-coverage 3 \
  --cpsat-solver-num-processes 6 \
  --mode cpsat_solver \
  --output-excel-file outputs/100peptides_5perpool_3x_noseqsim_sat-solver.xlsx \
  --assign-well-ids 1 \
  --num-plate-wells 96

echo "CP-SAT Solver Example 3: 120 Peptides, 12 Peptides per Pool, 3x Coverage."
ace generate \
  --num-peptides 120 \
  --num-peptides-per-pool 12 \
  --num-coverage 3 \
  --cpsat-solver-num-processes 6 \
  --mode cpsat_solver \
  --output-excel-file outputs/120peptides_12perpool_3x_noseqsim_sat-solver.xlsx \
  --assign-well-ids 1 \
  --num-plate-wells 96

echo "CP-SAT Solver Example 4: 120 Peptides, 8 Peptides per Pool, 3x Coverage."
ace generate \
  --num-peptides 120 \
  --num-peptides-per-pool 8 \
  --num-coverage 3 \
  --cpsat-solver-num-processes 6 \
  --cpsat-solver-max-peptides-per-block 64 \
  --mode cpsat_solver \
  --output-excel-file outputs/120peptides_8perpool_3x_noseqsim_sat-solver.xlsx \
  --assign-well-ids 1 \
  --num-plate-wells 96

echo "CP-SAT Solver Example 5: 220 Peptides, 11 Peptides per Pool, 3x Coverage."
ace generate \
  --num-peptides 220 \
  --num-peptides-per-pool 11 \
  --num-coverage 3 \
  --cpsat-solver-num-processes 6 \
  --mode cpsat_solver \
  --output-excel-file outputs/220peptides_11perpool_3x_noseqsim_sat-solver.xlsx \
  --assign-well-ids 1 \
  --num-plate-wells 96

echo "CP-SAT Solver Example 6: 240 Peptides, 12 Peptides per Pool, 3x Coverage."
ace generate \
  --num-peptides 240 \
  --num-peptides-per-pool 12 \
  --num-coverage 3 \
  --cpsat-solver-num-processes 6 \
  --mode cpsat_solver \
  --output-excel-file outputs/240peptides_12perpool_3x_noseqsim_sat-solver.xlsx \
  --assign-well-ids 1 \
  --num-plate-wells 96
echo""

echo "CP-SAT Solver Example 7: 25 Peptides, 5 Peptides per Pool, 3x Coverage (with Sequences)."
ace generate \
  --peptides-excel-file ../test/data/25peptide_sequences.xlsx \
  --num-peptides-per-pool 5 \
  --num-coverage 3 \
  --cpsat-solver-num-processes 6 \
  --mode cpsat_solver \
  --output-excel-file outputs/25peptides_5perpool_3x_seqsim_sat-solver.xlsx
