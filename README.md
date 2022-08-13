# ACE

## 01. Installation
```
python setup.py install
```

## 02. Usage
### 02-1. Generate an ELIspot configuration.
```
ace_generate_elispot_configuration.py [-h] 
    --num_peptides NUM_PEPTIDES 
    --num_peptides_per_pool NUM_PEPTIDES_PER_POOL 
    --num_coverage NUM_COVERAGE 
    --output_tsv_file OUTPUT_TSV_FILE

optional arguments:
  -h, --help            show this help message and exit
  --num_peptides NUM_PEPTIDES
                        Total number of peptides.
  --num_peptides_per_pool NUM_PEPTIDES_PER_POOL
                        Number of peptides per pool. Please make sure this integer is a factor of the total number of peptides.
  --num_coverage NUM_COVERAGE
                        Total coverage (i.e. number of peptide replicates).
  --output_tsv_file OUTPUT_TSV_FILE
                        Output TSV file.
```