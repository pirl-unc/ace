# ACE

## 01. Installation
```
python setup.py install
```

## 02. Usage
### 02-1. Generate an ELIspot configuration.
```
ace [-h] 
    --num_peptides NUM_PEPTIDES 
    --num_peptides_per_pool NUM_PEPTIDES_PER_POOL 
    --num_coverage NUM_COVERAGE 
    --num_threads NUM_THREADS 
    --output_tsv_file OUTPUT_TSV_FILE
    [--output_pdf_file OUTPUT_PDF_FILE]

Generates an ELIspot configuration.

required arguments:
  --num_peptides NUM_PEPTIDES
                        Total number of peptides.
  --num_peptides_per_pool NUM_PEPTIDES_PER_POOL
                        Number of peptides per pool. Please make sure this integer is a factor of the total number of peptides.
  --num_coverage NUM_COVERAGE
                        Total coverage (i.e. number of peptide replicates).
  --num_threads NUM_THREADS
                        Number of threads to parallelize the computation (default: 2). Recommended: 8.
  --output_tsv_file OUTPUT_TSV_FILE
                        Output TSV file.

optional arguments:
  --output_pdf_file OUTPUT_PDF_FILE
                        Output PDF file.
```

## 03. Generating a tar.gz release file
```
python setup.py sdist
```

## 04. Generating an eel program
pyinstaller generation
```
python -m eel acegui/ace_gui.py views --noconsole --onedir

python -m eel ace_gui.py views --onedir --hidden-import=pytorch --collect-data torch --copy-metadata torch --copy-metadata tqdm --copy-metadata regex --copy-metadata requests --copy-metadata packaging --copy-metadata filelock --copy-metadata numpy --copy-metadata tokenizers

```