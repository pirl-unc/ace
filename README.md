# ACE

## 01. Installation
```
pip install . --verbose
```

## 02. Usage
### 02-1. Generate an ELIspot configuration.
```
usage: ace [-h] [--version] {generate,identify,verify} ...

ACE: Assay Configurator for ELIspot.

positional arguments:
  {generate,identify,verify}
                        ACE sub-commands.
    generate            Generates an ELIspot experiment configuration.
    identify            Identify hit peptide IDs given read-outs from an ELIspot experiment.
    verify              Verifies whether an ELIspot configuration satisfies all ACE constraints.

optional arguments:
  -h, --help            show this help message and exit
  --version, -v         show program's version number and exit
```

## 03. Generating a tar.gz release file
```
python setup.py sdist
```

## 04. Generating an eel program
pyinstaller generation
```ls
python -m eel ace_gui.py views --onedir --noconsole --windowed --clean --hidden-import=pytorch --collect-data torch --copy-metadata torch --copy-metadata tqdm --copy-metadata regex --copy-metadata requests --copy-metadata packaging --copy-metadata filelock --copy-metadata numpy --copy-metadata tokenizers
```