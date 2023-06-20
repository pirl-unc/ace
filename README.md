# ACE

## 01. Installation
```
pip install . --verbose
```

## 02. Dependencies

## 03. Usage

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

## 04. Packaging
```
python -m build
```
