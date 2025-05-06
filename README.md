# ACE Configurator for ELISpot

ACE facilitates (1) generation of ELISpot configurations (peptide-pool assignments) 
using a deep learning approach to cluster similar peptides and (2) deconvolution 
of pool spot counts for identification of immunogenic peptides.

[![build](https://github.com/pirl-unc/ace/actions/workflows/main.yml/badge.svg?branch=main)](https://github.com/pirl-unc/ace/actions/workflows/main.yml)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

## 01. Installation

### 01-1. Standalone Graphical User Interface (Recommended)

Please note that the ACE GUI software will take a long time to load (~30 seconds :coffee:). 
We also recommend that you have [Google Chrome](https://www.google.com/chrome/) installed on your machine.

| Operating System | Link                                                                                                     | Version           |
|------------------|----------------------------------------------------------------------------------------------------------|-------------------|
| Mac              | [Download](https://github.com/pirl-unc/ace/releases/download/v0.1.2.0/ace-elispot-0.1.2.0-mac.zip)       | v0.1.2.0 (latest) | 
| Windows 10       | [Download](https://github.com/pirl-unc/ace/releases/download/v0.1.2.0/ace-elispot-0.1.2.0-windows10.zip) | v0.1.2.0 (latest) |

For Windows versions, first unzip the file and look for an application file called `ACE` inside the unzipped foler.
Previous versions of ACE are available [here](https://github.com/pirl-unc/ace/releases).

### 01-2. Python Package 

ACE is available on [PyPI](https://pypi.org/project/ace-elispot/)

```
pip install ace-elispot
```

You can also download a specific version of ACE from [here](https://github.com/pirl-unc/ace/releases).<br/>
Subsequently install the ACE package using pip:

```
pip install ace-elispot-<version>.tar.gz
```

#### Dependencies

- python==3.10.13
- pip==25.0
- numpy==1.26.4
- pandas==2.2.3
- scipy==1.15.2
- scikit-learn==1.6.1
- openpyxl=3.1.5
- setuptools==72.1.0
- pyinstaller==6.13.0
- transformers=4.30.2
- pillow==11.1.0
- eel==0.18.1
- golfy==2.5.0
- levenshtein==0.27.1
- networkx==3.4.2
- ortools==9.8.3296
- torch==2.2.2

## 02. Usage

ACE is available as a command-line interface after you install the python package:

```
usage: ace [-h] [--version] {generate,deconvolve,verify} ...

ACE Configurator for ELISpot.

positional arguments:
  {generate,deconvolve,verify}
                        ACE sub-commands.
    generate            Generates an ELISpot experiment configuration.
    deconvolve          Deconvolve hit peptide IDs given read-outs from an ELISpot experiment.
    verify              Verifies whether an ELISpot assignment satisfies all ACE constraints.

options:
  -h, --help            show this help message and exit
  -v, --version         show program version number and exit
```

Read the full documentation on the python package at https://pirl-unc.github.io/ace/

## 03. Citation

If you use ACE in a publication, please cite our [publication](https://doi.org/10.1093/bib/bbad495) describing ACE.
