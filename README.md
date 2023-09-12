# ACE Configurator for ELISpot

ACE facilitates (1) generation of ELISpot configurations (peptide-pool assignments) 
using a deep learning approach to cluster similar peptides and (2) deconvolution 
of pool spot counts for identification of immunogenic peptides.

[![Build Status](https://app.travis-ci.com/pirl-unc/ace.svg?branch=main)](https://app.travis-ci.com/pirl-unc/ace)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

## 01. Installation

### 01-1. Standalone Graphical User Interface (Recommended)

| Operating System | Link | Version           |
|------------------| ---- |-------------------|
| Mac              | [Download](https://github.com/pirl-unc/ace/releases/download/v0.1.0.4/ace-0.1.0.4.tar.gz) | v0.1.0.5 (latest) | 
| Windows          | [Download](https://github.com/pirl-unc/ace/releases/download/v0.1.0.4/ace-0.1.0.4.tar.gz) | v0.1.0.5 (latest) |

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

* python3 (>= 3.10)
* pandas (>=1.5.2)
* numpy (>=1.23.1)
* ortools (9.3.10497)
* torch
* transformers (==4.30.2)
* scikit-learn
* openpyxl
* golfy (>=2.5.0) 
* levenshtein

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
  --version, -v         show program's version number and exit
```

Read the full documentation on the python package at https://pirl-unc.github.io/ace/

## 03. Citation

If you use ACE in a publication, please cite our 
[preprint](https://www.biorxiv.org/content/10.1101/2023.09.02.554864v1) describing ACE.
