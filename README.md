[![Check code style](https://github.com/JetBrains-Research/Genegram/actions/workflows/check_code_style.yml/badge.svg)](https://github.com/JetBrains-Research/Genegram/actions/workflows/check_code_style.yml)
---
# Genegram

## Description

[comment]: <> (TODO)

## Install

We use [`Pipenv`](https://pipenv.pypa.io/en/latest/) to manage dependencies.

### [Install Pipenv](https://pipenv.pypa.io/en/latest/#install-pipenv-today)

```shell
pip install --user pipenv
```

### Install Genegram (from sources)

```shell
git clone https://github.com/JetBrains-Research/Genegram
cd Genegram
pipenv install --ignore-pipfile
```

## Usage

Run the following command with arguments.

```bash
python -m genegram
```

### **Arguments**

Argument | Required | Description
:--- | :---: | :---
-i, --inp | True | Path to the [`FASTA`](http://genetics.bwh.harvard.edu/pph/FASTA.html) file
-o, --out | True | Path to the folder where the [`Connectivity Tables`](http://rna.urmc.rochester.edu/Text/File_Formats.html#CT) will be saved
-m, --model | False | Type of the model to be used: </br> `main` -- The default model, the best on average </br> `mps` -- Multiplet prediction model </br> `pks` -- Pseudoknots prediction model
-l, --log | False | Type of the logging level to be used: </br> `INFO` -- Confirmation that things are working as expected </br> `WARNING` -- An indication that something unexpected happened, the software is still working as expected </br> `ERROR` -- Due to a more serious problem, the software has not been able to perform some function </br> `CRITICAL` -- A serious error, indicating that the program itself may be unable to continue running </br> `DEBUG` -- Detailed information, typically of interest only when diagnosing problems

## Code style

We recommend you use a [pre-commit](https://pre-commit.com/#install) hook, which runs [black](https://github.com/psf/black) when you type git commit.

### Install pre-commit

```shell
pipenv install --dev
pre-commit install
```

### Use pre-commit

```shell
pre-commit run --all-files --color always --verbose
```
