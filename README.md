[![Check code style](https://github.com/JetBrains-Research/Genegram/actions/workflows/check_code_style.yml/badge.svg)](https://github.com/JetBrains-Research/Genegram/actions/workflows/check_code_style.yml)
---
# Genegram

## Description

[comment]: <> (TODO)

## Installation

Simply clone the repository and run the following commands:

```bash
pip install -r requirements.txt
```

## Usage

Run the following command with arguments.

```bash
python -m genegram
```

### **Arguments**

Argument | Required | Description
:--- | :---: | :---
--seq_fasta | True | Path to the `seq.fasta` file
--out | True | Path to the folder where the predictions will be saved
--model | False | Type of the model to be used: </br> `main` -- The default model, the best on average </br> `mps` -- Multiplet prediction model </br> `pks` -- Pseudoknots prediction model

## Code style

We recommend you use a [pre-commit](https://pre-commit.com/#install) hook, which runs [black](https://github.com/psf/black) when you type git commit.

Run the following script at the root of the repository:

```bash
pip install pre-commit
pre-commit install
```
