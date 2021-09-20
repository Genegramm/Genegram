[![Check code style](https://github.com/JetBrains-Research/Genegram/actions/workflows/check_code_style.yml/badge.svg)](https://github.com/JetBrains-Research/Genegram/actions/workflows/check_code_style.yml)
---
# Genegram

## Description

[comment]: <> (TODO)

## Installation

Simply clone the repository and run the following commands:

```bash
pip install -r requirements.txt
pip install .
```

## Usage

Run the `python -m generam` with the arguments.

### **Required arguments**

Argument | Description
:--- | :---
--seq_fasta | Path to the `seq.fasta` file
--out | Path to the folder where the predictions will be saved

## Code style

We recommend you use a [pre-commit](https://github.com/pre-commit/pre-commit) hook, which runs [black](https://github.com/psf/black) when you type git commit:

```bash
pre-commit install
```
