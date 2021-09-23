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

### **Required arguments**

Argument | Description
:--- | :---
--seq_fasta | Path to the `seq.fasta` file
--out | Path to the folder where the predictions will be saved

### **NOT required arguments**

Argument | Description
:--- | :---
--model | Type of the model to be used

## Code style

We recommend you use a [pre-commit](https://pre-commit.com/#install) hook, which runs [black](https://github.com/psf/black) when you type git commit.

Run the following script at the root of the repository:

```bash
pip install pre-commit
pre-commit install
```
