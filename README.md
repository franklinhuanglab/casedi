# Cancer ASE driver identification (CASEDI)

<img src="casedi.png" width=70% height=70%>

CASEDI is an analysis framework that uses tumor allele-specific expression (ASE) data to prioritize cancer driver genes. CASEDI includes functions to identify genes enriched for ASE between datasets, conduct ASE-stratified survival analysis, and identify genes showing ASE outlier expression.

## Installation

The latest stable release can be installed from PyPI:

```python
pip install casedi
```
You may instead want to use the development version from Github:

```python
pip install git+https://https://github.com/maggietsui/casedi.git
```

## Usage

Examples are available in the jupyter notebook `examples/examples.ipynb`.

## Requirements

+ pandas >= 1.3.1
+ numpy >= 1.23.3
+ scikit-survival
+ scipy >= 1.9.1
+ statsmodels >= 0.13.2
+ matplotlib >= 3.4.2
