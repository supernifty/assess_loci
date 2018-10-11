# find significantly mutated loci when comparing groups

## Installation
This software requires Python3.

```
python -m venv assess-env
source ./assess-env/bin/activate
pip install -r requirements.txt
```

## Usage

```
python assess.py --panels exons.bed --groups 0 0 0 0 1 1 1 --names TUMOR --vcfs vcfs... vcfs... --verbose > assess.mlh1msh2.tsv
```
