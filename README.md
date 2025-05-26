# MotifLLM

To run the project, please follow the steps below.

## Setting up the Environment

To make your environment ready to run the project, we need to do two things:

1. Install python3.
2. Install the dependencies.

To install the dependencies required for this project, run the following command in your terminal or command prompt:

```bash
pip install -r requirements.txt

```

## Generate JSON tree from molecules (smiles)

To generate JSON tree from molecules (smiles), run the following command (You may replace the dataset with your own dataset):

```bash
python mol2json.py --data=qm9.csv

```
