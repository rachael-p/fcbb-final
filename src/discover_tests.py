# create new conda env according to github and switch py interpreter in vscode


import discover
from discover import pairwise_discover_test
import os
import pandas as pd

directory = "data/mutations"
matrix_list = []
sample_labels = {}

# load file data in to create unified df with all cohorts
for file in os.listdir(directory):
    if file.endswith(".txt"):
        pathname = os.path.join(directory, file)
        cancer_type = file.split("_")[0]
        data = pd.read_csv(pathname, sep="\t", index_col=0)
        data.columns = [f"{cancer_type}_{sample}" for sample in data.columns]
        data.index = data.index.str.strip().str.upper()  # standardize index names
        data = data[~data.index.duplicated(keep='first')]   # drop duplicate genes in that file
        for col in data.columns:
            sample_labels[col] = cancer_type
        matrix_list.append(data)

# data cleaning
matrix = pd.concat(matrix_list, axis=1)
matrix = matrix.groupby(matrix.index).max()  # removes duplicate genes across cohorts by taking max mutation call
sample_labels = pd.Series(sample_labels, name="Cancer Type")
matrix = matrix[(matrix.sum(axis=1) > 0) & (matrix.sum(axis=1) < matrix.shape[1])]  # drops genes that are never mutated or always mutated since not useful for co-occurrence/mutual exclusivity
matrix = matrix.dropna()  # remove NA vals
matrix = matrix.clip(0, 1).astype(int)  # cast non-binary values to 0 or 1

disc_matrix = discover.DiscoverMatrix(matrix)
result = pairwise_discover_test(disc_matrix)
print(result)


