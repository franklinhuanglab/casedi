import pandas as pd

# make matrix for different column values in main table
def make_matrix(ase_table, value_column):
    mat = (
        ase_table[["sampleID", "gene", value_column]]
        .set_index(["sampleID", "gene"])
        .unstack()
    )
    mat.columns = mat.columns.get_level_values(1)
    return mat