import scipy
import pandas as pd

from casedi_test.util import make_matrix

def ase_compare_expression(gene, ase_table, expression_table, alternative='two-sided'):
    """Wrapper function for scipy.stats.mannwhitneyu to test expression
    difference between samples showing ASE and BAE using the Mann-Whitney
    U test
    
    Parameters
    ----------
    gene: gene to compute expression difference for
    
    ase_table: DataFrame
        Dataframe where each row represents one gene in one sample, with
        a column 'ase' indicating that the gene shows ASE (1) or BAE (0). 
        Required columns: 'gene','ase','sampleID'
    
    expression_table: DataFrame
        A dataframe where rows are samples and columns are genes, with 
        expression values for each gene. Index must be the sampleID.
    
    alternative: {'two-sided', 'less', 'greater'}, Default: two-sided
        'greater' tests whether ASE distribution is greater than BAE
        distribution. 'less' tests whether ASE distribution is less 
        than BAE distribution. 
        
    Returns
    -------
    MannwhitneyuResult, including Mann Whitney U statistic and pvalue
    """
    ase_mat = make_matrix(ase_table,'ase')
    
    # Make sure tables have the same sampleIDs
    if len(ase_mat.index) != len(expression_table.index):
            raise Exception("Unequal number of Sample IDs in tables")
    elif set(ase_mat.index) != set(expression_table.index):
        raise Exception("Sample IDs of tables do not match")
            
    bae = ase_mat.loc[ase_mat[gene] == 0, gene].index
    ase = ase_mat.loc[ase_mat[gene] == 1, gene].index

    if len(bae) == 0:
        raise Exception("No samples showing BAE")
    elif len(ase) == 0:
        raise Exception("No samples showing ASE")
        
    return scipy.stats.mannwhitneyu(x=expression_table.loc[expression_table.index.isin(ase),gene].dropna(),
                             y=expression_table.loc[expression_table.index.isin(bae),gene].dropna(),
                             alternative=alternative)
