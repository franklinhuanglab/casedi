import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def ase_compare_groups(ase_table, group1_ids, group2_ids,testable_fraction=0.10):
    """Function to test for enrichment of ASE for all genes
    with alellic data in two different sample groups. Example: 
    tumor subtype A vs tumor subtype B, normal samples vs 
    tumor samples, etc
    
    Parameters
    ----------
    ase_table: DataFrame
        Dataframe where each row represents one gene in one sample, with
        a column 'ase' indicating that the gene shows ASE (1) or BAE (0). 
        Required columns: 'gene','ase','sampleID'
    
    testable_fraction: float, default: 0.10
        Minimum fraction of samples in each group that must be testable
        for each gene
        
    group1_ids: array-like
        SampleIDs in comparison group 1
        
    group2_ids: array-like
        SampleIDs in comparison group 2
        
    Returns
    -------
    fisher_table: DataFrame
        Table of counts and results of fisher's exact test for each gene tested
    """
    if (len(group1_ids) == 0) | (len(group2_ids) == 0):
        raise Exception("One or more groups is empty")
    
    mat1 = ase_table[ase_table['sampleID'].isin(group1_ids)]
    mat2 = ase_table[ase_table['sampleID'].isin(group2_ids)]
    
    if (len(set(mat1['sampleID']))!=len(group1_ids)) | (len(set(mat2['sampleID']))!=len(group2_ids)):
        raise Exception("Some group1 or group2 ids were not found in ASE table")

        
    # Get minimum # of samples testable in each group
    mat1_min = len(set(mat1['sampleID']))*testable_fraction
    mat2_min = len(set(mat2['sampleID']))*testable_fraction
    
    mat1_genes = mat1['gene'].value_counts()
    mat1_genes =mat1_genes[mat1_genes >= mat1_min].index
    mat1 = mat1[mat1['gene'].isin(mat1_genes)]
    
    mat2_genes = mat2['gene'].value_counts()
    mat2_genes =mat2_genes[mat2_genes >= mat2_min].index
    mat2 = mat2[mat2['gene'].isin(mat2_genes)]
        
    mat1_sig = mat1[(mat1['ase'] ==1)]
    mat2_sig = mat2[(mat2['ase'] ==1)]

    print("Group1 genes showing ASE: ",len(set(mat1_sig['gene'])))
    print("Group1 total genes: ",len(set(mat1['gene'])))
    print("Group2 genes showing ASE: ",len(set(mat2_sig['gene'])))
    print("Group2 total genes: ",len(set(mat2['gene'])))
    
    # common genes between group 1 and 2 that are testable in  
    # testable_fraction for each group. 
    common_ASE = list(set(mat1_genes).intersection(set(mat2_genes)))
    
    # Genes with at least 1 ASE in either group that are testable in both.
    # This will be the search space for association testing
    common_ASE_sig_mat1 = [x for x in set(mat1_sig['gene']) if x in common_ASE]
    common_ASE_sig_mat2 = [x for x in set(mat2_sig['gene']) if x in common_ASE]
    genes_to_test = set(common_ASE_sig_mat1 + common_ASE_sig_mat2)
    print("testing " + str(len(genes_to_test)) + " genes")
    
    
    # Contingency table will look like:
    #         ASE       BAE
    # group1     8        2
    # group2     1        5
    
    total_mat1 = mat1['gene'].value_counts().loc[genes_to_test]
    total_mat2 = mat2['gene'].value_counts().loc[genes_to_test]
    sig_mat1 = mat1_sig['gene'].value_counts().loc[common_ASE_sig_mat1]
    sig_mat2 = mat2_sig['gene'].value_counts().loc[common_ASE_sig_mat2]
    fisher_table = pd.concat([sig_mat1,total_mat1,sig_mat2,total_mat2],join= 'outer',axis=1)
    fisher_table.columns = ['group1_ASE','group1_total','group2_ASE', 'group2_total']
    fisher_table.fillna(0, inplace=True)
    fisher_table['group1_BAE'] = fisher_table['group1_total'] - fisher_table['group1_ASE']
    fisher_table['group2_BAE'] = fisher_table['group2_total'] - fisher_table['group2_ASE']

    # Run fisher exact test on each gene, getting p-value and odds ratio
    fisher_table['odds_ratio'] = fisher_table.apply(
        lambda row: fisher_exact([[row['group1_ASE'], row['group1_BAE']], [row['group2_ASE'], row['group2_BAE']]])[0],
        axis = 1)
    fisher_table['p_val'] = fisher_table.apply(
        lambda row: fisher_exact([[row['group1_ASE'], row['group1_BAE']], [row['group2_ASE'], row['group2_BAE']]])[1],
        axis = 1)

    # Correct p-values using benjamini-hochberg
    fisher_table['FDR'] = multipletests(fisher_table['p_val'], method='fdr_bh')[1]  # the 2nd element is the actual padjs
    
    return fisher_table
