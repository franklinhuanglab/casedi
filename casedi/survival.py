from sksurv.nonparametric import kaplan_meier_estimator
from sksurv.compare import compare_survival
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from casedi.util import make_matrix

def ase_compare_survival(gene, ase_table, survival_table,  return_stats=False,expression_table=None,plot=False,stratify_expression=None):
    """Wrapper function for sksurv kaplan_meier_estimator and 
    compare_survival that uses log-rank test to compares survival
    between samples showing ASE and BAE. Optionally, compare 
    survival between ASE + high expression and BAE + low expression
    or vice versa.
    
    Parameters
    ----------
    gene: gene to compute survival difference for
    
    ase_table: DataFrame
        Dataframe where each row represents one gene in one sample, with
        a column 'ase' indicating that the gene shows ASE (1) or BAE (0). 
        Required columns: 'gene','ase','sampleID'
    
    survival_table: DataFrame
        Dataframe of survival time and survival status for each sample.
        Index must be the sampleID. Survival status must be either 
        0 for no event or 1 for event. Required columns: 'survival_time',
        'survival_status'
        
    return_stats: bool, default: False
        Return summary from sksurv compare_survival function.
        
    stratify_expression: {'ase high', 'ase low', 'None'}, default: None
        Optional, set to further stratify samples into ASE + high expression
        (greater than median) and BAE + low expression (less than median)
        or vice versa.
    
    expression_table: DataFrame
        Required if stratify_expression is set. A dataframe where rows are
        samples and columns are genes, with expression values for each gene.
        Index must be the sampleID.
    
    plot: bool, default: False
        Set to plot Kaplan-Meier curve with p-value
    
    Returns
    -------
    chisq: float
        Test statistic.
    pvalue: float
        Two-sided p-value with respect to the null hypothesis
        that the hazard rates across all groups are equal.
    stats: pandas.DataFrame
        Summary statistics for each group:  number of samples,
        observed number of events, expected number of events,
        and test statistic.
        Only provided if `return_stats` is True.
    covariance: array, shape=(n_groups, n_groups)
        Covariance matrix of the test statistic.
        Only provided if `return_stats` is True.

    list of Line2D:
        Only if plot=True. A list of lines representing the plotted data
        from matplotlib
    """
    
    ase_mat = make_matrix(ase_table,'ase')
    if (ase_mat[gene] == 1).any() == False:
        raise Exception("No samples show ASE for input gene")
    if stratify_expression==None:
        if len(ase_mat.index) != len(survival_table.index):
            raise Exception("Unequal number of Sample IDs in tables")
        elif set(ase_mat.index) != set(survival_table.index):
            raise Exception("Sample IDs of tables do not match")
        elif gene not in set(ase_mat.columns):
            raise Exception("Input gene not found")
    else: # check if expression table meet criteria 
        if expression_table is None:
            raise ValueError("Expression table must be passed if stratify_expression is set")
        if len(ase_mat.index) != len(expression_table.index):
            raise Exception("Unequal number of Sample IDs in tables")
        elif set(ase_mat.index) != set(expression_table.index):
            raise Exception("Sample IDs of tables do not match")
        
    # Get samples showing ase and bae 
    ase_id = ase_mat[(ase_mat[gene] ==1)].index
    bae_id = ase_mat[ase_mat[gene]==0].index
    
    # Get survival time for ase and bae samples
    ase_vals = survival_table.loc[ase_id, "survival_time"].astype(float).dropna()
    bae_vals = survival_table.loc[bae_id, "survival_time"].astype(float).dropna()
    
    if stratify_expression==None:
        groups = pd.concat([pd.Series(index=bae_vals.index,data="BAE"), pd.Series(index=ase_vals.index,data="ASE")])
        label_ase='ASE'
        label_bae='BAE'
    else:
        # Get expr vals for pats with ase data
        expression_table = expression_table[(expression_table.index.isin(ase_vals.index)) | (expression_table.index.isin(bae_vals.index))]
    
        if stratify_expression=="ase high":
            # Get samples with high expression and ase
            # then samples with low expression and bae
            high = expression_table[(expression_table.index.isin(ase_id)) & (expression_table[gene] > expression_table[gene].median())].index
            low = expression_table[(expression_table.index.isin(bae_id)) & (expression_table[gene] <= expression_table[gene].median())].index

            # subset the survival vals
            ase_vals = ase_vals[ase_vals.index.isin(high)]
            bae_vals = bae_vals[bae_vals.index.isin(low)]

            label_ase ="ASE and high expression"
            label_bae ="BAE and low expression"
        elif stratify_expression=="ase low":
            # Get samples with high expression and bae
            # then samples with low expression and ase
            high = expression_table[(expression_table.index.isin(bae_id)) & (expression_table[gene] > expression_table[gene].median())].index
            low = expression_table[(expression_table.index.isin(ase_id)) & (expression_table[gene] <= expression_table[gene].median())].index

            ase_vals = ase_vals[ase_vals.index.isin(low)]
            bae_vals = bae_vals[bae_vals.index.isin(high)]

            label_ase ="ASE and low expression"
            label_bae ="BAE and high expression"

        high_group = pd.Series(index=high,data="high")
        low_group = pd.Series(index=low,data="low")   
        groups = pd.concat([low_group, high_group])

    time_treatment_ase, survival_prob_treatment_ase = kaplan_meier_estimator(
           survival_table.loc[ase_vals.index, "survival_status"] == 1,
            ase_vals)

    time_treatment_bae, survival_prob_treatment_bae = kaplan_meier_estimator(
            survival_table.loc[bae_vals.index,  "survival_status"] == 1,
            bae_vals)
    
    event_ind = survival_table.loc[list(bae_vals.index) + list(ase_vals.index), "survival_status"] == 1
    time_to = pd.concat([bae_vals,ase_vals])
    
    dt=np.dtype('bool,float')
    strarray = np.array(list(zip(event_ind,time_to)),dtype=dt)
    
    if plot: # Plot kaplan meier plot
        pvalue=compare_survival(y=strarray, group_indicator=groups)[1]
        plt.figure(figsize=(4,3.5),dpi=150)
        plt.step(time_treatment_ase, survival_prob_treatment_ase, where="post",
                     label= label_ase, color = "firebrick", linestyle='dashed')    
        plt.step(time_treatment_bae, survival_prob_treatment_bae, where="post",
                     label=label_bae, color = "slateblue")

        plt.ylabel("Survival")
        plt.xlabel("Time")
        plt.ylim(0,1.1)
        plt.legend(loc='lower left',title=gene,fontsize=8)
        
        # label pvalue
        plt.text(0, 1.12, "p = " + str(round(pvalue,4)))

    else:
        return compare_survival(
                    y=strarray, group_indicator=groups, return_stats=return_stats)
 