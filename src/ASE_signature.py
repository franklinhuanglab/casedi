import pandas as pd
import numpy as np

# Load imprinted genes
imp = pd.read_csv("/wynton/home/rotation/mtsui/ASE_project/GTEx/geneimprint_imprinted_genes.csv")
imp = imp[imp['Status'] == 'Imprinted']
imp.set_index("Gene",inplace=True)

loc = pd.read_csv("/wynton/home/rotation/mtsui/ASE_project/candidate_genes_tables/locvnorm_ase_genes.csv",index_col=0)
loc['n_ase_with_strongcnv'].fillna(0,inplace=True)

recurr_perc = 0.05
loc= loc[~(loc.index.isin(imp.index)) & # not imprinted
         (loc['ase_t'] >= 495*recurr_perc) & # Recurrent > 5%
         (loc['n_ase_with_strongcnv']/loc['ase_t'] < 0.2)]# CNV in <20%
