
# coding: utf-8
#!usr/bin/env python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


import glob

file = glob.glob("*.diff")




name = [f.rstrip(".diff") for f in file]
nam = [(a.split("_vs_")[0],a.split("_vs_")[1]) for a in ]

# In[26]:

for i,f in enumerate(file):
    dat = pd.read_table(f,skiprows=9)
    dat_sig = dat[(dat['GFOLD(0.01)'].abs() > 0.01) & (dat['log2fdc'].abs() >1)]
    dat_sig.columns= ['GeneSymbol', 'GeneName', 'GFOLD(0.01)', 'E-FDR', 'log2fdc', 
                      '1stRPKM.'+nam[i][0],'2ndRPKM.'+nam[i][1]]
    dat_sig['up_down'] = dat_sig['log2fdc'].apply(lambda x: 'up' if x > 0 else 'down')
    dat_sig.to_csv(f+".fc2.sig.csv", index=False)


diff_file = glob.glob("*.fc2.sig.csv")


#name = ['KO1_4d_vs_WT1_4', 'KO2_4d_vs_WT2_4', 'KO2_8d_vs_WT2_8', 'KO8_4d_vs_WT1_8']



file_dict = {}
for i,f in enumerate(diff_file):
    diff_merge = pd.read_csv(f)
    diff_list = diff_merge.GeneName.str.upper().squeeze().tolist()
    diff_list_up =  diff_merge[diff_merge['log2fdc'] > 0].GeneName.str.upper().squeeze().tolist()
    diff_list_dw =  diff_merge[diff_merge['log2fdc'] < 0].GeneName.str.upper().squeeze().tolist()
    file_dict.update({name[i] : [diff_list, diff_list_up, diff_list_dw]})



# In[34]:

import gseapy as gp


# In[35]:

genesets = ['GO_Biological_Process_2015','GO_Cellular_Component_2015','KEGG_2016','MGI_Mammalian_Phenotype_2013']


# In[36]:

#all_list = [diff_list, diff_list_dw, diff_list_up]
desc = ['all','up_regulate','down_regulate']


# In[40]:

plt.style.use('classic')

for gs in genesets:
    for k in file_dict.keys():
        alst = file_dict.get(k)
        for i in range(len(alst)):
            r = gp.enrichr(gene_list=alst[i], gene_sets=gs, 
                      description= desc[i], outdir='Enrichr_'+k+'_'+ desc[i]+'_'+ gs)




