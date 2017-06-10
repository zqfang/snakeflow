

import glob, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp


# Significant events are based on FDR < 5% and | deltaPSI | > 10%


as_rmats = glob.glob("./rMATS.HDE_vs_Ctrl/*.MATS.JCEC.txt")

outdir="./rMATS.HDE_vs_Ctrl_sig"
os.mkdir(outdir)

group_b1 = "HDE"
group_b2 = "Ctrl"



# In[14]:

as_type =[]
as_total = []
as_sig = []

for f in as_rmats:
    temp = f.split("/")
    ast =temp[-1].split(".")[0]
    outname= outdir+temp[-1].replace(".txt", ".sig.csv")
    s = pd.read_table(f)
    ss =  s[(s['IncLevelDifference'].abs() > 0.1) & (s['FDR'] < 0.05) ]
    ss.to_csv(outname,index=False)
    
    as_type.append(ast)
    as_total.append(len(s))
    as_sig.append(len(ss))

SE_sig = pd.read_csv("./rMATS.HDE_vs_Ctrl_sig/SE.MATS.JCEC.sig.csv", index_col='ID')


sampleInfo = pd.read_table("../snakeflow/example_sample_info.txt",header=None, sep=" ")
sampleInfo.columns="sample_name alias group time".split()
group_b1=sampleInfo.loc[sampleInfo.group==group_b1,'alias']
group_b2=sampleInfo.loc[sampleInfo.group==group_b2,'alias']
sampleInfo


#split psi values for each sample
data = []

for i, row in enumerate(group_b1):
    sample1 = SE_sig['IncLevel1'].str.split(",").str[i].astype('float')
    sample1.name="psi."+ row
    data.append(sample1)
for i, row in enumerate(group_b2):
    sample2 = SE_sig['IncLevel2'].str.split(",").str[i].astype('float')
    sample2.name="psi."+ row
    data.append(sample2)
dat = pd.concat(data,axis=1,)


#gsea data
data_ann = pd.concat([SE_sig[['GeneID','geneSymbol']],dat],axis=1)
data_ann.to_csv(outdir+"/Diff_skip_exons_table_for_gsea.txt",sep="\t")

# save psi to csv
data_ann2 = pd.concat([SE_sig, dat],axis=1)
data_ann2.to_csv(outdir+"/SE.MATS.JCEC.sig.annotated.csv")


#plotting
sns.set(font_scale=1.5, context='talk')
sg = sns.clustermap(dat,yticklabels=False,figsize=(6,6),z_score=0)
sg.fig.suptitle("differentially_skipped_exons")
sg.savefig("./differentially_skipped_exons.pdf")
sg.savefig("./differentially_skipped_exons.png",dpi=300)


#gene_expression_table
gene_exp=pd.read_excel("../differential_expression/diff_HDE_vs_Ctrl_results.annotated.xls",index_col='gene_id')

#remove .versions of each id
gene_exp.index = gene_exp.index.str.split(".").str[0]

#load RNA binding protein list
rbp = pd.read_csv("../temp/221RBPs.csv")
rbp = rbp.dropna(axis=1)

#save rbp expression profile
rbp_exp = gene_exp.loc[rbp.EnsemblGeneID]
rbp_exp.dropna(inplace=True)
rbp_exp.to_csv("RNA_Binding_Protein_gene_exp_table.csv")
rbp_exp.head()

#save significant changed RBPs
rbp_sig = rbp_exp[(rbp_exp.log2FoldChange.abs() > 1 ) & (rbp_exp['padj'] <=0.05)]
rbp_sig.to_csv("RNA_Binding_Protein_gene_exp_table.sig.fc2.csv")


#vacano plot
import matplotlib.transforms as trans
sns.set(style='whitegrid',context='talk',font_scale=1.5)


fig, ax = plt.subplots(figsize=(6,6))
sc = ax.scatter(x = rbp_exp.log2FoldChange, y=  - np.log10(rbp_exp.padj), c= np.log10(rbp_exp.baseMean),
           cmap=plt.cm.viridis_r, edgecolor='face')
ax.vlines(x=1,ymin=0, ymax=5,linestyles='dotted',linewidths=2)
ax.vlines(x=-1,ymin=0, ymax=5,linestyles='dotted',linewidths=2)
ax.set_ylim([-0.2,3])
ax.set_xlim([-4,4])
ax.set_xlabel("log$_2$FoldChange(%s/%s)"%(g_b1,g_b2))
ax.set_ylabel(" - log$_{10}$ padj")
#colorbar
cax=fig.add_axes([1.02,0.25,0.03,0.25])
cbar = fig.colorbar(sc, cax=cax,)
cbar.ax.tick_params(right='off',left='off',labelsize=14)
cbar.ax.set_title('log$_{10}$ baseMean',loc='left',fontsize=14)
#sns.despine()
fig.savefig(outdir+"RBP_vacano.png",bbox_inches='tight')
fig.savefig(outdir+"RBP_vacano.pdf",bbox_inches='tight')


#extract expression
b1_treat = [col for col in rbp_exp.columns if col.startswith("TPM."+g_b1)]
b2_treat = [col for col in rbp_exp.columns if col.startswith("TPM."+g_b2)]
g_b1_meanTPM = rbp_exp[b1_treat].mean(axis=1)
g_b2_meanTPM = rbp_exp[b2_treat].mean(axis=1)

#scatter plot
fig, ax = plt.subplots(figsize=(6,6))
sc = ax.scatter(x = np.log2(g_b1_meanTPM),
                y= np.log2(g_b2_meanTPM),
                c= rbp_exp.log2FoldChange,
                cmap=plt.cm.viridis_r, edgecolors='face',s=90)
#ax.plot(x=[-10,15],y=[-10,15],)
#colorbar
cax=fig.add_axes([1.02,0.25,0.03,0.25])
cbar = fig.colorbar(sc, cax=cax,)
cbar.ax.tick_params(right='off',left='off',labelsize=14)
cbar.ax.set_title('log$_2$FoldChange',loc='left',fontsize=14)

#for idx in rbp_top18.index:
#    ax.annotate(s=rbp_top18.geneSymbol.loc[idx],xy=(np.log2(rbp_exp.avgFPKM_DMSO.loc[idx]), np.log2(rbp_exp.avgFPKM_YPA.loc[idx])))
ax.set_xlabel("log$_2$(avgTPM %s)"%g_b1)
ax.set_ylabel("log$_2$(avgTPM %s)"%g_b2)
#sns.despine()
plt.show()
fig.savefig(os.path.join(outdir,"RBP_scatter.png"),bbox_inches='tight')
fig.savefig(os.path.join(outdir,"RBP_scatter.pdf"),bbox_inches='tight')



#bar plot of sig RBPs
rbp_sig = rbp_sig.sort_values('log2FoldChange',)

# the y coords of this transformation are data, and the
# x coord are axes
t = trans.blended_transform_factory(ax.transAxes, ax.transData)

fig, ax = plt.subplots(figsize=(6,4))
#ax.bar(left=np.arange(len(rbp_top18)),height=rbp_top18.log2FoldChange,color='#4C72B0', align='center',
#       tick_label=rbp_top18.geneSymbol)
bar = sns.barplot(x='gene_name',y='log2FoldChange',data=rbp_sig, ax=ax,palette="BuGn_d",edgecolor='none')
ax.set_xticklabels(rbp_sig.gene_name,rotation=90)
ax.set_ylabel("log$_2$(%s/%s)"%(g_b1,g_b2))
ax.set_title("Differential_RNA_Binding_Protein_Expression")
#ax.hlines(xmin=[0,0],xmax=[17,17],y=[1,-1], linestyles='dashed',linewidths=1,colors='k')
ax.set_xlabel("")
sns.despine()
fig.savefig(os.path.join(outdir,"RBP_Expression_sig.png"),bbox_inches='tight')
fig.savefig(os.path.join(outdir,"RBP_Expression_sig.pdf"),bbox_inches='tight')


# ## GO
SE_sig = SE_sig.sort_values(by='IncLevelDifference',ascending=False)
rank_list = SE_sig[['geneSymbol','IncLevelDifference']]
rank_list.head()


# In[145]:

rank_list = rank_list.reset_index()
rank_list = rank_list.drop('ID',axis=1)
rank_list_up = rank_list[rank_list.IncLevelDifference < 0]
rank_list_down = rank_list[rank_list.IncLevelDifference > 0]


outname = "GSEA_AS_HDE_vs_Ctrl/"
# go domain
GO_DOMAIN = ['GO_Cellular_Component_2015','GO_Molecular_Function_2015',
             'GO_Biological_Process_2015','Human_Phenotype_Ontology',
             'MSigDB_Oncogenic_Signatures','WikiPathways_2016',
             'KEGG_2016']

plt.style.use('classic')
for domain in GO_DOMAIN:
    try:
        prerank = gp.prerank(rnk=rank_list, gene_sets=domain, min_size=15, max_size=500,
                         pheno_pos=g_b1,pheno_neg=g_b2, outdir=outname+domain)
    except:
        print("something wrong, %s"%domain)


# In[149]:

GO_DOMAIN = ['MSigDB_Oncogenic_Signatures','WikiPathways_2016','KEGG_2016']
for domain in GO_DOMAIN:
    for glist, gl_type in zip([rank_list.geneSymbol.squeeze().tolist(),
                               rank_list_up.geneSymbol.squeeze().tolist(), 
                               rank_list_down.geneSymbol.squeeze().tolist()],['all','up','down']):
        enrichr = gp.enrichr(gene_list=glist, gene_sets=domain, description=domain, 
                             outdir='Enrichr_Skip_exons/%s_%s'%(domain, gl_type))





