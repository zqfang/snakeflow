def rmats_anno(indir, outdir, rbps, diff_exp, go):

    import glob, os
    import matplotlib 
    matplotlib.use('agg')

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import gseapy as gp

    #gene_expression_table
    gene_exp=pd.read_excel(diff_exp, index_col='gene_id')
    #remove .versions of each id
    gene_exp.index = gene_exp.index.str.split(".").str[0]


    # Significant events are based on FDR < 5% and | deltaPSI | > 10%
    as_rmats = glob.glob(os.path.join(indir, "*.MATS.JCEC.txt"))

    treat, ctrl = indir.split("/")[-1].lstrip("rMATS.").split("_vs_")

    as_type =[]
    as_total = []
    as_sig = []

    for f in as_rmats:
        temp = f.split("/")
        ast =temp[-1].split(".")[0]
        outname= os.path.join(outdir, "%s.MATS.JCEC.sig.csv"%(ast))
        s = pd.read_table(f)
        ss =  s[(s['IncLevelDifference'].abs() > 0.1) & (s['FDR'] < 0.05) ]
        ss.to_csv(outname,index=False)
        
        as_type.append(ast)
        as_total.append(len(s))
        as_sig.append(len(ss))

    #pie chart
    explode = (0, 0, 0, 0, 0.1)  # only "explode" the 2nd slice (i.e. SE)
    fig, ax=plt.subplots(ncols=2,nrows=1, figsize=(6,3))
    ax[0].pie(as_sig, explode=explode, labels=as_type, autopct='%1.1f%%',
              shadow=True, startangle=90)
    ax[0].axis('equal') # Equal aspect ratio ensures that pie is drawn as a circle.
    ax[0].set_title("Significant AS Events: %s vs %s"%(treat, ctrl))
    ax[1].pie(as_total, explode=explode, labels=as_type, autopct='%1.1f%%',
              shadow=True, startangle=90)
    ax[1].axis('equal')
    ax[1].set_title("Total AS Events: %s vs %s"%(treat, ctrl))
    fig.savefig(outdir+"/differential_AS_events.piechart.pdf",bbox_inches='tight')
    fig.savefig(outdir+"/differential_AS_events.piechart.png",bbox_inches='tight', dpi=300)


    SE_sig = pd.read_csv(os.path.join(outdir, "SE.MATS.JCEC.sig.csv"), index_col='ID')

    cols_ = [col for col in gene_exp.columns if col.startswith("TPM")]   
    cols_group = [col.split(".")[1] for col in cols_ ]

    group_b1  = [col for col, group in zip(cols_, cols_group) if treat == group] 
    group_b2  = [col for col, group in zip(cols_, cols_group) if ctrl == group]

    #split psi values for each sample
    data = []
    _b1 = [g.strip("TPM.") for g in group_b1] 
    _b2 = [g.strip("TPM.") for g in group_b2] 
    for i, row in enumerate(_b1):
        sample1 = SE_sig['IncLevel1'].str.split(",").str[i].astype('float')
        sample1.name="PSI."+ row
        data.append(sample1)
    for i, row in enumerate(_b2):
        sample2 = SE_sig['IncLevel2'].str.split(",").str[i].astype('float')
        sample2.name="PSI."+ row
        data.append(sample2)
    dat = pd.concat(data, axis=1,)
    dat = dat.dropna()
    
    
    outdir = outdir+"/Skip_Exons"
    #gsea data
    data_ann = pd.concat([SE_sig[['GeneID','geneSymbol']],dat],axis=1)
    data_ann.to_csv(outdir+"/Diff_skip_exons_table_for_gsea.txt",sep="\t")
    # save psi to csv
    data_ann2 = pd.concat([SE_sig, dat],axis=1)
    data_ann2.to_csv(outdir+"/SE.MATS.JCEC.sig.annotated.csv")

    #plotting
    sns.set(font_scale=1.5, context='talk')
    dat = dat[dat.std(axis=1) != 0 ]
    sg = sns.clustermap(dat,yticklabels=False,figsize=(6,6), z_score=0)

    sg.fig.suptitle("differentially_skipped_exons")
    sg.savefig(outdir+"/differentially_skipped_exons.pdf",bbox_inches='tight')
    sg.savefig(outdir+"/differentially_skipped_exons.png",bbox_inches='tight', dpi=300)



    #load RNA binding protein list
    rbp = pd.read_csv(rbps)
    rbp = rbp.dropna(axis=1)

    #save rbp expression profile
    rbp_exp = gene_exp.loc[rbp.EnsemblGeneID]
    rbp_exp.dropna(inplace=True)
    rbp_exp.to_csv(outdir+"/RNA_Binding_Protein_gene_exp_table.csv")
    rbp_exp.head()

    #save significant changed RBPs
    rbp_sig = rbp_exp[(rbp_exp.log2FoldChange.abs() > 1 ) & (rbp_exp['padj'] <=0.05)]
    rbp_sig.to_csv(outdir+"/RNA_Binding_Protein_gene_exp_table.sig.fc2.csv")


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
    ax.set_xlabel("log$_2$FoldChange(%s/%s)"%(treat, ctrl))
    ax.set_ylabel(" - log$_{10}$ padj")
    #colorbar
    cax=fig.add_axes([1.02,0.25,0.03,0.25])
    cbar = fig.colorbar(sc, cax=cax,)
    cbar.ax.tick_params(right='off',left='off',labelsize=14)
    cbar.ax.set_title('log$_{10}$ baseMean',loc='left',fontsize=14)
    #sns.despine()
    fig.savefig(outdir+"/RBP_vacano.png",bbox_inches='tight')
    fig.savefig(outdir+"/RBP_vacano.pdf",bbox_inches='tight')



    #select columns for gsea
    #b1_treat  = [col for col, group in zip(cols_, cols_group) if treat == group] 
    #b2_treat  = [col for col, group in zip(cols_, cols_group) if ctrl == group]

    #extract expression
    g_b1_meanTPM = rbp_exp[group_b1].mean(axis=1)
    g_b2_meanTPM = rbp_exp[group_b2].mean(axis=1)

    #scatter plot
    fig, ax = plt.subplots(figsize=(6,6))
    sc = ax.scatter(x = np.log2(g_b1_meanTPM),
                    y= np.log2(g_b2_meanTPM),
                    c= rbp_exp.log2FoldChange,
                    cmap=plt.cm.viridis_r, edgecolors='face', s=90)
    #ax.plot(x=[-10,15],y=[-10,15],)
    #colorbar
    cax=fig.add_axes([1.02,0.25,0.03,0.25])
    cbar = fig.colorbar(sc, cax=cax,)
    cbar.ax.tick_params(right='off',left='off',labelsize=14)
    cbar.ax.set_title('log$_2$FoldChange',loc='left',fontsize=14)

    #for idx in rbp_top18.index:
    #    ax.annotate(s=rbp_top18.geneSymbol.loc[idx],xy=(np.log2(rbp_exp.avgFPKM_DMSO.loc[idx]), np.log2(rbp_exp.avgFPKM_YPA.loc[idx])))
    ax.set_xlabel("log$_2$(avgTPM %s)"%treat)
    ax.set_ylabel("log$_2$(avgTPM %s)"%ctrl)
    #sns.despine()
    fig.savefig(os.path.join(outdir,"RBP_scatter.png"),bbox_inches='tight')
    fig.savefig(os.path.join(outdir,"RBP_scatter.pdf"),bbox_inches='tight')





    #bar plot of sig RBPs
    if len(rbp_sig) >= 1:
        rbp_sig = rbp_sig.sort_values('log2FoldChange',)

        fig, ax = plt.subplots(figsize=(6,4))
        # the y coords of this transformation are data, and the
        # x coord are axes
        t = trans.blended_transform_factory(ax.transAxes, ax.transData)

        bar = sns.barplot(x='gene_name',y='log2FoldChange',data=rbp_sig, ax=ax, palette="BuGn_d",edgecolor='none')
        ax.set_xticklabels(rbp_sig.gene_name,rotation=90)
        ax.set_ylabel("log$_2$(%s/%s)"%(treat, ctrl))
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


    
    # go domain
    GO_DOMAIN = go

    plt.style.use('classic')

    for domain in GO_DOMAIN:
        outname = os.path.join(outdir, "GSEA_AS_%s_vs_%s"%(treat, ctrl), domain)
        try:
            prerank = gp.prerank(rnk=rank_list, gene_sets=domain, min_size=15, max_size=500,
                             pheno_pos=treat,pheno_neg=ctrl, outdir=outname)
        except:
            print("something wrong, %s"%domain)

    for domain in GO_DOMAIN:

        for glist, gl_type in zip([rank_list.geneSymbol.squeeze().tolist(),
                                   rank_list_up.geneSymbol.squeeze().tolist(), 
                                   rank_list_down.geneSymbol.squeeze().tolist()],['all','up','down']):
            
            outname = os.path.join(outdir, "Enrichr_SkipExons_%s_vs_%s"%(treat, ctrl))  
            try:     
                enrichr = gp.enrichr(gene_list=glist, gene_sets=domain, description=domain, cutoff=0.2,
                                     outdir=outname+'/%s_%s'%(domain, gl_type))
            except Exception:
                print("Enrichr Server No Response. Try again later!!!")


rmats_anno(snakemake.params['indir'], snakemake.params['outdir'],
            snakemake.params['rbps'], snakemake.input[0], snakemake.params['go'])

