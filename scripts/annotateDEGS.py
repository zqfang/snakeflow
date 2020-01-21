def anno_genes(annotation, tpms, deseqs, diff_anno, samples, alias, group, treat, ctrl, log2fc, padj, threads):
    # python code
    import os
    import pandas as pd 

    anno = pd.read_csv(annotation, index_col='gene_id', sep="\t")
    anno.index = anno.index.str.split(".").str[0]

    tpm = pd.read_csv(tpms, index_col=0, sep="\t")
    tpm = tpm[samples]
    cols_ = ["TPM.%s.%s"%(g, a) for a, g, s in zip(alias, group, samples)]
    tpm.columns = cols_
    #parse treat and ctrl
    #diff_group = deseqs.split("/")[-1].lstrip("diff_").rpartition("_")[0]
    #treat, ctrl = diff_group.split("_vs_")

    #select columns for treat and ctrl group 
    #cols_ = [col for col in tpm.columns if col.startswith("TPM.")]    
    #cols_group = [col.split(".")[1] for col in cols_ ]
    cols  = [col for col, g in zip(cols_, group) if treat == g] +\
            [col for col, g in zip(cols_, group) if ctrl == g]
    tpm_diffs = tpm[cols]
   
    deseq = pd.read_csv(deseqs, index_col=0, sep="\t")
    #merge results
    merge = pd.concat([anno, deseq, tpm_diffs], axis=1, join='inner', sort=True)
    merge.index.name ='gene_id'
    
    #output a blacklist for group comparasion have not sinificant differential genes
    sig_deg = merge[(merge['log2FoldChange'].abs()> log2fc) & (merge['padj'] < padj) ]

    writer = pd.ExcelWriter(diff_anno)
    merge.to_excel(writer, sheet_name="gene_exp")

    if sig_deg.shape[0] > 0:
        sig_deg = sig_deg.sort_values('padj', axis=0)
        sig_deg.loc[:,'up_down'] = sig_deg.log2FoldChange.apply(lambda x : 'up' if x > 0 else 'down' )
        sig_deg_up = sig_deg[sig_deg.log2FoldChange>0]
        sig_deg_dw = sig_deg[sig_deg.log2FoldChange<0] 
        #write data to excel
        sig_deg.to_excel(writer, sheet_name="sig-all.log2fc%s-padj%s"%(log2fc, padj))        
        sig_deg_up.to_excel(writer, sheet_name="sig-up",)
        sig_deg_dw.to_excel(writer, sheet_name="sig-down",)
    else:
        with open("temp/blacklist.txt",'a') as black:
            black.write(diff_anno+'\n')
        os.system("touch %s"%diff_anno)
    writer.save()

anno_genes(snakemake.params['gene_anno'], snakemake.params['tpm'],
           snakemake.input[0], snakemake.output[0], 
             snakemake.params['samples'],snakemake.params['alias'],
             snakemake.params['group'], snakemake.params['treat'],
             snakemake.params['ctrl'], snakemake.params['log2fc'],
             snakemake.params['padj'], snakemake.threads)
