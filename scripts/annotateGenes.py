def anno_genes(annotation, tpms, deseq, sample_anno, diff_anno, log2fc, padj, threads):
    # python code
    from pandas import read_table, read_excel, concat, ExcelWriter

    anno = read_table(annotation, index_col='gene_id')

    tpm = read_table(tpms, index_col=0)
    tpm.columns = ['TPM.'+col for col in tpm.columns]
    deseq = read_table(deseq, index_col=0)
    #merge results
    merge = concat([anno, deseq, tpm], axis=1, join='inner')
    merge.index.name ='gene_id'
    merge.to_csv(sample_anno)

    sig_deg = merge[(merge['log2FoldChange'].abs()> params.log2fc) & (merge['padj'] < params.padj) ]
    sig_deg = sig_deg.sort_values('padj',axis=0)
    sig_deg.loc[:,'up_down'] = sig_deg.log2FoldChange.apply(lambda x : 'up' if x > 0 else 'down' )
     

    sig_deg_up = sig_deg[sig_deg['up_down'] == 'up']
    sig_deg_dw = sig_deg[sig_deg['up_down'] == 'down'] 

    writer = ExcelWriter(diff_anno)
    sig_deg.to_excel(writer,sheet_name="sig-log2fc%s-padj%s"%(log2fc, padj))
    sig_deg_up.to_excel(writer,sheet_name="sig-up",)
    sig_deg_dw.to_excel(writer,sheet_name="sig-down",)
    writer.save()

anno_genes(snakemake.input[0], snakemake.input[1], snakemake.input[2], 
             snakemake.output['sample_anno'],snakemake.output['diff_anno'],  
             snakemake.params['log2fc'],snakemake.params['padj'],
             snakemake.threads)
