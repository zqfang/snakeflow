def gsea_enrichr(diff, log2fc, padj, go, threads):
    # python code
    from pandas import read_table, read_excel, concat, ExcelWriter
    import gseapy as gp

    writer = diff
    sig_deg = read_excel(writer,sheet_name="sig-fc%s-padj%s"%(log2fc, padj))
    sig_deg_up= read_excel(writer,sheet_name="sig-up",)
    sig_deg_dw= read_excel(writer,sheet_name="sig-down",)
    
    degs_sig = [deg.gene_name.squeeze().tolist() for deg in[sig_deg, sig_deg_up,sig_deg_dw]]

    sig_deg_gsea = sig_deg[['gene_name','log2FoldChange']]
    sig_deg_gsea_sort = sig_deg_gsea.sort_values('log2FoldChange',ascending=False)
    sig_deg_gsea_sort = sig_deg_gsea_sort.reset_index(drop=True)


   
    for domain in go:
        prerank = gp.prerank(rnk=sig_deg_gsea_sort, gene_sets=domain, pheno_pos='', pheno_neg='', min_size=15, max_size=500, 
                             outdir='differential_expression/GSEA_prerank_'+domain)
    for domain in go:
        for glist, gl_type in zip(degs_sig, ['all','up','down']):
            enrichr = gp.enrichr(gene_list=glist, gene_sets=domain, 
                                 no_plot=True,  description=gl_type,
                                 outdir='differential_expression/Enrichr_%s_%s'%(domain, gl_type))

gsea_enrichr(snakemake.input[0], snakemake.params['log2fc'], snakemake.params['padj'],
            snakemake.params['go'], snakemake.threads)
