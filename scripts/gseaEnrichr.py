def gsea_enrichr(diff, log2fc, padj, go):
    # python code
    from pandas import read_excel
    import gseapy as gp

    sig_deg = read_excel(diff, sheet_name="sig-all.log2fc%s-padj%s"%(log2fc, padj))
    sig_deg_up= read_excel(diff, sheet_name="sig-up",)
    sig_deg_dw= read_excel(diff, sheet_name="sig-down",)
    
    degs_sig = [deg.gene_name.squeeze().tolist() for deg in[sig_deg, sig_deg_up,sig_deg_dw]]

    sig_deg_gsea = sig_deg[['gene_name','log2FoldChange']]
    sig_deg_gsea_sort = sig_deg_gsea.sort_values('log2FoldChange',ascending=False)
    sig_deg_gsea_sort = sig_deg_gsea_sort.reset_index(drop=True)


    outGSEAname = diff.split("/")[-1].lstrip("diff_").rpartition("_")[0]
    treat, ctrl =outGSEAname.split("_vs_")
    
    for domain in go:
        for glist, gl_type in zip(degs_sig, ['all','up','down']):
            try:
                outdir='GO/Enrichr_%s/%s_%s'%(outGSEAname, domain, gl_type)
                gp.enrichr(gene_list=glist, gene_sets=domain, description=gl_type,
                            outdir=outdir)
            except:
                print("connetion to the Enrichr Server is interupted by the host, retry again.")
    #run gseapy
    for domain in go:
        try:
            outdir="GO/GSEA_%s/%s"%(outGSEAname, domain)
            gp.prerank(rnk=sig_deg_gsea_sort, gene_sets=domain,
                        pheno_pos=treat, pheno_neg=ctrl, min_size=15, max_size=500, 
                        outdir=outdir)
        except:
            print("Oops...skip GSEA plotting for %s, please adjust paramters for GSEA input."%domain)


gsea_enrichr(snakemake.input[0], snakemake.params['log2fc'], snakemake.params['padj'],
            snakemake.params['go'])
