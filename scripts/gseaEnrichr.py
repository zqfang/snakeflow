def gsea_enrichr(diff, log2fc, padj, go):
    # python code
    import errno
    from pandas import read_excel
    import gseapy as gp
    al_res = read_excel(diff, sheetname=None)
    sig_deg = al_res["sig-all.log2fc%s-padj%s"%(log2fc, padj)]
    sig_deg_up = al_res['sig-up']
    sig_deg_dw = al_res['sig-down']
    
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
                           cutoff=0.2, outdir=outdir)
            except Exception:
                print("Server No response. Try again later!")
            #print("connetion to the Enrichr Server is interupted by the host, retry again.")
    #run prerank
    """
    for domain in go:
        try:
            outdir="GO/GSEA_prerank_%s/%s"%(outGSEAname, domain)
            gp.prerank(rnk=sig_deg_gsea_sort, gene_sets=domain,
                        pheno_pos=treat, pheno_neg=ctrl, min_size=15, max_size=500, 
                        outdir=outdir)
        except:
            print("Oops...%s_vs_%s: skip GSEA plotting for %s, please adjust paramters for GSEA input."%(treat, ctrl, domain))
    
    """
    #select columns for gsea
    cols_ = [col for col in sig_deg.columns if col.startswith("TPM")]
    
    cols_group = [col.split(".")[1] for col in cols_ ]
    cols  = [col for col, group in zip(cols_, cols_group) if treat == group] +\
            [col for col, group in zip(cols_, cols_group) if ctrl == group]

    col2 = ['gene_name']+ cols
    cls_vec = [g.split(".")[1] for g in cols]
    
    # run gsea
    for domain in go:
        try:
            outdir="GO/GSEA_%s/%s"%(outGSEAname, domain)
            gs = gp.gsea(data=sig_deg[col2], gene_sets=domain, cls=cls_vec,
                         min_size=10, max_size=500, outdir=outdir)
        except:
            print("Oops...%s_vs_%s: skip GSEA plotting for %s, please adjust paramters for GSEA input."%(treat, ctrl, domain))
        """
        # delete empty dirs
        try:
            os.rmdir(outdir)
        except OSError as ex:
            if ex.errno == errno.ENOTEMPTY:
                print("deleting empty directory."
        """

gsea_enrichr(snakemake.input[0], snakemake.params['log2fc'], snakemake.params['padj'],
            snakemake.params['go'])
