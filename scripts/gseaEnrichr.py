def gsea_enrichr(diff, treat, ctrl, log2fc, padj, go):
    # python code
    import os, errno
    from pandas import read_excel
    import gseapy as gp

    #outputfile name
    outGSEAname = diff.split("/")[-1].lstrip("diff_").rpartition("_")[0]
    #treat, ctrl =outGSEAname.split("_vs_")    

    #parse blacklist and skip no significant results
    if os.path.isfile("temp/blacklist.txt"):
        with open("temp/blacklist.txt") as black:
            blacklist = [ bla.strip('\n') for bla in black]
            # handle files with no significant genes
        if diff in blacklist:
            print("No significant degs in this comparasion group.")
            for domain in go:
                outfile1="GO/GSEA_%s/%s/gseapy.gsea.gene_sets.report.csv"%(outGSEAname, domain)
                os.system("touch %s"%outfile1)
                for gl_type in ['all','up','down']:
                    outfile2='GO/Enrichr_{n}/{d}_{t}/{d}.{t}.enrichr.reports.txt'.format(n=outGSEAname, d=domain, t=gl_type)
                    os.system("touch %s"%outfile2)
            return 

    #start to parse significant results
    al_res = read_excel(diff, sheetname=None)
    sig_deg = al_res["sig-all.log2fc%s-padj%s"%(log2fc, padj)]
    sig_deg_up = al_res['sig-up']
    sig_deg_dw = al_res['sig-down']
    
    degs_sig = [deg.gene_name.squeeze().tolist() for deg in[sig_deg, sig_deg_up,sig_deg_dw]]

    sig_deg_gsea = sig_deg[['gene_name','log2FoldChange']]
    sig_deg_gsea_sort = sig_deg_gsea.sort_values('log2FoldChange',ascending=False)
    sig_deg_gsea_sort = sig_deg_gsea_sort.reset_index(drop=True)

    
    for domain in go:
        for glist, gl_type in zip(degs_sig, ['all','up','down']):
            outdir='GO/Enrichr_%s/%s_%s'%(outGSEAname, domain, gl_type)
            outfile = "{o}/{d}.{t}.enrichr.reports.txt".format(o=outdir, d=domain, t=gl_type)
            #skip plotting while file exists
            if os.path.isfile(outfile): continue
            try:                
                gp.enrichr(gene_list=glist, gene_sets=domain, description=gl_type,
                           cutoff=0.1, outdir=outdir)
            except Exception as e:
                log1="Enrichr Server No response: %s vs %s, %s, %s \n"%(treat, ctrl, domain, gl_type,)
                log2="the lenght of input gene list = %s \n"%(len(glist))
                print(log1, log2)
                # touch file error exists
                os.system("touch  %s"%outfile)
                with open("temp/blacklist.enrichr.degs.%s_vs_%s.txt"%(treat, ctrl),'a') as black:
                    black.write(log1)
                    black.write(log2)            
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
    
    cols_group = [col.lstrip("TPM.") for col in cols_ ]
    cols  = [col for col, group in zip(cols_, cols_group) if group.startswith(treat)] +\
            [col for col, group in zip(cols_, cols_group) if group.startswith(ctrl)]

    col2 = ['gene_name']+ cols
    cls_vec = [g.split(".")[1] for g in cols]
    
    # run gsea
    for domain in go:
        outdir="GO/GSEA_%s/%s"%(outGSEAname, domain)
        outfile = "%s/gseapy.gsea.gene_sets.report.csv"%outdir
        #skip plotting while file exists
        if os.path.isfile(outfile): continue
        try:
            gs = gp.gsea(data=sig_deg[col2], gene_sets=domain, cls=cls_vec,
                         min_size=10, max_size=500, outdir=outdir)
        except:
            log1="Oops...%s_vs_%s: skip GSEA plotting for %s, please adjust paramters for GSEA input.\n"%(treat, ctrl, domain)
            log2="the lenght of input degs = %s \n"%sig_deg[col2].shape[0]  
            print(log1, log2)     
            os.system("touch %s/gseapy.gsea.gene_sets.report.csv"%outdir)
            with open("temp/blacklist.gsea.degs.%s_vs_%s.txt"%(treat, ctrl),'a') as black:
                black.write(log1)
                black.write(log2)


        """
        # delete empty dirs
        try:
            os.rmdir(outdir)
        except OSError as ex:
            if ex.errno == errno.ENOTEMPTY:
                print("deleting empty directory."
        """

gsea_enrichr(snakemake.input[0], snakemake.params['treat'], snakemake.params['ctrl'],
             snakemake.params['log2fc'], snakemake.params['padj'], snakemake.params['go'])
