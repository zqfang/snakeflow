def anno_genes(gene_anno, tx_anno, tpms, txtpms, tpm_out, txtpm_out, samples):
    # python code
    from pandas import read_table, concat

    annGE = read_table(gene_anno, index_col='gene_id')
    annGE.dropna(axis=1, how='all', inplace=True)
    annTX = read_table(tx_anno, index_col='tx_id')

    tpm = read_table(tpms, index_col=0)
    tpm = tpm[samples]

    txtpm = read_table(txtpms, index_col=0)
    txtpm = txtpm[samples]

    tpm.columns = ["TPM_%s"%a for a in samples]
    txtpm.columns = ["TPM_%s"%a for a in samples] # ["TPM_%s_%s"%(g,a) for a, g, s in zip(alias, group, samples)]

    mergeGene = concat([annGE, tpm], axis=1, join='inner')
    mergeGene.index.name ='gene_id'
    mergeGene.to_csv(tpm_out)

    mergeTransx = concat([annGE, tpm], axis=1, join='inner')
    mergeTransx.index.name ='transcript_id'
    mergeTransx.to_csv(txtpm_out)


anno_genes(snakemake.input[0], snakemake.input[1],
           snakemake.input[2], snakemake.input[3], 
           snakemake.output[0], snakemake.output[1], 
           snakemake.params['samples'])
