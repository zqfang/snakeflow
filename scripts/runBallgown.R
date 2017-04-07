ballgown <- function(ids, transx, genex, group) {
	    # R code

          library("ballgown")
          samples = unlist(strsplit(ids,","))
          bg = ballgown(samples = samples, meas='all')
          whole_tx_table = texpr(bg, 'all')
          gene_expression = gexpr(bg)
          write.csv(whole_tx_table,file=transx)
          write.csv(gene_expression,file=genex)

}

ballgown(snakemake@params[['ids']], snakemake@output[['transx']], snakemake@output[['genex']], snakemake@threads)


