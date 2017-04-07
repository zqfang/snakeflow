do_something <- function(tx2gene, out_tpm, out_counts, out_image, threads, sample_ids) {
  # R code
  library(tximport)
  library(readr)

  tx2gene <- read.table(tx2gene, header= T, sep="\t", stringsAsFactors = F)
  samples <- unlist(strsplit(sample_ids,","))

  salmon.files <- file.path('salmon',samples, "quant.sf")
  names(salmon.files) <- samples
  all(file.exists(salmon.files))

  #aggregate transcripts to genes 
  txi.salmon <- tximport(salmon.files, type = "salmon", 
                         tx2gene = tx2gene, reader = read_tsv,
                         countsFromAbundance = "lengthScaledTPM")

  salmon.counts<- txi.salmon$counts
  salmon.counts<- as.data.frame(salmon.counts)
  #salmon.counts$gene_name<- rownames(salmon.counts)
  write.table(salmon.counts, out_counts, sep="\t", quote=F)

  salmon.TPM<- txi.salmon$abundance
  salmon.TPM<- as.data.frame(salmon.TPM)
  #salmon.TPM$gene_name<- rownames(salmon.TPM)
  write.table(salmon.TPM, out_tpm, sep="\t", quote=F)
  save(txi.salmon, file=out_image)
}

do_something(snakemake@input[['tx2gene']], snakemake@output[['tpm']], snakemake@output[['counts']], snakemake@output[['image']], snakemake@threads, snakemake@params[['ids']])

