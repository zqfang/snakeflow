do_salmon <- function(tx2gene, out_tpm, outTrans_tpm, out_image, threads, sample_ids) {
  # R code
  library(tximport)
  library(readr)
  #suppressMessages(library('EnsDb.Hsapiens.v86'))

  #txdb <- EnsDb.Hsapiens.v86
  #k <- keys(txdb, keytype = "GENEID")
  #df <- select(txdb, keys = k, keytype = "GENEID", columns = c("TXID","GENEID"))
  #tx2gene <- df[, 2:1]  # tx ID, then gene ID

  tx2gene <- read.table(tx2gene, header= T, sep="\t", stringsAsFactors = F)
  samples <- unlist(strsplit(sample_ids,","))

  salmon.files <- file.path('salmon',samples, "quant.sf")
  names(salmon.files) <- samples
  all(file.exists(salmon.files))

  #aggregate transcripts to genes, and extract raw read counts  
  #add ignoreTxVersion to remove id versions
  txi.transcripts <- tximport(salmon.files, type = "salmon", 
                              txOut = TRUE, tx2gene = tx2gene,
                              ignoreTxVersion = TRUE) #countsFromAbundance="scaledTPM" if not use DESeq2friomtximport

  txi.salmon <- summarizeToGene(txi.transcripts, tx2gene ignoreTxVersion = TRUE) 
  # #countsFromAbundance="scaledTPM" if not use DESeq2friomtximport

  # #save raw counts after set tximport(..., countsFromAbundance="scaledTPM")
  # salmon.counts<- txi.salmon$counts
  # salmon.counts<- as.data.frame(salmon.counts)
  # write.table(salmon.counts, out_counts, sep="\t", quote=F)
  
  #save gene tpms
  salmon.TPM<- txi.salmon$abundance
  salmon.TPM<- as.data.frame(salmon.TPM)
  write.table(salmon.TPM, out_tpm, sep="\t", quote=F)
  #save transcripts tpms
  salmon.trans.TPM<- txi.transcripts$abundance
  salmon.trans.TPM<- as.data.frame(salmon.trans.TPM)
  write.table(salmon.trans.TPM, outTrans_tpm, sep="\t", quote=F)
   
  save(txi.salmon, file=out_image)
}

do_salmon(snakemake@input[['tx2gene']], snakemake@output[['tpm']], 
	  snakemake@output[['txtpm']],  
	  snakemake@output[['image']], 
    snakemake@threads, 
    snakemake@params[['ids']] #snakemake@output[['counts']],
    ) 


