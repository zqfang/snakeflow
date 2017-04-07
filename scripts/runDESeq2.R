deseq2 <- function(txi_image, out_file, threads, group) {
	    # R code

          library(DESeq2)
          load(txi_image)

          #assign each sample to differrent group.
          group=unlist(strsplit(group, " "))
          sampleTable <- data.frame(condition = factor(group))
          rownames(sampleTable) <- colnames(txi.salmon$counts)

          #run DESeq2
          dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~condition)
          dds <- DESeq(dds)
          res <- results(dds, addMLE=TRUE)
          resOrdered <- res[order(res$padj),]
          resOrdered = as.data.frame(resOrdered)
          write.table(resOrdered, file=out_file,sep="\t")

}

deseq2(snakemake@input[['image']], snakemake@output[['res']], snakemake@threads, snakemake@params[['group']])


