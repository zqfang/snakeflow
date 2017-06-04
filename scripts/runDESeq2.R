deseq2 <- function(txi_image, out_file, group, alias, threads) {
     # R code

     suppressMessages(library("ggplot2"))
     suppressMessages(library("RColorBrewer"))
     suppressMessages(library("gplots"))
     suppressMessages(library("pheatmap"))
     suppressMessages(library("DESeq2"))
     suppressMessages(library('BiocParallel'))
     register(MulticoreParam(threads))

     load(txi_image)
     #assign each sample to differrent group.
     group <- unlist(strsplit(group, " "))
     alias <- unlist(strsplit(alias, " "))
     sampleTable <- data.frame(condition = factor(group))
     #rownames(sampleTable) <- colnames(txi.salmon$counts)
     rownames(sampleTable) <- alias
     #ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=base_dir, design=~condition)
     #rownames(ddsHTSeq) <- gsub('\\.[0-9]+', '', rownames(ddsHTSeq))
     ## Filter genes with atleast 2 count
     #ddsHTSeq <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) > 1,  ]
     #colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c('control','knockdown'))   
 
     #run DESeq2
     dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~condition)
     dds <- DESeq(dds)

     ugr <- unique(group)
     group_num <- length(ugr)
     dds$condition <- relevel(dds$condition, ref=ugr[1])
     for (i in 2:group_num)
     {
         #res <- results(dds, contrast=c("condition","treated","control"))
         res <- results(dds, contrast=c("condition", ugr[i], ugr[i-1]))
         resOrdered <- res[order(res$padj),]
         resOrdered = as.data.frame(resOrdered)
         outRES=paste("differential_expression/diff", ugr[i], "vs",ugr[i-1],"results.txt",sep="_")
         write.table(resOrdered, file=outRES, quote=F, sep="\t")

          #MA plot
          outMA = paste("differential_expression/diff", ugr[i], "vs", ugr[i-1],"MAplot.pdf",sep="_")
          pdf(outMA, width = 5, height = 5)     
          plotMA(res, ylim=c(-5,5))
          dev.off()

          #TopGenes
          betas <- coef(dds)
          topGenes <- head(order(res$padj),20)
          mat <- betas[topGenes, -c(1,2)]
          thr <- 3 
          mat[mat < -thr] <- -thr
          mat[mat > thr] <- thr

          outGenes = paste("differential_expression/diff", ugr[i], "vs", ugr[i-1],"top20genes.pdf",sep="_")
          pdf(outGenes, width = 4,height = 3)
          pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101), cluster_col=FALSE)
          dev.off()

     } 





     rld <- rlog(dds)
     vsd <- varianceStabilizingTransformation(dds)
     rlogMat <- assay(rld)
     vstMat <- assay(vsd)

     #clustering plot
     hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
     distsRL <- dist(t(assay(rld)))
     mat <- as.matrix(distsRL)
     rownames(mat) <- colnames(mat) <- with(colData(dds), paste(rownames(colData(dds))))
    

     pdf("differential_expression/Samples.correlation.heatmap.pdf",width = 5, height = 5)
     hc <- hclust(distsRL)
     heatmap.2(mat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", 
               col = rev(hmcol), margin=c(13, 13))
     dev.off()

     #PCA plot.
     pdf("differential_expression/Samples.PCA.pdf",width = 4, height = 3)
     data <- plotPCA(rld, intgroup="condition", returnData=TRUE)
     percentVar <- round(100 * attr(data, "percentVar"))
     ggplot(data, aes(PC1, PC2, color=condition)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance"))
     dev.off()


     #save dds for further processing
     save(dds, file=out_file)

}

deseq2(snakemake@input[['image']], snakemake@output[['res']], snakemake@params[['group']],snakemake@threads)


