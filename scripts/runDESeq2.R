deseq2 <- function(txi_image, out_file, threads, group) {
     # R code

     library("DESeq2")
     library("ggplot2")
     library("RColorBrewer")
     library("gplots")
     library("pheatmap")

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

     #MA plot
     pdf("differential_expression/Sample.MAplot.pdf",width = 5, height = 5)     
     plotMA(res, ylim=c(-5,5))
     dev.off()


     rld <- rlog(dds)
     vsd <- varianceStabilizingTransformation(dds)
     rlogMat <- assay(rld)
     vstMat <- assay(vsd)

     #clustering plot
     hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
     distsRL <- dist(t(assay(rld)))
     mat <- as.matrix(distsRL)
     rownames(mat) <- colnames(mat) <- with(colData(dds), paste(rownames(colData(dds))))
    

     pdf("differential_expression/Sample.correlation.pdf",width = 5, height = 5)
     hc <- hclust(distsRL)
     heatmap.2(mat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", 
               col = rev(hmcol), margin=c(13, 13))
     dev.off()

     #PCA plot.
     pdf("differential_expression/Sample.PCA.pdf",width = 4, height = 3)
     data <- plotPCA(rld, intgroup="condition", returnData=TRUE)
     percentVar <- round(100 * attr(data, "percentVar"))
     ggplot(data, aes(PC1, PC2, color=condition)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance"))
     dev.off()


     #TopGenes
     betas <- coef(dds)
     topGenes <- head(order(res$padj),20)
     mat <- betas[topGenes, -c(1,2)]
     thr <- 3 
     mat[mat < -thr] <- -thr
     mat[mat > thr] <- thr
     pdf("differential_expression/Top20.genes.heatmap.pdf",width = 4,height = 3)
     pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101), cluster_col=FALSE)
     dev.off()

     #save dds for further processing
     save(dds, file="differential_expression/deseq2.dds.RData")

}

deseq2(snakemake@input[['image']], snakemake@output[['res']], snakemake@threads, snakemake@params[['group']])


