deseq2 <- function(txi_image, out_file, group, alias, threads) {
     # R code

     suppressMessages(library("ggplot2"))
     suppressMessages(library("RColorBrewer"))
     suppressMessages(library("gplots"))
     suppressMessages(library("pheatmap"))
     suppressMessages(library("DESeq2"))
     #suppressMessages(library('BiocParallel'))
     #register(MulticoreParam(threads))
     library(gtools)

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

     #set level
     ugr <- unique(group)
     ugr_len <- length(ugr)

     dds$condition <- relevel(dds$condition, ref=ugr[1])
     dds <- DESeq(dds)

     rld <- rlog(dds)
     vsd <- varianceStabilizingTransformation(dds)
     rlogMat <- assay(rld)
     vstMat <- assay(vsd)

     colnames(rld) <- alias
     
     comb <- combinations(ugr_len, 2, ugr)
     for (i in 1:dim(comb)[0])
     {
         #res <- results(dds, contrast=c("condition","treated","control"))
         res <- results(dds, contrast=c("condition", comb[i,2], comb[i,1]))
         resOrdered <- res[order(res$padj),]
         resOrdered = as.data.frame(resOrdered)
         outRES=paste("differential_expression/diff", comb[i,2], "vs", comb[i,1],"results.txt",sep="_")
         write.table(resOrdered, file=outRES, quote=F, sep="\t")

          #MA plot
          outMA = paste("differential_expression/diff", comb[i,2], "vs", comb[i,1],"MAplot.pdf",sep="_")
          pdf(outMA, width = 5, height = 5)     
          plotMA(res, ylim=c(-5,5))
          dev.off()

          #TopGenes
          #betas <- coef(dds)
          topGenes <- head(order(res$padj),20)
          df = data.frame(conditoin=group,row.names=colnames(rlogMat))
          outGenes = paste("differential_expression/diff", comb[i,2], "vs", comb[i,1],"top20genes.pdf",sep="_")
          pdf(outGenes, width = 4,height = 4)
          pheatmap(rlogMat[topGenes,], cluster_rows=T, show_rownames=T, cluster_cols=T, annotation_col = df)
          dev.off()
          
          #all Degs
          degs <- which(res$padj < 0.05)
          outDEGs = paste("differential_expression/diff", comb[i,2], "vs", comb[i,1], "all.degs.pdf",sep="_")
          pdf(outGenes, width = 5,height = 5)
          pheatmap(rlogMat[degs,], cluster_rows=T, show_rownames=F, cluster_cols=T,)


     } 


     #clustering plot
     hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
     distsRL <- dist(t(assay(rld)))
     mat <- as.matrix(distsRL)
     rownames(mat) <- colnames(mat) <- with(colData(dds), paste(alias, group, sep=":"))
    

     pdf("differential_expression/Samples.correlation.heatmap.pdf",width = 5, height = 5)
     hc <- hclust(distsRL)
     heatmap.2(mat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", 
               col = rev(hmcol), margin=c(13, 13))
     dev.off()

     #PCA plot.
     pdf("differential_expression/Samples.PCA.pdf",width = 4, height = 3)
     data <- plotPCA(rld, intgroup="condition", returnData=TRUE)
     percentVar <- round(100 * attr(data, "percentVar"))
     ggplot(data, aes(PC1, PC2, color=condition, label=rownames(data)))+
         geom_text(check_overlap = T,vjust = 0, nudge_y = 0.5) + geom_point(size=3) +
         xlab(paste0("PC1: ",percentVar[1],"% variance")) +
         ylab(paste0("PC2: ",percentVar[2],"% variance"))
     dev.off()


     #save dds for further processing
     save(dds, file=out_file)

}

deseq2(snakemake@input[['image']], snakemake@output[['res']], snakemake@params[['group']],snakemake@threads)


