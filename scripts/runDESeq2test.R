deseq2 <- function(txi_image, outdds, outntd, group, time, alias) {
  # R code
  
  suppressMessages(library("ggplot2"))
  suppressMessages(library("RColorBrewer"))
  suppressMessages(library("gplots"))
  suppressMessages(library("pheatmap"))
  suppressMessages(library("DESeq2"))
  suppressMessages(library('EnsDb.Hsapiens.v86'))
  #library("BiocParallel")
  #register(MulticoreParam(4))
  # and set DESeq() et.al with parallel=TRUE
  library(ggrepel)
  
  load(txi_image)
  #assign each sample to differrent group.
  group <- unlist(strsplit(group, " "))
  alias <- unlist(strsplit(alias, " "))
  sampleTable <- data.frame(condition = factor(group))
  rownames(sampleTable) <- colnames(txi.salmon$counts)
  
  #ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=base_dir, design=~condition)
  #rownames(ddsHTSeq) <- gsub('\\.[0-9]+', '', rownames(ddsHTSeq))
  ## Filter genes with atleast 2 count
  #ddsHTSeq <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) > 1,  ]
  
  #run DESeq2
  dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~condition)
  
  #set level
  ugr <- unique(group)
  ugr_len <- length(ugr)
  
  dds$condition <- relevel(dds$condition, ref=ugr[1])
  dds <- DESeq(dds)
  
  rld <- rlog(dds)
  colnames(rld) <- alias
  
  vsd <- varianceStabilizingTransformation(dds)
  #rlogMat <- assay(rld)
  #vstMat <- assay(vsd)
  
  # this gives log2(n + 1)
  ntd <- normTransform(dds)
  #ntd2 <- t(scale(t(as.matrix(assay(ntd)))))
  colnames(ntd) <- alias
  ntd <- assay(ntd)
  
  #remove tails(versions) of gene id
  rownames(ntd) <- gsub('\\.[0-9]+', '', rownames(ntd))
  
  #change the ensemble gene_id to gene_name for plotting
  edb <- EnsDb.Hsapiens.v86
  maps_names <- mapIds(edb, keys = rownames(ntd), column="GENENAME",
                       keytype =  "GENEID", multiVals = "first") 
  rownames(ntd) <- maps_names
  
  #annotate columns of heatmap
  df <- data.frame(treatment=group, row.names=colnames(ntd))
  
  #save dds for further processing
  save(dds, df,rld, vsd, ntd, maps_names, file=outdds)
  
  #save ntd
  save(ntd, df, group, file=outntd)  
  
  #clustering plot
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(255)
  distsRL <- dist(t(assay(rld)))
  mat <- as.matrix(distsRL)
  rownames(mat) <- colnames(mat) <- with(colData(dds), paste(alias, sep="")) #paste(alias, group, sep=":")
  
  pdf("differential_expression/Samples.correlation.heatmap.pdf", width = 8, height = 8)
  hc <- hclust(distsRL)
  heatmap.2(mat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", 
            col = rev(hmcol), margins=c(10, 10),
            main="Sample Correlation")
  dev.off()
  
  png("differential_expression/Samples.correlation.heatmap.png",res=300)
  heatmap.2(mat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", 
            col = rev(hmcol), margins=c(10, 10),
            main="Sample Correlation")
  dev.off()
  
  #PCA plot.
  pdf("differential_expression/Samples.PCA.pdf", width = 8, height = 8)
  data <- plotPCA(rld, intgroup="condition", returnData=TRUE)
  percentVar <- round(100 * attr(data, "percentVar"))
  #add geom_text(check_overlap = T, to remove overlap text)
  p <- ggplot(data, aes(PC1, PC2, color=condition, label=rownames(data)))
  p <- p+ geom_point(size=3) +
    ggtitle("Sampls PCA") +
    geom_text_repel(fontface = "bold")+ 
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))
  #you have to use print() when calling ggplot and save to pdf
  print(p)
  dev.off()
  
  png("differential_expression/Samples.PCA.png")
  print(p)
  dev.off()
  
  #save results for each group     
  comb <- t(combn(ugr,2))
  for (i in 1:dim(comb)[1])
  {
    #res <- results(dds, contrast=c("condition","treated","control"))
    res <- results(dds, contrast=c("condition", comb[i,2], comb[i,1]))
    resOrdered <- res[order(res$padj),]
    resOrdered = as.data.frame(resOrdered)
    
    #save results to outdir
    outDIR = paste0("differential_expression/diff_", comb[i,2], "_vs_", comb[i,1], "/diff")       
    #outRES=paste("differential_expression/diff", comb[i,2], "vs", comb[i,1],"results.txt",sep="_")
    outRES=paste(outDIR, comb[i,2], "vs", comb[i,1],"results.txt",sep="_")
    write.table(resOrdered, file=outRES, quote=F, sep="\t")
    
    #MAplot pdf
    outMA = paste(outDIR, comb[i,2], "vs", comb[i,1],"MAplot.pdf",sep="_")
    pdf(outMA, width = 5, height = 5)     
    plotMA(res, ylim=c(-5,5))
    dev.off()
    
    #MAplot png
    outMA = paste(outDIR, comb[i,2], "vs", comb[i,1],"MAplot.png",sep="_")
    png(outMA)     
    plotMA(res, ylim=c(-5,5))
    dev.off()
    
  } 
  
}

load("temp/txi.salmon.RData")

outdds="temp/deseq2.dds2.RData"
outntd="temp/deseq2.ntd.RData"
group=paste(sampleTab$V3,collapse=' ')
time=paste(sampleTab$V4,collapse=' ')
alias=paste(sampleTab$V2,collapse=' ')


pheatmap(mat,
         clustering_distance_rows=distsRL,
         clustering_distance_cols=distsRL,
         col=rev(hmcol))
