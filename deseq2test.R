suppressMessages(library("pheatmap"))
suppressMessages(library("RColorBrewer"))
library(DESeq2)
library(NMF)
#load ntd, group, df
setwd("./projects/neural_patterning/")
load("./temp/deseq2.dds.RData")
load("./temp/deseq2.ntd.RData")
#color palete for heatmap
heatcols <- colorRampPalette(brewer.pal(9, "RdBu"))(255)
#heatcols=colorRampPalette(c("royalblue4", "white", "red3"))(50)

treat= 'S21K_Untreatedd4'
ctrl = 'Ctrl_Untreatedd4'
outDIR <- paste0("differential_expression/diff_", treat, "_vs_", ctrl, "/diff")
#select groups for plotting
ntd_cols <- ifelse(group ==treat | group == ctrl, TRUE, FALSE)

#TopGenes
res <-results(dds, contrast=c('condition',treat, ctrl))
topGenes <- head(order(res$padj), 20)
#all Degs
degs <- which(res$padj < 0.05)

ntdMat = as.matrix(ntd[topGenes,ntd_cols])
aheatmap(ntdMat, color="-RdBu:100", scale="row", Rowv = "correlation")
#pdf
pheatmap(ntd[topGenes, ntd_cols], scale = "row", cluster_rows=T, show_rownames=T,
         cluster_cols=T, annotation_col = df, cellwidth = 15, fontsize = 12,
         color=heatcols, main = paste(treat, "vs", ctrl, sep="_"))
#png
pheatmap(ntd[topGenes,ntd_cols], scale = "row", cluster_rows=T, show_rownames=T,
         cluster_cols=T, annotation_col = df, cellwidth = 15, fontsize = 12,
         color=heatcols, main = paste(treat, "vs", ctrl, sep="_"),
         filename = paste(outDIR, treat, "vs", ctrl,"top20genes.png",sep="_"))
#save results ammong groups comparasion.
#pdf
pheatmap(ntd[topGenes,], scale = "row", cluster_rows=T, show_rownames=T,
         cluster_cols=T, annotation_col = df, cellwidth = 15, fontsize = 8,
         color=heatcols, main = paste(treat, "vs", ctrl, sep="_"),
         filename = paste(outDIR, treat, "vs", ctrl, "groups.top20genes.pdf",sep="_"))
#png
pheatmap(ntd[topGenes,], scale = "row", cluster_rows=T, show_rownames=T,
         cluster_cols=T, annotation_col = df, cellwidth = 15, fontsize = 8,
         color=heatcols, main = paste(treat, "vs", ctrl,sep="_"),
         filename = paste(outDIR, treat, "vs", ctrl, "groups.top20genes.png",sep="_"))

n <- 50; p <- 20
x <- abs(rmatrix(n, p, rnorm, mean=4, sd=1))
x[1:10, seq(1, 10, 2)] <- x[1:10, seq(1, 10, 2)] + 3
x[11:20, seq(2, 10, 2)] <- x[11:20, seq(2, 10, 2)] + 2
rownames(x) <- paste("ROW", 1:n)
colnames(x) <- paste("COL", 1:p)

aheatmap(x)

## Distance methods
aheatmap(x, Rowv = "correlation")
aheatmap(x, Rowv = "man") # partially matched to 'manhattan'
aheatmap(x, Rowv = "man", Colv="binary")

# Generate column annotations
annotation = data.frame(Var1 = factor(1:p %% 2 == 0, labels = c("Class1", "Class2")), Var2 = 1:10)
aheatmap(x, annCol = annotation)