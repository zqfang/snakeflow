snake_heatmap <- function(degstab, ntdimage, treat, ctrl, padj, top)
{
    suppressMessages(library("pheatmap"))
    suppressMessages(library("RColorBrewer"))
    suppressMessages(library("DESeq2"))
    suppressMessages(library("genefilter"))
    #load ntd, group, df
    load(ntdimage)
    #load("temp/deseq2.dds.RData")
    res <- results(dds, contrast=c("condition", treat, ctrl))
    #color palete for heatmap
    heatcols <- colorRampPalette(rev(brewer.pal(10, "RdBu")))(255)
    #heatcols=colorRampPalette(c("royalblue4", "white", "red3"))(50)

    outDIR <- paste0("differential_expression/diff_", treat, "_vs_", ctrl, "/diff")
    #select groups for plotting
    ntd_cols <- ifelse(group == treat | group == ctrl, TRUE, FALSE)

    #TopGenes
    topGenes <- head(order(res$padj), top)
    #all Degs
    degs <- which(res$padj < padj)


     #MAplot
     #outMA = paste(outDIR, comb[i,2], "vs", comb[i,1],"MAplot.pdf",sep="_")
     #pdf(outMA, width = 5, height = 5)
     #plotMA(res, ylim=c(-5,5))
     #dev.off()

    #pdf
    pheatmap(ntd[topGenes, ntd_cols], scale = "row", cluster_rows=T, show_rownames=T,
           cluster_cols=T, annotation_col = df, cellwidth = 15, fontsize = 8,
           color=heatcols, border_color = FALSE,
           main = paste(treat, "vs", ctrl, sep="_"),
           filename = paste(outDIR, treat, "vs", ctrl,"top20genes.pdf",sep="_"))
    #png
    pheatmap(ntd[topGenes,ntd_cols], scale = "row", cluster_rows=T, show_rownames=T,
           cluster_cols=T, annotation_col = df, cellwidth = 15, fontsize = 8,
           color=heatcols, border_color = FALSE,
           main = paste(treat, "vs", ctrl, sep="_"),
           filename = paste(outDIR, treat, "vs", ctrl,"top20genes.png",sep="_"))
    #save results ammong groups comparasion.
    #pdf
    pheatmap(ntd[topGenes,], scale = "row", cluster_rows=T, show_rownames=T,
       cluster_cols=T, annotation_col = df, cellwidth = 15, fontsize = 8,
       color=heatcols, border_color = FALSE,
       main = paste(treat, "vs", ctrl, sep="_"),
       filename = paste(outDIR, treat, "vs", ctrl, "groups.top20genes.pdf",sep="_"))
    #png
    pheatmap(ntd[topGenes,], scale = "row", cluster_rows=T, show_rownames=T,
       cluster_cols=T, annotation_col = df, cellwidth = 15, fontsize = 8,
       color=heatcols, border_color = FALSE,
       main = paste(treat, "vs", ctrl,sep="_"),
       filename = paste(outDIR, treat, "vs", ctrl, "groups.top20genes.png",sep="_"))

     # if no significant genens found
     if (length(degs) < 1) {

          print(paste(treat, "vs", ctrl,"has no significant degs when padj =",padj, sep=" "))
          cmd1 <-paste(outDIR, treat, "vs", ctrl, "all.degs.pdf",sep="_")
          cmd2 <-paste(outDIR, treat, "vs", ctrl, "all.degs.png",sep="_")
          cmd3 <-paste(outDIR, treat, "vs", ctrl, "groups.all.degs.pdf",sep="_")
          cmd4 <-paste(outDIR, treat, "vs", ctrl, "groups.all.degs.png",sep="_")

          system(paste("touch", cmd1, sep=" "))
          system(paste("touch", cmd2, sep=" "))
          system(paste("touch", cmd3, sep=" "))
          system(paste("touch", cmd4, sep=" "))

     } else {
          # bug: Pheatmap error Error in hclust(d, method = method) :
          # NA/NaN/Inf in foreign function call (arg 11)
          # remove rows with rowSds() == 0 fixed this bug
          dat1 <- ntd[degs, ntd_cols]
          dat1 <- dat1[rowSds(dat1) != 0,]
          dat2 <- ntd[degs,]
          dat2 <- dat2[rowSds(dat2) != 0,]

          #pdf
          #add scale = "row"
          pheatmap(dat1, scale = "row", cluster_rows=T, show_rownames=F,
              cluster_cols=T, annotation_col = df,
              cellwidth = 15, fontsize = 8,
              color=heatcols, border_color = FALSE,
              main = paste(treat, "vs", ctrl, sep="_"),
              filename = paste(outDIR, treat, "vs", ctrl, "all.degs.pdf",sep="_"))
          #png
          pheatmap(dat1, scale = "row", cluster_rows=T, show_rownames=F,
              cluster_cols=T, annotation_col = df,
              cellwidth = 15, fontsize = 8, color=heatcols,
              main = paste(treat, "vs", ctrl, sep="_"),
              filename = paste(outDIR, treat, "vs", ctrl, "all.degs.png",sep="_"))

          #pdf
          pheatmap(dat2, scale = "row", cluster_rows=T, show_rownames=F,
              cluster_cols=T, annotation_col = df,
              cellwidth = 15, fontsize = 8,
              color=heatcols, border_color = FALSE,
              main = paste(treat, "vs", ctrl, sep="_"),
              filename = paste(outDIR, treat, "vs", ctrl, "groups.all.degs.pdf",sep="_"))
          #png
          pheatmap(dat2, scale = "row", cluster_rows=T, show_rownames=F,
              cluster_cols=T, annotation_col = df,
              cellwidth = 15, fontsize = 8,
              color=heatcols, border_color = FALSE,
              main = paste(treat, "vs", ctrl, sep="_"),
              filename = paste(outDIR, treat, "vs", ctrl, "groups.all.degs.png",sep="_"))
  }

}
snake_heatmap(snakemake@input[['degstab']], snakemake@input[['image']], snakemake@params[['treat']],
              snakemake@params[['ctrl']],snakemake@params[['padj']],
              snakemake@params[['topgene']])
