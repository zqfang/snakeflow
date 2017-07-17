snake_heatmap <-  <- function(degstab, ntdimage, treat, ctrl, padj, top) 
{
    suppressMessages(library("pheatmap"))
    suppressMessages(library("RColorBrewer"))
    #load ntd, group, df
    load(ntdimage)

    #color palete for heatmap
    heatcols <- colorRampPalette(brewer.pal(10, "RdBu"))(255)
    #heatcols=colorRampPalette(c("royalblue4", "white", "red3"))(50)

    outDIR = paste0("differential_expression/diff_", treat, "_vs_", ctrl, "/diff")
    #select groups for plotting
    ntd_cols = ifelse(group ==treat | group == ctrl, TRUE, FALSE)

    #TopGenes
    res = read.table(degstab, stringsAsFactors=F )
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
           cluster_cols=T, annotation_col = df, cellwidth = 15, fontsize = 12,
           color=heatcols, main = paste(treat, "vs", ctrl, sep="_"),
           filename = paste(outDIR, treat, "vs", ctrl,"top20genes.pdf",sep="_"))
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

     # if no significant genens found    
     if (length(degs) < 1) {

          print(paste(treat, "vs", ctrl,"has no significant degs when padj =",padj, sep=" "))
          cmd1 <-paste(outDIR, treat, "vs", ctrl, "all.degs.pdf",sep="_")
          cmd2 <-paste(outDIR, treat, "vs", ctrl, "all.degs.png",sep="_")
          cmd3 <-paste(outDIR, treat, "vs", ctrl, "groups.all.degs.pdf",sep="_"))
          cmd4 <-paste(outDIR, treat, "vs", ctrl, "groups.all.degs.png",sep="_"))

          system(paste("touch", cmd1, sep=" ")
          system(paste("touch", cmd2, sep=" ")
          system(paste("touch", cmd3, sep=" ")
          system(paste("touch", cmd4, sep=" ")

     } else {
          #pdf
          pheatmap(ntd[degs, ntd_cols], scale = "row", cluster_rows=T, show_rownames=F, 
              cluster_cols=T, annotation_col = df, 
              cellwidth = 15, fontsize = 12, color=heatcols,
              main = paste(treat, "vs", ctrl, sep="_"),
              filename = paste(outDIR, treat, "vs", ctrl, "all.degs.pdf",sep="_"))
          #png
          pheatmap(ntd[degs, ntd_cols], scale = "row", cluster_rows=T, show_rownames=F, 
              cluster_cols=T, annotation_col = df, 
              cellwidth = 15, fontsize = 12, color=heatcols,
              main = paste(treat, "vs", ctrl, sep="_"),
              filename = paste(outDIR, treat, "vs", ctrl, "all.degs.png",sep="_"))

          #pdf
          pheatmap(ntd[degs,], scale = "row", cluster_rows=T, show_rownames=F, 
              cluster_cols=T, annotation_col = df, 
              cellwidth = 15, fontsize = 8, color=heatcols,
              main = paste(treat, "vs", ctrl, sep="_"),
              filename = paste(outDIR, treat, "vs", ctrl, "groups.all.degs.pdf",sep="_"))
          #png
          pheatmap(ntd[degs,], scale = "row", cluster_rows=T, show_rownames=F, 
              cluster_cols=T, annotation_col = df, 
              cellwidth = 15, fontsize = 8, color=heatcols,
              main = paste(treat, "vs", ctrl, sep="_"),
              filename = paste(outDIR, treat, "vs", ctrl, "groups.all.degs.png",sep="_"))
  }
          
}
snake_heatmap(snakemake@input[['degstab']], snakemake@input[['image']], snakemake@params[['treat']],
              snakemake@params[['ctrl']],snakemake@params[['padj']],
              snakemake@params[['topgene']])
