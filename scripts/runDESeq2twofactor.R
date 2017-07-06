suppressMessages(library("ggplot2"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("gplots"))
suppressMessages(library("pheatmap"))
suppressMessages(library("DESeq2"))
suppressMessages(library("tximport"))
suppressMessages(library('clusterProfiler'))
library('EnsDb.Hsapiens.v86')
library('org.Hs.eg.db')
library('topGO')
library('Rgraphviz')

setwd("./public-seq/H170013-P001/")

sampleInfo = read.table("./sample_info_single.txt", sep=" ")

condition=c("DMSO","DMSO","DMSO","RA","RA","RA","DMSO","DMSO","DMSO","RA","RA","RA")
treatment=c("PBS","PBS","PBS","PBS","PBS","PBS","RK","RK","RK","RK","RK","RK")
sampTab = data.frame(condition=condition, treatment=treatment,row.names = sampleInfo[,1])
load("./diff_RKA_vs_RA_specific/deseq2.dds.RKA_vs_RA_specific.RData")


#run DESeq2
#two factor design
dds <- DESeqDataSetFromTximport(txi.salmon, sampTab, ~condition+treatment+condition:treatment)
dds <- DESeq(dds)
#see ?results for details
RA_RKvsPBS <- results(dds, name="conditionRA.treatmentRK")
res = RA_RKvsPBS
resOrdered <- res[order(res$padj),]
resOrdered = as.data.frame(resOrdered)
outRES="differential_expression/diff_RA_RKvsPBS_specific_changed.txt"
write.table(resOrdered, file=outRES, quote=F, sep="\t")

#remove tails(versions) of gene id
rownames(resOrdered) <- gsub('\\.[0-9]+', '', rownames(resOrdered))

#annotate genes
edb <- EnsDb.Hsapiens.v86
anno_cols=c( "TXBIOTYPE","ENTREZID","SYMBOL")

resOrdered_new = resOrdered
for (i in 1:length(anno_cols))
{
  maps_anno <- mapIds(edb, keys = rownames(resOrdered), column=anno_cols[i],
                 keytype =  "GENEID", multiVals = "first")  
  resOrdered_new <- cbind(maps_anno, resOrdered_new)
  names(resOrdered_new)[1] = anno_cols[i]
}

#sig genes
res_sig =subset( resOrdered_new, padj < 0.05 )
res_sig2 =subset( res_sig, log2FoldChange > 1 | log2FoldChange < -1) 

sort_uniq=res_sig2[order(res_sig2$padj),]#按照矫正p值排序
sort_uniq=res_sig[order(res_sig$padj),]#按照矫正p值排序

#标记上下调基因
sort_uniq$up_down=ifelse(sort_uniq$log2FoldChange > 0,'up','down')
write.csv(res_sig, file="./diff_RKA_vs_RA_specific/diff_RARK_vs_RA_specific_changed_annotated_sig.csv",quote = F)

#GO analysis
onto = c("MF", "BP", "CC")
ontoDB = 'org.Hs.eg.db'
outname="./GO_RKA_vs_RA_specific/diff_genes_RKA_vs_RA_specific.sig-"

#down regulated genes
diff_ENTREZID = sort_uniq[sort_uniq$log2FoldChange < 0 ,'ENTREZID']
diff_ENTREZID=na.omit(diff_ENTREZID)
for (ot in onto){
  ego <- enrichGO(gene=diff_ENTREZID,ontoDB,ont=ot, pvalueCutoff=0.05,readable=TRUE)
  write.csv(as.data.frame(ego),paste0(outname,ot,".down.csv"),row.names =F)
  
  pdf(paste0(outname,ot,".down.pdf"))
  print(dotplot(ego,showCategory=10,font.size=12))
  dev.off() 
  
  pdf(paste0(outname,ot,".topGO.down.pdf"))
  print(plotGOgraph(ego))
  dev.off()
}
#KEGG
ekk <- enrichKEGG(gene=diff_ENTREZID,organism = "hsa",pvalueCutoff=0.05)
write.csv(as.data.frame(ekk),paste0(outname,"KEGG.down.csv"),row.names =F)
pdf(paste0(outname,"KEGG.down.pdf"))
print(dotplot(ekk,showCategory=20,font.size=12))
dev.off()

#up-regulated genes
diff_ENTREZID = sort_uniq[sort_uniq$log2FoldChange > 0 ,'ENTREZID']
diff_ENTREZID=na.omit(diff_ENTREZID)
for (ot in onto) {
  ego <- enrichGO(gene=diff_ENTREZID,ontoDB,ont=ot, pvalueCutoff=0.05,readable=TRUE)
  write.csv(as.data.frame(ego),paste0(outname,ot,".up.csv"),row.names =F)
  pdf(paste0(outname,ot,".up.pdf"))
  print(dotplot(ego,showCategory=10,font.size=12))
  dev.off() 
  
  pdf(paste0(outname,ot,".topGO.up.pdf"))
  print(plotGOgraph(ego))
  dev.off()
}
#KEGG
ekk <- enrichKEGG(gene=diff_ENTREZID,organism = "hsa",pvalueCutoff=0.05)
write.csv(as.data.frame(ego),paste0(outname,"KEGG.up.csv"),row.names =F)

pdf(paste0(outname,ot,"KEGG.up.pdf"))
print(dotplot(ekk,showCategory=20,font.size=12))
dev.off()

dotplot(ekk,showCategory=20,font.size=12)
clusterProfiler::plotGOgraph(ekk)


## categorySize can be scaled by 'pvalue' or 'geneNum'
cnetplot(ego, categorySize="pvalue", foldChange=sort_uniq[sort_uniq$log2FoldChange < 0 ,'log2FoldChange'])

enrichMap(ego)


rld <- rlog(dds)
# this gives log2(n + 1)
ntd <- normTransform(dds)
#ntd2 <- t(scale(t(as.matrix(assay(ntd)))))
colnames(ntd) = colnames(rld)
ntd = assay(ntd)
rownames(ntd) <- gsub('\\.[0-9]+', '', rownames(ntd))
#colnames(ntd2) = colnames(rld)
edb <- EnsDb.Hsapiens.v86
maps_names <- mapIds(edb, keys = rownames(ntd), column="GENENAME",
                     keytype =  "GENEID", multiVals = "first") 
rownames(ntd) <- maps_names



#heatmap
df = data.frame(condition=condition, treatment=treatment, row.names=colnames(ntd))
degs <- which(res$padj < 0.05)


pheatmap(ntd[degs,], scale = "row", cluster_rows=T, show_rownames=F,
         cluster_cols=T, annotation_col = df,cellwidth = 15, fontsize = 12, 
         main ="RKA_vs_RK",
         filename = paste0(outname,"all.degs.pdf"))


topGenes <- head(order(res$padj),20)
pheatmap(ntd[topGenes,], scale = "row", cluster_rows=T, show_rownames=T,
         cluster_cols=T, annotation_col = df,cellwidth = 15, fontsize = 12,
         main ="RKA_vs_RK",
         filename = paste0(outname,"top20genes.pdf"))
