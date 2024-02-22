## see tutorial here: 
## https://www.jianshu.com/p/5ddef11ddbe9 # for ATAC
## https://www.jianshu.com/p/f849bd55ac27
## https://www.jieandze1314.com/post/cnposts/214/

library(DiffBind)
dbObj <- dba(sampleSheet="1.csv")
dbObj=dba.count(dbObj)
pdf(file="PCA_plot.pdf",width = 7,height = 7)
dba.plotPCA(dbObj, attributes=DBA_FACTOR, label=DBA_ID)
dev.off()
pdf(file="heatmap_plot.pdf",width = 7,height = 7)
plot(dbObj)
dev.off()


# Establishing a contrast 
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR,minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
#  summary of results
dba.show(dbObj, bContrasts=T)
#  overlapping peaks identified by the two different tools (DESeq2 and edgeR)
pdf(file="overlap_DESeq2_edgeR.pdf",width = 7,height = 7)
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS) # method 
dev.off()

# extract results
comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)


# EdgeR
out <- as.data.frame(comp1.edgeR)
write.table(out, file="results_edgeR.txt", sep="\t", quote=F, col.names = NA)
# DESeq2
out <- as.data.frame(comp1.deseq)
write.table(out, file="results_deseq2.txt", sep="\t", quote=F, col.names = NA)



# Create bed files for each keeping only significant peaks (p < 0.05)
# EdgeR
out <- as.data.frame(comp1.edgeR)
edge.bed <- out[ which(out$FDR < 0.05), 
                 c("seqnames", "start", "end", "strand", "Fold")]
write.table(edge.bed, file="results_edgeR_sig.bed", sep="\t", quote=F, row.names=F, col.names=F)

# DESeq2
out <- as.data.frame(comp1.deseq)
deseq.bed <- out[ which(out$FDR < 0.05), 
                  c("seqnames", "start", "end", "strand", "Fold")]
write.table(deseq.bed, file="results_deseq2_sig.bed", sep="\t", quote=F, row.names=F, col.names=F)
