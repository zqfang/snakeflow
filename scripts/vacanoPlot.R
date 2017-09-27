s1k=c("S1K_ESd0","S1K_IWP2d4","S1K_IWP2d8",
      "S1K_CT0.4d4","S1K_CT0.4d8",
      "S1K_CT1.0d4","S1K_CT1.0d8")
wt=c("Ctrl_ESd0","Ctrl_IWP2d4","Ctrl_IWP2d8",
      "Ctrl_CT0.4d4","Ctrl_CT0.4d8",
      "Ctrl_CT1.0d4","Ctrl_CT1.0d8")


outDIR="differential_expression_20170924"
dir.create(outDIR, recursive = TRUE)

for (i in 1:length(wt)) {
    treat=s1k[i]
    ctrl=wt[i]
    outRES=paste0("differential_expression_20170924/diff_",treat,"_vs_",ctrl,"_results.txt")
    res <- results(dds, contrast=c("condition", treat, ctrl))
    resOrdered <- res[order(res$padj),]
    resOrdered = as.data.frame(resOrdered)
    write.table(resOrdered, file=outRES, quote=F, sep="\t")
    
    outDIR="differential_expression_20170924/diff"
    #MAplot pdf
    outMA = paste(outDIR, treat, "vs", ctrl,"MAplot.pdf",sep="_")
    pdf(outMA, width = 5, height = 5)
    plotMA(res, ylim=c(-5,5), main=paste(treat, "vs", ctrl,sep="_"))
    dev.off()
    
    #MAplot png
    outMA = paste(outDIR, treat, "vs", ctrl,"MAplot.png",sep="_")
    png(outMA, width = 5, height = 5, units = 'in', res = 600)
    plotMA(res, ylim=c(-5,5), main=paste(treat, "vs", ctrl,sep="_"))
    dev.off()
}


# vacanoplot
# Download the data from github (click the "raw" button, save as a text file called "results.txt").
# https://gist.github.com/stephenturner/806e31fce55a8b7175af
res <- read.table("results.txt", header=TRUE)
head(res)

# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(res, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=.8))