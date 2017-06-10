source("https://bioconductor.org/biocLite.R")
biocLite("EnsDb.Hsapiens.v75")

tx2gene = "~/genome/Human_GRCh38/gencode.v26.annotation.extracted.transx2gene.txt"
gene_anno="~/genome/Human_GRCh38/gencode.v26.annotation.extracted.genes.annotation.txt"

library(EnsDb.Hsapiens.v75)
txdb <- EnsDb.Hsapiens.v75
k <- keys(txdb, keytype = "TXID")
df_trans <- select(txdb, keys = k, keytype = "TXID", columns = c("TXID","TXNAME","GENEID","ENTREZID","GENENAME"))
df_genes <- select(txdb, keys = row.names(), keytype = "TXID", columns = c("TXID","TXNAME","GENEID","ENTREZID","GENENAME"))
setwd("~/public-seq/H170012-P001/trim_results/")
gene_exp = read.table("gene_expression/gene_expression.TPM.txt")
trans_exp = read.table("gene_expression/transcripts_expression.TPM.txt")

sampleInfo <-read.table("./snakeflow/example_sample_info.txt",header = F,sep = " ")
colnames(gene_exp) = sampleInfo[,2]
colnames(trans_exp) = sampleInfo[,2]
rownames(gene_exp) <- gsub('\\.[0-9]+', '', rownames(gene_exp))
rownames(trans_exp) <- gsub('\\.[0-9]+', '', rownames(trans_exp))
gene_exp.anno<-merge(gene_exp, df, by.x=rownames(gene_exp), by.y=df$GENEID, all.x=T,)

anno = mapIds(txdb, keys = rownames(gene_exp), 
              keytype = "GENEID", columns=c("GENEID","ENTREZID","GENENAME"), multiVals = "first")
df_genes <- select(txdb, keys = row.names(gene_exp), keytype = "GENEID", columns = c("GENEID","ENTREZID","GENENAME"))
df_trans <- select(txdb, keys = row.names(trans_exp), keytype = "TXID", columns = c("TXID","GENEID","ENTREZID","GENENAME"))

geneMerge = merge(df_genes, gene_exp, by.x='GENEID',by.y="row.names",all.y=T)
transxMerge = merge(df_trans, trans_exp, by.x='TXID',by.y="row.names",all.y=T)
write.csv(geneMerge, file="gene_expression/gene_expression.TPM.annotated.csv", quote = F, row.names = F)
write.csv(transxMerge, file="gene_expression/transcripts_expression.TPM.annotated.csv", quote = F, row.names = F)
