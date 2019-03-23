library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ape)
setwd("./projects/SOX1_low_high/")
load("./temp/deseq2.dds.RData")
#rename
s1EGFPgroup = group
s1EGFPdds = dds
s1EGFPdf = df
s1EGFPntd = ntd

# load wt samples deseq2 results
load("~/projects/neural_patterning/temp/deseq2.dds.RData")

## extract low high degs
res <- results(s1EGFPdds, contrast=c("condition",'High' ,'Low'))
#res <- results(s1EGFPdds)
resOrdered <- res[order(res$padj),]
resOrdered = as.data.frame(resOrdered)
pheno = read.csv("~/projects/neural_patterning/wgcna/sampleinfo.wt.csv")
pheno2 = read.table("./sample_info_single.txt", sep = " ", comment.char = '#')
patt_pheno = pheno[pheno$treat != 'ES',]
patt_pheno_d8 = patt_pheno[patt_pheno$time == 'D8',]
indx = which(abs(res$log2FoldChange)> 1 & res$padj < 0.05)
# extract log normalized counts 
s1EGFP_degs_ntd = s1EGFPntd[indx,]
s1EGFP_degs_names = rownames(s1EGFP_degs_ntd)
## exact wt day 8 ntd
wt_ntd = ntd[s1EGFP_degs_names, c(as.character(patt_pheno_d8$alias))]
# combine two dataframe
wt_s1EGFP_ntd = cbind(wt_ntd, s1EGFP_degs_ntd)

# format colnames
coln = colnames(wt_s1EGFP_ntd)
# sufix = substring(coln, 7)
# prefix = substr(coln,1,5)
# prefix = gsub("CTRL1", "_A", prefix)
# prefix = gsub("CTRL2", "_B", prefix)
coln = gsub("CTRL", "WT", coln)
coln = gsub("Untreated", "CT0.0", coln)
colnames(wt_s1EGFP_ntd) = coln
group_names = c("IWP2","IWP2","CT0.0","CT0.0","CT0.4","CT0.4",
                "CT0.8","CT0.8","CT1.0","CT1.0","CT2.0","CT2.0",
                "CT3.0","CT3.0","CT4.0","CT4.0","CT4.0RA","CT4.0RA",
                "SOX1_low","SOX1_low","SOX1_low",
                "SOX1_high","SOX1_high","SOX1_high",)


### averaging replicates ntds
group_wt_ntd=data.frame(row.names = rownames(wt_s1EGFP_ntd))
for (i in seq(1,18,2)){
  print(coln[c(i,i+1)])
  temp = rowMeans(wt_s1EGFP_ntd[,coln[c(i,i+1)]])
  names(temp) = group_names[i]
  group_wt_ntd = cbind(group_wt_ntd,temp)
  
}
names(group_wt_ntd) = group_names[seq(1,18,2)]

group_s1_ntd=data.frame(row.names = rownames(wt_s1EGFP_ntd))
temp1 = rowMeans(s1EGFP_degs_ntd[,1:3])
temp2 = rowMeans(s1EGFP_degs_ntd[,4:6])
group_s1_ntd = cbind(group_s1_ntd,temp1,temp2)
names(group_s1_ntd) = c("Low","High")
group_wt_s1_ntd = cbind(group_wt_ntd, group_s1_ntd)

### plotting correlation


# save data
write.csv(wt_s1EGFP_ntd, "S1EGFP.HighLow.DEGS.ntd.csv",quote = F)
# correlation to wt samples
corrs_pearson = cor(wt_s1EGFP_ntd)
corrs_spearman = cor(wt_s1EGFP_ntd, method='spearman')

corrs_pearson_mean = cor(group_wt_s1_ntd)
corrs_spearman_mean  = cor(group_wt_s1_ntd, method='spearman')


pheatmap(corrs_pearson, main="Pearson")
pheatmap(corrs_spearman, main="Spearman")
pheatmap(corrs_pearson_mean, main="Pearson")
pheatmap(corrs_spearman_mean, main="Spearman")

cmap =  colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100)
pheatmap(corrs_pearson, color = cmap, main="Pearson", filename = "S1EGFP.HighLow.and.WT.Day8.pearson.pdf" )
pheatmap(corrs_spearman, color = cmap,main="Spearman", filename="S1EGFP.HighLow.and.WT.Day8.spearman.pdf")
pheatmap(corrs_pearson_mean, color = cmap, main="Pearson", filename = "S1EGFP.HighLow.and.WT.Day8.pearson.average.pdf" )
pheatmap(corrs_spearman_mean, color = cmap,main="Spearman", filename="S1EGFP.HighLow.and.WT.Day8.spearman.average.pdf")

# hclust
s21.dist=dist(t(wt_s1EGFP_ntd),method="euclidean")
out.hclust=hclust(s21.dist,method="complete") #根据距离聚类
plot(out.hclust)

y = wt_s1EGFP_ntd
y = y[rowSums(y)>0,]
mydatascale=t(scale(t(y)))
hc=hclust(as.dist(1-cor(mydatascale, method="spearman")), method="complete")
plot(hc)

# convert to dendrogram
hcd = as.dendrogram(hc)
plot(hcd, horiz = TRUE, type ="triangle", xlab='Height') # type="rectangle"
# Unrooted
plot(as.phylo(hc), type = "unrooted", cex = 0.6,
     no.margin = TRUE)
plot(as.phylo(hc), type = "fan")
# Radial
plot(as.phylo(hc), type = "radial")


## scatter plot of degs
indx = which(abs(res$log2FoldChange)> 1 & res$padj < 0.05)
load("./temp/txi.salmon.RData")
degs_tpm = txi.salmon$abundance[indx,]


## vocano plot
#1.导入数据包：
library(ggplot2)
# 2. 准备数据
res2 = as.data.frame(res)
res3 = res2[complete.cases(res2),]

## anotate genenames
library('EnsDb.Hsapiens.v86')
#remove tails(versions) of gene id
gene_id <- gsub('\\.[0-9a-zA-Z]+', '', rownames(res3))

#change the ensemble gene_id to gene_name for plotting
edb <- EnsDb.Hsapiens.v86
maps_names <- mapIds(edb, keys = gene_id, column="GENENAME",
                     keytype =  "GENEID", multiVals = "first")
res3$gene_name <- maps_names
res3$change = as.factor(ifelse(abs(res3$log2FoldChange)>1 & res3$padj < 0.0507,
                               ifelse(res3$log2FoldChange > 0,'Up','Down'),'Not'))
res3 = res3[complete.cases(res3),]
this_tile = paste0("The number of up genes are", nrow(res3[res3$change == 'Up',]),
                   "\nThe number of down genes are", nrow(res3[res3$change == 'Down',]))
## add gene name to vacano plot
shown=strsplit("SOX1 OTX1 OTX2 LMX1A LMX1B EN1 WNT3A WNT1 HOXA2 HOXA1 GBX2 PHOX2B EPHA3"," ")[[1]]
#shown=strsplit("SOX1 OTX1 OTX2 LMX1A LMX1B EN1 WNT3A WNT1 HOXA2 HOXA1 PHOX2B EPHA3"," ")[[1]]
res3$shown_name = as.factor(ifelse(res3$gene_name %in% shown, 'yes','no'))

res4 = res3[order(res3$padj),]
# 3.设置横轴和纵轴
library(ggrepel)
r04=ggplot(data = res3,aes(x=log2FoldChange,y=-1*log10(padj), color=change))
# 4.显示火山图
r04 = r04+geom_point(alpha=0.5, size=1.75,stroke = 0)+
         theme_set(theme_set(theme_bw(base_size=16)))+
         geom_text_repel(data=subset(res3, shown_name == 'yes'),
                         aes(label=gene_name),
                         segment.size  = 0.2,
                         segment.color = "grey50")+
         xlab(expression(log[2]('Fold Change')))+
         ylab(expression(-log[10]('adjusted p-value')))+
         scale_color_manual(values=c('blue','black','red'))+
         ggtitle(this_tile) + 
         theme(plot.title = element_text(size=16,hjust=0.5, face="bold"),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank())
         #geom_hline(yintercept=1.3)+geom_vline(xintercept=c(-1,1))

print(r04)


### move label to left and right
r06=ggplot(data = res3,aes(x=log2FoldChange,y=-1*log10(padj), color=change))
# 4.显示火山图
r06 = r06+geom_point(alpha=0.5, size=1.75,stroke = 0)+
  theme_set(theme_set(theme_bw(base_size=16)))+
  geom_text_repel(data=subset(res3, shown_name == 'yes'& log2FoldChange < 0),
                  aes(label=gene_name),
                  xlim  = c(NA, -5),
                  segment.size  = 0.2,
                  segment.color = "grey50")+
  geom_text_repel(data=subset(res3, shown_name == 'yes'& log2FoldChange > 0),
                  aes(label=gene_name),
                  xlim  = c(5, NA),
                  segment.size  = 0.2,
                  segment.color = "grey50")+
  xlab(expression(log[2]('Fold Change')))+
  ylab(expression(-log[10]('adjusted p-value')))+
  scale_color_manual(values=c('blue','black','red'))+
  ggtitle(this_tile) + 
  theme(plot.title = element_text(size=16,hjust=0.5, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
#geom_hline(yintercept=1.3)+geom_vline(xintercept=c(-1,1))

print(r06)


# 6.设置坐标轴范围和标题 #横坐标范围：xlim()，纵坐标范围：ylim()函数，添加标签：labs(title=“..”,x=“..”,y=“..”)函数r03xy=addcolor+xlim(-4,4)+ ylim(0,30)+ labs(title="Volcanoplot",x=expression(log2(log2FoldChange)),y=expression(-log10(pvalue)))

# 8.添加阈值线（y轴截距，横坐标范围）
addline=volcano+geom_hline(yintercept=1.3)+geom_vline(xintercept=c(-1,1))
addline

# 9.保存图片(名称，图，宽，高)：
ggsave("volcano8.png",volcano,width=8,height=8)


##res3$sign
res3$sign <- ifelse(res3$padj < 0.05 & abs(res3$log2FoldChange) > 2,res3$gene_name,'')
for (i in 1:17){
  print(coln[c(i,i+1)])
  print("\n")
}
