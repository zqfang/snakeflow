# This script is for ploting a heatmap 
dat_fil = data[rowMeans(data)>1,]
z = log2(dat_fil+1)
y=data.matrix(z)

library(RColorBrewer)
library(gplots)
#heatmap color range
hcol= colorRampPalette(brewer.pal(9,"GnBu"))(100)
#this is for colside lable color, "7" used is to define the numbers of colors(e.g. clusters, tissues, cell types)
color=colorpanel(5,low="red",mid="purple",high="green") 

# sacle your expression data fisrt, calculate the dist and clustering
mydatascale=t(scale(t(y)))
hc=hclust(as.dist(1-cor(mydatascale, method="spearman")), method="complete")
plot(hc)

#define the numbers of clusters you wanted.
rect.hclust(hc, k=2,border="red")
myc2=cutree(hc, k=2)#k=N N is the number is clusters you want to generate

#assign colors to your clusters,tissues,cell types......
mycolhr=sample(color)
mycolhr=mycolhr[as.vector(myc2)]



# use heatmap.2 to plot the same dengrogram with "hc",row_dendrogram using row_cluster(see below)
heatmap.2(mydatascale,Colv=as.dendrogram(hc),Rowv=as.dendrogram(row_cluster),col=hcol,scale="none",trace="none",labCol=colnames(y),labRow=F,density.info="none",symkey=F,margins=c(5,5),ColSideColors=mycolhr)
# define the clustering and dist method your own

hclustfunc <- function(x) hclust(x, method="average")
distfunc <- function(x) dist(x,method="euclidean")
heatmap.2(mydatascale,col=hcol,scale="none",trace="none",hclust=hclustfunc,distfun=distfunc,labCol=colnames(y),labRow=F,density.info="none",symkey=F,margins=c(5,5),ColSideColors=mycolhr)

# or
row_distance = dist(mydatascale, method = "euclidean")
row_cluster = hclust(row_distance, method = "average")
col_distance = dist(t(mydatascale), method = "manhattan")
col_cluster = hclust(col_distance, method = "ward.D2")
heatmap.2(mydatascale,Colv=as.dendrogram(col_cluster),Rowv=as.dendrogram(row_cluster),col=hcol,scale="none",trace="none",labCol=colnames(y),labRow=F,density.info="none",symkey=F,margins=c(8,16),ColSideColors=mycolhr)

# If you don't want to reorder your columns
# If you want to plot redgreen plot, use col=redgreen(10) or col=greenred(10)

heatmap.2(mydatascale,col=hcol,Colv=FALSE,Rowv=as.dendrogram(row_cluster),scale="none",trace="none",labCol=colnames(y),labRow=F,density.info="none",symkey=F,margins=c(5,5),ColSideColors=mycolhr)



pr = prcomp(scale(t(y)))
plot(pr$x,col="white",main="PCA Plot")
text(pr$x[,1],pr$x[,2],labels = colnames(y))

biplot(pr,cex=c(1,0.5),main="Biplot",col=c("black","grey"))