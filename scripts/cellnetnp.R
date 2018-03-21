#CellNet protocol
setwd("~/projects/CellNet")
library(CellNet)
cn_setup(local = TRUE)
iFileHuman = ""
fetchIndexHandler(destination = "~/projects/CellNet", species = "human",iFile=iFileHuman)
# load metadata
# download.file("https://s3.amazonaws.com/CellNet/rna_seq/human/examples/SRP043684/st_SRP043684_example.rda","st_SRP043684_example.rda")
# stQuery<-utils_loadObject("st_SRP043684_example.rda")

####### Bugs ######
## if you use your own fastq file in your computor, you should note::
## metadata fnameCol must endswith ".fastq", or the source code report bugs (source code use "fastq" to split filename and suffix)
## metadata must contained column sra_id, use your own unique id if use local files.

stQuery<-read.csv("./CellNetMeta2.csv", header = T, stringsAsFactors=F)
# pathToSalmonIndex
iFileHuman <-"salmon.index.human.052617"
pathToSalmon <-"~/miniconda/bin"
# refDir is where salmonIndex located
expList<-cn_salmon(stQuery,
                   salmonIndex=iFileHuman,
                   refDir="ref",
                   fnameCol = "fname",
                   geneTabfname="geneToTrans_Homo_sapiens.GRCh38.80.exo_Jul_04_2015.R",
                   salmonPath=pathToSalmon)
# save expression table
fname<-paste0("expList_NP.03162018.rda")
save(expList,file=fname)
# load(fname)
# analysis query data

# fetch CellNet kernel
download.file("https://s3.amazonaws.com/CellNet/rna_seq/human/cnProc_RS_hs_Oct_25_2016.rda",
              dest="./cnProc_RS_hs_Oct_25_2016.rda")
cnProc<-utils_loadObject("cnProc_RS_hs_Oct_25_2016.rda")

# Apply CellNet to query data and save results.

cnRes1<-cn_apply(expList[['normalized']], stQuery, cnProc)
fname<-paste0("cnRes_NP.03162018.rda")
save(cnRes1, file=fname)
# load(fname)

# Plot C/T classification results.
pdf(file='hmclass_NPattern.pdf', width=7, height=5)
cn_HmClass(cnRes1)
dev.off()


# plot GRN statues
fname<-'grnstats_esc_subset_SRP043684.pdf'
bOrder<-c("esc_train", unique(as.vector(stQuery$description2)), "neuron_train")
cn_barplot_grnSing(cnRes1,cnProc,"esc", c("esc", "neuron"), bOrder,
                   sidCol="sra_id", dlevel="description2")
ggplot2::ggsave(fname, width=5.5, height=5)
dev.off()
fname<-'grnstats_neuron_subset_SRP043684.pdf'
bOrder<-c("esc_train", unique(as.vector(stQuery$description2)), "neuron_train")
cn_barplot_grnSing(cnRes1,cnProc,"neuron", c("esc", "neuron"), bOrder, sidCol=
                     "sra_id", dlevel='description2')
ggplot2::ggsave(fname, width=5.5, height=5)
dev.off()


# compute NISs
rownames(stQuery)<-as.vector(stQuery$sra_id)
tfScores<-cn_nis_all(cnRes1, cnProc, "neuron")
fname='nis_neuron_subset_example_ctrlipsNeurons.pdf'
plot_nis(tfScores, "neuron", stQuery, "Control iPS neurons",
         dLevel="description2", limitTo=0)
ggplot2::ggsave(fname, width=4, height=12)
dev.off()
save(expList, file=fname)
