o#!bin/bash

#this script used for initialize working directory and 
#install requirements files for runing snakemake workflows.
#the indpendencies for the tools will be automatically installed by conda.
set -e

log () {
    echo
    echo "[`date`] Step: $1 "
    echo
}

PY_VERSION=3.7
#setup bioconda channel.
conda config --set always_yes yes 
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge


log "install snakemake requirements"

name="snakeflow"
#clone a environment
#conda list --export > snakemake-env-packages.txt
#conda create -n snakeflow-clone --file snakemake-env-packages.txt

conda env list | grep -q $name && conda env remove -y -n $name
conda create -y -n $name  python=${PY_VERSION}  #--file "snakemake-env-packages.txt" 

    
source activate $name

log "install snakemake-env-packages"
#python packages
conda install ipython cython numpy scipy pandas matplotlib seaborn snakemake  xlrd xlwt 

#other commandline tools
conda install hisat2 stringtie salmon star samtools bedtools fastqc graphviz multiqc deeptools gseapy rseqc macs2
#R packages
#conda install bioconductor-deseq2 bioconductor-tximport bioconductor-readr bioconductor-ballgwon r-pheatmap r-ggrepel

#log "install gene annotation requirements"
wget http://bioconductor.org/biocLite.R
Rscript -e "BiocManager::install(c('DESeq2','readr','pheatmap','tximport','ballgown','ggrepel','topGO','clusterProfiler','org.Hs.eg.db','Rgraphviz','EnsDb.Hsapiens.v86'))"

log "all files are ready. Proceed to next step now."
