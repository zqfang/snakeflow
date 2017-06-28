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

PY_VERSION=3.6


#if [[ "$PY_VERSION" == "2.7" ]]; then
#        wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh;
#else
#        wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
#fi

#bash miniconda.sh -b -p /miniconda
#export PATH="/miniconda/bin:$PATH"

#setup bioconda channel.
conda config --set always_yes yes 
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels bioconda


log "install snakemake requirements"

#python2 packages
#conda install macs2 rseqc

source deactivate
conda install rseqc
name="snakeflow"
#clone a environment
#conda list --export > snakemake-env-packages.txt
#conda create -n snakeflow-clone --file snakemake-env-packages.txt

conda env list | grep -q $name && conda env remove -y -n $name
conda create -y -n $name  python=${PY_VERSION}
               # --file "snakemake-env-packages.txt" 

    
source activate $name

log "install snakemake-env-packages"
#python packages
conda install ipython cython numpy scipy pandas matplotlib seaborn snakemake  xlrd xlwt 
conda install -c bioninja gseapy

#other commandline tools
conda install hisat2 stringtie salmon star samtools bedtools fastqc graphviz multiqc
#R packages
#conda install bioconductor-deseq2 bioconductor-tximport bioconductor-readr bioconductor-ballgwon r-pheatmap r-ggrepel

#log "install gene annotation requirements"
wget http://bioconductor.org/biocLite.R
Rscript -e "source('biocLite.R');options(BioC_mirror='http://mirrors.ustc.edu.cn/bioc/');biocLite();biocLite(c('DESeq2','readr','pheatmap','tximport','ballgown','ggrepel'))"

log "all files are ready. Proceed to next step now."
