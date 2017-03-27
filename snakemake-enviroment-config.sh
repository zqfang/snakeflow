#!bin/bash

#this script used for initialize working directory and 
#install requirements files for runing snakemake workflows.

set -e

log () {
    echo
    echo "[`date`] Step: $1 "
    echo
}


#if [[ "$TRAVIS_PYTHON_VERSION" == "2.7" ]]; then
#        wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh;
#else
#        wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
#    fi

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p /miniconda
export PATH="/miniconda/bin:$PATH"

#setup bioconda channel.
conda config --set always_yes yes 
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda


log "install snakemake requirements"
source deactivate
name="snakemake"
PY_VERSION=3.5

#clone a environment
#conda list --export > snakemake-env-packages.txt
#conda create -n snakeflow-clone --file snakemake-env-packages.txt

conda env list | grep -q $name && conda env remove -y -n $name
conda create -y -n $name  python=${PY_VERSION} --file "snakemake-env-packages.txt" 

    
source activate $name



log "copy index files for salmon"

#cp -r /reference/Indices/Salmon/  $HOME/genome
#cp -r /reference/gtf/             $HOME/genome
log "install snakemake-env-packages"
conda install ipython numpy scipy pandas matplotlib snakmake rpy2 gseapy

conda install hisat2 stringtie salmon star samtools bedtools fastqc rseqc multiqc graphviz

conda isntall bioconductor-deseq2 bioconductor-tximport bioconductor-readr bioconductor-ballgwon

log "install salmon requirements"
wget http://bioconductor.org/biocLite.R
Rscript -e "source('biocLite.R');options(BioC_mirror='http://mirrors.ustc.edu.cn/bioc/');biocLite();biocLite('EnsDb.Hsapiens.v86')"

log "all files are ready. Proceed to next step now."