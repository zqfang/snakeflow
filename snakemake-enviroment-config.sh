#!bin/bash

#this script used for initialize working directory and 
#install requirements files for runing snakemake workflows.
#all tools'indpendencies will be automatically installed by conda.
set -e

log () {
    echo
    echo "[`date`] Step: $1 "
    echo
}

PY_VERSION=2.7

if [[ "$PY_VERSION" == "2.7" ]]; then
        wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh;
else
        wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
fi

bash miniconda.sh -b -p /miniconda
export PATH="/miniconda/bin:$PATH"

#setup bioconda channel.
conda config --set always_yes yes 
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda


log "install snakemake requirements"

#python2 packages
conda install rseqc

source deactivate
name="snakeflow"
#clone a environment
#conda list --export > snakemake-env-packages.txt
#conda create -n snakeflow-clone --file snakemake-env-packages.txt

conda env list | grep -q $name && conda env remove -y -n $name
conda create -y -n $name  python=${PY_VERSION}
               # --file "snakemake-env-packages.txt" 

    
source activate $name



#log "copy index files for salmon"
#cp -r /reference/Indices/Salmon/  $HOME/genome
#cp -r /reference/gtf/             $HOME/genome

log "install snakemake-env-packages"
#python packages
conda install ipython cython numpy scipy pandas matplotlib snakemake gseapy xlrd xlwt multiqc 

#other commandline tools
conda install hisat2 stringtie salmon star samtools bedtools fastqc graphviz
#R packages
conda isntall bioconductor-deseq2 bioconductor-tximport bioconductor-readr bioconductor-ballgwon

#log "install salmon requirements"
#wget http://bioconductor.org/biocLite.R
#Rscript -e "source('biocLite.R');options(BioC_mirror='http://mirrors.ustc.edu.cn/bioc/');biocLite();biocLite('EnsDb.Hsapiens.v86')"

log "all files are ready. Proceed to next step now."
