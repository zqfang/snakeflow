# snakeflow

### RNA-Seq workflow based on snakemake

## Dependency
* python 3.5
  - numpy
  - pandas
  - snakemake

* hisat2, stringtie, ballgwon
* salmon, tximport, deseq2
* samtools
* fastqc, rseqc, multiqc
* graphviz


## Installation

### Set up running environment. This config file will create a python 3.5 env.

    bash snakemake-enviroment-config.sh


### Set up environment mannually if the internet connection to AWS is blocked

#### Step1: Create a python 3.5 enviroment named ``snakeflow``. 

    conda create -n snakflow python=3.5

#### Step2: Install necessary packages using **pip**

    source activate snakeflow
    pip install snakemake
    pip install pandas

#### Step3: install other required bioinfo softwares above... 

    
## usage
    
    #Step1: activate snakemake
    source activate snakeflow

    #Step2: clone this repo
    
    #Step3: copy all your fastq files into fastq dir
    find . -name "*fastq.gz" | while read id; do cp $id fastq/; done;
    
    #Step4: modify config.yml with your own paramter
    #Note: put config.yml in the same dir with your snakefile.
    vim  config.yml

    #Step5: run snakemake with -np option. this test your workflow runs without errors.
    snakemake -s isoRNA-salmon-tximport-deseq2-v0.1.snakefile -np

    #Step6: export workflow charts
    snakemake -s isoRNA-salmon-tximport-deseq2-v0.1.snakefile --dag | dot -Tpdf > dag.pdf

    #Step7: run your workflow dependent on an isolated py2 conda env for rmats with 8 threads.
    snakemake -s isoRNA-hisat2-rmats-v0.1.snakefile --use-conda -p -j 8

    #Step7: or using the default snakemake environment you've created above.
    snakemake -s isoRNA-salmon-tximport-deseq2-v0.1.snakefile -p -j 8
