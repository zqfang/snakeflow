# snakeflow
My NGS workflows based on snakemake

**workflows include**:
- RNA-seq (Salmon,hisat2)
- ChIP-seq (MACS2)
- ATAC-seq (MACS2)
- CITE-seq (Antibody captured, 10X genomics)
- Germline SNV calling (GATK, BCFtools)
- Germline Structural Variant calling 
    - short-read: Speedseq + svtools
    - long-read:
        - Pacbio: ngmlr/minimap2 + sniffiles


## Dependency

* General
  - samtools, deeptools, bedtools
  - fastqc, rseqc, multiqc, fastp
  - graphviz

* python 3
  - numpy
  - pandas
  - snakemake
  - matplotlib
  - seaborn
  - gseapy
  - macs2
  - rseqc

* R
  - DESeq2
  - tximport
  - readr
  - pheatmap
  - ggplot2
  - ggrepel
  - clusterProfiler
  - ChIPSeeker
  - EnsDb.Hsapiens.v86


* Variant calling
  - GATK (> 4.0)
  - BCFtools
  - Speedseq + svtools
  - minimap2
  - ngmlr
  - sniffiles

* RNA-seq
  - hisat2, salmon
  - rMATS-turbo, rmats2sashimiplot

* Single cell genomics
 - cellranger

## Installation

### Set up running environment. This config file will create a python 3.x env.

    bash snakeflow-enviroment-setup.sh
    
## usage
    
    # Step1: activate snakemake
    source activate snakeflow

    # Step2: clone this repo
    
    # Step3: copy all your fastq files into fastq dir
    find . -name "*fastq.gz" | while read id; do cp $id fastq/; done;
    
    # Step4: modify config.yml with your own paramter
    # Note: put config.yml in the same dir with your snakefile.
    vim  config.yml

    # Step5: run snakemake with -np option. this test your ``commands`` runs without any errors.
    snakemake -s salmon-tximport-deseq2-v0.2.snakefile -np

    # Step6: export workflow charts
    snakemake -s salmon-tximport-deseq2-v0.2.snakefile --dag | dot -Tpdf > dag.pdf

    # Step7: or using the default snakemake environment you've created above.
    snakemake -s salmon-tximport-deseq2-v0.1.snakefile -p -j 8
