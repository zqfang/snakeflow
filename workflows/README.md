# Snakemake pipeline for SNV calling

## Usage
### 1. Edit the `config.yaml` file for required files:
Install `snakemake` first

### 2 Run on local computer
```shell
# modify the file path in haplomap and run with 12 cores
snakemake -s bcftools.call.smk  --configfile config.yaml \
          -k -p -j 12   
```

## Two pipeline developed:
1. bcftools call 
  - prefer pipeline for inbred mouse and HBCGM input
  - more accuracy for inbred mouse ?


2. GATK best practice
  - designed for human genetics 
  - have to play with ``VQSR`` or ``hardfiltering`` parameters if use non-human data


**Caution !**: Both pipelines take a long time to run.


## Ensemble VEP

### Install VEP

[Download and Install](https://uswest.ensembl.org/info/docs/tools/vep/script/vep_download.html)

or use conda to install `vep`

### Install Genome
Install data for offline mode
```shell
INSTALL.pl -a cfp -s mus_musculus -y GRCm38 --CONVERT --PLUGINS CADD,GO,TSSDistance,LoF,SpliceAI
```

or Huam
```shell
INSTALL.pl -a cfp -s homo_sapiens -y GRCh38 --CONVERT
```

### Notes


## Why not GATK ?

One of my colleague who studies mouse genetics, said, 

> I tried the haplotype caller from GATK. But it seems that the haplotype caller is designed for heterogeneous genome like human than for mice. Therefore, the result coming out of HC is worse than samtools, as I manually inspected a few regions that HC calls didn't make sense.

> In addition, in one of their mouse genomic paper that we reviewed, they even skipped the second recalibration step. We asked them why and they said it was because of the same reason: good for human but not that good for the homogeneous inbred mouse.

From my point of view, I think GATK works OK.