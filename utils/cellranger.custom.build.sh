## see tutorial here: 
## https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_mr

##### summaries the code for handy use

## make copies 
mkdir custom-GRCh38-2020-A
cd custom-GRCh38-2020-A
cp ../refdata-gex-GRCh38-2020-A/genes/genes.gtf customref-GRCh38-2020-A.gtf
cp ../refdata-gex-GRCh38-2020-A/fasta/genome.fa customref-GRCh38-2020-A.fa

## count base numbers
cat GFP.fa | grep -v "^>" | tr -d "\n" | wc -c

## The results of above command shows there are 922 bases.
## add record  
echo -e 'GFP\tunknown\texon\t1\t922\t.\t+\t.\tgene_id "GFP"; transcript_id "GFP"; gene_name "GFP"; gene_biotype "protein_coding";' > GFP.gtf

## append to genome fasta
cat GFP.fa >> customref-GRCh38-2020-A.fa

## append to gtf
cat GFP.gtf >> customref-GRCh38-2020-A.gtf

## mkref
cellranger mkref --genome=csutomref_build \
  --fasta=customref-GRCh38-2020-A.fa \
  --genes=customref-GRCh38-2020-A.gtf