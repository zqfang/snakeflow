
""" see a short tutorial here for metagenome analysis
http://www.int-res.com/articles/suppl/m648p169_supp3.txt

#What: Computational pipeline for the V3/V4 region of the 16S rRNA gene 
(i.e., phylogenetic marker for bacteria) QIIME2 (v 2019.1), SILVA (v 132), and Amplicon Sequence Variants (ASVs)
"""
import os, sys 


workdir: "/data/bases/fangzq/Amplicon_seq"


POOL_SAMPLES = ['CRISPaint_293']
SPLIT_SAMPLES = ['S501-N701','S502-N702']
FINAL_SAMPLES = ['TUBB-mNeon','TUBB-HOT','H4C3-mNeon','H4C3-HOT']
# OUT = expand("{sample}.deblur.rep.align.mask.qza", sample=SAMPLES)
OUT = expand("{sample}.trimmed.join.filtered.qzv", sample=POOL_SAMPLES) 
FASTQ = ["%s_demux/%s_%s_L001_R1_001.fasta"%(pool, split, i) for i, split in enumerate(SPLIT_SAMPLES) for pool in POOL_SAMPLES]
# ALIGNMENT = expand("{pool}_demux/{split}_{i}.msa.fasta", split=SPLIT_SAMPLES, pool=POOL_SAMPLES, i=range(len(SPLIT_SAMPLES)))

record = []
with open("barcodes-2nd.txt", 'r') as p:
    for line in p:
        if line.startswith("#"): continue
        record.append(line.strip().split())

OUT2= []
ALIGNMENT = []
for i, split in enumerate(SPLIT_SAMPLES):
   for pool in POOL_SAMPLES:
       for fx in FINAL_SAMPLES:
           OUT2.append(f"{pool}_demux/{split}_{i}.{fx}.fasta")
           ALIGNMENT.append(f"{pool}_demux/{split}_{i}.{fx}.msa.txt")


rule all:
    input: OUT, OUT2, ALIGNMENT


# rule rename_r1:
#     input: "{sample}_R1.fastq.gz"
#     output: "{sample}/forward.fastq.gz"
#     shell:
#         "mv {input} {output}"

# rule rename_r2:
#     input: "{sample}_R2.fastq.gz"
#     output: "{sample}/reverse.fastq.gz"
#     shell:
#         "mv {input} {output}"



rule data_import:
    """
    The import format for paired-end reads with the barcodes still in the sequence 
    is MultiplexedPairedEndBarcodeInSequence 

    this format expects two files in a directory (forward.fastq.gz and reverse.fastq.gz).
    """
    input: 
        R1="{sample}/forward.fastq.gz", 
        R2="{sample}/reverse.fastq.gz",
    output: "{sample}.multiplexed.qza"
    shell:
        "qiime tools import "
        "--type MultiplexedPairedEndBarcodeInSequence "
        "--input-path {wildcards.sample} "
        "--output-path {output} "


rule demulex_cutadapt:
    """
    cut barcode and demuxled
    """
    input: 
        qza = "{sample}.multiplexed.qza",
        barcode = "barcodes.txt"
    output: 
        demulexed="{sample}.demulexed.qza",
        unmatched="{sample}.unmatched.qza"
    shell:
        "qiime cutadapt demux-paired "
        "--i-seqs {input.qza}  "
        "--m-forward-barcodes-file {input.barcode} "
        "--m-reverse-barcodes-column barcode-rev "
        "--m-forward-barcodes-column barcode-fwd "
        "--m-reverse-barcodes-file {input.barcode}  "
        "--o-per-sample-sequences {output.demulexed}  "
        "--o-untrimmed-sequences {output.unmatched} "

rule trim_primer:
    """
    If there are sequencing adapters or PCR primers in the reads which you’d like to remove, you can do that next as follows.
    """
    input: "{sample}.demulexed.qza"
    output: "{sample}.trimmed.qza"
    params:
        primer_f= "GCACTCACTAGCCGTGACTATCGACTAATAGAT",#"GCACTCACTAGCCGTGACT",
        primer_r= "CGACCGCACGTCTATTTAGTGGAGCCGAGATGT",#CGACCGCACGTCTATTTAGT",
    threads: 8
    shell:
        "qiime cutadapt trim-paired "
        "--i-demultiplexed-sequences {input} "
        "--p-front-f {params.primer_f} "
        "--p-front-r {params.primer_r} "
        "--p-error-rate 0.1 "
        "--o-trimmed-sequences {output} "
        "--verbose "
        "--p-cores {threads} "

rule join_paired:
    input: "{sample}.trimmed.qza"
    output: "{sample}.demulexed.join.qza"
    shell:
        "qiime vsearch join-pairs "
        "--i-demultiplexed-seqs {input} "
        "--o-joined-sequences {output} "


rule qc_filter:
    input: "{sample}.demulexed.join.qza"
    output: 
        filt="{sample}.trimmed.join.filtered.qza",
        stat="{sample}.trimmed.join.filtered.stat.qza",
    shell:
        "qiime quality-filter q-score "
        "--i-demux {input} "
        "--p-min-quality 20 "
        "--o-filtered-sequences {output.filt} "
        "--o-filter-stats {output.stat} "

rule summarise:
    input: "{sample}.trimmed.join.filtered.qza"
    output: "{sample}.trimmed.join.filtered.qzv"
    shell:
        "qiime demux summarize "
        "--i-data {input} "
        "--o-visualization {output} "


rule export:
    input: "{sample}.trimmed.join.filtered.qza"
    output: ["{sample}_demux/%s_%s_L001_R1_001.fastq.gz"%(split, i) for i, split in enumerate(SPLIT_SAMPLES)]
    shell:
        "qiime tools export --input-path {input} "
        "--output-path {wildcards.sample}_demux"


rule amplicon_demux:
    input: 
        fq="{sample}_demux/{split}_{i}_L001_R1_001.fastq.gz", 
        primers="barcodes-2nd.txt"
    output: ["{sample}_demux/{split}_{i}.%s.fasta"%(f) for f in FINAL_SAMPLES],
    run:
        # mamba install -c bioconda seqkit 
        with open(input.primers, 'r') as p:
            for line in p:
                if line.startswith("#"): continue
                record = line.strip().split()
                start, end = len(record[1]), len(record[2])
                cmd = f"seqkit amplicon -F {record[1]} -R {record[2]} -r {start+1}:-{end+1} > "        
                cmd = "seqkit fq2fa {input.fq} |  " + cmd 
                cmd += "{wildcards.sample}_demux/{wildcards.split}_{wildcards.i}.%s.fasta"%(record[0])
                shell(cmd)


rule clustalo:
    input: "{sample}_demux/{split}_{i}.{final}.fasta"
    output: "{sample}_demux/{split}_{i}.{final}.msa.fasta"
    threads: 8
    shell:
        # mamba install -c bioconda clustalo
        "clustalo -i {input} -o {output} --force --outfmt fasta --threads {threads}"

rule fx2tabular:
    input: "{sample}_demux/{split}_{i}.{final}.msa.fasta"
    output: "{sample}_demux/{split}_{i}.{final}.msa.txt"
    shell:
        "seqkit fx2tab {input} > {output}"


# # ### DENOISE PROCESSED READS ###        
# rule deblur:
#     input: 
#         filt="{sample}.trimmed.join.filtered.qza",
#         ref = ""
#     output:
#         rep = "{sample}.deblur.rep.qza",
#         tab = "{sample}.deblur.table.qza",
#         stat = "{sample}.debulr.stat.qza"
#     params:
#         length = 0 # -1, disable trim
#     threads: 8
#     shell:
#         "qiime deblur denoise-other " # amplicon data not 16S
#         "--p-jobs-to-start {threads} "
#         "--i-demultiplexed-seqs {input.filt} "
#         "--i-reference-seqs {input.ref} "
#         "--p-trim-length {params.length} "
#         "--o-representative-sequences {output.rep} "
#         "--o-table {output.tab} "
#         "--p-sample-stats "
#         "--o-stats {output.stat} "


# # ### STEP 7: CREATE PHYLOGENY ###
# # #ALIGNMENT OF REPRESENTATIVE SEQUENCES
# rule mafft:
#     input: "{sample}.trimmed.join.filtered.qza" #"{sample}.deblur.rep.qza"
#     output: "{sample}.deblur.rep.align.qza"
#     threads: 8
#     shell:
#         "qiime alignment mafft --p-n-threads {threads} "
#         "--i-sequences {input} "
#         "--o-alignment {output} "

# # #MASK HIGHLY VARIABLE NOISY POSITIONS IN ALIGNMENT
# rule mask:
#     input: "{sample}.deblur.rep.align.qza"
#     output: "{sample}.deblur.rep.align.mask.qza"
#     shell:
#         "qiime alignment mask "
#         "--i-alignment {input} "
#         "--o-masked-alignment {output}"



