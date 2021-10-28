
"""
This is a pipeline for CRISPR Editing Amplicon analysis

NGS sequencing with 250 bp in each end of amplicon

R1: 
    barcode| ----- primer ---- | ----- amplicon ----
    >>>>>>>>===================AATCGCGTTCCGGAAAATCCTG......

R2: 

    barcode| ----- primer ---- | ----- amplicon ----
    >>>>>>>>===================CCTGAACCTCGATCGGATAGCT......

pipeline:
>> demux barcode (split sample) >> trim primer >> join R1-R2 >> amplicon >> qc >> CRISPR Resso 
-> cutadapt -> CRISPResso2 -> NGmerge or fastp (for join pair emd reads)

NOTE: run demultiplexing first, then run CRISPResso2
"""


import os, sys, csv 
import pandas as pd

workdir: "/data/bases/fangzq/Amplicon_seq/30-585200003"


POOL_SAMPLES = ['C3-ALB-Clone', 'PNPLA3-C3-clone','PNPLA3-Cellpool']
SPLIT_SAMPLES = {}
## NOTE: 
# barcodes metadata file has following columns
# ['sample-id', 'barcode-fwd', 'adaptor-fwd', 'barcode-rev', 'adaptor-rev', 'sgRNA','PCR-fwd','PCR-rev','reference']
for s in POOL_SAMPLES:
    _SAMPLES = []
    with open(f"{s}.barcodes.txt", newline='') as f:
        reader = csv.DictReader(f, delimiter='\t',)
        for row in reader:
            _SAMPLES.append(row['sample-id'])
            #SPLIT_TEMPLATES[row['sample-id']] = row['template']
    SPLIT_SAMPLES [s] = _SAMPLES




#### OUTPUT FILES ##########
OUTDIR = directory(expand("{sample}_Results", sample=POOL_SAMPLES))
CRISPResso_OUTPUT = ["%s_Results/CRISPResso_on_%s/%s.Alleles_frequency_table.zip"%(sample, split, split) for sample in POOL_SAMPLES for split in SPLIT_SAMPLES[sample]]
AMPLICONS = ["%s_Results/%s_amplicon.fastq.gz"%(sample, split) for sample in POOL_SAMPLES for split in SPLIT_SAMPLES[sample]]


rule all:
    input: OUTDIR, CRISPResso_OUTPUT, AMPLICONS# ALIGNMENT


rule fastp:
    input: 
        R1="00_fastq/{sample}_R1_001.fastq.gz",
        R2="00_fastq/{sample}_R2_001.fastq.gz",
    output:
        R1="{sample}/forward.fastq.gz",
        R2="{sample}/reverse.fastq.gz",
        j = "temp/{sample}.fastp.json",
        h = "temp/{sample}.fastp.html"
    shell:
        "fastp -i {input.R1} -I {input.R2} "
        "-o {output.R1} -O {output.R2} "
        "--json {output.j} --html {output.h} "



checkpoint demulex_barcodes:
    input:
        barcode = "{sample}.barcodes.txt",
        R1="{sample}/forward.fastq.gz",
        R2="{sample}/reverse.fastq.gz",
    output:
        t = temp(directory("temp/{sample}")),
        d = directory("{sample}_Results")
    run:
        bfwd = []
        brev = []
        os.makedirs(output.t, exist_ok=True)
        os.makedirs(output.d, exist_ok=True)
        with open(input.barcode, newline='') as f:
            reader = csv.DictReader(f, delimiter='\t',)
            for row in reader:
                bfwd.append(row['barcode-fwd'])
                brev.append(row['barcode-rev']) 
                s = row['sample-id']
                # split sample by barcode
                cmd1 = "cutadapt -e 0.15 --no-indels --pair-adapters -m 1 --trimmed-only " +\
                    "-g ^%s -G ^%s "%(bfwd[-1], brev[-1]) +\
                    "-o temp/%s/%s_R1.fastq.gz -p temp/%s/%s_R2.fastq.gz "%(wildcards.sample, s, wildcards.sample, s) +\
                    "%s/forward.fastq.gz %s/reverse.fastq.gz"%(wildcards.sample, wildcards.sample)

                shell(cmd1)
                # remove adaptor per sample
                cmd2 = "cutadapt -e 0.1 --no-indels --pair-adapters -m 1 --trimmed-only " +\
                    "-g ^%s -G ^%s "%(row['adaptor-fwd'], row['adaptor-rev']) +\
                    "-o %s_Results/%s_R1.fastq.gz -p %s_Results/%s_R2.fastq.gz "%(wildcards.sample, s, wildcards.sample, s) +\
                    "temp/%s/%s_R1.fastq.gz temp/%s/%s_R2.fastq.gz "%(wildcards.sample, s, wildcards.sample, s)
                shell(cmd2)



rule CRISPResso_Analysis:
    """
    CRISPR Editing analysis
    """
    input: 
        #R1=lambda wildcards: "%s/{split}_R1.fastq.gz"%(checkpoints.demulex_barcodes.get(sample=wildcards.sample)),
        #R2=lambda wildcards: "%s/{split}_R2.fastq.gz"%(checkpoints.demulex_barcodes.get(sample=wildcards.sample)),
        #R2=lambda wildcards: "%s_Results/{split}_R2.fastq.gz"%(wildcards.sample),
        R1="{sample}_Results/{split}_R1.fastq.gz",
        R2="{sample}_Results/{split}_R2.fastq.gz",
        barcode = "{sample}.barcodes.txt",
    output: 
        msa="{sample}_Results/CRISPResso_on_{split}/{split}.Alleles_frequency_table.zip",
        #outdir = directory("{sample}_demux/CRISPResso/{split}")
    threads: 1
    params: 
        outdir = "{sample}_Results/",
    run:
        barcodes = pd.read_table(input.barcode, index_col=0)
        amplicon = barcodes.loc[wildcards.split, 'reference']
        sgRNA = barcodes.loc[wildcards.split, 'sgRNA']
        cmd = "CRISPResso --amplicon_seq " + amplicon # this could be WT allel sequences.
        # cmd = " --expected_hdr_amplicon_seq " + hdr_amplicon  # if you know the exepected sequence
        cmd += " --fastq_r1 {input.R1} --fastq_r2 {input.R2} "
        cmd += "--name {wildcards.split}  --output_folder {params.outdir} --file_prefix {wildcards.split} "
        cmd += " --amplicon_name {wildcards.sample} " 
        #cmd += " --guide_seq %s "%sgRNA   # exclude PAM
        #cmd += " --guide_name %s "%name
        cmd += "--max_paired_end_reads_overlap 250 " # 250 most for amplicon seq
        cmd += "--min_paired_end_reads_overlap 10 "
        cmd += " --quantification_window_size 10 " #--quantification_window_center -10 "
        cmd += " --write_cleaned_report "# --place_report_in_output_folder"
        shell(cmd)




rule merge_paired_end:
    """
    merge paired end reads into one amplicons.
    """
    input:
        #R1="{sample}_demux/{split}_R1.fastq.gz",
        #R2="{sample}_demux/{split}_R2.fastq.gz",
        R1=lambda wildcards: "%s/{split}_R1.fastq.gz"%(checkpoints.demulex_barcodes.get(**wildcards).output.d),
        R2=lambda wildcards: "%s/{split}_R2.fastq.gz"%(checkpoints.demulex_barcodes.get(**wildcards).output.d),
    output:
        "{sample}_demux/{split}_amplicon.fastq.gz"
    shell:
        "NGmerge -1 {input.R1} -2 {input.R2} -m 10 -o {output}"
