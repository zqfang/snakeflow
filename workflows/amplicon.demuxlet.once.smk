
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
-> cutadapt -> CRISPResso2 -> NGmerge or fastp (for join pair end reads)

NOTE: run demultiplexing first, then run CRISPResso2
"""


import os, sys, csv 
import pandas as pd

########### INPUTS ########
# change working diretory here
workdir: "/data/bases/shared/Amplicon_seq/30-585200003"

# raw fastq path
FASTQ = "/data/bases/shared/Amplicon_seq/30-585200003/00_fastq"

# sequencing sub-libraries name, this pipline support multiple libraries !!!
LIBRARIES = ['C3-ALB-Clone','PNPLA3-C3-clone','PNPLA3-Cellpool'] 


## NOTE: all you need is to prepared a "{library}.barcodes.txt" file for each library like this 
## the barcode metadata file looks like this:
# sample-id       barcode-fwd     adaptor-fwd     barcode-rev     adaptor-rev     sgRNA   PCR-fwd PCR-rev ref_wt ref_hdr 
# S501-N701       TAGATCGC        GCACTCACTAGCCGTGACT     TCGCCTTA        CGACCGCACGTCTATTTAG     GCTGCAAGTCAAGCTGCCTT
# S501-N702       TAGATCGC        GCACTCACTAGCCGTGACT     CTAGTACG        CGACCGCACGTCTATTTAG     GCTGCAAGTCAAGCTGCCTT
# S501-N703       TAGATCGC        GCACTCACTAGCCGTGACT     TTCTGCCT        CGACCGCACGTCTATTTAG     GCTGCAAGTCAAGCTGCCTT
  

########### OUTPUTS ############################
CRISPRESSO = []
AMPLICON = []
SUB_LIBRARIES_SPLIT = {}
# parse each sub-libraries
for pool in LIBRARIES:
    ## read samples, barcodes
    SPLIT_SAMPLES = []
    BARCODE_1 = pool+".barcodes.txt"
    if not os.path.exists(BARCODE_1):
        print("File not exists: %s"%BARCODE_1)
        sys.exit()
    with open(BARCODE_1, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t',)
        for row in reader:
            s = row['sample-id']
            SPLIT_SAMPLES.append(s)
            #splits[s] = [row['adaptor-fwd'], row['adaptor-rev']]
    SUB_LIBRARIES_SPLIT[pool] = SPLIT_SAMPLES
    ## defined output
    for i, split in enumerate(SPLIT_SAMPLES):
        CRISPRESSO.append(f"{pool}_demux/CRISPResso/CRISPResso_on_{split}/CRISPResso_mapping_statistics.txt")
        AMPLICON.append(f"{pool}_demux/{split}.amplicon.fastq.gz")

############# WORKFLOWS ########################

rule all:
    input: AMPLICON, CRISPRESSO,

rule fastp:
    """
    qc and rename file
    """
    input:
        R1=os.path.join(FASTQ, "{library}_R1_001.fastq.gz"),
        R2=os.path.join(FASTQ, "{library}_R2_001.fastq.gz"),
        #barcodes = "{library}.barcodes.txt"
    output:
        R1="temp/{library}_R1.fastq.gz",
        R2="temp/{library}_R2.fastq.gz",
        j = "temp/{library}.fastp.json",
        h = "temp/{library}.fastp.html"
    shell:
        "fastp -i {input.R1} -I {input.R2} "
        "-o {output.R1} -O {output.R2} "
        "--json {output.j} --html {output.h} "


rule get_barcodes:
    input: "{library}.barcodes.txt",
    output:
        bf="temp/{library}.forward.barcodes.fasta",
        br="temp/{library}.reverse.barcodes.fasta",
    run:
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio import SeqIO
        bfwd = set()
        brev = set()
        r1_out = open(output.bf, 'w')
        r2_out = open(output.br, 'w')
        # os.makedirs(output.d, exist_ok=True)
        with open(input[0], newline='') as f:
            reader = csv.DictReader(f, delimiter='\t',)
            for row in reader:
                s = row['sample-id']
                s = s.split("-") ## TODO for paired adaptors ? 
                if row['barcode-fwd'] not in bfwd: 
                    sr1 = SeqRecord(Seq(row['barcode-fwd']), s[0], '', '')
                    SeqIO.write(sr1, r1_out, 'fasta')
                if row['barcode-rev'] not in brev:
                    sr2 = SeqRecord(Seq(row['barcode-rev']), s[1], '', '')
                    SeqIO.write(sr2, r2_out, 'fasta')
                bfwd.add(row['barcode-fwd'])
                brev.add(row['barcode-rev']) 
        r1_out.close()
        r2_out.close()


rule cutadapt_demulex_barcodes:
    """
    To support multi-library input, I have to iterate each library here in a for loop
    """
    input:
        b1=expand("temp/{library}.forward.barcodes.fasta",library=LIBRARIES),
        b2=expand("temp/{library}.reverse.barcodes.fasta",library=LIBRARIES),
        R1=expand("temp/{library}_R1.fastq.gz",library=LIBRARIES),
        R2=expand("temp/{library}_R2.fastq.gz", library=LIBRARIES),
    output:
        r1 = ["temp/%s_demux_1st/%s_R1.fastq.gz"%(lib, ss) for lib, split in SUB_LIBRARIES_SPLIT.items() for ss in split ],
        r2 = ["temp/%s_demux_1st/%s_R2.fastq.gz"%(lib, ss) for lib, split in SUB_LIBRARIES_SPLIT.items() for ss in split ],
    threads: 6
    params: 
        libraries = LIBRARIES
    run:
        os.makedirs("temp/logs", exist_ok=True)
        for lib in params.libraries:
            cmd = "cutadapt -j {threads} -e 0.15 --no-indels --discard-untrimmed "#--pair-adapters " ## paired adaptors input
            cmd += " -g ^file:temp/%s.forward.barcodes.fasta"%lib
            cmd += " -G ^file:temp/%s.reverse.barcodes.fasta"%lib
            cmd += " -o 'temp/%s_demux_1st/{{name1}}-{{name2}}_R1.fastq.gz' "%lib
            cmd += " -p 'temp/%s_demux_1st/{{name1}}-{{name2}}_R2.fastq.gz' "%lib
            cmd += " temp/%s_R1.fastq.gz temp/%s_R2.fastq.gz "%(lib, lib)
            cmd += " &> temp/logs/%s.cutadapt.demuxlet.log"%lib
            shell(cmd)


# re-evalue input and output since demulex_barcodes generate scatter files
rule cutadapt_primers_remove:
    input:
        barcode = "{library}.barcodes.txt",
        R1="temp/{library}_demux_1st/{split}_R1.fastq.gz",
        R2="temp/{library}_demux_1st/{split}_R2.fastq.gz",
    output:
        R1="temp/{library}_demux_1st/{split}_trim_R1.fastq.gz",
        R2="temp/{library}_demux_1st/{split}_trim_R2.fastq.gz",
    threads: 1
    log: "temp/logs/{library}.cutadapt.rm.primers.{split}.log"
    run:
        splits = {}
        with open(input.barcode, newline='') as f:
            reader = csv.DictReader(f, delimiter='\t',)
            for row in reader:
                s = row['sample-id']
                splits[s] = [row['adaptor-fwd'], row['adaptor-rev']]
        cmd2 = "cutadapt -j {threads} -e 0.1 --no-indels --pair-adapters -m 1 --discard-untrimmed " +\
            "-g ^%s -G ^%s "%(splits[wildcards.split][0], splits[wildcards.split][1] ) +\
            "-o {output.R1} -p {output.R2} " +\
            "{input.R1} {input.R2} &> {log} "
        shell(cmd2)

rule merge_paired_end:
    input:
        r1 = "temp/{library}_demux_1st/{b1st}_trim_R1.fastq.gz",
        r2 = "temp/{library}_demux_1st/{b1st}_trim_R2.fastq.gz",
    output: 
        amp="{library}_demux/{b1st}.amplicon.fastq.gz",
        j = temp("temp/{library}.{b1st}.fastp.json"),
        h = temp("temp/{library}.{b1st}.fastp.html")
    shell:
        "fastp --merge --merged_out {output.amp} "
        "-i {input.r1} -I {input.r2} "
        "--json {output.j} --html {output.h}"


rule CRISPResso_Analysis:
    """
    CRISPR Editing analysis
    """
    input: 
        #R1=lambda wildcards: "%s/{split}_R1.fastq.gz"%(checkpoints.demulex_barcodes.get(sample=wildcards.sample)),
        #R2=lambda wildcards: "%s/{split}_R2.fastq.gz"%(checkpoints.demulex_barcodes.get(sample=wildcards.sample)),
        #R2=lambda wildcards: "%s_Results/{split}_R2.fastq.gz"%(wildcards.sample),
        R1="temp/{library}_demux_1st/{split}_trim_R1.fastq.gz",
        R2="temp/{library}_demux_1st/{split}_trim_R2.fastq.gz",
        barcode = "{library}.barcodes.txt",
    output: 
        msa="{library}_demux/CRISPResso/CRISPResso_on_{split}/CRISPResso_mapping_statistics.txt",
        #outdir = directory("{sample}_demux/CRISPResso/{split}")
    threads: 1
    params: 
        outdir = "{library}_demux/CRISPResso",
    run:
        barcodes = pd.read_table(input.barcode, index_col=0)
        amplicon = barcodes.loc[wildcards.split, 'ref_wt']
        hdr_amplicon = barcodes.loc[wildcards.split, 'ref_hdr']
        # hdr_amplicon = barcodes.loc[wildcards.split, 'hdr']
        sgRNA = barcodes.loc[wildcards.split, 'sgRNA']
        cmd = "CRISPResso --amplicon_seq " + amplicon #+ ","+ amplicon2 # this should be the WT allel sequences.
        #cmd += " --amplicon_name Ref1,Ref2 " 
        cmd += " --amplicon_name Ref1 " 
        cmd += " --expected_hdr_amplicon_seq " + hdr_amplicon  # if you know the exepected sequence
        cmd += " --fastq_r1 {input.R1} --fastq_r2 {input.R2} "
        cmd += "--name {wildcards.split}  --output_folder {params.outdir} " #--file_prefix {wildcards.split} "
        cmd += " --guide_seq %s "%sgRNA   # exclude PAM only need 20 bp
        cmd += " --guide_name %s "%sgRNA
        cmd += "--max_paired_end_reads_overlap 250 " # 250 most for amplicon seq
        cmd += "--min_paired_end_reads_overlap 10 "
        cmd += " --quantification_window_size 10 " #--quantification_window_center -10 "
        cmd += " --write_cleaned_report "# --place_report_in_output_folder"
        shell(cmd)


# rule amplicon_demux_P2A:
#     """
#     second round of demuxleting
#     """
#     input: 
#         fq="{library}_demux/{b1st}.amplicon.fastq.gz", 
#         primers="{library}.barcodes-2nd.txt",
#     output: 
#         first = ["{library}_demux/{split}.%s.fastq"%(f) for f in FINAL_SAMPLES],
#         p2a = ["{library}_demux/{split}.%sP2A.fastq"%(f) for f in P2A_SAMPLES],
#     params:
#         P2A = P2A_SEQUENCE,
#         P2A_SAMPLES = P2A_SAMPLES
#     run:
#         # mamba install -c bioconda seqkit 
#         with open(input.primers, 'r') as p:
#             reader = csv.DictReader(p, delimiter='\t',)
#             for row in reader:
#                 s = row['sample-id']
#                 start = len(row['adaptor-fwd'])
#                 end =  len(row['adaptor-rev'])
#                 # start, end = len(record[1]), len(record[2])
#                 outfile = "{wildcards.sample}_demux/{wildcards.split}.%s.fastq"%(s)
#                 # command     
#                 cmd = "seqkit amplicon -F %s -R %s {input.fq}"%(row['adaptor-fwd'], row['adaptor-rev']) # -r {start+1}:-{end+1} > 
#                 if s in params.P2A_SAMPLES:
#                     op2a = "{wildcards.sample}_demux/{wildcards.split}.%sP2A.fastq"%(s)
#                     shell(cmd + " | seqkit grep -s -i -p {params.P2A} > %s"%op2a)
#                     shell(cmd + " | seqkit grep -s -i -v -p {params.P2A} > %s"%outfile)
#                 else:
#                     shell(cmd + " > " +  outfile)
