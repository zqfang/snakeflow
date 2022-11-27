
"""
This is a pipeline for CRISPR Editing Amplicon analysis

This pipeline runs a two rounds of demuxleting

NGS sequencing with 250 bp in each end of amplicon

R1: 
    barcode| ----- primer ---- | ----- amplicon ----
    >>>>>>>>===================AATCGCGTTCCGGAAAATCCTG......

R2: 

    barcode| ----- primer ---- | ----- amplicon ----
    >>>>>>>>===================CCTGAACCTCGATCGGATAGCT......

pipeline:
>> qc (fastp) >> demux barcode (split library into sample, cutadapt) 
>> trim primer (cutadpat)>> join R1-R2 (fastp) >> amplicon  >> CRISPR Resso 

"""
import os, sys, csv
import pandas as pd 

########### INPUTS ########
# change working diretory here
workdir: "/data/bases/shared/Amplicon_seq/30-554549436_CRIPaint_293"

# change raw fastq path
FASTQ = "raw"

# change library name. use the file prefix of raw the fastq file
#  e.g. "raw/CRISPaint_R1_001.fastq.gz" => "{LIBRARY}_R1_001.fastq.gz"
LIBRARY = "CRISPaint"

# if need to split P2A samples, please specify sample-ids from "{library}.barcodes-2nd.txt"
P2A_SAMPLES = [] # ["TUBB-HOT", "H4C3-HOT"]
P2A_SEQUENCE = "gctactaacttcagcctgctgaagcaggctggagacgtggaggagaaccctggacct"

# prepare the two barcodes files 
# first round barcodes to split libraries
BARCODE_1ST = "{library}.barcodes.txt".format(libary=LIBRARY) 
# second round barcode to split amplicons into final samples
BARCODE_2ND = "{library}.barcodes-2nd.txt".format(libary=LIBRARY) 

## NOTE: all you need is to prepared a "{library}.barcodes.txt" file for each library like this 
## the barcodes-2nd.txt file looks like this:
# sample-id       adaptor-fwd     adaptor-rev     sgRNA   ref_wt  ref_hdr
# TUBB-mNeon      agtccacaccgttgatggagcca CGGAGCAGTTCACTGCCATGTTC CGGTGAGGAGGCCGAAGAGG    CGGAGCAGTTCACTGCCATGTTCCGCCGGA
# TUBB-HOT        cttcagggtcagcttgccg     CGGAGCAGTTCACTGCCATGTTC CGGTGAGGAGGCCGAAGAGG    CGGAGCAGTTCACTGCCATGTTCCGCCGGA
# H4C3-mNeon      agtccacaccgttgatggagcca TCAAGCGCATTTCCGGTCTTATC GGGCGCACTCTGTATGGCTT    TCAAGCGCATTTCCGGTCTTATCTATGAGG
# H4C3-HOT        cttcagggtcagcttgccg     TCAAGCGCATTTCCGGTCTTATC GGGCGCACTCTGTATGGCTT    TCAAGCGCATTTCCGGTCTTATCTATGAGG


########### OUTPUTS ############################

## read library, barcodes
SPLIT_SAMPLES = []
if not os.path.exists(BARCODE_1ST):
    print("%s Not Found ! Is input path correct?"%BARCODE_1ST)
    sys.exit(0)

with open(BARCODE_1ST, newline='') as f:
    reader = csv.DictReader(f, delimiter='\t',)
    for row in reader:
        s = row['sample-id']
        SPLIT_SAMPLES.append(s)
        #splits[s] = [row['adaptor-fwd'], row['adaptor-rev']]

FINAL_SAMPLES = [] # ['TUBB-mNeon','TUBB-HOT','H4C3-mNeon','H4C3-HOT']
if not os.path.exists(BARCODE_2ND):
    print("%s Not Found ! Is input path correct?"%BARCODE_2ND)
    sys.exit(0)
with open(BARCODE_2ND, newline='') as f:
    reader = csv.DictReader(f, delimiter='\t',)
    for row in reader:
        s = row['sample-id']
        FINAL_SAMPLES.append(s)
## defined output
CRISPRESSO = []
#LIBRARY = BARCODE_1ST.replace(".barcodes.txt", "")
for i, split in enumerate(SPLIT_SAMPLES):
    for fx in FINAL_SAMPLES:
        CRISPRESSO.append(f"{LIBRARY}_demux/CRISPResso/CRISPResso_on_{split}.{fx}/CRISPResso_mapping_statistics.txt")
        if fx in P2A_SAMPLES:
            CRISPRESSO.append(f"{LIBRARY}_demux/CRISPResso/CRISPResso_on_{split}.{fx}P2A/CRISPResso_mapping_statistics.txt")

############# WORKFLOWS ########################
rule all:
    input: CRISPRESSO

rule fastp:
    """
    qc and rename file
    """
    input:
        R1=os.path.join(FASTQ, "{library}_R1_001.fastq.gz"),
        R2=os.path.join(FASTQ, "{library}_R2_001.fastq.gz"),
    output:
        R1="temp/{library}_R1.fastq.gz",
        R2="temp/{library}_R2.fastq.gz",
        j = "temp/{library}.fastp.json",
        h = "temp/{library}.fastp.html"
    shell:
        "fastp -i {input.R1} -I {input.R2} "
        "-o {output.R1} -O {output.R2} "
        "--json {output.j} --html {output.h} "


rule get_barcodes_1st:
    input: "{library}.barcodes.txt",
    output:
        bf="temp/{library}.forward.barcodes.fasta",
        br="temp/{library}.reverse.barcodes.fasta",
    run:
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio import SeqIO
        bfwd = []
        brev = []
        r1_out = open(output.bf, 'w')
        r2_out = open(output.br, 'w')
        # os.makedirs(output.d, exist_ok=True)
        with open(input[0], newline='') as f:
            reader = csv.DictReader(f, delimiter='\t',)
            for row in reader:
                bfwd.append(row['barcode-fwd'])
                brev.append(row['barcode-rev']) 
                s = row['sample-id']
                s = s.split("-") ## TODO
                sr1 = SeqRecord(Seq(bfwd[-1]), s[0], '', '')
                sr2 = SeqRecord(Seq(brev[-1]), s[1], '', '')
                SeqIO.write(sr1, r1_out, 'fasta')
                SeqIO.write(sr2, r2_out, 'fasta')
                # split sample by barcode
        r1_out.close()
        r2_out.close()


rule cutadapt_demulex_barcodes_1st:
    input:
        b1="temp/{library}.forward.barcodes.fasta",
        b2="temp/{library}.reverse.barcodes.fasta",
        R1="temp/{library}_R1.fastq.gz",
        R2="temp/{library}_R2.fastq.gz",
    output:
        #t = temp(directory("temp/{library}")),
        #d = directory("temp/{wildcards.library}_demux_1st"),
        r1 = ["temp/{library}_demux_1st/%s_R1.fastq.gz"%b1st for b1st in SPLIT_SAMPLES],
        r2 = ["temp/{library}_demux_1st/%s_R2.fastq.gz"%b1st for b1st in SPLIT_SAMPLES],
    threads: 6
    log: "temp/log/cutadapt.demuxlet_1st.{library}.log"
    run:
        #os.makedirs("temp/%s_demux_1st"%wildcards.library, exist_ok=True)
        shell("cutadapt -j {threads} -e 0 --no-indels "
            "--discard-untrimmed " # --pair-adapters
            "-g ^file:{input.b1} -G ^file:{input.b2} "
            "-o 'temp/{wildcards.library}_demux_1st/{{name1}}-{{name2}}_R1.fastq.gz' "
            "-p 'temp/{wildcards.library}_demux_1st/{{name1}}-{{name2}}_R2.fastq.gz' "
            "{input.R1} {input.R2} &> {log}")


rule cutadapt_primers_remove:
    input:
        barcode = "{library}.barcodes.txt",
        R1="temp/{library}_demux_1st/{split}_R1.fastq.gz",
        R2="temp/{library}_demux_1st/{split}_R2.fastq.gz",
    output:
        R1="temp/{library}_demux_1st/{split}_trim_R1.fastq.gz",
        R2="temp/{library}_demux_1st/{split}_trim_R2.fastq.gz",
    threads: 1
    log: "temp/log/cutadapt.demuxlet_1st.rm.primer.{library}.{split}.log"
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
            "{input.R1} {input.R2} &> {log}"
        shell(cmd2)


rule merge_paired_end:
    input:
        r1 = "temp/{library}_demux_1st/{b1st}_trim_R1.fastq.gz",
        r2 = "temp/{library}_demux_1st/{b1st}_trim_R2.fastq.gz",
    output: 
        amp = "temp/{library}_demux_1st/{b1st}.amplicon.fastq.gz",
        j = temp("temp/{library}_demux_1st/{b1st}.fastp.json"),
        h = temp("temp/{library}_demux_1st/{b1st}.fastp.html")
    shell:
        "fastp --merge --merged_out {output.amp} "
        "-i {input.r1} -I {input.r2} "
        "--json {output.j} --html {output.h} "

rule amplicon_demuxlet_2nd:
    input: 
        fq="temp/{library}_demux_1st/{split}.amplicon.fastq.gz", 
        primers="{library}.barcodes-2nd.txt",
    output: 
        second = ["{library}_demux/{split}.%s.fastq"%(f) for f in FINAL_SAMPLES],
    params:
        P2A = P2A_SEQUENCE,
        P2A_SAMPLES = P2A_SAMPLES
    run:
        # mamba install -c bioconda seqkit 
        with open(input.primers, 'r') as p:
            reader = csv.DictReader(p, delimiter='\t',)
            for row in reader:
                s = row['sample-id']
                start = len(row['adaptor-fwd'])
                end =  len(row['adaptor-rev'])
                outfile = "{wildcards.library}_demux/{wildcards.split}.%s.fastq"%(s)
                # command     
                cmd = "seqkit amplicon -F %s -R %s {input.fq}"%(row['adaptor-fwd'], row['adaptor-rev']) 
                #cmd += f"-r {start+1}:-{end+1} " # -r will cut primer.   
                if s in params.P2A_SAMPLES:
                    continue
                    # op2a = "{wildcards.library}_demux/{wildcards.split}.%sP2A.fastq"%(s)
                    # shell(cmd + " | seqkit grep -s -i -p {params.P2A} > %s"%op2a)
                    # shell(cmd + " | seqkit grep -s -i -v -p {params.P2A} > %s"%outfile)
                else:
                    shell(cmd + " > " +  outfile)


if P2A_SAMPLES:
    rule amplicon_demuxlet_2nd_p2a:
        input: 
            fq="temp/{library}_demux_1st/{split}.amplicon.fastq.gz", 
            primers="{library}.barcodes-2nd.txt",
        output: 
            #second = ["{library}_demux/{split}.%s.fastq"%(f) for f in FINAL_SAMPLES],
            p2a = ["{library}_demux/{split}.%sP2A.fastq"%(f) for f in P2A_SAMPLES],
        params:
            P2A = P2A_SEQUENCE,
            P2A_SAMPLES = P2A_SAMPLES
        run:
            # mamba install -c bioconda seqkit 
            with open(input.primers, 'r') as p:
                reader = csv.DictReader(p, delimiter='\t',)
                for row in reader:
                    s = row['sample-id']
                    start = len(row['adaptor-fwd'])
                    end =  len(row['adaptor-rev'])
                    outfile = "{wildcards.library}_demux/{wildcards.split}.%s.fastq"%(s)
                    # command     
                    cmd = "seqkit amplicon -F %s -R %s {input.fq}"%(row['adaptor-fwd'], row['adaptor-rev']) 
                    #cmd += f"-r {start+1}:-{end+1} " # -r will cut primer.   
                    if s in params.P2A_SAMPLES:
                        op2a = "{wildcards.library}_demux/{wildcards.split}.%sP2A.fastq"%(s)
                        shell(cmd + " | seqkit grep -s -i -p {params.P2A} > %s"%op2a)
                        shell(cmd + " | seqkit grep -s -i -v -p {params.P2A} > %s"%outfile)
                    # else:
                    # shell(cmd + " > " +  outfile)

# re-eval outputs from amplicon_demuxlet_2nd output, so not files are missing
checkpoint CRISPResso_alignment:
    input: 
        fq="{libary}_demux/{split}.{final}.fastq",
        barcode="{libary}.barcodes-2nd.txt",
    output: 
        msa="{libary}_demux/CRISPResso/CRISPResso_on_{split}.{final}/CRISPResso_mapping_statistics.txt",
    threads: 1
    params: 
        outdir = "{libary}_demux/CRISPResso",
        final = "{final}"
    message: "run this command need to activate CRISPResso env!!!"
    run:
        final = wildcards.final
        # since not prior information specified P2A samples in the barcode-2nd.txt
        # need to handle errors for missing keys. just use the existed keys
        if final.endswith("P2A"): final = final[:-3]
        barcodes = pd.read_table(input.barcode, index_col=0)
        amplicon = barcodes.loc[final, 'ref_wt']
        hdr_amplicon = barcodes.loc[final, 'ref_hdr']
        sgRNA = barcodes.loc[final, 'sgRNA']
        name = final.split("-")[0]
        cmd = "CRISPResso --amplicon_seq " + amplicon # this could be WT allel sequences.
        cmd += " --expected_hdr_amplicon_seq " + hdr_amplicon  # if you know the exepected sequence
        cmd += " --fastq_r1 {input.fq} --output_folder {params.outdir} "
        cmd += " --amplicon_name {wildcards.final} " 
        cmd += " --guide_seq %s "%sgRNA   # exclude PAM only need 20 bp
        cmd += " --guide_name sgRNA_%s "%name
        cmd += " --quantification_window_size 20 " #--quantification_window_center -10 "
        cmd += " --write_cleaned_report --exclude_bp_from_right 0  --exclude_bp_from_left 0 "# --place_report_in_output_folder"
        cmd += " --max_paired_end_reads_overlap 200 " 
        #cmd += " --file_prefix {wildcards.final} "
        shell(cmd)