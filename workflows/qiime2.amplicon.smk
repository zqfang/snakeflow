
"""
This is a pipeline for CRISPR Editing Amplicon analysis using QIIME2

NGS sequencing with 250 bp in each end of amplicon

R1: 
    barcode| ----- primer ---- | ----- amplicon ----
    >>>>>>>>===================AATCGCGTTCCGGAAAATCCTG......

R2: 

    barcode| ----- primer ---- | ----- amplicon ----
    >>>>>>>>===================CCTGAACCTCGATCGGATAGCT......

pipeline:
>> demux barcode (split sample) >> trim primer >> join R1-R2 >> amplicon >> qc >> CRISPR Resso 


two pipline ->  qiime2 +  CRISPRResso2 


see a short tutorial here for metagenome analysis
http://www.int-res.com/articles/suppl/m648p169_supp3.txt

#What: Computational pipeline for the V3/V4 region of the 16S rRNA gene 
(i.e., phylogenetic marker for bacteria) QIIME2 (v 2019.1), SILVA (v 132), and Amplicon Sequence Variants (ASVs)
"""
import os, sys 


workdir: "/data/bases/fangzq/Amplicon_seq"


POOL_SAMPLES = ['CRISPaint_293']
SPLIT_SAMPLES = ['S501-N701','S502-N702']
FINAL_SAMPLES = ['TUBB-mNeon','TUBB-HOT','H4C3-mNeon','H4C3-HOT']
P2A_SAMPLES = ['TUBB-HOT', 'H4C3-HOT']
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
CRISPRESSO = []
for i, split in enumerate(SPLIT_SAMPLES):
   for pool in POOL_SAMPLES:
       for fx in FINAL_SAMPLES:
           OUT2.append(f"{pool}_demux/{split}_{i}.{fx}.fastq")
           ALIGNMENT.append(f"{pool}_demux/{split}_{i}.{fx}.msa.txt")
           CRISPRESSO.append(f"{pool}_demux/CRISPResso/CRISPResso_on_{split}_{i}.{fx}/Alleles_frequency_table.zip")
           if fx in P2A_SAMPLES:
               OUT2.append(f"{pool}_demux/{split}_{i}.{fx}P2A.fastq") 
               ALIGNMENT.append(f"{pool}_demux/{split}_{i}.{fx}P2A.msa.txt")
               CRISPRESSO.append(f"{pool}_demux/CRISPResso/CRISPResso_on_{split}_{i}.{fx}P2A/Alleles_frequency_table.zip")


rule all:
    input: OUT, OUT2, CRISPRESSO# ALIGNMENT


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
        #"--p-mixed-orientation "
        "--p-error-rate 0.20 "


rule trim_index_primer:
    """
    If there are sequencing adapters or PCR primers in the reads which you’d like to remove, you can do that next as follows.
    """
    input: 
        demux="{sample}.demulexed.qza",
        barcode = "{sample}.barcodes.txt"
    output: "{sample}.trimmed.qza"
    threads: 8
    run:
        fwd = set()
        rev = set()
        with open(input.barcode, newline='') as f:
            reader = csv.DictReader(f, delimiter='\t',)
            for row in reader:
                fwd.add(row['adaptor-fwd'])
                rev.add(row['adaptor-rev'])    

        primer_f = fwd.pop()
        primer_r = rev.pop()    

        cmd =  "qiime cutadapt trim-paired " +\
                "--i-demultiplexed-sequences {input.demux} " +\
                f"--p-front-f {primer_f} " +\
                f"--p-front-r {primer_r} " +\
                "--p-error-rate 0.1 " +\
                "--o-trimmed-sequences {output} " +\
                "--verbose " +\
                "--p-cores {threads} "
        shell(cmd)

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
    output: 
        first = ["{sample}_demux/{split}_{i}.%s.fastq"%(f) for f in FINAL_SAMPLES],
        second = ["{sample}_demux/{split}_{i}.%sP2A.fastq"%(f) for f in P2A_SAMPLES],
    params:
        P2A = "gctactaacttcagcctgctgaagcaggctggagacgtggaggagaaccctggacct",
        P2A_SAMPLES = P2A_SAMPLES
    run:
        # mamba install -c bioconda seqkit 
        with open(input.primers, 'r') as p:
            for line in p:
                if line.startswith("#"): continue
                record = line.strip().split()
                start, end = len(record[1]), len(record[2])
                outfile = "{wildcards.sample}_demux/{wildcards.split}_{wildcards.i}.%s.fastq"%(record[0])
                # command     
                cmd = "seqkit amplicon -F %s -R %s {input.fq}"%(record[1], record[2]) # -r {start+1}:-{end+1} > 
                if record[0] in params.P2A_SAMPLES:
                    op2a = "{wildcards.sample}_demux/{wildcards.split}_{wildcards.i}.%sP2A.fastq"%(record[0])
                    shell(cmd + " | seqkit grep -s -i -p {params.P2A} > %s"%op2a)
                    shell(cmd + " | seqkit grep -s -i -v -p {params.P2A} > %s"%outfile)
                else:
                    shell(cmd + " > " +  outfile)


rule CRISPResso_alignment:
    input: 
        fq="{sample}_demux/{split}_{i}.{final}.fastq",
        ref = "template.ref.fa", 
    output: 
        msa="{sample}_demux/CRISPResso/CRISPResso_on_{split}_{i}.{final}/Alleles_frequency_table.zip",
    threads: 1
    params: 
        outdir = "{sample}_demux/CRISPResso",
        final = "{final}",
        tubb = "GAGGCCGAAGAGGAGGCCTA", # sgRNA
        h4c3 = "GGGCGCACTCTGTATGGCTT",
    run:
        from Bio import SeqIO
        amplicons = list(SeqIO.parse("template.ref.fa", "fasta"))
        ids = { t.id : i for i, t in enumerate(amplicons) }
        idx = ids[params.final]
        name = params.final.split("-")[0].lower()
        sgRNA = params.tubb if name == 'tubb' else params.h4c3
        amplicon = str(amplicons[idx].seq)       
        cmd = "CRISPResso --amplicon_seq " + amplicon # this could be WT allel sequences.
        # cmd = " --expected_hdr_amplicon_seq " + hdr_amplicon  # if you know the exepected sequence
        cmd += " --fastq_r1 {input.fq} --output_folder {params.outdir} "
        cmd += " --amplicon_name {wildcards.final} " 
        #cmd += " --guide_seq %s "%sgRNA   # exclude PAM
        #cmd += " --guide_name %s "%name
        cmd += " --quantification_window_size 20 " #--quantification_window_center -10 "
        cmd += " --write_cleaned_report "# --place_report_in_output_folder"
        shell(cmd)


# rule amplicon_demux:
#     input: 
#         fq="{sample}_demux/{split}_{i}_L001_R1_001.fastq.gz", 
#         primers="barcodes-2nd.txt"
#     output: 
#         first = ["{sample}_demux/{split}_{i}.%s.fasta"%(f) for f in FINAL_SAMPLES],
#         second = ["{sample}_demux/{split}_{i}.%sP2A.fasta"%(f) for f in P2A_SAMPLES],
#     params:
#         P2A = "gctactaacttcagcctgctgaagcaggctggagacgtggaggagaaccctggacct",
#         P2A_SAMPLES = P2A_SAMPLES
#     run:
#         # mamba install -c bioconda seqkit 
#         with open(input.primers, 'r') as p:
#             for line in p:
#                 if line.startswith("#"): continue
#                 record = line.strip().split()
#                 start, end = len(record[1]), len(record[2])
#                 outfile = "{wildcards.sample}_demux/{wildcards.split}_{wildcards.i}.%s.fasta"%(record[0])
#                 # command     
#                 cmd = "seqkit fq2fa {input.fq} |  " 
#                 cmd = cmd + "seqkit amplicon -F %s -R %s  "%(record[1], record[2]) # -r {start+1}:-{end+1} > 

#                 if record[0] in params.P2A_SAMPLES:
#                     op2a = "{wildcards.sample}_demux/{wildcards.split}_{wildcards.i}.%sP2A.fasta"%(record[0])
#                     shell(cmd + " | seqkit grep -s -i -p {params.P2A} > %s"%op2a)
#                     shell(cmd + " | seqkit grep -s -i -v -p {params.P2A} > %s"%outfile)
#                 else:
#                     shell(cmd + " > " +  outfile)
               


# rule clustalo:
#     input: 
#         fasta="{sample}_demux/{split}_{i}.{final}.fasta",
#         ref = "template.ref.fa"
#     output: 
#         msa="{sample}_demux/{split}_{i}.{final}.msa.fasta",
#         ref = temp("tmp.{sample}.{split}_{i}.{final}.ref.fa"),
#         ref_query =temp("tmp.{sample}.{split}_{i}.{final}.fa")
#     threads: 8
#     run:
#         shell("seqkit grep -p {wildcards.final} {input.ref} > {output.ref}") # extract fasta record by ID
#         shell("cat {output.ref} {input.fasta} > {output.ref_query}")
#         # mamba install -c bioconda clustalo
#         shell("clustalo -i {output.ref_query} -o {output.msa} --force --outfmt fasta --threads {threads}")

# rule fx2tabular:
#     input: "{sample}_demux/{split}_{i}.{final}.msa.fasta"
#     output: "{sample}_demux/{split}_{i}.{final}.msa.txt"
#     shell:
#         "seqkit fx2tab {input} > {output}"