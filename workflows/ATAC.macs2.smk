
"""
This pipline is modified for mouse from here:
https://github.com/crazyhottommy/pyflow-ATACseq
http://barcwiki.wi.mit.edu/wiki/SOPs/atac_Seq

"""

import os

# shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
# shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

workdir: "/data/bases/fangzq/20220711_CC/ATAC"

ALL_SAMPLES = ['Fadu_control', 'Fadu_NSD1-sh_KDM2A-sh', 'Fadu_NSD1-sh', 
               'PCI-13_control', 'PCI-13_NSD1-sh_KDM2A-sh', 'PCI-13_NSD1-sh']

BOWTIE2_INDEX = "/home/fangzq/genome/human/bowtie2Indices_GRCh38_noalt_as/GRCh38_noalt_as"
## download from https://sites.google.com/site/anshulkundaje/projects/blacklists
BLACKLIST = "/home/fangzq/genome/human/hg38-blacklist.v3.bed" 


TARGET_READS =  25000000 ## number of reads down sample to (25million default)
# genome size for bedtools slop
GENOME_SIZE =  "/home/fangzq/genome/human/GRCh38.p13.genome.size"
## genome fasta for nucleoATAC
GENOME_FASTA = "/home/fangzq/genome/human/GRCh38.p13.genome.fa"

MACS2_GENOME = "hs"
MACS2_PVAL = 1e-5
MACS2_PVAL_BROAD =  1e-5




############### OUTPUTS ##########
ALL_BIGWIG = expand("bigwig/{sample}.bw", sample = ALL_SAMPLES)
ALL_PEAKS = expand("peaks/{sample}.macs2_peaks.clean.broadPeak", sample = ALL_SAMPLES)
ALL_NUCLEO = expand("peaks/{sample}_nucleoATAC.occpeaks.bed.gz",  sample = ALL_SAMPLES)

ALL_QC = ["multiQC/multiQC_log.html"]
ALL_ATAQV = expand("BAM/{sample}.sorted.bam.ataqv.json", sample = ALL_SAMPLES)

rule all:
	input: ALL_PEAKS, ALL_BIGWIG, ALL_QC, ALL_ATAQV


rule fastqc:
    input:  "FASTQ/{sample}_R1.fastq.gz", "FASTQ/{sample}_R2.fastq.gz"
    output: "fqc/{sample}_R1_fastqc.zip", "fqc/{sample}_R2_fastqc.zip"
    log:    "log/{sample}_fastqc"
    threads: 1
    params : jobname = "{sample}"
    message: "fastqc {input}: {threads}"
    shell:
        """
        # fastqc works fine on .gz file as well
        fastqc -o fqc -f fastq --noextract {input[0]} {input[1]} 2> {log}
        """


rule trim_adapter:
    input: "FASTQ/{sample}_R1.fastq.gz", "FASTQ/{sample}_R2.fastq.gz"
    output: 
        "trim/{sample}_R1_val_1.fq.gz", 
        "trim/{sample}_R2_val_2.fq.gz",
        "trim/{sample}_R1_val_1_fastqc.zip", 
        "trim/{sample}_R2_val_2_fastqc.zip",
    log: "log/{sample}_trim_adaptor.log"
    threads: 6
    params: jobname = "{sample}"
    message: "trim_adaptor {input}: {threads}"
    shell:
        "trim_galore --cores {threads} --fastqc --gzip -nextera "
        "--paired --length 30 -o trim/ {input[0]} {input[1]} 2> {log}"

## the later step will remove chrM from the bam file and coordinate sort the bam
## so I did not cooridnate sort the bam at this step to save some time.
rule align:
    input: "trim/{sample}_R1_val_1.fq.gz", "trim/{sample}_R2_val_2.fq.gz"
    output: "BAM/{sample}.sorted.bam"
    threads: 12
    params:
        jobname = "{sample}",
        bt2 = BOWTIE2_INDEX, # bowtiwe2 index
    message: "aligning {input}: {threads} threads"
    log:
        bowtie2 = "log/{sample}.align",
        markdup = "log/{sample}.markdup",
    shell:
        ## samblaster mark duplicates for read id grouped reads. I do not coordinate sort the bam
        "bowtie2 --very-sensitive --no-discordant --threads 5  "
        "-X 2000 -x {params.bt2} -1 {input[0]} -2 {input[1]} 2> {log.bowtie2} "
        "| samblaster 2> {log.markdup} "
        "| samtools view -Sb - > {output[0]} "




# check number of reads mapped by samtools flagstat
rule flagstat_bam:
    input:  "BAM/{sample}.sorted.bam"
    output: "BAM/{sample}.sorted.bam.flagstat"
    log:    "BAM/{sample}.flagstat_bam"
    threads: 1
    params: jobname = "{sample}"
    message: "flagstat_bam {input}: {threads} threads"
    shell:
        """
        samtools flagstat {input} > {output} 2> {log}
        """


rule ataqv:
    input: "BAM/{sample}.sorted.bam"
    output: "BAM/{sample}.sorted.bam.ataqv.json"
    log: "log/{sample}_ataqv.log"
    threads: 1
    params:
        jobname = "{sample}",
        genome = 'human'
    message: "ataqv quality control for {input}"
    shell:
        """
        ataqv {params.genome} {input} --metrics-file {output} 2> {log}
        """


rule json_to_html:
	input: expand("BAM/{sample}.sorted.bam.ataqv.json", sample = ALL_SAMPLES)
	output: "ATAC_qc_html"
	log: "log/ATAC_qc_html.log"
	threads: 1
	message: "compiling json files to html ATAC-seq QC"
	shell:
		"""	
		mkarv 11ATAC_qc_html {input}
		"""


## shifting the reads are only critical for TF footprint, for peak calling and making bigwigs, it should be fine using the bams without shifting
# https://sites.google.com/site/atacseqpublic/atac-seq-analysis-methods/offsetmethods
rule remove_chrM_bam:
	input: "BAM/{sample}.sorted.bam"
	output: "BAM/{sample}_exclude_chrM.sorted.bam", "BAM/{sample}_exclude_chrM.sorted.bam.bai"
	log: "log/{sample}_exclude_chrM_bam.log"
	threads: 12
	message: "excluding chrM from bam {input} : {threads} threads"
	params: jobname = "{sample}"
	shell:
		# remove duplicates and reads on chrM, coordinate sort the bam
		# samblaster expects name sorted bamq
		"samtools view -h {input} | samblaster --removeDups "
		"| grep -v -P '\tchrM\t' "
		"| samtools view -Sb -F 4 - "
		"| samtools sort -m 10G -@ {threads} -T {input}.tmp -o {output[0]} "
		"samtools index {output[0]}" 

rule phantom_peak_qual:
    input: "BAM/{sample}_exclude_chrM.sorted.bam"
    output: "phantompeakqual/{sample}_phantom.txt"
    log: "log/{sample}_phantompeakqual.log"
    threads: 4
    params: jobname = "{sample}"
    message: "phantompeakqual for {input} : {threads} threads"
    shell:
        """
        Rscript  /home/fangzq/program/phantompeakqualtools/run_spp_nodups.R \
                -c={input} -savp -rf  -p={threads} -odir=phantompeakqual  \
                -out={output} -tmpdir=phantompeakqual 2> {log}
        """

rule down_sample:
    input: "BAM/{sample}_exclude_chrM.sorted.bam", "BAM/{sample}_exclude_chrM.sorted.bam.bai",
	 	   "BAM/{sample}_exclude_chrM.sorted.bam.flagstat"
    output: "BAM/{sample}-downsample.sorted.bam", "BAM/{sample}-downsample.sorted.bam.bai"
    log: "log/{sample}_downsample.log"
    threads: 12
    params: 
        jobname = "{sample}",
        target_reads = TARGET_READS,
    message: "downsampling for {input}"
    run:
        import re
        import subprocess
        with open (input[2], "r") as f:
            # fifth line contains the number of mapped reads
            line = f.readlines()[5]
            match_number = re.match(r'(\d.+) \+.+', line)
			## how many paired reads, roughly total #reads/2
            total_reads = float(match_number.group(1))/2

        target_reads = params.target_reads # 15million reads  by default, set up in the config.yaml file
        if total_reads > target_reads:
            down_rate = target_reads/total_reads
        else:
            down_rate = 1

        shell("sambamba view -f bam -t {threads} --subsampling-seed=3 -s %s {input[0]} | "%down_rate +\
              "samtools sort -m 2G -@ {threads} -T {output[0]}.tmp > {output[0]} 2> {log}")
        shell("samtools index {outbam}".format(outbam = output[0]))


rule make_bigwigs:
    input : "BAM/{sample}-downsample.sorted.bam", "BAM/{sample}-downsample.sorted.bam.bai"
    output: "bigwig/{sample}.bw"
    log: "log/{sample}.makebw"
    threads: 8
    params: jobname = "{sample}"
    message: "making bigwig for {input} : {threads} threads"
    shell:
        """
    	# no window smoothing is done, for paired-end, bamCoverage will extend the length to the fragement length of the paired reads
        bamCoverage -b {input[0]} --ignoreDuplicates \
                    --skipNonCoveredRegions --normalizeUsing RPKM \
                    -p {threads} --extendReads -o {output} 2> {log}
        """



# https://github.com/taoliu/MACS/issues/145
rule call_peaks_macs2:
    input: "BAM/{sample}-downsample.sorted.bam", "BAM/{sample}-downsample.sorted.bam.bai"
    output: bed = "peaks/{sample}_macs2_peaks.broadPeak"
    log: "log/{sample}_call_peaks_macs2.log"
    params:
        name = "{sample}_macs2",
        jobname = "{sample}",
        g = MACS2_GENOME,
        pval = MACS2_PVAL,
        pval_broad = MACS2_PVAL_BROAD,
    message: "call_peaks macs2 {input}: {threads} threads"
    shell:
        """
        ## for macs2, when nomodel is set, --extsize is default to 200bp, this is the same as 2 * shift-size in macs14.
        macs2 callpeak -t {input[0]} \
            --keep-dup all -f BAMPE -g {params.g} \
            --outdir peaks -n {params.name} -p {params.pval} \
            --broad --broad-cutoff {params.pval_broad} &> {log}
        """

rule multiQC:
    input :
        expand("log/{sample}.align", sample = ALL_SAMPLES),
        expand("BAM/{sample}.sorted.bam.flagstat", sample = ALL_SAMPLES),
        expand("trim/{sample}_{read}_fastqc.zip", sample = ALL_SAMPLES, read = ["R1_val_1", "R2_val_2"]),
        expand("fqc/{sample}_{read}_fastqc.zip", sample = ALL_SAMPLES, read = ["R1", "R2"])
    output: "multiQC/multiQC_log.html"
    log: "log/multiqc.log"
    message: "multiqc for all logs"
    shell:
        """
        multiqc fqc trim BAM log -o multiQC -d -f -v -n multiQC_log 2> {log}
        """

## extend the broad peak a bit for nucelosome analysis by nuceloATAC
rule make_bed_nucleoATAC:
    input: "peaks/{sample}_macs2_peaks.broadPeak"
    output: "peaks/{sample}_nucleo.bed"
    log: "log/{sample}_make_nucleoATAC_bed.log"
    threads: 1
    message: "making bed for nucleoATAC from {input}"
    params: 
        jobname= "{sample}",
        gs = GENOME_SIZE
    shell:
        """
        cat {input} | bedtools slop -b 200 -g {params.gs]} | sort -k1,1 -k2,2n | bedtools merge > {output} 2> {log}
        """

## nucleoATAC works on the non-shifted bam, and shift the reads internally!
# https://github.com/GreenleafLab/NucleoATAC/issues/58
rule nucleoATAC:
    input: "BAM/{sample}-downsample.sorted.bam", "BAM/{sample}-downsample.sorted.bam.bai", "peaks/{sample}_nucleo.bed"
    output: "peaks/{sample}_nucleoATAC.occpeaks.bed.gz"
    log: "log/{sample}_nucleoATAC.log"
    threads: 5
    message: "calling nucleosome by nucleoATAC for {input} : {threads} threads"
    params:
        jobname = "{sample}",
        outputdir = os.path.dirname(srcdir("log")),
        fa = GENOME_FASTA
    shell:
        """
        cd nucleoATAC
        nucleoatac run --bed {params.outputdir}/{input[2]} --bam {params.outputdir}/{input[0]} --cores {threads} \
                        --fasta {params.fa} --out {wildcards.sample} 2> {params.outputdir}/{log}
        """


rule subtract_blacklist:
    input:
        bed = "peaks/{sample}_macs2_peaks.broadPeak",
        blacklist= BLACKLIST,
    output:
        "peaks/{sample}.macs2_peaks.clean.broadPeak"
    shell:
        "sort -k1,1 -k2,2n {input.bed} | bedtools subtract -a stdin -b {input.blacklist} -A > {output}"