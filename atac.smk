
"""
This pipline is modified for mouse from here:
https://github.com/crazyhottommy/pyflow-ATACseq

"""

import os

# shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
# shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

workdir: "/data/bases/fangzq/MouseEpigenomeAtlas"

with open("atac.accession.txt") as atac:
    ALL_SAMPLES = atac.read().strip().split() 

#this is the name sorted bam, not coordinate sorted bam, samtools index only for coordiante sorted bam.
BAM_DIR = "/data/bases/fangzq/MouseEpigenomeAtlas/bam"
BLACKLIST = "/home/fangzq/genome/mouse/mm10-blacklist.v2.bed" 
## downsample from the chrM excluded bam files
ALL_DOWNSAMPLE_BAM = expand("bam/{sample}.downsample.sorted.bam", sample = ALL_SAMPLES)
ALL_DOWNSAMPLE_INDEX = expand("bam/{sample}.downsample.sorted.bam.bai", sample = ALL_SAMPLES)

ALL_BIGWIG = expand("bigwig/{sample}.bw", sample = ALL_SAMPLES)
ALL_PEAKS = expand("beds/{sample}.macs2_peaks.blacklist_removed.broadPeak", sample = ALL_SAMPLES)

# ALL_QC = ["10multiQC/multiQC_log.html"]
# ALL_ATAQV = expand("04aln/{sample}.sorted.bam.ataqv.json", sample = ALL_SAMPLES)

rule all:
	input: ALL_PEAKS

rule namesort_bam:
    input: "bam/{sample}.bam"
    output: bam = temp("bam/{sample}.sorted.bam")
    threads: 6
    shell:
        "samtools sort -n -m 2G -@ {threads} -T {input}.tmp -o {output.bam} {input}"

# rule bam_index:
#     input: "bam/{sample}.sorted.bam"
#     output: temp("bam/{sample}.sorted.bam.bai")
#     shell: 
#         "samtools index {input}"

## shifting the reads are only critical for TF footprint, for peak calling and making bigwigs, it should be fine using the bams without shifting
# https://sites.google.com/site/atacseqpublic/atac-seq-analysis-methods/offsetmethods
rule remove_chrM_bam:
    input: 
        bam="bam/{sample}.sorted.bam",
    output:
        bam="bam/{sample}.exclude_chrM.sorted.bam",
    log: "logs/{sample}.exclude_chrM_bam.log"
    threads: 6
    message: "excluding chrM from bam {input.bam} : {threads} threads"
    params: jobname = "{sample}"
	shell:
		# remove duplicates and reads on chrM, coordinate sort the bam
		# samblaster expects name sorted bamq
		"samtools view -h {input.bam} | samblaster --removeDups "
		"| grep -v -P '\tchrM\t' | samtools view -Sb -F 4 - " 
		"| samtools sort -m 2G -@ {threads} -T {output.bam}.tmp -o {output.bam}"

rule bam_index:
    input: 
        "bam/{sample}.exclude_chrM.sorted.bam"
    output: 
        "bam/{sample}.exclude_chrM.sorted.bam.bai"
    shell: 
        "samtools index {input}"

## consider how to reuse the rules.
rule flagstat_bam:
    input: 
        "bam/{sample}.exclude_chrM.sorted.bam"
    output: 
        "bam/{sample}.exclude_chrM.sorted.bam.flagstat"
    log:    "logs/{sample}.exclude_chrM_flagstat_bam"
    threads: 1
    params: jobname = "{sample}"
    message: "flagstat_bam {input}: {threads} threads"
    shell:
        "samtools flagstat {input} > {output} 2> {log}"

rule down_sample:
    input: 
        "bam/{sample}.exclude_chrM.sorted.bam", 
        "bam/{sample}.exclude_chrM.sorted.bam.bai",
        rules.flagstat_bam.output
	 	#"bam/{sample}.exclude_chrM.sorted.bam.flagstat"
    output: 
        "bam/{sample}.downsample.sorted.bam", 
        "bam/{sample}.downsample.sorted.bam.bai"
    log: "logs/{sample}.downsample.log"
    threads: 5
    params: jobname = "{sample}"
    message: "downsampling for {input}"
    run:
        import re
        with open (input[2], "r") as f:
            # fifth line contains the number of mapped reads
            line = f.readlines()[5]
            match_number = re.match(r'(\d.+) \+.+', line)
			## how many paired reads, roughly total #reads/2
            total_reads = float(match_number.group(1))/2

        # target_reads = config["target_reads"] # 25million reads  by default, set up in the config.yaml file
        target_reads = 25000000
        if total_reads > target_reads:
            down_rate = target_reads/total_reads
        else:
            down_rate = 1
        
        cmd = "sambamba view -f bam -t {threads} --subsampling-seed=3 -s {rate} {inbam} | "
        cmd += "samtools sort -m 2G -@ {threads} -T {outbam}.tmp > {outbam} 2> {log}"
        cmd = cmd.format(rate = down_rate, inbam = input[0], outbam = output[0], log = log, threads=threads)
        shell(cmd)
        shell("samtools index {output[0]}")

# https://github.com/taoliu/MACS/issues/145
rule callpeaks_macs2:
    input: 
        "bam/{sample}.downsample.sorted.bam", 
        "bam/{sample}.downsample.sorted.bam.bai"
    output: 
        bed = "macs2/{sample}_macs2_peaks.broadPeak"
    log: "logs/{sample}.call_peaks_macs2.log"
    params:
        name = "{sample}_macs2",
        jobname = "{sample}",
        g = 'mm', ## hs, mm, ce, dm
        pvalue = 1e-5,
        pvalue_broad = 1e-5

    message: "call_peaks macs2 {input}: {threads} threads"
    shell: 
        ## for macs2, when nomodel is set, --extsize is default to 200bp, 
        ## this is the same as 2 * shift-size in macs14.
        "macs2 callpeak -t {input[0]} "
        "--keep-dup all -f BAMPE -g {params.g} "
        "--outdir macs2 -n {params.name} -p {params.pvalue} "
        "--broad --broad-cutoff {params.pvalue_broad} &> {log}"



rule make_bigwigs:
    input: 
        "bam/{sample}.downsample.sorted.bam", 
        "bam/{sample}.downsample.sorted.bam.bai"
    output: 
        "bigwig/{sample}.bw"
    log: "logs/{sample}.makebw"
    threads: 5
    params: jobname = "{sample}"
    message: "making bigwig for {input} : {threads} threads"
    shell:
    	# no window smoothing is done, for paired-end, bamCoverage will extend the length to the fragement length of the paired reads
        "bamCoverage -b {input[0]} --ignoreDuplicates --skipNonCoveredRegions "
        "--normalizeUsingRPKM -p {threads} --extendReads -o {output} 2> {log}"



# ENCFF606PRC
# ENCFF643SJB
# ENCFF471PVN
# ENCFF819UNG
rule subtract_blacklist:
    input:
        bed = "macs2/{sample}_macs2_peaks.broadPeak",
        blacklist= BLACKLIST,
    output:
        "beds/{sample}.macs2_peaks.blacklist_removed.broadPeak"
    shell:
        "sort -k1,1 -k2,2n {input} | bedtools subtract -a stdin -b {input.blacklist} -A > {output}"