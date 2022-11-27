
"""
This pipline is modified for mouse from here:
http://barcwiki.wi.mit.edu/wiki/SOPs/atac_Seq
https://github.com/crazyhottommy/pyflow-ATACseq


Some tips:
ENCODE ATAC-seq standards: https://www.encodeproject.org/atac-seq/#standards

1) Post-alignment filtering: remove mito/ChrM reads, remove duplicates, remove mapping artifacts: < 38bp and fragments > 2kb, and discordant reads

2) Unless you have high coverage, use all reads for downstream analysis and not remove mono-, di-, tri-, etc. nucleosome fragments or reads. 
   Otherwise, you may miss many open chromatin regions during peak calling.

3) The nucleosome free region (NFR) should be > 38bp to < 147 bases (one nucleosome); mono nucleosome fragments are in the range of 147-200 bp
   Shorter reads lengths, e.g. 50x50 or 75x75, should be used instead of longer reads, e.g. 100x100 or longer, to ensure NFR/fragments are sequenced


4) ATAC-Seq depth or coverage:
   calling open/accessible regions, at least ~50M reads are recommended
   transcription factor footprinting, at least ~200M reads are recommended or optimal to ensure high coverage of the NFR


5) Tn5 produces 5' overhangs of 9 bases long: pos. strand +4 and neg strand -5 (see shiftGAlignmentsList and shiftReads functions in ATACseqQC package)
splitGAlignmentsByCut (in ATACseqQC package): creates different bins of reads, e.g. NFR, mono, di, etc. 
Shifted reads that do not fit into any of the bins should be discarded.

Use housekeeping genes to check QC: signal enrichment is expected in the regulatory regions of housekeeping genes in good ATAC-seq experiments. 
Use IGVSnapshot function with geneNames param.

"""

import os
import pandas as pd
# shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
# shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

workdir: "/data/bases/fangzq/20220711_CC/ATAC"

ALL_SAMPLES = ['Fadu_control', 'Fadu_NSD1-sh', 'Fadu_NSD1-sh_KDM2A-sh', 
               'PCI-13_control', 'PCI-13_NSD1-sh', 'PCI-13_NSD1-sh_KDM2A-sh']

BOWTIE2_INDEX = "/home/fangzq/genome/human/bowtie2Indices_GRCh38_noalt_as/GRCh38_noalt_as"
## download from https://sites.google.com/site/anshulkundaje/projects/blacklists
BLACKLIST = "/home/fangzq/genome/human/hg38-blacklist.v3.bed" 


TARGET_READS =  25000000 ## number of reads down sample to (25million default)
# genome size for bedtools slop
GENOME_SIZE =  "/home/fangzq/genome/human/GRCh38.p13.genome.size"
## genome fasta for nucleoATAC
GENOME_FASTA = "/home/fangzq/genome/human/hg38_masked.fa" # masked genome for meme-chip
JASPAR_MOTIF = "/home/fangzq/genome/human/Jaspar_hs_core_homer.motifs"
MEME_MOTIF = "/home/fangzq/genome/human/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme"
TSS_BED =  "/home/fangzq/genome/human/GRCh38.refGene.TSS.chr.bed"
MACS3_GENOME = "hs"
MACS3_PVAL = 1e-2
MACS3_PVAL_BROAD =  1e-2




############### OUTPUTS ##########
ALL_BIGWIG = expand("bigwig/{sample}.bw", sample = ALL_SAMPLES)
ALL_PEAKS = expand("peaks/{sample}_macs2_peaks.clean.narrowPeak", sample = ALL_SAMPLES)
ALL_NUCLEO = expand("peaks/{sample}_nucleoATAC.occpeaks.bed.gz",  sample = ALL_SAMPLES)
ALL_ANNOT = expand("peaks_annot/{sample}.peaksAnnotate.txt", sample = ALL_SAMPLES)
ALL_QC = ["multiQC/multiQC_log.html"]
ALL_ATAQV = "ATAC_qc/index.html"
ALL_MOTIF = expand("motif_enrichment/{sample}_motif/motifFindingParameters.txt", sample=ALL_SAMPLES)
PHANTOM = expand("phantompeakqual/{sample}_phantom.txt", sample = ALL_SAMPLES)
MEME = "meme_motif/meme.html"
HEATMAP = "figures/matrix.tss.gz"
rule all:
	input: ALL_PEAKS, ALL_BIGWIG,  ALL_MOTIF, ALL_QC, ALL_ATAQV, HEATMAP, ALL_ANNOT, MEME


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
        # using the following settings:
        # --very-sensitive 
        # --no-discordant    suppress discordant alignments for paired reads 
        # -X/--maxins <int>  maximum fragment length (default=500). Increase to 2000 to get a better nucleosome distribution.
        "bowtie2 --very-sensitive --threads {threads} "
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
    input: 
        bam = "BAM/{sample}.sorted.bam",
        #bai = "BAM/{sample}.sorted.bam.bai",
        tss = "refGene.hg38.TSS.bed"
    output: "BAM/{sample}.sorted.bam.ataqv.json"
    log: "log/{sample}_ataqv.log"
    threads: 1
    params:
        jobname = "{sample}",
        genome = 'human'
    message: "ataqv quality control for {input}"
    shell:
        """
        ataqv {params.genome} {input.bam} --metrics-file {output} 2> {log}
        """


rule json_to_html:
	input: expand("BAM/{sample}.sorted.bam.ataqv.json", sample = ALL_SAMPLES)
	output: "ATAC_qc/index.html"
	log: "log/ATAC_qc_html.log"
	threads: 1
	message: "compiling json files to html ATAC-seq QC"
	shell:
		"""	
		mkarv --force ATAC_qc  {input}
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
	run:
		# remove duplicates and reads on chrM, coordinate sort the bam
		# samblaster expects name sorted bamq
		shell("samtools view -h {input} | samblaster --removeDups "
              "| grep -v -P '\tchrM\t' "
              "| samtools view -Sb -F 4 - "
              "| samtools sort -m 10G -@ {threads} -T {input}.tmp -o {output[0]} ")
		shell("samtools index {output[0]}") 

# rule phantom_peak_qual:
#     input: "BAM/{sample}_exclude_chrM.sorted.bam"
#     output: "phantompeakqual/{sample}_phantom.txt"
#     log: "log/{sample}_phantompeakqual.log"
#     threads: 4
#     params: jobname = "{sample}"
#     message: "phantompeakqual for {input} : {threads} threads"
#     shell:
#         "/usr/bin/Rscript  /home/fangzq/program/phantompeakqualtools/run_spp.R "
#         "-c='{input}' -savp -rf  -p={threads} -odir=phantompeakqual  "
#         "-out='{output}' -tmpdir=phantompeakqual 2> {log} "


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


rule make_bigwigs_downsample:
    input : "BAM/{sample}-downsample.sorted.bam", "BAM/{sample}-downsample.sorted.bam.bai"
    output: "bigwig/{sample}.dowmsample.bw"
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

rule make_bigwigs:
    input : "BAM/{sample}_exclude_chrM.sorted.bam", "BAM/{sample}_exclude_chrM.sorted.bam.bai",
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



rule call_peaks_macs3:
    input: "BAM/{sample}_exclude_chrM.sorted.bam", "BAM/{sample}_exclude_chrM.sorted.bam.bai",
    output: 
        bed = "peaks/{sample}_macs2_peaks.narrowPeak",
        summit = "peaks/{sample}_macs2_summits.bed"
    log: "log/{sample}_call_peaks_macs2.log"
    params:
        name = "{sample}_macs2",
        jobname = "{sample}",
        g = MACS3_GENOME,
        qval = MACS3_PVAL,
    message: "call_peaks macs2 {input}: {threads} threads"
    shell:
        "/home/fangzq/miniconda/envs/fastai/bin/macs3 callpeak "
        "-t {input[0]} -f BAMPE -g {params.g} --call-summits "
        "--outdir peaks -n {params.name} --qvalue {params.qval} &> {log}"

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


rule genome_size:
    input: GENOME_FASTA
    output: GENOME_SIZE
    shell:
        # samtools faidx {input}
        # cut -f1,2 {input}.fai > {input}
        # pip install pyfaidx
        """
        faidx {input} -i chromsizes > {output}
        """


## extend the broad peak a bit for nucelosome analysis by nuceloATAC
rule make_bed_nucleoATAC:
    input: "peaks/{sample}_macs2_peaks.narrowPeak"
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
    input: 
        "BAM/{sample}-downsample.sorted.bam", 
        "BAM/{sample}-downsample.sorted.bam.bai", 
        "peaks/{sample}_nucleo.bed",
        GENOME_SIZE,
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
        bed = "peaks/{sample}_macs2_peaks.narrowPeak",
        blacklist= BLACKLIST,
    output:
        "peaks/{sample}_macs2_peaks.clean.narrowPeak"
    shell:
        "grep -P 'chr[\dXY]+[ \\t]' {input.bed} | "
        "sort -k1,1 -k2,2n | "
        "bedtools subtract -a stdin -b {input.blacklist} -A > {output}"
# bedtools intersect -v -a ${PEAK} -b ${BLACKLIST} \
#    		 | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
#    		 | grep -P 'chr[\dXY]+[ \t]'  | gzip -nc > ${FILTERED_PEAK}

## HOMER
# To download species specific Jaspar motifs, convert to homer motif format, and save to motif file.
# R code
# library(TFBSTools)
# # library(JASPAR2022)
# opts <- list()
# opts["collection"] <- "CORE"
# opts["species"] = 9606
# Jaspar_hs_core <- getMatrixSet(JASPAR2022::JASPAR2022, opts)
# # convert to homer motif format:
# library(universalmotif)
# write_homer (Jaspar_hs_core, file="Jaspar_hs_core_homer.motifs")


rule annotate_peaks:
    input: 
        peak = "peaks/{sample}_macs2_peaks.clean.narrowPeak",
        motif = JASPAR_MOTIF,
    output: 
        annot = "peaks_annot/{sample}.peaksAnnotate.txt",
        go = "peaks_annot/{sample}.GO/biological_process.txt",
        motif = "peaks_annot/{sample}.motif.bed"
    log: "log/{sample}.annotate_peaks.log"
    shell:
        # -m: motifs can be combined first and save as a file. 
        # In the output file, this will link motifs associated with a peak together. 
        # -mbed <filename> (Output motif positions to a BED file to load at genome browser)
        "annotatePeaks.pl {input.peak} hg38 "
        "-go peaks_annot/{wildcards.sample}.GO "
        "-m {input.motif} -mbed {output.motif} "
        "> {output.annot} 2> {log}"



rule moitf_enrichment:
    input:
        peak = "peaks/{sample}_macs2_peaks.clean.narrowPeak",
        motif = JASPAR_MOTIF,   
    output: "motif_enrichment/{sample}_motif/motifFindingParameters.txt"
    log: "log/{sample}.find_motif.log"
    params:
        genome = "hg38",
    message: "motif enrichment for {wildcards.sample}"
    threads: 6
    shell:  
        # snakemake use bash -c mode 
        # Escape certain characters, such as \t by \\t, $ by $, and { by {{.
        # Use triple quotation marks to surround the command line call. 
        """ 
        awk '{{print $4"\\t"$1"\\t"$2"\\t"$3"\\t+"}}' {input.peak} | \
        findMotifsGenome.pl - {params.genome} motif_enrichment/{wildcards.sample}_motif \
        -size 300 -S 2 -p {threads} -cache 100 -fdr 5 -mask \
        -mknown {input.motif} \
        -mcheck {input.motif} 2> {log}
        """
        # input parameters:
        # -mask: use the repeat-masked sequence
        # -size: (default 200). Explanation from homer website: 
        #        "If analyzing ChIP-Seq peaks from a transcription factor, 
        #        Chuck would recommend 50 bp for establishing the primary motif bound by a given transcription factor and 
        #        200 bp for finding both primary and "co-enriched" motifs for a transcription factor. 
        #        When looking at histone marked regions, 500-1000 bp is probably a good idea (i.e. H3K4me or H3/H4 acetylated regions). 
        # -mknown <motif file> (known motifs to check for enrichment.
        # -mcheck <motif file> (known motifs to check against de novo motifs,
        # -S: Number of motifs to find (default 25)




### promoters

# regions contains -a only
# bedtools window -a {tss} \
#                 -b macs_out_q0.01/SOX21_peaks.narrowPeak \
#                 -w 2000 -u > enhancers/SOX21.peaks.tss.bed


## get promoters that overlaped with promoter

#library(rtracklayer)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# promoter <- promoters(genes(txdb), upstream = 3000, downstream = 3000)
# peak <- import(peak_path, format = "narrowPeak")
# overlap_index <- findOverlaps(promoter, peak)
# keep_promoter <- unique(promoter[queryHits(overlap_index)])
# export.bed(keep_promoter, bed_path)

rule tss_overlap:
    input:
        peak = "peaks/{sample}_macs2_peaks.clean.narrowPeak",
        tss = TSS_BED,
    output: temp("peaks/{sample}.tss.bed")
    message: "get tss bed that overlapped with {wildcards.sample}"
    shell:
        # get promoters (-3k, +3k ) overlap with peaks, only keep autochromosome using grep -P 
        "grep -P 'chr[\dXY]+[ \\t]' {input.tss} | bedtools window -a stdin -b {input.peak} " 
        "-w 3000 -u > {output}"
        # # output regions contains -a -b, need cut to get -b regions only
        # cut -f 7- ${sample}.tss.bed > out
rule consensus_peak:
    input: 
        narrow = expand("peaks/{sample}_macs2_peaks.clean.narrowPeak", sample=ALL_SAMPLES),
        peak = expand("peaks/{sample}_macs2_summits.bed", sample=ALL_SAMPLES),
        tss = expand("peaks/{sample}.tss.bed", sample=ALL_SAMPLES)
    output:
        tss = "peaks/consensus.tss.bed",
        peak = "peaks/consensus.summits.bed",
        narrow = "peaks/consensus.narrowPeak"
    run:
        shell("cat {input.peak} | grep -P 'chr[\dXY]+[ \\t]' | sort -k1,1 -k2,2n | bedtools merge > {output.peak} ")
        shell("cat {input.narrow} | grep -P 'chr[\dXY]+[ \\t]' | sort -k1,1 -k2,2n | bedtools merge > {output.narrow} ")
        shell("cat {input.tss} | grep -P 'chr[\dXY]+[ \t]' | sort -k1,1 -k2,2n | bedtools merge | "
              "awk 'OFS=\"\\t\" {{if ($2 == $3) $3+=1; print}}' > {output.tss} ") # end add 1 if start == end

        
rule tss_enrichment:
    input:
        tss = "peaks/consensus.tss.bed",
        bw = expand("bigwig/{sample}.bw", sample=ALL_SAMPLES),
    output:
        mat = "figures/matrix.tss.gz",
        bed = "figures/genes.tss.bed",
        heatmap = "figures/heatmap.tss.pdf",
        profile = "figures/profile.tss.pdf",
    message: "tss enrichment heatmap and profile"
    threads: 8
    run:
        ## NOTE: the INPUT TSS file, each record have range length 1. so use reference-point cmd
        shell("computeMatrix reference-point --referencePoint TSS  -p {threads} "
              "-b 3000 -a 3000 -R {input.tss} -S {input.bw}  "
              "--skipZeros -o {output.mat}  "
              "--outFileSortedRegions {output.bed}")
        ## both plotHeatmap and plotProfile will use the output from   computeMatrix
        shell("plotHeatmap -m {output.mat} -out {output.heatmap} --plotFileFormat pdf  --dpi 720")  
        shell("plotProfile -m {output.mat} -out {output.profile} --plotFileFormat pdf --perGroup --dpi 720")



# # motif 的分析往往受假阳性困扰，这称为无效定理（Futility Theorem）。
# # 比如说往往在基因组序列中观察到大量的潜在转录因子结合位点，其中很少是真正起作用的，
# # 大部分预测的转录因子结合位点是无效的。所以 motif 分析后，要想办法进行人工筛选。

# # MEME-ChIP 主要执行以下步骤：

# # 在输入序列的中间区域（默认 100bp）进行 motif 发现（MEME, STREME）
# # CentriMo 分析哪些 motif 是富集在区域中心的
# # Tomtom 分析他们与已知 motif 相似性
# # 根据 motif 相似性对显著结果归类
# # motif spacing analysis (SpaMo) (分析 motif 与结合在相邻位置 motif 的物理作用) (-spamo-skip 参数跳过这步分析)
# # 分析 motif 结合位置（FIMO），默认选取每一类别最显著的 motif 进行分析
# # MEME-ChIP 默认输入序列是约 500bp 长且中间 100bp 是 motif 区域。
# # ATAC-Seq 分析得到 peak 是长度不一的，因此取 summit 向两边各延申 250bp 得到 500bp 序列。

## see tutorial here https://github.com/roonysgalbi/memeMotifs
## the ATAC-seq part
rule get_peak_fasta:
    input: 
        peak=  "peaks/consensus.summits.bed",
        genome = GENOME_FASTA,
    output: "meme_motif/consensus.peaks.fasta"
    shell:
        # get  mid length
        "awk -v FS='\\t' -v OFS='\\t' '{{midpos=$2;if ((midpos-250) > 0) print $1,midpos-250,midpos+250;}}' {input.peak}  | "
        "bedtools getfasta -fo {output} -fi {input.genome} -bed stdin "
        # cut -f 1-5 ${beddir}/${bed} | bedtools slop -i stdin -g ${genome_size} -b 100 | \
        # bedtools getfasta -fi ${genome_seq} -bed stdin -bedOut > peaks_fasta/${bed}.slop100.fasta.bed

# # sort by q-value
# sort -k9nr sample.narrowPeak >sample.sorted.narrowPeak
# # select the top 1000 peaks
# head -1000 sample.sorted.narrowPeak >sample.top1000.narrowPeaks
# # create a bed file of 500bp regions centered on the peak summits
# awk 'BEGIN{ OFS="\t";}{ midPos=$2+$10; print $1, midPos-250, midPos+250; }' sample.top1000.narrowPeaks >sample.regions.bed
# # create fasta file
# fastaFromBed -fi mm10_masked.fa -bed sample.regions.bed -fo sample.sequences.fa

# rule get_background_model:
#     input: "meme_motif/consensus.peaks.fasta"
#     output: "meme_motif/background.model"
#     shell:
#         "fasta-get-markov -m 2 -dna -nostatus -nosummary {input} {output}"

# # motif 数据库在 https://meme-suite.org/meme/db/motifs -> 下载 https://meme-suite.org/meme/doc/download.html。
# # MEME-ChIP 输出结果在 meme.html 查看；结果的汇总在 summary.tsv 文件；combined.meme 文件包含所有被鉴定出的 Motif.
# # motif E value 表示该 motif 多大概率是统计错误出现，而不是真的 motif. E 值越小说明发现的 motif 越大概率为真。

## see tutorial here https://github.com/roonysgalbi/memeMotifs
## the ATAC-seq part

rule meme_chip:
    input: 
        narrow = "peaks/consensus.narrowPeak",
        meme_db = MEME_MOTIF,
        #bg_model = "meme_motif/background.model",
        fasta = "meme_motif/consensus.peaks.fasta"
    output:
        "meme_motif/summary.tsv",
        "meme_motif/meme.html"
    params: 
        outdir = "meme_motif",
    threads: 36
    run:
        # peaks = pd.read_table(input.narrow, header=None)
        # mean_peak_len = (peaks.iloc[:,2] - peaks.iloc[:,1]).abs().mean() # int number
        # max_size = int(600 * mean_peak_len)
        # -spamo-skip -> skip, takes too long
        # -db: can set multiple times
        # # -meme-maxsize -> see docs in meme -h
        # # -ccut 0 for ATAC-seq to search whole input sequence, not the center 100
        shell("meme-chip -maxw 30 -minw 6 -oc {params.outdir} -spamo-skip "
              "-dna -ccut 0 -db {input.meme_db} "
              "-meme-p {threads} -meme-mod zoops -meme-nmotifs 10 -meme-nrand {input.fasta}") 

# for MEME-CHIP < v5.0
# For each peak pf 500bp the default is to use the centre 100 for motif searching and the remaining 400 for background. 
# but for ATAC-seq i want to search the whole 500 bases. So I must provide a background. #
# The max data size limit is 100000, calculated by ccut x nmeme (default: 100x600=60000). 
# To increase max size use -meme-maxsize [600 x mean peak length] and -ccut 0, 
# so no cutting occurs and all of each region is used. Assuming MACS2 was used to call peaks:

# calculate the mean peak length from narrow_peak.xls (column 4)
# get peak summit from narrow_peak.xls (column 5) or narrow_summits.bed (column 3)
# bed file region length = midPos+/- (1/2 mean peak length)
# create background file (fasta-get-markov -m 2)
# add to meme-chip command: -ccut = 0 --meme-maxsize = [600 x mean peak length]


# # # MEME-ChIP 的结果如果发现不方便批量提取结果，也可以自己单独运行它的子分析。比方说它 FIMO 只分析了部分 motif 而你需要分析全部的，那么可以单独进行 FIMO 分析。
rule meme_fimo:
    input: 
        meme_db = MEME_MOTIF,
        fasta = "meme_motif/consensus.peaks.fasta"
    output: "meme_motif/fimo_out/fimo.tsv"
    params: 
        outdir = "meme_motif/fimo_out/"
    shell:
        "fimo --parse-genomic-coord -oc {params.outdir} "
        "--motif SP1_HUMAN.H11MO.1.A "
        "--motif SMAD1_HUMAN.H11MO.0.D "
        "--motif SMAD2_HUMAN.H11MO.0.A "
        "--motif SP1_HUMAN.H11MO.1.A " 
        "{input.meme_db} {input.fasta}"