from os.path import join



#Workding directory
workdir: config['workdir']

# Full path to a folder that holds all of your FASTQ files.
FASTQ_DIR = config['fastq_dir']

PAIRED = 'paired' if config['paired'] else 'single'
# Full path to a Genome.
GENOME = config['genome']
#CDNA =           join(GENOME,"gencode.v25.transcripts.fa")
# genome sequence
FASTA_REF =     config['dna']
# index_dir
BOWTIE_REFDIR= config['bowtie2_index']
# index basename
INDEX_PREFIX = config['index_prefix']
_SAMPLES = config['samples']

############ Samples ##################
# A Snakemake regular expression matching the forward mate FASTQ files.
# the part in curly brackets {} will be saved, so the variable SAMPLES
# is a list of strings #['Sample1','Sample2'].

#notice that SAMPLES, has a trailing comma.
#you must include this trailing comma, or else the code won’t work correctly.

#SAMPLES, = glob_wildcards(join(FASTQ_DIR, '{sample, SRR[^/]+}_R1.fastq.gz'))

# Patterns for the 1st mate and the 2nd mate using the 'sample' wildcard.
#PATTERN_R1 = '{sample}_R1.fastq.gz'
#PATTERN_R2 = '{sample}_R2.fastq.gz'
PATTERN_R1 = config['read_pattern']['r1']
PATTERN_R2 = config['read_pattern']['r2']

with open(_SAMPLES, 'r') as s:
    SAMPLES = [ l.strip() for l in s]
print(SAMPLES)
#PERRETY_SAMPLE = expand("mapped/{sample}.{replicate}.bam", sample=SAMPLES, replicate=[0, 1])
#SAMPLES,GROUP,IDS,= glob_wildcards(join(FASTQ_DIR, '{sample}_{group}_{id}_R1.fastq.gz'))
#######output######
_BW = expand("mapped/{sample}.q25.bw", sample=SAMPLES)
# _MACS = "macs_out"
# _HOMER 
# _MOTIF ="motif"
# _GO = "peaks.GO.HOMMER"
# _GREAT = "peaks.GO.GREAT"


rule target:
    input: _BW

# rule bowtie_index:
#     input:
#         fasta = FASTA_REF,
#         gtf = GTF_FILE,
#     output: expand(join(BOWTIE_REFDIR,INDEX_PREFIX)+".{ids}.bt2",ids=range(1,4))
#     params:
#         basename=join(BOWTIE_REFDIR, INDEX_PREFIX)
#     log: "logs/bowtie2/bowtie2.index.build.log"
#     threads: 12
#     shell: "bowtie2-build -f {input.fasta}  -p {threads}  {params.basename} &> {log}"


# rule botiwe_align_single:
#     input:
#         index=expand(join(BOWTIE_REFDIR, INDEX_PREFIX)+".{ids}.bt2", ids=range(1,4)),
#         r1 = join(FASTQ_DIR, "{sample}.fastq")
#     output:
#         temp('mapped/{sample}.bam')
#     log:
#         "logs/bowtie2/{sample}.align.log"
#     threads: 12
#     params:
#         ref = join(BOWTIE_REFDIR, INDEX_PREFIX),
#     shell:
#         "(bowtie2 -p {threads} -x {params.ref} -U {input.r1} "
#         " | samtools view -Sbh  -q 25 -@ {threads}  -o {output} - ) 2> {log}"


rule botiwe_align_paired:
    input:
        index=expand(join(BOWTIE_REFDIR, INDEX_PREFIX)+".{ids}.bt2", ids=range(1,4)),
        r1 = join(FASTQ_DIR, PATTERN_R1),
        r2 = join(FASTQ_DIR, PATTERN_R2)
    output:
        temp('mapped/{sample}.bam')
    log:
        "logs/bowtie2/{sample}.align.log"
    threads: 12
    params:
        ref = join(BOWTIE_REFDIR, INDEX_PREFIX),
    shell:
        "(bowtie2 -p {threads} -x {params.ref} -1 {input.r1} -2 {input.r2} "
        " | samtools view -Sbh  -q 25 -@ {threads}  -o {output} - ) 2> {log}"

rule bam_sort:
    input: "mapped/{sample}.bam"
    output: protected("mapped/{sample}.q25.sorted.bam")
    threads: 12
    shell: 
        "samtools sort -@ {threads} {input} > {output}"


rule bam_index:
    input: "mapped/{sample}.q25.sorted.bam"
    output: "mapped/{sample}.q25.sorted.bam.bai"
    shell: 
        "samtools index {input}"

rule bam2bw:
    """deeptools"""
    input: 
        bam = "mapped/{sample}.q25.sorted.bam",
        bai =   "mapped/{sample}.q25.sorted.bam.bai"
    output:
        "mapped/{sample}.q25.bw"
    log: "logs/deeptools/{sample}.bam2bw.log"
    threads: 8
    params:
        fsize=200,
        extra=" --normalizeUsing CPM --centerReads",
    shell:
        "bamCoverage -b {input.bam} -p {threads} --extendReads {params.fsize} -o {output} {params.extra} &> {log}"


# rule prepare_macs:
#     input: config['groups']
#     output: "temp/{sample}.{treat}.{ctrl}.txt"
#     params: 
#     run:
#         groups=[]
#         with open(input[0],'r') as f:
#             for line in f:
#                 groups.append(line.strip("\n").split("\t"))
       
#         "touch temp/{sample}.{treat}.{ctrl}.txt"

# rule macs_narrow:
#     """point source factors peaks calling.
#        Transprition factors
#        H3K4me3
#     """
#     input: 
#         treat="mapped/{treat}.highQuality.q25.sorted.bam",
#         ctrl="mapped/{ctrl}.highQuality.q25.sorted.bam",
#         temp="temp/{sample}.{treat}.{ctrl}.txt"
#     output: 
#         "macs_out/{sample}_summits.bed",
#         "macs_out/{sample}_peaks.narrowPeak",
#     conda=MACS2_ENV,
#     log: "logs/macs/{sample}.macs2.log"
#     prams:
#         extra=" -f BAM -g hs -B --SPMR -q 0.01 ",
#         extra2=" --SPMR --nomodel --extsize 200",
#     shell:
#         "macs2 callpeak -t {input.treat} -c {input.ctrl}"
#         " --outdir macs2_highQuality_results -n {wildcards.sample} "
#         " {params.extra} 2> {log}"

# rule macs_broad:
#     """use for braod domains calling, H3K27ac, H3K4me1.et. al"""
#     input: 
#         treat="mapped/{treat}.highQuality.q25.sorted.bam",
#         ctrl="mapped/{ctrl}.highQuality.q25.sorted.bam",
#         temp="temp/{sample}.{treat}.{ctrl}.txt"
#     output: 
#         "macs_out/{sample}_summits.bed",
#         "macs_out/{sample}_peaks.broadPeak",
#     conda=MACS2_ENV,
#     log: "logs/macs/{sample}.macs2.log"
#     params:
#         extra=" -f BAM -g hs -B --SPMR --broad --call-summits",
#         extra2=" --nomodel --extsize 147  -q 0.1 --fe-cutoff 1.5",
#     shell:
#        "macs2 callpeak -t {input.treat} -c {input.ctrl} "
#        "--outdir macs2_highQuality_results -n {wildcards.sample}  --broad"
#        "{params.extra} {params.extra2} 2> {log}"

# rule annotatepeaks:
#     input:
#         bed="{sample}.bed"
#     output:
#         "{sample}.peaksAnnotate.txt"
#     params:
#         go_outdir=".",
#         genome = "hg19"
#     shell:
#         "annotatePeaks.pl {input.bed} {params.genome} -go {params.go_outdir} > {output}"

# rule findmoitf:
#     input:
#         bed="{sample}.bed"
#     output:
#     params:
#         go_outdir="",
#         genome = "hg19"
#     shell:   
#         """ 
#         awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' {input.bed} | \
#         findMotifsGenome.pl - hg19 ./findMotif/NANOG_motif -len 8,10,12 -size given 
#         """
# rule beta:
#     """differential motif finding"""

# rule region_profile:
#     input:
#     output:
#     shell:

# rule region_heatmap:
#     input:
#     output:
#     shell:

# rule screen_shoot:
#     input:
#     output:
#     shell:


