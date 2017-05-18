from os.path import join



configfile: 'config.yml'
#Workding directory
workdir: config['workdir']

# Full path to a folder that holds all of your FASTQ files.
FASTQ_DIR = config['fastq_dir']

PAIRED = 'paired' if config['paired'] else 'single'
# Full path to a Genome.
GENOME = config['genome']
#CDNA =           join(GENOME,"gencode.v25.transcripts.fa")
# genome sequence
FASTA_REF =     config['fasta']
# index_dir
BOWTIE_REFDIR= config['index_dir']
# index basename
INDEX_PREFIX = config['index_prefix']

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

#PERRETY_SAMPLE = expand("mapped/{sample}.{replicate}.bam", sample=SAMPLES, replicate=[0, 1])
#SAMPLES,GROUP,IDS,= glob_wildcards(join(FASTQ_DIR, '{sample}_{group}_{id}_R1.fastq.gz'))
#######output######
_MACS
_HOMER
_MOTIF
_GO
_GREAT


rule target:
    input: _MACS,_MACS, _HOMER, _MOTIF, _GO, _GREAT


rule fastqc:
    input:
        join(FASTQ_DIR,"{prefix}.fastq.gz"),
    output:
        "qc/fastqc/{prefix}_fastqc.html",
        "qc/fastqc/{prefix}_fastqc.zip",

    shell: "fastqc -o qc/fastqc {input}"

rule bowtie2_index:
    input:
        fasta = FASTA_REF,
        gtf = GTF_FILE,
    output: expand(join(BOWTIE_REFDIR,INDEX_PREFIX)+".{ids}.ht2",ids=range(1,9))
    params:
        basename=join(BOWTIE_REFDIR, INDEX_PREFIX)
    log: "logs/hisat2/hisat2.index.build.log"
    threads: 12
    shell: "bowtie2-build -f {input.fasta}  -p {threads}  {params.basename} &> {log}"


rule hisat2_align:
    input:
        index=expand(join(HISAT2_REFDIR,INDEX_PREFIX)+".{ids}.ht2", ids=range(1,9)),
        site = join(BOWTIE_REFDIR, "splicesites.txt"),
        r1 = join(FASTQ_DIR, PATTERN_R1),
        r2 = join(FASTQ_DIR, PATTERN_R2)
    output:
        temp('mapped/{sample}.highQuality.q25.bam')
    log:
        "logs/hisat2/{sample}.align.log"
    threads: 12
    params:
        ref = join(BOWTIE_REFDIR, INDEX_PREFIX),
    shell:
        "(bowtie2 {params.extra} -p {threads} -x {params.ref} -U {input.r1} "
        " | samtools view -Sbh  -q 25 -@ {threads}  -o {output} - ) 2> {log}"


rule bam_sort:
    input: "mapped/{sample}.bam"
    output: protected("mapped/{sample}.highQuality.q25.sorted.bam")
    threads: 12
    shell: 
        "samtools sort -T mapped/{wildcards.sample} -@ {threads} -O bam  {input} > {output}"


rule bam_index:
    input: "mapped/{sample}.highQuality.q25.sorted.bam"
    output: "mapped/{sample}.highQuality.q25.sorted.bam.bai"
    shell: 
        "samtools index {input}"

rule bam2bw:
    input: "mapped/{sample}.highQuality.q25.sorted.bam",
           "mapped/{sample}.highQuality.q25.sorted.bam.bai"
    output:
        "mapped/{sample}.highQuality.q25.bigwig"
    log: "logs/deeptools/{sample}.bam2bw.log"
    params:
        fragement_size=200,
        extra=" --normalizeUsingRPKM",
    shell:
        "bamCoverage -b {input} -e {params.fragement_size} -o {output} &> {log}"
rule macs_narrow:
        input: 
            treat="mapped/{sample}.highQuality.q25.sorted.bam",
        output: _MACS
        conda=MACS2_ENV,
        log: ""
        prams:
            extra=" ",
        shell:
            " macs2 callpeak -t ./bowtie_out/SOX1-CT10-ChIP.highQuality.sorted.bam "
             "-c ./bowtie_out/SOX1-CT10-Input.highQuality.sorted.bam -f BAM -g hs "
             "--outdir macs_out2 -n SOX1_CT1.0_ChIP -B --SPMR --call-summits --nomodel --extsize 200 "

rule macs_broad:

           "macs2 callpeak -t ./bowtie_out/KDme2ChIP.sam  -c ./bowtie_out/KDinput.sam"
           "-f SAM -g mm --outdir macs_out -n KDme2ChIP -B --SPMR --nomodel --extsize 147 --broad"



