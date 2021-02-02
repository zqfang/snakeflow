
import os,re
############### Globals ########################
#configfile: "config.yaml"
workdir: config['WORKSPACE']

GENOME = config['GENOME']
dbSNP = config['dbSNP']
dbINDEL = config['dbINDEL']
BAM_DIR = config['BAM_DIR']
TMPDIR = config['TMPDIR']
STRAINS = config['STRAINS']

#CHROMSOME = [ str(c) for c in range(1,20)] + ["X", "Y", "MT"]
CHROMSOME = ['1'] + [ str(c) for c in range(10,20)] + [ str(c) for c in range(2,10)]+ ["X", "Y"]

#outputs
BAMS = expand("BAM/{sample}.marked.fixed.BQSR.bam", sample=STRAINS)

### refined each time running !!!
# automatic search sample names using regex
pat = re.compile(r"([A-Z0-9-]+)-([A-Z0-9a-z]{2,4})_[0-9A-Z]{2}_(L[0-9]{3})_([RI][12])_[0-9]{3}.fastq.gz")
# name, capture, lane, reads = pat.search(f).groups()
files = glob.glob(os.path.join(config['FASTQS']['PATH'], "*fastq.gz"))
files = [pat.search(f).groups() for f in files]
# remove duplicate entries
files = {s[0]: { x[1]: x[0] for x in files } for s in files}

#################### rules #######################
rule target:
    input: BAMS

rule bwa_index:
    input: GENOME,
    output: GENOME +".bwt"
    params: 
        prefix="GRCm38"
    log: "logs/bwa.index.log"
    shell:
        "bwa index -a bwtsw -p {params.prefix} {input} 2> {log}"

rule bwa_men:
    input:
        r1="fastq/{sample}_R1.fastq.gz",
        r2="fastq/{sample}_R2.fastq.gz",
        index= GENOME + ".bwt"
    ouput: temp("{sample}.sam")
    threads: 8
    log: "logs/bwa.mem.log"
    params:
        # PL has to be one of ILLUMINA，SLX，SOLEXA，SOLID，454，LS454，COMPLETE，
        # PACBIO，IONTORRENT，CAPILLARY，HELICOS, or UNKNOWN
        # ID：输入reads集的ID号； LB： reads集的文库名； SM：样本名称； PL：测序平台
        RG="@RG\tID:SAMPLE_ID\tPL:illumina\tSM:{sample}"
    shell:
        "bwa mem -t {threads} -R '{param.RG}' {input.index} "
        "{input.r1} {input.r2} 1> {output} 2> {log}"
        
rule sam_sort:
    input: "BAM/{sample}.sam"
    output: temp("BAM/{sample}.bam")
    params:
        tmpdir=TMPDIR,
        java_ops= "-Xmx32G -Djava.io.tmpdir=%s"%TMPDIR
    shell:
        "gatk --java-options '{params.java_ops}' SortSam -I {input} -O {output} "
        "-SO coordinate --CREATE_INDEX true"

rule mark_dups:
    input:  "BAM/{sample}.bam"
    output: 
        bam = temp("{sample}.marked.bam")
        metric = temp("{sample}.marked.bam.metrics")
    params:
        extra = "--CREATE_INDEX" # if skip fixMateInformation, you need this
    shell:
        "gatk MarkDuplicates -I {input} -O {output.bam} -M {output.metrics}"

rule fix_mate:
    input: "BAM/sample}.marked.bam"
    output: temp("BAM/{sample}.marked.fixed.bam")
    shell:
        "gatk FixMateInformation -I {input} -O {output.bam} "
        "-SO coordinate --CREATE_INDEX true"

rule genome_dict:
    input: GENOME
    output: GENOME + ".dict"
    shell:
        "gatk CreateSequenceDictionary -R {input} -O {output}"

## dbSNP, dbINDEL need to be indexed
# rule vcf_index:
#     input: dbINDEL, dbSNP
#     output: dbINDEL +".idx", dbSNP +".idx"
#     shell:
#         """gatk IndexFeatureFile -F {input[0]}
#            gatk IndexFeatureFile -F {input[1]}
#         """

## GATK HaplotypeCaller will assemble haplotype itself, you skip Recalibrate if use GATK 
rule baseRecalibrate:
    input: 
        genome=GENOME,
        gdict = GENOME + ".dict",
        bam = "BAM/{sample}.marked.fixed.bam",
        bai = "BAM/{sample}.marked.fixed.bai",
    output: 
        "BAM/{sample}.marked.fixed.table"
    params:
        dbsnp = dbSNP,
        dbindel = dbINDEL,
    shell:
        "gatk BaseRecalibrator -R {input.genome} -I {input.bam} -O {output} "
        "--known-sites {params.dbsnp} --known-sites {params.dbindel}"

rule applyBQSR:
    input: 
        genome=GENOME,
        gdict = GENOME + ".dict",
        bam = "BAM/{sample}.marked.fixed.bam",
        bai = "BAM/{sample}.marked.fixed.bai",
        table = "BAM/{sample}.marked.fixed.table",
    output: 
        protected("BAM/{sample}.marked.fixed.BQSR.bam")
    params:
        dbsnp=dbSNP,
        dbindel=dbINDEL
    shell:
        "gatk ApplyBQSR -R {input.genome} -I {input.bam} "
        "--bqsr-recal-file {input.table} -O {output}"
