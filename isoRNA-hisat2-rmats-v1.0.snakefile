from os.path import join
from snakemake.utils import R
from snakemake.shell import shell

configfile: 'config.yml'
#Workding directory
workdir: config['workdir']

# utils function
def unique(seq):
    """Remove duplicates from a list in Python while preserving order.
    :param seq: a python list object.
    :return: a list without duplicates while preserving order.
    """    
    seen = set()
    seen_add = seen.add

    return [x for x in seq if x not in seen and not seen_add(x)]

################### globals #############################################

# Full path to an uncompressed FASTA file with all chromosome sequences.
CDNA = config['cdna']

# Full path to a folder that holds all of your FASTQ files.
FASTQ_DIR = config['fastq_dir']
READ_LEN = config['read_length']
PAIRED = 'paired' if config['paired'] else 'single'
# Full path to a Genome.
GENOME = config['genome']
#CDNA =           join(GENOME,"gencode.v25.transcripts.fa")
# genome sequence
FASTA_REF =     config['fasta']
# index_dir
HISAT2_REFDIR= config['index_dir']
# index basename
INDEX_PREFIX = config['index_prefix']
# gtf
GTF_FILE =       config['gtf']
GTF_Genes =      GTF_FILE.rstrip(".gtf")+".extracted.genes.annotation.txt"
GTF_Trans =      GTF_FILE.rstrip(".gtf")+".extracted.transx2gene.txt"
############ Samples ##################
# A Snakemake regular expression matching the forward mate FASTQ files.
# the part in curly brackets {} will be saved, so the variable SAMPLES
# is a list of strings #['Sample1','Sample2'].

#notice that SAMPLES, has a trailing comma.
#you must include this trailing comma, or else the code won’t work correctly.

#SAMPLES, = glob_wildcards(join(FASTQ_DIR, '{sample, SRR[^/]+}_R1.fastq.gz'))
SAMPLES = config['samples']['name'].split()
GROUP=config['samples']['group'].split()
TIME=config['samples']['time'].split()

group_unique = unique(GROUP)
grpname1=group_unique[0]
grpname2=group_unique[1]
GROUP1 = [s for s, g in zip(SAMPLES, GROUP) if g == grpname1]
GROUP2 = [s for s, g in zip(SAMPLES, GROUP) if g == grpname2]
# Patterns for the 1st mate and the 2nd mate using the 'sample' wildcard.
#PATTERN_R1 = '{sample}_R1.fastq.gz'
#PATTERN_R2 = '{sample}_R2.fastq.gz'
PATTERN_R1 = config['read_pattern']['r1']
PATTERN_R2 = config['read_pattern']['r2']



####### Tools dir#################
RMATS_DIR=config['rmats']['dir']

#rmats runing envs
RMATS_ENV = config['rmats']['conda_env']

# dirs 
DIRS = ['qc','mapped','counts','alternative_splicing', 'gene_expression',
        'differential_expression','logs','temp']


########### Target output files #################

FASTQC = expand("qc/fastqc/{prefix}_fastqc.{suf}", prefix=SAMPLES,suf=['html','zip'])

RMATS = expand("alternative_splicing/{sample1}_vs_{sample2}/MATS_output/{type}.MATS.ReadsOnTargetAndJunctionCounts.txt",
                sample1= grpname1, sample2= grpname2, type=['A3SS','A5SS','MXE','RI','SE'])
################## Rules #######################################

rule target:
    input: RMATS

rule hisat2_index:
    input:
        fasta = FASTA_REF,
        gtf = GTF_FILE,
    output: expand(join(HISAT2_REFDIR,INDEX_PREFIX)+".{ids}.ht2",ids=range(1,9))
    params:
        basename=join(HISAT2_REFDIR, INDEX_PREFIX)
    log: "logs/hisat2/hisat2.index.build.log"
    threads: 12
    shell: "hisat2-build -f {input.fasta}  -p {threads}  {params.basename} &> {log}"


rule hisat2_extract_splicesites:
    input: GTF_FILE
    output:
        splice = join(HISAT2_REFDIR, 'splicesites.txt'),
        exon =   join(HISAT2_REFDIR, 'exon.txt')
    threads: 12
    shell:
        """
        hisat2_extract_splice_sites.py {input} > {output.splice}
        hisat2_extract_exons.py {input} > {output.exon}
        """

rule hisat2_align:
    input:
        index=expand(join(HISAT2_REFDIR,INDEX_PREFIX)+".{ids}.ht2", ids=range(1,9)),
        site = join(HISAT2_REFDIR, "splicesites.txt"),
        r1 = join(FASTQ_DIR, PATTERN_R1),
        r2 = join(FASTQ_DIR, PATTERN_R2)
    output:
        temp('mapped/{sample}.bam')
    log:
        "logs/hisat2/{sample}.align.log"
    threads: 12
    params:
        ref = join(HISAT2_REFDIR, INDEX_PREFIX),
        extra="--min-intronlen 1000 --dta -t"
    shell:
        "(hisat2 {params.extra} --threads {threads} -x {params.ref}"
        " -1 {input.r1} -2 {input.r2}  --known-splicesite-infile {input.site}"
        " | samtools view -Sbh -@ {threads}  -o {output} - ) 2> {log}"

rule bam_sort:
    input: "mapped/{sample}.bam"
    output: protected("mapped/{sample}.sorted.bam")
    threads: 12
    shell: "samtools sort -T mapped/{wildcards.sample} -@ {threads} -O bam  {input} > {output}"

rule bam_index:
    input: "mapped/{sample}.sorted.bam"
    output: "mapped/{sample}.sorted.bam.bai"
    shell: "samtools index {input}"


# Alternative splicing
rule rmats:
    input:
        bam=expand("mapped/{sample}.sorted.bam", sample=SAMPLES),
        bai=expand("mapped/{sample}.sorted.bam.bai", sample=SAMPLES),
        gtf=GTF_FILE
    output: RMATS
    conda: RMATS_ENV
    params:
        b1=",".join(expand("mapped/{sample}.sorted.bam", sample=GROUP1)),
        b2=",".join(expand("mapped/{sample}.sorted.bam", sample=GROUP2)),
        rmats=RMATS_DIR,
        prefix="alternative_splicing/{sample1}_vs_{sample2}".format(sample1=grpname1, sample2=grpname2),
        extra="-t %s -len %s -a 1 -c 0.0001 -analysis U -novelSS 0"%(PAIRED, READ_LEN)

    shell:
        "python2 {params.rmats} -b1 {params.b1} -b2 {params.b2} "
        "-gtf {input.gtf}  -o {params.prefix} {params.extra}"

