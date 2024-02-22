
import os,re
############### Globals ########################
#configfile: "config.yaml"
workdir: "/data/bases/fangzq/UDN/Hearing_loss/test"
BUNDLE="/home/fangzq/genome/gatk_bundle/hg38"
GENOME = os.path.join(BUNDLE, "Homo_sapiens_assembly38.fasta")
dbSNP =  os.path.join(BUNDLE, "dbsnp_146.hg38.vcf.gz")
db1K =  os.path.join(BUNDLE, "1000G_phase1.snps.high_confidence.hg38.vcf.gz")
dbINDEL =  os.path.join(BUNDLE,"Mills_and_1000G_gold_standard.indels.hg38.vcf.gz")

SAMPLES = ["RQ82599"]
#CHROMSOME = [ str(c) for c in range(1,20)] + ["X", "Y", "MT"]
CHROMSOME = ['1'] + [ str(c) for c in range(10,20)] + [ str(c) for c in range(2,10)]+ ["X", "Y"]

#outputs
BAMS = expand("BAM/{sample}.sorted.bam", sample=SAMPLES)
BQSR = expand("BAM/{sample}.markdup.fixed.BQSR.bam", sample=SAMPLES)
VCFS = expand("VCF/{sample}.call.raw.vcf.gz", sample=SAMPLES)
### refined each time running !!!
# automatic search sample names using regex
# pat = re.compile(r"([A-Z0-9-]+)-([A-Z0-9a-z]{2,4})_[0-9A-Z]{2}_(L[0-9]{3})_([RI][12])_[0-9]{3}.fastq.gz")
# # name, capture, lane, reads = pat.search(f).groups()
# files = glob.glob(os.path.join(config['FASTQS']['PATH'], "*fastq.gz"))
# files = [pat.search(f).groups() for f in files]
# # remove duplicate entries
# files = {s[0]: { x[1]: x[0] for x in files } for s in files}

#################### rules #######################
rule target:
    input: BAMS, BQSR, VCFS

# rule fastp_pe:
#     input:
#         r1="fastq_path/{sample}_R1.fastq.gz",
#         r2="fastq_path/{sample}_R2.fastq.gz",
#     output:
#         r1="trim/{sample}_clean_R1.fastq.gz", 
#         r2="trim/{sample}_clean_R2.fastq.gz",
#         # html="report/{sample}.fastp.html",
#         # json="report/{sample}.fastp.json",
#     # log:
#     #     "logs/fastp.{sample}.log"
#     threads: 8
#     shell:
#         "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2}"

rule bwa_index:
    input: GENOME,
    output: GENOME +".bwt"
    params: 
        prefix="Homo_sapiens_assembly38"
    log: "logs/bwa.index.log"
    shell:
        # for genome more than 2G , add bwtsw
        "bwa index -a bwtsw {input} 2> {log}"

rule bwa_men:
    input:
        r1="trim/{sample}_clean_R1.fastq.gz", 
        r2="trim/{sample}_clean_R2.fastq.gz",
        index= GENOME + ".bwt",
        genome=GENOME
    output: temp("BAM/{sample}.temp.bam")
    threads: 8
    log: "logs/bwa.mem.{sample}.log"
    params:
        # PL has to be one of ILLUMINA，SLX，SOLEXA，SOLID，454，LS454，COMPLETE，
        # PACBIO，IONTORRENT，CAPILLARY，HELICOS, or UNKNOWN
        # ID：input read group's ID； LB： read group's library name； SM：sample name； PL：sequencing platform
        RG="@RG\tID:{sample}\tPL:illumina\tSM:{sample}"
    shell:
        "bwa mem -t {threads} -R '{params.RG}' {input.genome} "
        "{input.r1} {input.r2} | samtools view -S -b - > {output}"
        
rule bam_sort:
    input: "BAM/{sample}.temp.bam"
    output: "BAM/{sample}.sorted.bam"
    threads: 8
    shell:
        "samtools sort -@ {threads} -O bam -o {output} {input}"

rule mark_dups:
    input:  "BAM/{sample}.sorted.bam"
    output: 
        bam = temp("BAM/{sample}.markdup.bam"),
        metric = "BAM/{sample}.markdup.metrics"
    params:
        extra = "--CREATE_INDEX" # if skip fixMateInformation, you need this
    shell:
        "gatk MarkDuplicates -I {input} -O {output.bam} -M {output.metric}"

rule fix_mate:
    input: "BAM/{sample}.markdup.bam"
    output: 
        bam=temp("BAM/{sample}.markdup.fixed.bam"),
        bai=temp("BAM/{sample}.markdup.fixed.bai")
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
#     output: dbINDEL +".tbi", dbSNP +".tbi"
#     shell:
#         """gatk IndexFeatureFile -I {input[0]}
#            gatk IndexFeatureFile -I {input[1]}
#         """

## GATK HaplotypeCaller will assemble haplotype itself, you skip Recalibrate if use GATK 
rule baseRecalibrate:
    input: 
        genome= GENOME,
        gdict = GENOME + ".dict",
        bam = "BAM/{sample}.markdup.fixed.bam",
        bai = "BAM/{sample}.markdup.fixed.bai",
    output: 
        "BAM/{sample}.markdup.fixed.table"
    params:
        dbsnp = dbSNP,
        dbsnp1k= db1K,
        dbindel = dbINDEL,
    shell:
        "gatk BaseRecalibrator -R {input.genome} -I {input.bam} -O {output} "
        "--known-sites {params.dbsnp} --known-sites {params.dbindel} "
        "--known-sites {params.dbsnp1k}"

# --known-sites b37_1000G_phase1.indels.b37.vcf.gz \
# --known-sites b37_Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
# --known-sites b37_dbsnp_138.b37.vcf.gz -O ERR250971.BQSR.table

rule applyBQSR:
    input: 
        genome=GENOME,
        gdict = GENOME + ".dict",
        bam = "BAM/{sample}.markdup.fixed.bam",
        bai = "BAM/{sample}.markdup.fixed.bai",
        table = "BAM/{sample}.markdup.fixed.table",
    output: 
        protected("BAM/{sample}.markdup.fixed.BQSR.bam")
    params:
        dbsnp=dbSNP,
        dbindel=dbINDEL
    shell:
        """
        gatk ApplyBQSR -R {input.genome} -I {input.bam} --bqsr-recal-file {input.table} -O {output} --create-output-bam-index true
        # samtools index {output}
        """

rule call:
    input:
        bam="BAM/{sample}.markdup.fixed.BQSR.bam",
        genome=GENOME,
    output: protected("VCF/{sample}.call.raw.vcf.gz")
    shell:
        "gatk HaplotypeCaller -R {input.genome} -I {input.bam} -O {output}"