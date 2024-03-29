
import os,re
############### Globals ########################
#configfile: "config.yaml"
workdir: "/data/bases/fangzq/UDN/Hearing_loss/test"
BUNDLE="/home/fangzq/genome/gatk_bundle/hg38"
GENOME = os.path.join(BUNDLE, "Homo_sapiens_assembly38.fasta")
dbSNP =  os.path.join(BUNDLE, "dbsnp_146.hg38.vcf.gz")
db1K =  os.path.join(BUNDLE, "1000G_phase1.snps.high_confidence.hg38.vcf.gz")
dbINDEL =  os.path.join(BUNDLE,"Mills_and_1000G_gold_standard.indels.hg38.vcf.gz")
VEPBIN = "/home/fangzq/github/ensembl-vep"
SAMPLES = ["RQ82599"]
#CHROMSOME = [ str(c) for c in range(1,20)] + ["X", "Y", "MT"]
CHROMSOME = ['1'] + [ str(c) for c in range(10,20)] + [ str(c) for c in range(2,10)]+ ["X", "Y"]

#outputs
BAMS = expand("BAM/{sample}.sorted.bam", sample=SAMPLES)
BQSR = expand("BAM/{sample}.markdup.fixed.BQSR.bam", sample=SAMPLES)
VCFS = expand("VCFs/{sample}.call.raw.vcf.gz", sample=SAMPLES)
SNPS = expand("VCFs/{sample}.snp.filter.vcf.gz",sample=SAMPLES)
INDELS = expand("VCFs/{sample}.indel.filter.vcf.gz",sample=SAMPLES)
VEP = expand("VEP/{sample}.{snv}.vep.txt.gz",sample=SAMPLES, snv=["snp", "indel"])
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
    input: BAMS, BQSR, VCFS, SNPS, INDELS, VEP

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
        RG="@RG\\tID:{sample}\\tPL:illumina\\tLB:library\\tSM:{sample}"
    shell:
        "bwa mem -t {threads} -R '{params.RG}' {input.genome} "
        "{input.r1} {input.r2} | samtools view -Sb - > {output}"
        
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
    output: 
        vcf = protected("VCFs/{sample}.call.raw.vcf.gz"),
        vcfi = "VCFs/{sample}.call.raw.vcf.gz.tbi"
    shell:
        "gatk HaplotypeCaller -R {input.genome} -I {input.bam} -O {output}"

################# Hard filering ######################
rule selectSNPs:
    input: 
        vcf="VCFs/{sample}.call.raw.vcf.gz",
        vcfi="VCFs/{sample}.call.raw.vcf.gz.tbi",
    output: 
        vcf=temp("VCFs/{sample}.snp.vcf.gz"), # snp only, before filter
        vcfi=temp("VCFs/{sample}.snp.vcf.gz.tbi"), # 
    shell:
        "gatk SelectVariants -select-type SNP " 
        "-V {input.vcf} -O {output.vcf} 2>/dev/null"

rule hardFilterSNPs:
    input: 
        vcf="VCFs/{sample}.snp.vcf.gz",
        vcfi="VCFs/{sample}.snp.vcf.gz.tbi",
    output: 
        vcf="VCFs/{sample}.snp.filter.vcf.gz",
        vcfi="VCFs/{sample}.snp.filter.vcf.gz.tbi",
    log: "logs/{sample}.snp.filter.log"
    shell:
        "gatk VariantFiltration "
        "-filter 'QD < 2.0' --filter-name 'QD2' " 
        "-filter 'QUAL < 30.0' --filter-name 'QUAL30' " 
        "-filter 'SOR > 3.0' --filter-name 'SOR3' " 
        "-filter 'FS > 60.0' --filter-name 'FS60' " 
        "-filter 'MQ < 40.0' --filter-name 'MQ40' " 
        "-filter 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' " 
        "-filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8' " 
        "-V {input.vcf} -O {output.vcf} 2> {log}"

rule selectINDELs:
    input: 
        vcf="VCFs/{sample}.call.raw.vcf.gz",
        vcfi="VCFs/{sample}.call.raw.vcf.gz.tbi",
    output: 
        vcf=temp("VCFs/{sample}.indel.vcf.gz"),
        vcfi=temp("VCFs/{sample}.indel.vcf.gz.tbi"),
    shell:
        "gatk SelectVariants -select-type INDEL " 
        "-V {input.vcf} -O {output.vcf} 2>/dev/null "   

rule hardFilterINDELs:
    input: 
        vcf="VCFs/{sample}.indel.vcf.gz",
        vcfi="VCFs/{sample}.indel.vcf.gz.tbi",
    output: 
        vcf="VCFs/{sample}.indel.filter.vcf.gz",
        vcfi="VCFs/{sample}.indel.filter.vcf.gz.tbi",
    log: "logs/{sample}.indel.filter.log"
    shell:
        "gatk VariantFiltration " 
        "-filter 'QD < 2.0' --filter-name 'QD2' " 
        "-filter 'QUAL < 30.0' --filter-name 'QUAL30' " 
        "-filter 'FS > 200.0' --filter-name 'FS200' " 
        "-filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20' " 
        "-V {input.vcf} -O {output.vcf} 2> {log}"     

# rule mergeHardVCFs:
#     input:
#         snp= "VCFs/{sample}.snp.filter.vcf.gz",
#         snpi= "VCFs/{sample}.snp.filter.vcf.gz.tbi",
#         indel= "VCFs/{sample}.indel.filter.vcf.gz",
#         indeli= "VCFs/{sample}.indel.filter.vcf.gz.tbi",
#     output: 
#         vcf=protected("VCFs/{sample}.hardfilter.vcf.gz"),
#         vcfi=protected("VCFs/{sample}.hardfilter.vcf.gz.tbi")
#     shell:
#         "gatk MergeVcfs -I {input.snp} -I {input.indel} "
#         "-O {output.vcf} 2>/dev/null "



rule VEP4SNPs:
    """emsemble-vep"""
    input: 
        vcf="VCFs/{sample}.snp.filter.vcf.gz",
        reference=GENOME,
    output: "VEP/{sample}.snp.vep.txt.gz"
    params:
        genome_build = " -a GRCh38 --species homo_sapiens ",
        #genome_build = "GRCh38",
        VEPBIN = VEPBIN,
        #extra=" --dir_cache "  + config['VEP']['CACHE_DIR']
    threads: 4
    shell:
        ## emsemble-vep
        # https://github.com/Ensembl/ensembl-vep
        "bcftools view {input.vcf} | {params.VEPBIN}/vep --fasta {input.reference} {params.genome_build} "
        "--format vcf --fork {threads} --force_overwrite "
        "--everything --sift --polyphen "
        "-o {output} --tab --compress_output gzip --database"


rule VEP4Indels:
    """emsemble-vep"""
    input: 
        vcf="VCFs/{sample}.indel.filter.vcf.gz",
        reference=GENOME,
    output: "VEP/{sample}.indel.vep.txt.gz"
    params:
        genome_build = " -a GRCh38 --species homo_sapiens ",
        #genome_build = "GRCh38",
        VEPBIN = VEPBIN,
        # extra=" --dir_cache "  + config['VEP']['CACHE_DIR']
    threads: 4
    shell:
        ## emsemble-vep
        # https://github.com/Ensembl/ensembl-vep
        "bcftools view {input.vcf} | {params.VEPBIN}/vep --fasta {input.reference} {params.genome_build} "
        "--format vcf --fork {threads} --force_overwrite "
        "--everything --sift --polyphen "
        "-o {output} --tab --compress_output gzip --database"
