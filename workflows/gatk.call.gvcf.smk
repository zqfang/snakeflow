import os
############### Globals ########################

#configfile: "config.yaml"
workdir: config['GATK']['WORKSPACE']

GENOME = config['GENOME']
dbSNP = config['dbSNP']
dbINDEL = config['dbINDEL']
BAM_DIR = config['BAM_DIR']
TMPDIR = config['TMPDIR']
STRAINS = sorted(config['STRAINS'])
HAPLOMAP = config['HBCGM']['BIN']# path to haplomap binary
VEP = config['VEP']

#CHROMOSOME = [ str(c) for c in range(1,20)] + ["X", "Y", "MT"]
CHROMOSOME = ['1'] + [ str(c) for c in range(10,20)] + [ str(c) for c in range(2,10)]+ ["X", "Y"]
# OUTPUT
VCF_VQSR = expand("VCFs/combined.chr{i}.VQSR.vcf.gz", i=CHROMOSOME)
VCF_HFILTER_PASS = expand("VCFs/combined.chr{i}.hardfilter.pass.vcf.gz", i=CHROMOSOME)
VEP_ANNO = expand("VEP/combined.chr{i}.hardfilter.pass.vep.txt.gz", i=CHROMOSOME)
VCF_RAW = expand("VCFs/combined.chr{i}.raw.vcf.gz", i=CHROMOSOME)
GVCF = expand("GVCF/{sample}.chr{i}.raw.g.vcf.gz", sample=STRAINS, i=CHROMOSOME)
SNPDB = expand("SNPs/chr{i}.txt", i=CHROMOSOME)

############## Rules ##########################
rule all:
    input: GVCF

## temp output for combineGVCFs.
rule chroms: 
    output: expand("GVCF/chr{i}.combine.tmp", i = CHROMOSOME)
    run:
        for out in output:
            shell("touch %s"%out)

## split into chromosomes to speedup
## because Haplotypecaller is not a multi-threading module
rule singleCalling:
    input:
        chrom="GVCF/chr{i}.combine.tmp",
        dbSNP = dbSNP,
        genome = GENOME,
        genomedict = GENOME.replace(".fa", ".dict"),
        #bam = "BAM/{sample}.marked.fixed.BQSR.bam", # if read from align.recal.pipeline
        bam = os.path.join(BAM_DIR, "{sample}/output.GATKrealigned.Recal.bam")
    output: 
        # save for future use when more strains are added
        gvcf= protected("GVCF/{sample}.chr{i}.raw.g.vcf.gz"), 
        gvcfi= protected("GVCF/{sample}.chr{i}.raw.g.vcf.gz.tbi") 
    threads: 4 # haplotypecaller
    log: "logs/{sample}.chr{i}.haplotypecaller.log"
    params:
        java_ops="-Xmx16G -Djava.io.tmpdir=%s"%TMPDIR, # increase memory if stack overflow
    shell:
        # --dbsnp <- annotation ID column to tell whether is known or not
        "gatk --java-options '{params.java_ops}' HaplotypeCaller "
        "-ERC GVCF "
        "--native-pair-hmm-threads {threads} "
        "--dbsnp {input.dbSNP} "
        "-L {wildcards.i} "
        "-R {input.genome} "
        "-I {input.bam} "
        "-O {output.gvcf} 2> {log} "
