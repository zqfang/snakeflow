import os
############### Globals ########################

#configfile: "config.yaml"
workdir: config['GATK']['WORKSPACE']

GENOME = config['GENOME']
dbSNP = config['dbSNP']
BAM_DIR = config['BAM_DIR']
STRAINS = sorted(config['STRAINS'])


#CHROMOSOME = [ str(c) for c in range(1,20)] + ["X", "Y", "MT"]
CHROMOSOME = ['1'] + [ str(c) for c in range(10,20)] + [ str(c) for c in range(2,10)]+ ["X", "Y"]
# OUTPUT
VCF_VQSR = expand("VCFs/combined.chr{i}.vqsr.vcf.gz", i=CHROMOSOME)
VCF_HFILTER_PASS = expand("VCFs/combined.chr{i}.hardfilter.vcf.gz", i=CHROMOSOME)
VEP_ANNO = expand("VEP/combined.chr{i}.hardfilter.pass.vep.txt.gz", i=CHROMOSOME)
VCF_RAW = expand("VCFs/combined.chr{i}.raw.vcf.gz", i=CHROMOSOME)

SNPDB = expand("SNPs/chr{i}.txt", i=CHROMOSOME)

############## Rules ##########################
rule all:
    input: VCF_HFILTER_PASS, SNPDB, #VEP_ANNO#VCF_VQSR

# include: "rules/align.recal.smk"

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
        dbsnp = dbSNP,
        genome = GENOME,
        genomedict = GENOME.replace(".fa", ".dict"),
        #bam = "BAM/{sample}.marked.fixed.BQSR.bam", # if read from align.recal.pipeline
        bam = os.path.join(BAM_DIR, "{sample}/output.GATKrealigned.Recal.bam")
    output: 
        # save for future use when more strains are added
        gvcf= protected("GVCF/{sample}.chr{i}.raw.g.vcf.gz"), 
        gvcfi= protected("GVCF/{sample}.chr{i}.raw.g.vcf.gz.tbi") 
    threads: 2 # haplotypecaller
    log: "logs/{sample}.chr{i}.haplotypecaller.log"
    params:
        java_ops="-Xmx16G -Djava.io.tmpdir=%s"%config['GATK']['TMPDIR'], # increase memory if stack overflow
    shell:
        # --dbsnp <- annotation ID column to tell whether is known or not
        "gatk --java-options '{params.java_ops}' HaplotypeCaller "
        "-ERC GVCF "
        "--native-pair-hmm-threads {threads} "
        "--dbsnp {input.dbsnp} "
        "-L {wildcards.i} "
        "-R {input.genome} "
        "-I {input.bam} "
        "-O {output.gvcf} 2> {log} "

# CombinedGVCFs, has to be an interval 
## parallel computing for joint_calling       
rule concatGVCFs:
    input:
        genome=GENOME,
        chrom="GVCF/chr{i}.combine.tmp",
        gvcf=["GVCF/%s.chr{i}.raw.g.vcf.gz"%s for s in STRAINS],
        gvcfi=["GVCF/%s.chr{i}.raw.g.vcf.gz.tbi"%s for s in STRAINS]
    output:
        gvcf=temp("GVCF/combined.chr{i}.g.vcf.gz"),
        gvcfi=temp("GVCF/combined.chr{i}.g.vcf.gz.tbi")
    params:
        gvcf=" --variant ".join(["GVCF/%s.chr{i}.raw.g.vcf.gz"%s for s in STRAINS]),
        java_ops="-Xmx8G -Djava.io.tmpdir=%s"%config['GATK']['TMPDIR']
    log: "logs/chr{i}.combineGVCFs.log"
    shell:
        "gatk --java-options '{params.java_ops}' "
        "CombineGVCFs -L {wildcards.i} -R {input.genome} "
        "--variant {params.gvcf} -O {output.gvcf} 2> {log}"

# NOTE：GATK 相比于samtools 更容易找到杂合突变
rule jointCalling:
    input: 
        gvcf="GVCF/combined.chr{i}.g.vcf.gz",
        gvcfi="GVCF/combined.chr{i}.g.vcf.gz.tbi",
        genome=GENOME
    output: 
        vcf=protected("VCFs/combined.chr{i}.raw.vcf.gz"),
        vcfi=protected("VCFs/combined.chr{i}.raw.vcf.gz.tbi"),
    params:
        java_ops= "-Xmx16G -Djava.io.tmpdir=%s"%config['GATK']['TMPDIR']
    log: "logs/chr{i}.GenotypeGVCFs.log"
    shell:
        "gatk --java-options '{params.java_ops}' "
        "GenotypeGVCFs -R {input.genome} -V {input.gvcf} -O {output.vcf} 2> {log}"


################# VQSR ###############################
# for human, first choice is VQSR
# for non-human, better to do is hardfilering
rule variantRecalSNPs:
    input: 
        vcf="VCFs/combined.chr{i}.raw.vcf.gz",
        vcfi="VCFs/combined.chr{i}.raw.vcf.gz.tbi",
        genome=GENOME,
    output: 
        recal="VCFs/combined.chr{i}.snp.Recal",
        tranches= "VCFs/combined.chr{i}.snp.tranches",
        rscript= "VCFs/combined.chr{i}.snp.R", 
    params:
        dbsnp="snps,known=true,training=true,truth=true,prior=2.0:"+config['GATK']['dbSNP'],  
    log: "logs/combined.chr{i}.snp.Recal.log"   
    shell:
        "gatk VariantRecalibrator -mode SNP "
        "-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP "
        "--resource {params.dbsnp} "
        "-R {input.genome} -V {input.vcf} -O {output.recal} "
        "--tranches-file {output.tranches} --rscript-file {output.rscript} 2> {log}"

rule applyVQSR4SNPs:
    input:
        genome=GENOME,
        vcf="VCFs/combined.chr{i}.raw.vcf.gz",
        vcfi="VCFs/combined.chr{i}.raw.vcf.gz.tbi",
        recal="VCFs/combined.chr{i}.snp.Recal",
        tranches= "VCFs/combined.chr{i}.snp.tranches",
        rscript= "VCFs/combined.chr{i}.snp.R",        
    output:
        vcf=temp("VCFs/combined.chr{i}.snp.VQSR.vcf.gz"),
        vcfi=temp("VCFs/combined.chr{i}.snp.VQSR.vcf.gz.tbi"),
    log: "logs/combined.chr{i}.snp.VQSR.log"
    shell:
        "gatk ApplyVQSR -mode SNP "
        "-R {input.genome} -V {input.vcf} -O {output} "
        "--truth-sensitivity-filter-level 99.5 "
        "--tranches-file {input.tranches} "
        "--rscript-file {input.rscript} "
        "--recal-file {input.recal} 2> {log} "

rule variantRecalINDELs:
    input: 
        vcf="VCFs/combined.chr{i}.raw.vcf.gz",
        vcfi="VCFs/combined.chr{i}.raw.vcf.gz.tbi",
        genome=GENOME,
    output: 
        recal ="VCFs/combined.chr{i}.indel.Recal",
        tranches = "VCFs/combined.chr{i}.indel.tranches",
        rscript = "VCFs/combined.chr{i}.indel.plots.R",
    params:
        #format: "{},known={},training={},truth={},prior={}:{}"
        dbindel="indels,known=true,training=true,truth=true,prior=12.0:"+config['GATK']['dbINDEL'],
        dbsnp="snps,known=true,training=true,truth=true,prior=2.0:"+config['GATK']['dbSNP'], 
    log: "logs/combined.chr{i}.indel.Recal.log"
    shell:
        "gatk VariantRecalibrator -mode INDEL "
        "--max-gaussians 4"
        "-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum "
        "--resource {params.dbindel} "
        "--resource {params.dbsnp} "
        "-R {input.genome} -V {input.vcf} -O {output.recal} "
        "--tranches-file {output.tranches} --rscript-file {output.rscript} 2> {log}"

rule applyVQSR4INDEL:
    input:
        genome = GENOME,
        vcf = "VCFs/combined.chr{i}.raw.vcf.gz",
        vcfi= "VCFs/combined.chr{i}.raw.vcf.gz.tbi",
        recal = "VCFs/combined.chr{i}.indel.Recal",
        tranches = "VCFs/combined.chr{i}.indel.tranches",
        rscript = "VCFs/combined.chr{i}.indel.plots.R", 
    output:
        vcf=temp("VCFs/combined.chr{i}.indel.VQSR.vcf.gz"),
        vcfi=temp("VCFs/combined.chr{i}.indel.VQSR.vcf.gz.tbi"),
    log: "logs/combined.chr{i}.indel.VQSR.log"
    shell:
        "gatk ApplyVQSR -mode INDEL "
        "-R {input.genome} -V {input.vcf} -O {output.vcf} "
        "--tranches-file {input.tranches} "
        "--rscript-file {input.rscript} "
        "--recal-file {input.recal} "
        "--truth-sensitivity-filter-level 99.0 2> {log} "

rule mergeVQSRVCFs:
    input:
        snp = "VCFs/combined.{chrom}.snp.VQSR.vcf.gz",
        snpi = "VCFs/combined.{chrom}.snp.VQSR.vcf.gz.tbi",
        indel = "VCFs/combined.{chrom}.indel.VQSR.vcf.gz",
        indeli = "VCFs/combined.{chrom}.indel.VQSR.vcf.gz.tbi",
    output:
        vcf=protected("VCFs/combined.{chrom}.vqsr.vcf.gz"),
        vcfi=protected("VCFs/combined.{chrom}.vqsr.vcf.gz.tbi")
    shell:
        "gatk MergeVcfs -I {input.snp} -I {input.indel} "
        "-O {output.vcf} 2>/dev/null"


################# Hard filering ######################
rule selectSNPs:
    input: 
        vcf="VCFs/combined.{chrom}.raw.vcf.gz",
        vcfi="VCFs/combined.{chrom}.raw.vcf.gz.tbi",
    output: 
        vcf=temp("VCFs/combined.{chrom}.snp.vcf.gz"), # snp only, before filter
        vcfi=temp("VCFs/combined.{chrom}.snp.vcf.gz.tbi"), # 
    shell:
        "gatk SelectVariants -select-type SNP " 
        "-V {input.vcf} -O {output.vcf} 2>/dev/null"

rule hardFilterSNPs:
    input: 
        vcf="VCFs/combined.{chrom}.snp.vcf.gz",
        vcfi="VCFs/combined.{chrom}.snp.vcf.gz.tbi"
    output: 
        vcf=temp("VCFs/combined.{chrom}.snp.filter.vcf.gz"),
        vcfi=temp("VCFs/combined.{chrom}.snp.filter.vcf.gz.tbi"),
    log: "logs/combined.{chrom}.snp.filter.log"
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
        vcf="VCFs/combined.{chrom}.raw.vcf.gz",
        vcfi="VCFs/combined.{chrom}.raw.vcf.gz.tbi",
    output: 
        vcf=temp("VCFs/combined.{chrom}.indel.vcf.gz"),
        vcfi=temp("VCFs/combined.{chrom}.indel.vcf.gz.tbi"),
    shell:
        "gatk SelectVariants -select-type INDEL " 
        "-V {input.vcf} -O {output.vcf} 2>/dev/null "   

rule hardFilterINDELs:
    input: 
        vcf="VCFs/combined.{chrom}.indel.vcf.gz",
        vcfi="VCFs/combined.{chrom}.indel.vcf.gz.tbi",
    output: 
        vcf=temp("VCFs/combined.{chrom}.indel.filter.vcf.gz"),
        vcfi=temp("VCFs/combined.{chrom}.indel.filter.vcf.gz.tbi"),
    log: "logs/combined.{chrom}.indel.filter.log"
    shell:
        "gatk VariantFiltration " 
        "-filter 'QD < 2.0' --filter-name 'QD2' " 
        "-filter 'QUAL < 30.0' --filter-name 'QUAL30' " 
        "-filter 'FS > 200.0' --filter-name 'FS200' " 
        "-filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20' " 
        "-V {input.vcf} -O {output.vcf} 2> {log}"     

rule mergeHardVCFs:
    input:
        snp= "VCFs/combined.{chrom}.snp.filter.vcf.gz",
        snpi= "VCFs/combined.{chrom}.snp.filter.vcf.gz.tbi",
        indel= "VCFs/combined.{chrom}.indel.filter.vcf.gz",
        indeli= "VCFs/combined.{chrom}.indel.filter.vcf.gz.tbi",
    output: 
        vcf=protected("VCFs/combined.{chrom}.hardfilter.vcf.gz"),
        vcfi=protected("VCFs/combined.{chrom}.hardfilter.vcf.gz.tbi")
    shell:
        "gatk MergeVcfs -I {input.snp} -I {input.indel} "
        "-O {output.vcf} 2>/dev/null "


rule strainOrder:
    output: "strain.order.snpdb.txt"
    run:
        with open(output[0], 'w') as s:
            s.write("\n".join(STRAINS) +"\n")
        
rule snp2NIEHS:
    input:  
        strain = "strain.order.snpdb.txt",
        vcf = "VCFs/combined.chr{i}.hardfilter.vcf.gz",
    output: 
        protected("SNPs/chr{i}.txt")
    params:
        qual=config['GATK']['qual'], 
        het = config['GATK']['phred_likelihood_diff'],
        ad = config['GATK']['allele_depth'],
        ratio = config['GATK']['allele_mindepth_ratio'],
        mq = config['GATK']['mapping_quality'],
        sb = config['GATK']['strand_bias_pvalue'], 
        BIN = config['HBCGM']['BIN']# path to haplomap binary
    log: "logs/combined.chr{i}.snp2niehs.log"
    shell:
        # MARK: bcftools view -v snps won't work for GATK VCFs, 
        # haplomap niehs will handle indels
        "bcftools view -v snps {input.vcf} | "
        "{params.BIN}/haplomap niehs -o {output} -a {params.ad} -r {params.ratio} "
        "-q {params.qual} -p {params.het} -m {params.mq} -b {params.sb} "
        "-s {input.strain} -v > {log}"


## only do this for VEP input
# rule selectPASS:
#     input: 
#         vcf="VCFs/combined.{chrom}.hardfilter.vcf.gz",
#         vcfi="VCFs/combined.{chrom}.hardfilter.vcf.gz.tbi",
#         genome=GENOME
#     output: 
#         temp("VCFs/combined.{chrom}.hardfilter.pass.vcf.gz"),
#         temp("VCFs/combined.{chrom}.hardfilter.pass.vcf.gz.tbi")
#     shell:
#         "gatk SelectVariants -R {input.genome} -V {input.vcf} -O {output[0]} "
#         " -select 'vc.isNotFiltered()' 2>/dev/null"

rule variantEeffectPrediction:
    """emsemble-vep"""
    input: 
        vcf="VCFs/combined.{chrom}.hardfilter.vcf.gz",
        reference=GENOME,
    output: "VEP/combined.{chrom}.hardfilter.pass.vep.txt.gz"
    params:
        #genome_build = " -a GRCm38 --species mus_musculus ",
        genome_build = config['VEP']['GENOME_BUILD'],
        VEPBIN = config['VEP']['BIN'],
        extra=" --dir_cache "  + config['VEP']['CACHE_DIR']
    threads: 1
    shell:
        ## emsemble-vep
        # https://github.com/Ensembl/ensembl-vep
        "bcftools view -f .,PASS {input.vcf} | "
        "{params.VEPBIN}/vep --fasta {input.reference} {params.genome_build} "
        "--format vcf --fork {threads} --hgvs --force_overwrite "
        "--uniprot --domains --symbol --regulatory --distance 1000 --biotype "
        "--gene_phenotype MGI --check_existing  --pubmed --numbers "
        "--offline --cache --variant_class "
        "--gencode_basic --no_intergenic --individual all "
        "-o {output} --tab --compress_output gzip"