import os
from snakemake.shell import shell
############### Globals ########################
WORKSPACE = "/home/fangzq/data/20200327"
workdir: WORKSPACE

GENOME="/home/fangzq/genome/mouse/GRCm38_68.fa"
dbSNP="/home/fangzq/genome/mouse/mgp.v5.merged.snps_all.dbSNP142.sorted.vcf"
dbINDEL="/home/fangzq/genome/mouse/mgp.v5.merged.indels.dbSNP142.normed.vcf"
STRAINS_FILE = "/data/bases/fangzq/strain"
BAM_DIR = "/data/bases/fangzq/strains"
TMPDIR = "/home/fangzq/TMPDATA"

#CHROMSOME = [ str(c) for c in range(1,20)] + ["X", "Y", "MT"]
CHROMSOME = ['1'] + [ str(c) for c in range(10,20)] + [ str(c) for c in range(2,10)]+ ["MT", "X", "Y"]

# with open(STRAINS_FILE, 'r') as s:
#     STRAINS = s.read().strip().split()
# STRAINS_ALL = \
# "B10 BTBR BUB CEJ DBA1J FVB KK NON NUJ NZB NZW RFJ RHJ RIIIS SJL " +\ 
# "129P2 129S1 129S5 A_J AKR B_C C3H C57BL10J C57BL6NJ C57BRcd C57LJ " +\
# "C58 CBA DBA ILNJ LGJ LPJ MAMy MRL NOD NZO PJ PLJ SEA SMJ ST SWR " +\
# "CAST MOLF PWD PWK SPRET WSB"  # <- wild derived
STRAINS = "129S1 129S5 B10 B_C C3H C57BL10J C57BRcd C57LJ "+\
          "C58 CBA CEJ DBA ILNJ LGJ LPJ MAMy NOD NON NUJ NZO PJ PLJ RHJ RIIIS SEA SJL SMJ SWR " +\
          "CAST MOLF PWD PWK SPRET WSB"
STRAINS = STRAINS.split(" ")
# OUTPUT
VCF_VQSR = expand("VCFs/combined.chr{i}.VQSR.vcf.gz", i=CHROMSOME)
VCF_HFILTER = expand("VCFs/combined.chr{i}.hardfilter.vcf.gz", i=CHROMSOME)
VCF_HFILTER_PASS = expand("VCFs/combined.chr{i}.hardfilter.pass.vcf.gz", i=CHROMSOME)
VCF_RAW = expand("VCFs/combined.chr{i}.raw.vcf", i=CHROMSOME)
# SNPs = expand("SNPs/combined.chr{i}.txt", i=CHROMSOME)


############## Rules ##########################
rule all:
    input: VCF_RAW, VCF_HFILTER, VCF_VQSR

# include: "rules/gatk.getbam.smk"

rule sampleCalling:
    input:
        dbSNP = dbSNP,
        genome = GENOME,
        genomedict = GENOME.replace(".fa", ".dict"),
        #bam = "BAM/{sample}.marked.fixed.BQSR.bam",
        bam = os.path.join(BAM_DIR, "{sample}/output.GATKrealigned.Recal.bam")
    output: 
        gvcf= temp("GVCF/{sample}.raw.g.vcf"),
        gvcfi= temp("GVCF/{sample}.raw.g.vcf.idx")
    threads: 12
    log: "logs/{sample}.haplotypecaller.log"
    params:
        java_ops="-Xmx32G -Djava.io.tmpdir=%s"%TMPDIR
    shell:
        "gatk --java-options '{params.java_ops}' HaplotypeCaller "
        "-ERC GVCF "
        "--native-pair-hmm-threads {threads} "
        "--dbsnp {input.dbSNP} "
        "-R {input.genome} "
        "-I {input.bam} "
        "-O {output.gvcf} 2> {log} "

# rule gatherVCFs:
#     input:
#         gvcf=expand("/data/bases/fangzq/strains/GATK_TMP/{strain}.chr{i}.raw.g.vcf", i=CHROMSOME),
#         gvcfi=expand("/data/bases/fangzq/strains/GATK_TMP/{strain}.chr{i}.raw.g.vcf.idx", i=CHROMSOME)
#     output:
#         gvcf="/data/bases/fangzq/GVCF/{strain}.raw.g.vcf",
#         gvcfi="/data/bases/fangzq/GVCF/{strain}.raw.g.vcf.idx"
#     log: "/data/bases/fangzq/strains/{strain}.gatherVCFs.log"
#     run:
#        g = " -I ".join(input.gvcf)
#        shell("gatk GatherVcfs -I {gvcf} -O {output.gvcf}".format(gvcf=g))

# 1 CombinedGVCFs, has to be an interval, so seperate each chromosome
rule combineGVCFs:
    input:
        genome=GENOME,
        gvcf=expand("GVCF/{sample}.raw.g.vcf", sample=STRAINS),
        gvcfi=expand("GVCF/{sample}.raw.g.vcf.idx", sample=STRAINS)
    output:
        temp(expand("GVCF/combined.chr{i}.g.vcf", i=CHROMSOME)),
        temp(expand("GVCF/combined.chr{i}.g.vcf.idx", i=CHROMSOME))
    params:
        chrom=CHROMSOME
    log: "logs/combineGVCFs.log"
    run:
        variant = " --variant ".join(input.gvcf) 
        for i in params['chrom']:
            shell("gatk CombineGVCFs -L {chr} -R {input.genome} --variant {var} -O combined.chr${chr}.g.vcf 2>> {log}".format(chr=i, var=variant))

# NOTE：GATK 相比于samtools 更容易找到杂合突变
rule jointCalling:
    input: 
        gvcf="GVCF/combined.chr{i}.g.vcf",
        gvcfi="GVCF/combined.chr{i}.g.vcf.idx",
        genome=GENOME
    output: 
        vcf=protected("VCFs/combined.chr{i}.raw.vcf"),
        vcfi=protected("VCFs/combined.chr{i}.raw.vcf.idx"),
    params:
        tmpdir=TMPDIR,
        java_ops= "-Xmx32G -Djava.io.tmpdir=%s"%TMPDIR
    log: "logs/chr{i}.GenotypeGVCFs.log"
    shell:
        "gatk --java-options '{params.java_ops}' "
        "GenotypeGVCFs -R {input.genome} -V {input.gvcf} -O {output.vcf} 2> {log}"


################# VQSR ###############################
# for human, first choice is VQSR
# for non-human, better to do is hardfilering
rule variantRecalSNPs:
    input: 
        vcf="VCFs/combined.chr{i}.raw.vcf",
        genome=GENOME,
    output: 
        recal="VCFs/combined.chr{i}.snp.Recal",
        tranches= "VCFs/combined.chr{i}.snp.tranches",
        rscript= "VCFs/combined.chr{i}.snp.R", 
    params:
        dbsnp=f"dbsnp,known=true,training=true,truth=true,prior=2.0:{dbSNP}",     
    shell:
        "gatk VariantRecalibrator -mode SNP "
        "-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP "
        "--resource {params.dbsnp} "
        "-R {input.genome} -V {input.vcf} -O {output.recal} "
        "--tranches-file {output.tranches} --rscript-file {output.rscript}"

rule applyVQSR4SNPs:
    input:
        genome=GENOME,
        vcf="VCFs/combined.chr{i}.raw.vcf",
        recal="VCFs/combined.chr{i}.snp.Recal",
        tranches= "VCFs/combined.chr{i}.snp.tranches",
        rscript= "VCFs/combined.chr{i}.snp.R",        
    output:
        temp("VCFs/combined.chr{i}.snp.VQSR.vcf"),
    shell:
        "gatk ApplyVQSR -mode SNP "
        "-R {input.genome} -V {input.vcf} -O {output} "
        "--truth-sensitivity-filter-level 99.5 "
        "--tranches-file {input.tranches} "
        "--rscript-file {input.rscript} "
        "--recal-file {input.recal} "

rule variantRecalINDELs:
    input: 
        vcf="VCFs/combined.chr{i}.raw.vcf",
        genome=GENOME,
    output: 
        recal ="VCFs/combined.chr{i}.indel.Recal",
        tranches = "VCFs/combined.chr{i}.indel.tranches",
        rscript = "VCFs/combined.chr{i}.indel.plots.R",
    params:
        #format: "{},known={},training={},truth={},prior={}:{}"
        dbindel=f"snps_all,known=true,training=true,truth=true,prior=12.0:{dbINDEL}",
        dbsnp=f"indels,known=true,training=true,truth=true,prior=2.0:{dbSNP}", 
    shell:
        "gatk VariantRecalibrator -mode INDEL "
        "--max-gaussians 4"
        "-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum "
        "--resource {params.dbindel} "
        "--resource {params.dbsnp} "
        "-R {input.genome} -V {input.vcf} -O {output.recal} "
        "--tranches-file {output.tranches} --rscript-file {output.rscript}"

rule applyVQSR4INDEL:
    input:
        genome = GENOME,
        vcf = "VCFs/combined.chr{i}.raw.vcf",
        recal ="VCFs/combined.chr{i}.indel.Recal",
        tranches = "VCFs/combined.chr{i}.indel.tranches",
        rscript = "VCFs/combined.chr{i}.indel.plots.R", 
    output:
        temp("VCFs/combined.chr{i}.indel.VQSR.vcf"),
    shell:
        "gatk ApplyVQSR -mode INDEL "
        "-R {input.genome} -V {input.vcf} -O {output} "
        "--tranches-file {input.tranches} "
        "--rscript-file {input.rscript} "
        "--recal-file {input.recal} "
        "--truth-sensitivity-filter-level 99.0"

rule mergeVQSRVCFs:
    input:
        snp = "VCFs/combined.{chrom}.snp.VQSR.vcf",
        indel = "VCFs/combined.{chrom}.indel.VQSR.vcf",
    output:
        protected("VCFs/combined.{chrom}.VQSR.vcf")
    shell:
        "gatk MergeVcfs -I {input.snp} -I {input.indel} "
        "-O {output}"

rule compressVQSRVCFs:
    input: "VCFs/combined.{chrom}.VQSR.vcf"
    output: protected("VCFs/combined.{chrom}.VQSR.vcf.gz")
    shell:
        """bgzip -f {input} 
           tabix -p vcf {output}
        """

################# Hard filering ######################
rule selectSNPs:
    input: "VCFs/combined.{chrom}.raw.vcf"
    output: temp("VCFs/combined.{chrom}.snp.vcf")
    shell:
        "gatk SelectVariants -select-type SNP " 
        "-V {input} -O {output} "

rule hardFilterSNPs:
    input: "VCFs/combined.{chrom}.snp.vcf"
    output: temp("VCFs/combined.{chrom}.snp.filter.vcf")
    params:
        snp_hard= "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
    shell:
        "gatk VariantFiltration -select-type SNP "
        "--filter-expression '{params.snp_hard}' " 
        "--filter-name 'SNP_FILTER' -V {input} -O {output} "    

rule selectINDELs:
    input: "VCFs/combined.{chrom}.raw.vcf"
    output: temp("VCFs/combined.{chrom}.indel.vcf"),
    shell:
        "gatk SelectVariants -select-type INDEL " 
        "-V {input} -O {output} "   

rule hardFilterINDELs:
    input: "VCFs/combined.{chrom}.indel.vcf"
    output: temp("VCFs/combined.{chrom}.indel.filter.vcf"),
    params:
        indel_hard="QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -20.0"
    shell:
        "gatk VariantFiltration -select-type SNP " 
        "--filter-expression '{params.indel_hard}' " 
        "--filter-name 'INDEL_FILTER' -V {input} -O {output} "     

rule mergeHardVCFs:
    input:
        snp= "VCFs/combined.{chrom}.snp.filter.vcf",
        indel= "VCFs/combined.{chrom}.indel.filter.vcf",
    output: protected("VCFs/combined.{chrom}.hardfilter.vcf")
    shell:
        "gatk MergeVcfs -I {input.snp} -I {input.indel} "
        "-O {output}"

rule compressHardVCF:
    input: "VCFs/combined.{chrom}.hardfilter.vcf"
    output: protected("VCFs/combined.{chrom}.hardfilter.vcf.gz")
    shell:
        """bgzip -f {input} 
           tabix -p vcf {output}
        """

rule extractPASS:
    input: "VCFs/combined.{chrom}.hardfilter.vcf"
    output: "VCFs/combined.{chrom}.hardfilter.pass.vcf.gz"
    shell:
        "gatk SelectVariants -R {input.genome} -V {input.vcf} -O {output} "
        " -select 'vc.isNotFiltered()'"
