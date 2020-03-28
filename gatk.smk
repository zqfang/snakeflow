import os
from snakemake.shell import shell
############### Globals ########################
WORKSPACE = "/home/fangzq/data/20200327"
workdir: WORKSPACE

GENOME="/home/fangzq/genome/mouse/GRCm38_68.fa"
dbSNP="/home/fangzq/genome/mouse/mgp.v5.merged.snps_all.dbSNP142.sorted.vcf"
STRAINS_FILE = "/data/bases/fangzq/strain"
BAM_DIR = "/data/bases/fangzq/strains"
TMPDIR = "/home/fangzq/TMPDATA"

#CHROMSOME = [ str(c) for c in range(1,20)] + ["X", "Y", "MT"]
CHROMSOME = ['1'] + [ str(c) for c in range(10,20)] + [ str(c) for c in range(2,10)]+ ["MT", "X", "Y"]

with open(STRAINS_FILE, 'r') as s:
    STRAINS = s.read().strip().split()

# STRAINS = "129P2 129S1 129S5 B10 B_C BPL BPN BUB C3H C57BL10J C57BL6NJ C57BRcd C57LJ "+\
#           "C58 CBA CEJ DBA ILNJ LGJ LPJ MAMy NOD NON NUJ NZB NZO PJ PLJ RHJ RIIIS SEA SJL SMJ ST SWR"

# OUTPUT
VCFs = expand("VCFs/combined.chr{i}.vcf", i=CHROMSOME)
SNPs = expand("SNPs/combined.chr{i}.txt", i=CHROMSOME)


############## Rules ##########################
rule all:
    input: VCFs, SNPs

    
rule sample_calling:
    input:
        dbSNP=dbSNP,
        genome=GENOME,
        bam=os.path.join(BAM_DIR, "{strain}/output.GATKrealigned.Recal.bam"),
        bai=os.path.join(BAM_DIR, "{strain}/output.GATKrealigned.Recal.bai"),
    output: 
        gvcf="GVCF/{strain}.raw.g.vcf",
        gvcfi="GVCF/{strain}.raw.g.vcf.idx"
        # gvcf=expand("/data/bases/fangzq/strains/GATK_TMP/{strain}.chr{i}.raw.g.vcf", i=CHROMSOME),
        # gvcfi=expand("/data/bases/fangzq/strains/GATK_TMP/{strain}.chr{i}.raw.g.vcf.idx", i=CHROMSOME)
    threads: 12
    log: "logs/{strain}.haplotypecaller.log"
    params:
        java_ops="-Xmx32G -Djava.io.tmpdir=%s"%TMPDIR,
        chrs=CHROMSOME,
        tmpdir=TMPDIR,
        strain="{strain}"
    shell:
        "gatk --java-options '{params.java_ops}' HaplotypeCaller "
        "-ERC GVCF "
        "--native-pair-hmm-threads {threads} "
        "--dbsnp {input.dbSNP} "
        "-R {input.genome} "
        "-I {input.bam} "
        "-O {output.gvcf} 2> {log} "
    ## split run into chromosomes
    # run:
    #     for ch in params['chrs']:
    #         shell("""gatk HaplotypeCaller  \
    #                 -ERC GVCF --tmp-dir /home/fangzq/TMPDATA \
    #                 --native-pair-hmm-threads 2 \
    #                 --dbsnp {input.dbSNP} \
    #                 -R {input.genome} \
    #                 -I {input.bam} \
    #                 -L {i} \
    #                 -O /data/bases/fangzq/strains/GATK_TMP/{params.strain}.chr{i}.raw.g.vcf 
    #                """.format(i=ch))

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
        gvcf=expand("GVCF/{strain}.raw.g.vcf", strain=STRAINS),
        gvcfi=expand("GVCF/{strain}.raw.g.vcf.idx", strain=STRAINS)
    output:
        expand("GVCF/combined.chr{i}.g.vcf", i=CHROMSOME),
        expand("GVCF/combined.chr{i}.g.vcf.idx", i=CHROMSOME)
    params:
        chrs=CHROMSOME
    log: "logs/combineGVCFs.log"
    run:
        variant = " --variant ".join(input.gvcf) 
        for i in params['chrs']:
            shell("gatk CombineGVCFs -L {chr} -R {input.genome} --variant {var} -O combined.chr${chr}.g.vcf 2>> {log}".format(chr=i, var=variant))

rule joint_calling:
    input: 
        gvcf="GVCF/combined.chr{i}.g.vcf",
        gvcfi="GVCF/combined.chr{i}.g.vcf.idx",
        genome=GENOME
    output: 
        vcf=protected("VCFs/combined.chr{i}.vcf"),
        vcfi=protected("VCFs/combined.chr{i}.vcf.idx"),
    params:
        tmpdir=TMPDIR,
        java_ops= "-Xmx32G -Djava.io.tmpdir=%s"%TMPDIR
    log: "logs/chr{i}.GenotypeGVCFs.log"
    shell:
        "gatk --java-options '{params.java_ops}' "
        "GenotypeGVCFs -R {input.genome} -V {input.gvcf} -O {output.vcf} 2> {log}"
        # """
        # gatk GenotypeGVCFs --tmp-dir {params.tmpdir} -R {input.genome} \
        #     -V {input.gvcf} \
        #     -O {output} 2> {log}
        # """


# g.vcf文件用CombineGVCFs方式或者GenomicsDBImport方式合并成一个文件，前者（比较传统）是一个总的g.vcf文件，后者是一个GenomicsDB（XX.db）
# 2 GenomicsDBImport方法（这里需要注意的是其必须要输入一个区间One or more genomic intervals，所以可以选择分染色体进行）：
#for i in $(seq 1 19) X Y MT;
# do
# gatk GenomicsDBImport $(fori in $(ls *.vcf);do echo "-V $i";done) \
#    -L $i \
#    --genomicsdb-workspace-path DB.chr${i}
#

#gatk --java-options "-Xmx32G -Djava.io.tmpdir=/home/fangzq/TMPDATA" GenotypeGVCFs \
#    -R $GENOME \
#    -V gendb://DB.chr${i} \
#    -O chr${i}.combined.vcf
# done

rule vcf2strains:
    input:  
        vcf = "VCFs/combined.chr{i}.vcf", 
        vcfi = "VCFs/combined.chr{i}.vcf.idx",
    output: 
        temp("VCFs/chr{chrom}.strains.txt")
    shell:
        # NOTE: '\t' is default delim for cut
        "head -n 1000 {input.vcf} | grep '^#CHROM' | "
        "cut -f10-  > {output}"
   
    
rule vcf2niehs:
    input:  
        vcf = "VCFs/chr{chrom}.vcf", 
        vcfi =  "VCFs/chr{chrom}.vcf.idx",
        stains = "VCFs/chr{chrom}.strains.txt"
    output: 
        protected("SNPs/chr{chrom}.txt")
    params:
        outdir= "SNPs",
        chrom="{chrom}",
        qual_samtools=50, 
        heterzygote_cutoff = 20
    script:
        "scripts/vcf2NIEHS.py"
