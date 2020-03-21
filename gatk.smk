from snakemake.shell import shell

GENOME="/home/fangzq/genome/mouse/GRCm38_68.fa"
dbSNP="/home/fangzq/genome/mouse/mgp.v5.merged.snps_all.dbSNP142.sorted.vcf"
STRAINS = "129P2 129S1 129S5 A_J AKR B10 B_C BPL BPN BTBR BUB C3H C57BL10J C57BL6NJ C57BRcd C57LJ C58 CBA CEJ DBA1J DBA FVB ILNJ KK LGJ LPJ MAMy NOD NON NOR NOR NUJ NZB NZO NZW PJ PLJ RBF RFJ RHJ RIIIS SEA SJL SMJ ST SWR TALLYHO"
STRAINS = STRAINS.split(" ")
TMPDIR = "/home/fangzq/TMPDATA"
#CHROMSOME = [ str(c) for c in range(1,20)] + ["X", "Y", "MT"]
CHROMSOME = ['1'] + [ str(c) for c in range(10,20)] + [ str(c) for c in range(2,10)]+ ["MT", "X", "Y"]
OUTPUT = expand("combined.chr{i}.vcf", i=CHROMSOME)

rule all:
    input: OUTPUT

    
rule sample_calling:
    input:
        dbSNP=dbSNP,
        genome=GENOME,
        bam="/data/bases/fangzq/strains/{strain}/output.GATKrealigned.Recal.bam",
        bai="/data/bases/fangzq/strains/{strain}/output.GATKrealigned.Recal.bai",
    output: 
        gvcf="/data/bases/fangzq/GVCF/{strain}.raw.g.vcf",
        gvcfi="/data/bases/fangzq/GVCF/{strain}.raw.g.vcf.idx"
        # gvcf=expand("/data/bases/fangzq/strains/GATK_TMP/{strain}.chr{i}.raw.g.vcf", i=CHROMSOME),
        # gvcfi=expand("/data/bases/fangzq/strains/GATK_TMP/{strain}.chr{i}.raw.g.vcf.idx", i=CHROMSOME)
    threads: 12
    log: "/data/bases/fangzq/strains/{strain}.haplotypecaller.log"
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
        gvcf=expand("/data/bases/fangzq/GVCF/{strain}.raw.g.vcf", strain=STRAINS),
        gvcfi=expand("/data/bases/fangzq/GVCF/{strain}.raw.g.vcf.idx", strain=STRAINS)
    output:
        expand("combined.chr{i}.g.vcf", i=CHROMSOME),
        expand("combined.chr{i}.g.vcf.idx", i=CHROMSOME)
    params:
        chrs=CHROMSOME
    log: "/data/bases/fangzq/strains/combineGVCFs.log"
    run:
        variant = " --variant ".join(input.gvcf) 
        for i in params['chrs']:
            shell("gatk CombineGVCFs -L {chr} -R {input.genome} --variant {var} -O combined.chr${chr}.g.vcf >> {log}".format(chr=i, var=variant))

rule joint_calling:
    input: 
        gvcf="combined.chr{i}.g.vcf",
        gvcfi="combined.chr{i}.g.vcf.idx",
        genome=GENOME
    output: "combined.chr{i}.vcf",
    params:
        tmpdir=TMPDIR,
        java_ops= "-Xmx32G -Djava.io.tmpdir=%s"%TMPDIR
    log: "/data/bases/fangzq/strains/chr{i}.GenotypeGVCFs.log"
    shell:
        "gatk --java-options '{params.java_ops}' "
        "GenotypeGVCFs -R {input.genome} -V {input.gvcf} -O {output} 2> {log}"
        # """
        # gatk GenotypeGVCFs \
        #     --tmp-dir {params.tmpdir} \
        #     -R {input.genome} \
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