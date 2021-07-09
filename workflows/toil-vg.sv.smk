
## this is pipeline for genotyping SVs given known VCFs.
## we use toil-vg insteal of vg to simply the commands

## see the tutorial here: 
## https://github.com/vgteam/toil-vg/wiki/Genotyping-Structural-Variants  
## https://github.com/vgteam/sv-genotyping-paper
# however, this tutorial is outdated.
# I've updated the commands and run toil-vg locally without container

# require package:
# 
# R package ->  BiocManager::install('jmonlong/sveval') 
# conda install vg=1.33.0 rtg-tools
# pip install toil==3.24.0
# pip install toil-vg==1.6.0
# apt-get install jq

## note: change toil-vg/src/vcf_eval.py "min.cov" -> "min.ol" will run sveval sucessfully.
import os, glob

## Commands
WKDIR = "/data/bases/fangzq/20200815_SV/PacBio_project"
workdir: WKDIR
THREADS = 48
GENOME = "/home/fangzq/genome/mouse/GRCm38_68.fa"
VCF_5_STRAIN = "/data/bases/fangzq/20200815_SV/PacBio_project/VCFs/merged_gt_SURVIVOR_1kbpdist.filter.vcf"
VCF_DELETION =  "merged_deletion.vcf"# "BTBR.deletion.vcf"
SAMPLES =  ['129P2', '129S1', '129S5', 'AKR', 'A_J', 'B10', 
        'BPL', 'BPN', 'BTBR', 'BUB', 'BALB', 'C3H', 'C57BL10J',
        'C57BL6NJ', 'C57BRcd', 'C57LJ', 'C58', 'CBA', 'CEJ', 
        'DBA', 'DBA1J', 'FVB', 'ILNJ', 'KK', 'LGJ', 'LPJ', 
        'MAMy', 'MRL','NOD', 'NON', 'NOR', 'NUJ', 'NZB', 'NZO', 'NZW', 
        'PJ', 'PLJ', 'RFJ', 'RHJ', 'RIIIS', 'SEA', 'SJL', 'SMJ', 'ST', 'SWR', 'TALLYHO', 'RBF',
        'CAST', 'MOLF', 'PWD','PWK', 'SPRET', 'WSB']# ["SJL","129S1","BALB"] #"C57BL6J",
#CHROMSOME = [ str(c) for c in range(1,20)] + ["X", "Y", "MT"]
CHROMOSOME = ['1'] + [ str(c) for c in range(10,20)] + [ str(c) for c in range(2,10)]+ ["X"]
VG = expand("vg/GRCm38_chr{i}.{g}", i=CHROMOSOME, g=['xg','gcsa'])
#GT = expand("vg/{sample}.genotypes.vcf", sample=SAMPLES)
VCALL = expand("gams/{sample}_chr{i}/GRCm38_chr{i}_{sample}.vcf.gz", i=CHROMOSOME,  sample=SAMPLES)
VEVAL = expand("gams/{sample}_chr{i}/sv_accuracy.tsv", i=CHROMOSOME,  sample=SAMPLES)


rule target:
    input: VG, VEVAL, 

rule vg_chrom:
    output:  temp(expand("vgraph.chr{i}.tmp", i=CHROMOSOME))
    run:
        for v in output:
            shell("touch %s"%v)

# rule toilvg_construct:
#     input:
#         genome = GENOME,
#         tmp = "vgraph.chr{i}.tmp",
#         vcf = VCF_DELETION + ".gz",
#         vcfi =  VCF_DELETION + ".gz.tbi",
#     output:  
#         "vg/GRCm38_chr{i}.gcsa",
#         "vg/GRCm38_chr{i}.xg",
#         "vg/GRCm38_chr{i}.vg",
#         "vg/GRCm38_chr{i}_alts.gam",
#     params:
#         wkdir = WKDIR,
#     threads: 1
#     log: "logs/vg_construct_chr{i}.log"
#     run:
#         shell("rm -rf {params.wkdir}/jobStore_construct_chr{wildcards.i}")
#         shell("mkdir -p {params.wkdir}/vg ")
#         shell("toil-vg construct "
#         "--vcf {input.vcf} --fasta {input.genome} --regions {wildcards.i} "
#         "--out_name GRCm38_chr{wildcards.i}  --container None " # --container None 
#         "--xg_index --gcsa_index --realTimeLogging --pangenome "
#         "--flat_alts --handle_svs --alt_path_gam_index --validate --normalize "
#         "--merge_graphs --gcsa_index_cores {threads} --workDir {params.wkdir} "
#         "{params.wkdir}/jobStore{wildcards.i} {params.wkdir}/vg "
#         "2>&1 | tee {log} ")

rule split_bams:
    input: "/data/bases/fangzq/strains/{sample}/output.GATKrealigned.Recal.bam"
    output: temp("tmp_bams/{sample}.chr{i}.bam")
    shell:
        "samtools view -h -O BAM {input} {wildcards.i} > {output}"

rule vg_map:
    input:
        bam= "tmp_bams/{sample}.chr{i}.bam",
        gcsa = "vg/GRCm38_chr{i}.gcsa",
        xg = "vg/GRCm38_chr{i}.xg",
    output:
        gam = "gams/{sample}_chr{i}/{sample}_chr{i}_default.gam",
    threads: 1
    params:
        wkdir = WKDIR,
    log: "logs/vg_map_{sample}_chr{i}.log"
    run:
        shell("rm -rf {params.wkdir}/jobStore_{wildcards.sample}_{wildcards.i}")
        shell("mkdir -p {params.wkdir}/gams/{wildcards.sample}_chr{wildcards.i} ")
        shell("toil-vg map {params.wkdir}/jobStore_{wildcards.sample}_{wildcards.i} "
              "{wildcards.sample}_chr{wildcards.i} {params.wkdir}/gams/{wildcards.sample}_chr{wildcards.i} "
              "--xg_index {params.wkdir}/{input.xg} "
              "--gcsa_index {params.wkdir}/{input.gcsa} "
              "--bam_input_reads {params.wkdir}/{input.bam} "
              "--interleaved --alignment_cores {threads} "
              "--single_reads_chunk --realTimeLogging  --container None --workDir {params.wkdir} "
              "2>&1 | tee {log} ")

rule split_vcf:
    input:
        tmp = "vgraph.chr{i}.tmp",
        vcf = VCF_DELETION + ".gz",
        vcfi = VCF_DELETION + ".gz.tbi"
    output:
        vcf = VCF_DELETION + ".chr{i}.vcf.gz",
        vcfi = VCF_DELETION + ".chr{i}.vcf.gz.tbi"
    run:
        shell("bcftools view {input.vcf} {wildcards.i} | bcftools sort -Oz > {output.vcf}")
        shell("tabix -f -p vcf {output.vcf}")


rule vg_call:
    input:
        vcf = VCF_DELETION + ".chr{i}.vcf.gz",
        vcfi = VCF_DELETION + ".chr{i}.vcf.gz.tbi",
        bam= "tmp_bams/{sample}.chr{i}.bam",
        #altgam = "vg/GRCm38_chr{i}_alts.gam",
        gam = "gams/{sample}_chr{i}/{sample}_chr{i}_default.gam",
        vg = "vg/GRCm38_chr{i}.vg",
    output:
        vcf = "gams/{sample}_chr{i}/GRCm38_chr{i}_{sample}.vcf.gz"
    params:
        wkdir = WKDIR,
    threads: 1
    log: "logs/vg_call_{sample}_chr{i}.log"
    run:
        shell("rm -rf {params.wkdir}/jobStore_call_{wildcards.sample}_{wildcards.i}")
        shell("mkdir -p {params.wkdir}/gams/{wildcards.sample}_chr{wildcards.i} ")
        shell("toil-vg call {params.wkdir}/jobStore_call_{wildcards.sample}_{wildcards.i} "
              "{params.wkdir}/gams/{wildcards.sample}_chr{wildcards.i} --graph {input.vg} "
              "--gam {input.gam} "# --alt_path_gam {input.altgam} --chroms {wildcards.i}
              "--genotype_vcf {input.vcf} --sample {wildcards.sample} " # makesure you know the sample name in vcf file 
              "--realTimeLogging --calling_cores {threads} "
              "--container None --workDir {params.wkdir} "
              "2>&1 | tee {log} ")

              
## compare back to original vcfs
rule vg_eval:
    input: 
        vcf_truth = VCF_DELETION + ".chr{i}.vcf.gz",
        vcfi = VCF_DELETION + ".chr{i}.vcf.gz.tbi",
        vcf_call = "gams/{sample}_chr{i}/GRCm38_chr{i}_{sample}.vcf.gz"
    output:
       acc = "gams/{sample}_chr{i}/sv_accuracy.tsv",
       evals = "gams/{sample}_chr{i}/sv_evaluation.tar.gz"
    params:
        wkdir = WKDIR,
    threads: 1
    log: "logs/vg_eval_{sample}_chr{i}.log"
    run:
        # need to install R package BiocManager::install('jmonlong/sveval') 
        # need to install rtg-tools, conda install rtg-tools
        shell("rm -rf {params.wkdir}/jobStore_eval_{wildcards.sample}_{wildcards.i}")
        shell("mkdir -p {params.wkdir}/gams/{wildcards.sample}_chr{wildcards.i} ")
        shell("toil-vg vcfeval {params.wkdir}/jobStore_{wildcards.sample}_{wildcards.i} "
              "{params.wkdir}/gams/{wildcards.sample}_chr{wildcards.i} "
              "--vcfeval_baseline {input.vcf_truth} --call_vcf {input.vcf_call} "
              "--genotype_eval --sveval --vcfeval_sample {wildcards.sample} "
              "--realTimeLogging --container None --workDir {params.wkdir} "
              "--vcfeval_cores {threads} "
              "2>&1 | tee {log} " ) # --vcfeval_fasta {input.genome}

