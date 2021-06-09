




import os, glob
import pysam

## see the tutorial here: 
## https://github.com/vgteam/vg/wiki/Working-with-a-whole-genome-variation-graph
## https://gtpb.github.io/CPANG18/pages/toy_examples
## https://github.com/vgteam/sv-genotyping-paper

## Commands
WKDIR = "/data/bases/fangzq/20200815_SV/PacBio_project"
workdir: WKDIR
THREADS = 48
GENOME = "/home/fangzq/genome/mouse/GRCm38_68.fa"
VCF_5_STRAIN = "/data/bases/fangzq/20200815_SV/PacBio_project/VCFs/merged_gt_SURVIVOR_1kbpdist.filter.vcf"
VCF_DELETION = "BTBR.deletion.vcf" # "merged_deletion.vcf"
SAMPLES =  [ "BTBR"] #, "BALB","SJL","129S1","AJ", "BTBR"] # ["SJL","129S1","BALB"] #"C57BL6J",
#CHROMSOME = [ str(c) for c in range(1,20)] + ["X", "Y", "MT"]
CHROMOSOME = ['1'] + [ str(c) for c in range(10,20)] + [ str(c) for c in range(2,10)]+ ["X"]
VG = expand("vg/GRCm38_chr{i}.{g}", i=CHROMOSOME, g=['xg','gcsa'])
#GT = expand("vg/{sample}.genotypes.vcf", sample=SAMPLES)
VMAP = expand("vcfs/{sample}_chr{i}.vcf.gz", i=CHROMOSOME,  sample=SAMPLES)



rule target:
    input: VG, VMAP

rule vg_chrom:
    output:  temp(expand("vgraph.chr{i}.tmp", i=CHROMOSOME))
    run:
        for v in output:
            os.system("touch %s"%v)

rule toilvg_construct:
    input:
        genome = GENOME,
        tmp = "vgraph.chr{i}.tmp",
        vcf = VCF_DELETION + ".gz",
        vcfi =  VCF_DELETION + ".gz.tbi",
    output:  
        "vg/GRCm38_chr{i}.gcsa",
        "vg/GRCm38_chr{i}.xg",
        "vg/GRCm38_chr{i}.vg",
        "vg/GRCm38_chr{i}_alts.gam",
    params:
        wkdir = WKDIR,
    threads: 1
    log: "logs/vg_construct_chr{i}.log"
    run:
        shell("rm -rf {params.wkdir}/jobStore{wildcards.i}")
        shell("mkdir -p {params.wkdir}/vg ")
        shell("toil-vg construct "
        "--vcf {input.vcf} --fasta {input.genome} --regions {wildcards.i} "
        "--out_name GRCm38_chr{wildcards.i}  --container None " # --container None 
        "--xg_index --gcsa_index --realTimeLogging --pangenome "
        "--flat_alts --handle_svs --alt_path_gam_index --validate --normalize "
        "--merge_graphs --gcsa_index_cores {threads} --workDir {params.wkdir} "
        "{params.wkdir}/jobStore{wildcards.i} {params.wkdir}/vg "
        "2>&1 | tee {log} ")

rule split_bams:
    input: "../bams/{sample}.realign.bam"
    output: "tmp_bams/{sample}.chr{i}.bam"
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
        shell("mkdir -p {params.wkdir}/gams ")
        shell("toil-vg map {params.wkdir}/jobStore_{wildcards.sample}_{wildcards.i} "
              "{wildcards.sample}_chr{wildcards.i} {params.wkdir}/gams/{wildcards.sample}_chr{wildcards.i} "
              "--xg_index {params.wkdir}/{input.xg} "
              "--gcsa_index {params.wkdir}/{input.gcsa} "
              "--bam_input_reads {params.wkdir}/{input.bam} "
              "--interleaved --alignment_cores {threads} "
              "--single_reads_chunk --realTimeLogging  --container None --workDir {params.wkdir} "
              "2>&1 | tee {log} ")

rule vg_call:
    input:
        vcf = VCF_DELETION + ".gz",
        vcfi = VCF_DELETION + ".gz.tbi",
        bam= "tmp_bams/{sample}.chr{i}.bam",
        altgam = "vg/GRCm38_chr{i}_alts.gam",
        gam = "gams/{sample}_chr{i}/{sample}_chr{i}_default.gam",
        gcsa = "vg/GRCm38_chr{i}.gcsa",
        xg = "vg/GRCm38_chr{i}.xg",
    output:
        vcf = "vcfs/{sample}_chr{i}.vcf.gz"
    params:
        wkdir = WKDIR,
    threads: 1
    log: "logs/vg_call_{sample}_chr{i}.log"
    run:
        shell("rm -rf {params.wkdir}/jobStore_{wildcards.sample}_{wildcards.i}")
        shell("mkdir -p {params.wkdir}/vcfs ")
        shell("toil-vg call {params.wkdir}/jobStore_{wildcards.sample}_{wildcards.i} {input.xg} "
              "{wildcards.sample}_chr{wildcards.i} {params.wkdir}/vcfs "
              "--chroms {wildcards.i} --gams {input.gam} --alt_path_gam {input.altgam} "
              "--genotype_vcf {input.vcf} " 
              "--realTimeLogging --call_chunk_cores {threads} "
              "--container None --workDir {params.wkdir} "
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
        shell("bcftools view {input.vcf} -t {whildcards.i} -O z > {output.vcf}")
        shell("tabix -f -p vcf {output.vcf}")

## compare back to original vcfs
rule vg_eval:
    input: 
        vcf = VCF_DELETION + ".chr{i}.vcf.gz",
        vcfi = VCF_DELETION + ".chr{i}.vcf.gz.tbi",
        vg_vcf = "vcfs/{sample}_chr{i}.vcf.gz"
    output:
       acc = "vcfs/{sample}_chr{i}/sv_accuracy.tsv",
       eval = "vcfs/{sample}_chr{i}/sv_evaluation.tar.gz"
    params:
        wkdir = WKDIR,
    threads: 1
    log: "logs/vg_eval_{sample}_chr{i}.log"
    run:
        shell("rm -rf {params.wkdir}/jobStore_{wildcards.sample}_{wildcards.i}")
        shell("mkdir -p {params.wkdir}/vcfs/{sample}_chr{wildcards.i} ")
        shell("toil-vg vcfeval {params.wkdir}/jobStore_{wildcards.sample}_{wildcards.i} "
              "{params.wkdir}/vcfs/{sample}_chr{wildcards.i} "
              "--vcfeval_baseline {input.vcf} --call_vcf {input.vg_vcf} "
              "--genotype_eval --sveval --vcfeval_sample {wildcard.sample}_chr{wildcards.i} "
              "--realTimeLogging --workDir {params.wkdir} "
              "2>&1 | tee {log} " )

