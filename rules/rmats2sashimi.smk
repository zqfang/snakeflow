
from os.path import join, isfile
from snakemake.shell import shell

SASHIMIPLOT="logs/sashimiplot/test.txt"
SAMPLES, = glob_wildcards(join('fastq_clean', '{sample, [^/]+}_R1.fastq.gz'))

rule target:
    input: SASHIMIPLOT


rule sashimiplot:
    input:
        evt="alternative_splicing/rMATS.{treat}_vs_{ctrl}_sig/{event}.MATS.JCEC.sig.txt",
        bam=expand("mapped/{sample}.sorted.bam", sample=SAMPLES),
        bai=expand("mapped/{sample}.sorted.bam.bai", sample=SAMPLES),
        b1="temp/rmats/b_{treat}.txt",
        b2="temp/rmats/b_{ctrl}.txt",
    output: 
        "alternative_splicing/rMATS.{treat}_vs_{ctrl}_sig/{event}.sashimi",
        "logs/sashimiplot/test.txt"
    params:
        extra="--exon_s 1 --intron_s 5",
        outdir="alternative_splicing/rMATS.{treat}_vs_{ctrl}_sig/{event}.sashimi",
    run:
        # blacklist to skip
        if isfile("temp/blacklist.txt"):
            with open("temp/blacklist.txt") as black:
                blacklist = [ bla.strip().split("/")[-1] for bla in black]
            # groups you want to skip 
            bk = "diff_%s_vs_%s_results.annotated.xls"%(wildcards.treat, wildcards.ctrl)
            if bk in blacklist:
                return
        t = open(input.b1).read().strip("\n").lstrip("/data/")
        c = open(input.b2).read().strip("\n").lstrip("/data/")
        shell("""
              rmats2sashimiplot --b1 %s --b2 %s -t %s -e %s \
              --l1 %s --l2 %s %s -o %s """.format(t, c, wildcards.event, 
                                                  input.evt, wildcards.treat, wildcards.ctrl,
                                                  params.extra, params.outdir))
        shell("touch logs/sashimiplot/test.txt")

