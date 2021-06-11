
### pacbio long-reads 
### mapp reads using minimap2
### denovo assembly with flye
### align assembly to refenrece with unimap

## there'er difficulties in installing svviz2, since it's out of maintainance
# install svviz2 in this way will work
# conda create -n svviz python=3.7
# pip install numpy pandas scipy rpy2 pysam pyfaidx
# pip install cython
# pip install seqlib
# pip install -U git+https://github.com/nspies/genomeview.git
# pip install -U git+git://github.com/nspies/svviz2.git



## Commands
workdir: "/data/bases/fangzq/20200815_SV/PacBio_project"
GENOME = "/home/fangzq/genome/mouse/GRCm38_68.fa"
SAMPLES =  ["BALB","SJL","129S1","AJ", "C57BL6J", "BTBR"] # ["SJL","129S1","BALB"] #
BAM = expand("flye_sr/{sample}.unimap.sorted.bam", sample=SAMPLES)
BAM_chr6 = expand("BAMs/{sample}.chr6.sorted.bam", sample=SAMPLES)


rule target:
    input: BAM, BAM_chr6, "ly49_sv/results.tsv"

# ### 1 Bam2Fq


# rule bam2fastq:
#     input: "{sample}.bam"
#     output: "{sample}.raw.fastq"
#     shell:
#         "samtools bam2fq {input}  > {output}

# ### 2 Map with minimap2 or ngmlr


# rule minimap2:
#     input: 
#         ref = "ref.fa",
#         fastq = "{sample}.raw.fastq"
#     output: 
#         "{sample}.sam"
#     shell:
#         "minimap2 -ax map-pb {input.ref} {input.fastq} > {output}"

#    #  ngmlr -r human_g1k_v37.fasta -q raw.fq -x pacbio >ngmlr.sam

### 3 Extract Split-Reads 

rule splitRead2fastq:
    input: "BAMs/{sample}.sorted.MD.bam"
    output: 
        fastq="FASTQs/{sample}.fastq",
        #sam=temp("{sample}.sam")
    threads: 6
    shell:
        "samtools bam2fq --threads {threads} {input} > {output.fastq}"

### 4 Assemble Contigs with Flye


# -g is the genome size
## --meta is for metagenomes / uneven coverage. Since we are assembling split-reads 
## the coverage is very uneven. Might want to experiment with this option, though
rule flye:
    input: "FASTQs/{sample}.fastq"
    output: "flye_sr/{sample}/assembly.fasta"
    threads: 36
    shell:
        "flye -t {threads} --pacbio-raw {input} -o flye_sr/{wildcards.sample} " # --meta 

# 
rule unimap:
    input:
        genome=GENOME,
        assembly="flye_sr/{sample}/assembly.fasta"
    output:
        temp("flye_sr/{sample}.unimap.sam")
    threads: 12
    shell:
        "unimap -a -x asm5 --cs -t {threads} {input.genome} {input.assembly} "
        "> {output}"
        
rule samtools_sort:
    input: "BAMs/{sample}.unimap.sam"
    output: 
        bam = protected("flye_sr/{sample}.unimap.sorted.bam"),
        bai = protected("flye_sr/{sample}.unimap.sorted.bam.bai"),
    # threads: 8
    run:
        shell("samtools view -Sb {input} | samtools sort -m 2G  > {output.bam}")
        shell("samtools index {output.bam}")


rule samtools_chr6:
    input: "BAMs/{sample}.unimap.sorted.bam"
    output: 
        bam = protected("BAMs/{sample}.chr6.sorted.bam"),
        bai = protected("BAMs/{sample}.chr6.sorted.bam.bai"),
    # threads: 8
    run:
        shell("samtools view -O BAM {input} 6 > {output.bam}")
        shell("samtools index {output.bam}")


rule svviz2:
    input:
        genome= GENOME,
        bams = expand("BAMs/{sample}.sorted.MD.bam", sample=SAMPLES),
        vcf = "SJL.ly49.vcf",
    output:
        # directory("ly49_sv"),
        "ly49_sv/results.tsv",
    params:
        outdir="ly49_sv",
    shell:
        "/home/fangzq/miniconda/envs/svviz2/bin/svviz2 --ref {input.genome} "
        "--savereads --outdir {params.outdir} "
        "--variants {input.vcf} {input.bams} "

