
# pipelines from here
# https://github.com/dantaki/Local-De-Novo-Assembly-with-Long-Reads


# Local De Novo Assembly with Long Reads
# A quick pipeline on local de novo assembly of long reads for SV validation
## Tutorial

# ## Methodology

# 1. Convert BAM -> FASTQ 
# 2. Map with Minimap2/NGMLR
# 3. Extract Split-Reads and Convert to FASTQ
# 4. Local De Novo Assembly of Split-Reads with Flye
# 5. Align Contigs with LAST


## Commands
workdir: "/data/bases/fangzq/20200815_SV/PacBio_project"
GENOME = "/home/fangzq/genome/mouse/129S1/Mus_musculus_129s1svimj.129S1_SvImJ_v1.dna_sm.chromosome.6.fa"
LASTDB = expand("GRCm38.chr6.last.{suf}", suf=["bck","des","prj","sds","ssp","suf", "tis"])
SAMPLES =  ["BALB","SJL","129S1","AJ", "C57BL6J", "BTBR"] # ["SJL","129S1","BALB"] #
SR_SAM = expand("flye_sr/{sample}/minimap2_flye_ref129s1.sorted.bam.bai", sample=SAMPLES)

rule target:
    input: SR_SAM

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
        # """
        # # # get the header
        # # samtools view -H {input} > {output.sam}

        # # # grep split-reads
        # # samtools view {input} | grep "SA:Z"  >> {output.sam}

        # # convert to FASTQ
        # samtools bam2fq --threads {threads} {input} > {output.fastq}
        # """


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


### 5 Align with LAST

# make a reference, just from chromosome 10
rule makeReference:
    input:
       ref = GENOME,
    output:
        "GRCm38.chr6.fa",
        "GRCm38.chr6.dict",
        LASTDB,
    params:
        chrom = 6    
    shell:
        """
        samtools faidx {input.ref} {params.chrom} > GRCm38.chr{params.chrom}.fa

        gatk CreateSequenceDictionary -R GRCm38.chr{params.chrom}.fa -O GRCm38.chr{params.chrom}.dict 
        
        lastdb -uNEAR -R01 GRCm38.chr{params.chrom}.last GRCm38.chr{params.chrom}.fa
        """

# optional, clean up contig names for LAST
rule CleanContigName:
    input: "flye_sr/{sample}/assembly.fasta",
    output: "flye_sr/{sample}/minimap2_sr_flye.fa"
    shell:
        """
        awk '/>/ {{$0 = ">" ++n}} 1' {input} > {output}
        """

rule LastParames:
    input: 
        last = LASTDB,
        flye = "flye_sr/{sample}/minimap2_sr_flye.fa"
    output:
        "flye_sr/{sample}/minimap2_ref129s1.flye.par"
    params:
        lastdb = "GRCm38.chr6.last",
    threads: 24
    shell:
        # create paramter file for LAST
        "last-train -P{threads} -Q0 {params.lastdb} {input.flye} > {output}"

# align with LAST and convert MAF file to SAM format
rule last_align:
    input: 
        flye = "flye_sr/{sample}/minimap2_sr_flye.fa",
        param = "flye_sr/{sample}/minimap2_ref129s1.flye.par",
    output: 
        "flye_sr/{sample}/minimap2_ref129s1.maf"
    params:
       lastdb = "GRCm38.chr6.last",
    threads: 24
    shell:
        "lastal -P{threads} -p {input.param} {params.lastdb} {input.flye} | last-split -m1 > {output}"


rule maf2bam:
    input: 
        maf="flye_sr/{sample}/minimap2_ref129s1.maf",
        fadict = "GRCm38.chr6.dict"
    output: "flye_sr/{sample}/minimap2_flye_ref129s1.sorted.bam"
    shell:
        "maf-convert -f {input.fadict} sam {input.maf} | "
        "samtools view -Sb | "
        "samtools sort > {output}"

rule bamIndex:
    input: "flye_sr/{sample}/minimap2_flye_ref129s1.sorted.bam"
    output: "flye_sr/{sample}/minimap2_flye_ref129s1.sorted.bam.bai"
    shell:
        "samtools index {input}"


