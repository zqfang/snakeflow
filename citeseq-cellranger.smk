# snakemake pipeline for 10x genomics single cell CITE-seq
# @author: Zhuoqing Fang
# @email: maxzqfang@stanford.edu
# @date: May 7th, 2020
import glob, re, os

WORKSPACE = "/data/bases/fangzq/20200505_HCC"
workdir: WORKSPACE
#### globals ########
#configfile: "config.yaml"
#workdir: config["workdir"]

#SAMPLES = config["samples"]
#READ1 = condfig["r1"]
#READ2 = config["r2"]
#READi = config["ri"]
#TRANSCRIPTOM10X = config["transcriptom"]

########### INPUTS ###########################################

# fastqs and hastag sequence
HASHTAG = os.path.join(WORKSPACE, "hashtag.totalseqB.V3.feature_ref.csv")
FASTQ_DIR = "/data/bases/shared/HCC/30-363886939/00_fastq"

# automatic search sample names using regex
pat = re.compile(r"([A-Z0-9-]+)-([A-Z0-9a-z]{2,4})_[0-9A-Z]{2}_(L[0-9]{3})_([RI][12])_[0-9]{3}.fastq.gz")
# name, capture, lane, reads = pat.search(f).groups()
files = glob.glob(os.path.join(FASTQ_DIR, "*fastq.gz"))
files = [pat.search(f).groups() for f in files]
# remove duplicate entries
files = {s[0]: { x[1]: x[0] for x in files } for s in files}

lane = "S1_L001"
PROJECT = "HCC"
SAMPLES = ["HCC1-GE", "HCC2-GE", "HCC3-2-GE"]
SAMPLE_HASH=["HCC1-HASH", "HCC3-HASH", "HCC3-HASH"]

# cellranger
TRANSCRIPTOM10X = "/home/fangzq/genome/10XRefdata/refdata-cellranger-GRCh38-3.0.0" 

# vireo
SNV = "/home/fangzq/genome"
KNOWN_SNV = os.path.join(SNV, "genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf")
DONOR_NUM_PER_SAMPLE = 3 # per 10x run samples/patients ...


#### Target output files #####
MTX =  expand("{sample}_cnt/outs/filtered_feature_bc_matrix/matrix.mtx.gz",sample=SAMPLES)
BARCODE = expand("{sample}_cnt/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", sample=SAMPLES)
FEATURES = expand("{sample}_cnt/outs/filtered_feature_bc_matrix/features.tsv.gz",sample=SAMPLES)
AGGR="{project}_aggr/outs/filtered_feature_bc_matrix.h5".format(project = PROJECT)
BAM = expand("{sample}_cnt/outs/possorted_genome_bam.bam", sample=SAMPLES)
VIREO = expand("{sample}_vireo/results.demux/donor_ids.tsv", sample=SAMPLES)

#### Rules ########
rule target: 
    input: MTX, BARCODE, FEATURES, BAM, AGGR, VIREO


rule prepare_count:
    output: 
        expand("{sample}.libraries.csv", sample=SAMPLES)
    params:
        samples=SAMPLES,
        fastqs=FASTQ_DIR,
        capture="HASH", # tag name represent which file is hashtag fastq
        files=files,
    run:
        for sample in params.samples:
            hashtag = params.files[sample]
            out = sample+".libraries.csv"
            with open(out, 'w') as w:
                w.write("fastqs,sample,library_type\n")
                for k, v in hashtag.items():
                    if k == params.capture:
                        w.write("%s,%s,Antibody Capture\n"%(params.fastqs, sample))
                    else:
                        w.write("%s,%s,Gene Expression\n"%(params.fastqs, sample))



# snakemake creates the folder if it does not exisit. but cellranger fails if the folder is 
# not created by itself. 

# trick 1: rm {sample}_cnt/_lock file and re-run pipeline, it will be OK.
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/troubleshooting
# trick 2: rm output folder before cellranger run
rule citeseq_count:
    input:
        read1 = os.path.join(FASTQ_DIR, "{sample}_%s_R1_001.fastq.gz"%lane),
        read2 = os.path.join(FASTQ_DIR, "{sample}_%s_R2_001.fastq.gz"%lane),
        readi = os.path.join(FASTQ_DIR, "{sample}_%s_I1_001.fastq.gz"%lane), 
        transcriptom = os.path.join(TRANSCRIPTOM10X, "fasta/genome.fa"),
        antibody = HASHTAG,
        libraries = "{sample}.libraries.csv",
        #libraries = rules.pre_count.output
    output:
        bam = "{sample}_cnt/outs/possorted_genome_bam.bam",
        barcode = "{sample}_cnt/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        features = "{sample}_cnt/outs/filtered_feature_bc_matrix/features.tsv.gz",
        mtx = "{sample}_cnt/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
        mh5 =  "{sample}_cnt/outs/molecule_info.h5",
    threads: 12
    log: "logs/{sample}.count.log"
    version: "0.1"
    params:
        extra=" --expect-cells=2000 --nosecondary", # --expect-cells=12000
        fastqs=FASTQ_DIR,
        transcriptom=TRANSCRIPTOM10X,
    run:
        # non antibody captured code
        # "cellranger count --id={wildcards.sample}.cnt"
        # "--fastqs={params.fastqs} "
        # "--sample={wildcards.sample} "
        # rm outputdir to solve "pipestance directory" error from cellranger
        shell("rm -rf {wildcards.sample}_cnt")
        shell("cellranger count --id={wildcards.sample}_cnt " +\ 
                "--transcriptome={params.transcriptom} " +\
                "--libraries {input.libraries} " +\
                "--feature-ref {input.antibody} " +\
                "--localcores={threads} {params.extra} > {log}") 
        
rule prepare_aggr:
    input: expand("{sample}_cnt/outs/molecule_info.h5", sample=SAMPLES)
    output: "%s.aggr.csv"%PROJECT
    params: 
        samples = SAMPLES,
        ids = [s[:4] for s in SAMPLES], # only used for loupe browser
        condition = [lane]*3 #CONDITION, ### ->>>>
    run:
        run = [s[:4] for s in SAMPLES]
        with open(output[0],'w') as out:
            # only add 'batch' when you need Chemistry Batch Correction
            # batch is optional, see docs if you have bath info
            # https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate
            out.write("library_id,molecule_h5,run\n")
            for s, b in zip(params.samples, params.ids):
                out.write(f"{s},{s}_cnt/outs/molecule_info.h5,{b}\n")


rule citeseq_aggr:
    input:
        mh5 = expand("{sample}_cnt/outs/molecule_info.h5", sample=SAMPLES),
        aggr = "{project}.aggr.csv"
    output:
        barcode = "{project}_aggr/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        features = "{project}_aggr/outs/filtered_feature_bc_matrix/features.tsv.gz",
        mtx = "{project}_aggr/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
        mh5 =  "{project}_aggr/outs/filtered_feature_bc_matrix.h5" 
    log: "logs/{project}.aggr.log"
    threads: 12
    run:
        # rm outputdir to solve "pipestance directory" error from cellranger
        shell("rm -rf {wildcards.project}_aggr")
        shell("cellranger aggr --id={wildcards.project}_aggr " +\ 
              "--csv={input.aggr} " +\
              "--normalize=mapped > {log}")

# demuxlet --sam ${sample}citeseq/outs/possorted_genome_bam.bam  \
#          --vcf /home/fangzq/projects/merged.sorted.vcf \
#          --out ${sample}.demuxlet \
#          --min-mac 1 --field GT --alpha 0  \
#          --group-list /home/fangzq/projects/SFGF-Peltz-YG-13709/${sample}citeseq/outs/filtered_feature_bc_matrix/barcodes.tsv.gz

# call single cell SNP

# a single BAM/SAM file, e.g., from cellranger, 
# a list of cell barcodes, a VCF file for common SNPs. 
# This mode is recommended comparing to mode 2, if a list of common SNP is known, e.g., human (see Candidate SNPs below)
rule cell_snp:
    input:
        bam="{sample}_cnt/outs/possorted_genome_bam.bam",
        barcode="{sample}_cnt/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        vcf=KNOWN_SNV,
    output:
        "{sample}_vireo/CELL_DATA/cellSNP.cells.vcf.gz",
        "{sample}_vireo/CELL_DATA/cellSNP.base.vcf.gz"
    threads: 12
    log: "logs/{sample}.cellsnp.log"
    shell:
        "cellSNP -s {input.bam} -b {input.barcode} -R {input.vcf} "
        "-p {threads} --minMAF 0.1 --minCOUNT 20 "
        "-O {wildcards.sample}_vireo/CELL_DATA > {log}"

# # deconvolution without known genotype
# You have to know how many individals, e.g. n_donor=3
rule vireo_demux:
    input: 
        "{sample}_vireo/CELL_DATA/cellSNP.cells.vcf.gz",
        "{sample}_vireo/CELL_DATA/cellSNP.base.vcf.gz",
    output: 
        "{sample}_vireo/results.demux/donor_ids.tsv",
        "{sample}_vireo/results.demux/GT_donors.vireo.vcf.gz",
    params:
        donor_num = DONOR_NUM_PER_SAMPLE,
        CELL_DATA = "{sample}_vireo/CELL_DATA"
        # with known genotype, add below argument
        # known_gt = "-d {DONOR_GT_FILE} -t GT" # DONOR_GT_FILE(vcf) obtained from input samples
    log: "logs/{sample}.vireo.log"
    shell:
        "vireo -c {params.CELL_DATA} -N {params.donor_num} -o {wildcards.sample}_vireo/results.demux > {log}"

# Mode 3: pileup a list of SNPs for one or multiple BAM/SAM files
# rule cell_SNP3:
#     input:
#         bam=expand("{sample}_cnt/outs/possorted_genome_bam.bam", sample=SAMPLES)
#         barcode=expand("{sample}_cnt/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", sample=SAMPLES)
#         vcf=
#     output:
#     params:
#         extra="--UMItag None" # only include this if smart-seq or bulk rna-seq
#         bams=",".join(expand("{sample}_cnt/outs/possorted_genome_bam.bam", sample=SAMPLES)),
#         ids=",".join(SAMPLES)
#     log:
#     shell:
#         "cellSNP -s {params.bams} -I {params.ids} -o $OUT_FILE -R {input.vcf} -p 20"

# https://github.com/10XGenomics/vartrix
rule vatrix:
    input:
        bam="{sample}_cnt/outs/possorted_genome_bam.bam",
        barcode="{sample}_cnt/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        vcf=KNOWN_SNV,
        genome=os.path.join(TRANSCRIPTOM10X, "fasta/genome.fa")
    output:
        "{sample}.vatrix/matrix.mtx"
    threads: 12
    shell:
        "vartrix -v {input.vcf} -b {input.bam} -f {input.genome} -c {input.barcode} "
        "--threads {threads} -o {wildcards.sample}.vatrix"

rule vcf2SNVloci:
    input: KNOWN_SNV
    output: "SNV.loci.txt"
    shell:
        "awk '{{print $1,$2}}' {input} | sed -i 's/\s/:/g' > {output}"
