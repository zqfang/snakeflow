# snakemake pipeline for 10x genomics single cell RNA-seq 
# @author: Zhuoqing Fang
# @email: fangzq@stanford.edu
# @date: May 31th, 2022
import glob, re, os

WORKSPACE = "/data/bases/fangzq/20220531_JJ_scRNAseq/"
workdir: WORKSPACE

########### INPUTS ###########################################

FASTQ_DIR = os.path.join(WORKSPACE, "fastq_path")
# split-seq reference file
REF = "/data/bases/shared/ParseBio/hg38_gfp_ref"
### GENOME and annotation
GENOME = "/data/bases/shared/ParseBio/GRCh38.gencode.v33.annotation.gfp.gtf"
GTF = "/data/bases/shared/ParseBio/GRCh38.p13.genome.gfp.fa"
CMD = "/data/bases/shared/ParseBio/ParseBiosciences-Pipeline.1.0.3p/split-pipe"


# automatic search sample names using regex
pat = re.compile(r"(.+)_(S\d{1})_(L\d{3})_(\w{1}\d{1})_\d{3}.fastq.gz")
# name, capture, lane, reads = pat.search(f).groups()
files = glob.glob(os.path.join(FASTQ_DIR, "*fastq.gz"))
names = [pat.search(f).groups() for f in files]
# remove duplicate entries
samples = sorted(list(set([f[0].split("/")[-1] for f in names])))

PROJECT = "Liver_organoid"
SAMPLES = ['NC', "HB", "T", "B"] #samples 
SUBLIBRARIES = ["L1", "L2"]
###
SAMPLE_LIST = "sample-list.txt"
# List file contents:
# cat sample-list.txt
# NC A1-A6
# T C1-C6
# P D1-D6
# HB A10-A12


### get fastqs for Gene Expression
FASTQ_GE = {s:[os.path.join(FASTQ_DIR, f) for f in os.listdir(FASTQ_DIR) if f.startswith(s)] for s in SAMPLES}



#### Target output files #####
MTX =  expand("{sublibray}_cnt/all-well/DGE_filtered/DGE.mtx",sample=SUBLIBRARIES)
BARCODE = expand("{sublibray}_cnt/all-well/DGE_filtered/DGE.mtx", sample=SUBLIBRARIES)
FEATURES = expand("{sublibray}_cnt/all-well/DGE_filtered/DGE.mtx",sample=SUBLIBRARIES)
AGGR="{project}_aggr/all-well/DGE_filtered/DGE.mtx".format(project = PROJECT)
BAM = expand("{sublibray}_cnt/outs/possorted_genome_bam.bam", sample=SUBLIBRARIES)
#### Rules ########
rule target: 
    input: MTX, BARCODE, FEATURES, BAM, AGGR

rule splitpipe_ref:
    input:
        genome = GENOME,
        gtf = GTF
    output:
        os.path.join(REF, "SAindex"),
        os.path.join(REF,"SA"),
        os.path.join(REF, "genome.fas.gz")
    params:
        genome_name = "hg38",
        outdir = REF
    threads: 16
    shell:
        "./split-pipe --mode mkref "
        "--genome_name {params.genome_name} "
        "--fasta {input.genome} "
        "--genes {input.gtf} "
        "--output_dir {params.outdir} "
        "--nthreads {threads}"



# snakemake creates the folder if it does not exisit. but cellranger fails if the folder is 
# not created by itself. 

# trick 1: rm {sample}_cnt/_lock file and re-run pipeline, it will be OK.
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/troubleshooting
# trick 2: rm output folder before cellranger run
rule splitpipe_count:
    input:
        fastq_ge = lambda wildcards: sorted(FASTQ_GE[wildcards.sublib]),
        ref= os.path.join(REF, "genome.fas.gz"),
        wells = SAMPLE_LIST # sample to well metadata
    output:
        bam = "{sublib}_cnt/process/barcode_headAligned_anno.bam",
        barcode = "{sublib}_cnt/all-well/DGE_filtered/cell_metadata.csv",
        features = "{sublib}_cnt/all-well/DGE_filtered/all_genes.csv",
        mtx = "{sublib}_cnt/all-well/DGE_filtered/DGE.mtx",
        h5ad = "{sublib}_cnt/all-well/DGE_filtered/anndata.h5ad"
    threads: 16
    log: "logs/{sublib}.count.log"
    version: "0.1"
    params:
        extra=" --expect-cells=6000 --nosecondary", # loading cells 9000
        fastqs=FASTQ_DIR,
        transcriptom=TRANSCRIPTOM10X,
    run:
        # non antibody captured code
        # "cellranger count --id={wildcards.sample}.cnt"
        # "--fastqs={params.fastqs} "
        # "--sample={wildcards.sample} "
        # rm outputdir to solve "pipestance directory" error from cellranger
        shell("rm -rf {wildcards.sublib}_cnt")
        shell("cellranger count --nthreads 36 "
              "--mode all --chemistry v2 
              "--run_name {wildcards.sublib} " 
              "--genome_dir {input.ref} "
              "--fq1 {input.fastq_ge[0]} "
              "--fq2 {input.fastq_ge[1]} "
              "--samp_list {input.wells} "
              "--output_dir {wildcards.sublib}_cnt"
              " > {log}") 
        
# rule aggregate_list:
#     input: expand("{sublib}_cnt/all-well/DGE_filtered/anndata.h5ad", sublib=SUBLIBRARIES)
#     output: "%s.aggr.csv"%PROJECT
#     params: 
#         samples = SAMPLES,
#         ids = SAMPLES, # only used for loupe browser
#         # condition = [lane]*3 #CONDITION, ### ->>>>
#     run:
#         with open(output[0],'w') as out:
#             # only add 'batch' when you need Chemistry Batch Correction
#             # batch is optional, see docs if you have bath info
#             # https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate
#             out.write("sample_id,molecule_h5,run\n")
#             for s, b in zip(params.samples, params.ids):
#                 out.write(f"{s},{s}_cnt/outs/molecule_info.h5,{b}\n")


rule splitseq_aggr:
    """
    combine sublibraries
    """
    input:
        mh5 = expand("{sublib}_cnt/all-well/DGE_filtered/anndata.h5ad", sublib=SUBLIBRARIES),
    output:
        barcode = "%s_aggr/all-well/DGE_filtered/cell_metadata.csv"%PROJECT,
        features = "%s_aggr/all-well/DGE_filtered/all_genes.csv"%PROJECT,
        mtx = "%s_aggr/all-well/DGE_filtered/DGE.mtx"%PROJECT,
        mh5 =  "%s_aggr/all-well/DGE_filtered/anndata.h5ad"%PROJECT,
    params:
        project = PROJECT,
        libs = SUBLIBRARIES
    log: "logs/%s.aggr.log"%PROJECT 
    threads: 16
    run:
        shell("rm -rf {params.project}_aggr")
        shell("cellranger aggr --mode comb " 
              "--sublibraries {params.libs}" 
              "--nthreads {threads} "
              "> {log}")
