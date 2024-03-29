# snakemake pipeline for 10x genomics single cell RNA-seq 
# @author: Zhuoqing Fang
# @email: fangzq@stanford.edu
# @date: May 31th, 2022
import glob, re, os

WORKSPACE = "/data/bases/fangzq/NK/20221220_iPSC-nk"
workdir: WORKSPACE

########### INPUTS ###########################################

FASTQ_DIR = os.path.join(WORKSPACE, "FASTQ")
# cellranger reference file
TRANSCRIPTOM10X = "/home/fangzq/genome/10XRefdata/refdata-gex-GRCh38-2020-A"

# automatic search sample names using regex
pat = re.compile(r"(.+)_(S\d{1})_(L\d{3})_(\w{1}\d{1})_\d{3}.fastq.gz")
# name, capture, lane, reads = pat.search(f).groups()
files = glob.glob(os.path.join(FASTQ_DIR, "*fastq.gz"))
names = [pat.search(f).groups() for f in files]
# remove duplicate entries
samples = sorted(list(set([f[0].split("/")[-1] for f in names])))



WORKSPACE = "/oak/stanford/groups/bpulend/lab.members/zhuoqing/X202SC23010264-Z01-F002-Pulendran-ShengyangWu-031023/"
workdir: WORKSPACE

########### INPUTS ###########################################

FASTQ_DIR = os.path.join(WORKSPACE, "01.RawData")
# cellranger reference file
TRANSCRIPTOM10X = "/oak/stanford/projects/genomics/references/refdata-gex-Mmul_10.109"

# automatic search sample names using regex
pat = re.compile(r"(.+)_(S\d{1})_(L\d{3})_(\w{1}\d{1})_001.fastq.gz")
# name, capture, lane, reads = pat.search(f).groups()
files = glob.glob(os.path.join(FASTQ_DIR, "**/*fastq.gz"), recursive = True)
names = [pat.search(f).groups() for f in files]
# remove duplicate entries
samples = sorted(list(set([f[0].split("/")[-1] for f in names])))



PROJECT = "iPSC-NK"
SAMPLES = ["CiNK", "FiNK", "PBNK", "PBNKIL15"] #samples 
### get fastqs for Gene Expression
# SAMPLES = ['RNA_%s'%i for i in range(1,7) ] #samples
#FASTQ_SAMPLES = ['RNA_1_CKDL230004165-1A_HT3KFDSX5',
# 'RNA_2_CKDL230004166-1A_HT3KFDSX5',
# 'RNA_3_CKDL230004167-1A_HT3KFDSX5',
# 'RNA_4_CKDL230004168-1A_HT3KFDSX5',
# 'RNA_5_CKDL230004169-1A_HT3KFDSX5',
# 'RNA_6_CKDL230004170-1A_HT3KFDSX5']
FASTQ_SAMPLES = { s : s for s in samples} # in case in have a shorted name in SAMPLES than the FASTQ_name prefiex (samples)
SAMPLES = sorted(list(FASTQ_SAMPLES.keys()))

### get fastqs for Gene Expression
FASTQ_GE = { k:[f for f in files if os.path.basename(f).startswith(s+"_")] for k,s in FASTQ_SAMPL
ES.items() }
#### Target output files #####
MTX =  expand("{sample}_cnt/outs/filtered_feature_bc_matrix/matrix.mtx.gz",sample=SAMPLES)
BARCODE = expand("{sample}_cnt/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", sample=SAMPLES)
FEATURES = expand("{sample}_cnt/outs/filtered_feature_bc_matrix/features.tsv.gz",sample=SAMPLES)
AGGR="{project}_aggr/outs/count/filtered_feature_bc_matrix.h5".format(project = PROJECT)
BAM = expand("{sample}_cnt/outs/possorted_genome_bam.bam", sample=SAMPLES)
#### Rules ########
rule target: 
    input: MTX, BARCODE, FEATURES, BAM, AGGR


# snakemake creates the folder if it does not exisit. but cellranger fails if the folder is 
# not created by itself. 

# trick 1: rm {sample}_cnt/_lock file and re-run pipeline, it will be OK.
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/troubleshooting
# trick 2: rm output folder before cellranger run
rule cellranger_count:
    input:
        fastq_ge = lambda wildcards: FASTQ_GE[wildcards.sample],
        transcriptom = os.path.join(TRANSCRIPTOM10X, "fasta/genome.fa"),
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
        extra=" --expect-cells=6000 --nosecondary", # loading cells 9000
        fastqs=  lambda wildcards: os.path.dirname(FASTQ_GE[wildcards.sample][0]),
        sample = lambda wildcards: FASTQ_SAMPLES[wildcards.sample],
        transcriptom=TRANSCRIPTOM10X,
    run:
        # non antibody captured code
        # "cellranger count --id={wildcards.sample}.cnt"
        # "--fastqs={params.fastqs} "
        # "--sample={wildcards.sample} "
        # rm outputdir to solve "pipestance directory" error from cellranger
        shell("rm -rf {wildcards.sample}_cnt")
        shell("cellranger count --id={wildcards.sample}_cnt "
                "--transcriptome={params.transcriptom} "
                "--fastqs={params.fastqs} "
                "--sample={params.sample} "
                "--localcores={threads} {params.extra} > {log}")
        
rule aggregate_list:
    input: expand("{sample}_cnt/outs/molecule_info.h5", sample=SAMPLES)
    output: "%s.aggr.csv"%PROJECT
    params: 
        samples = SAMPLES,
        ids = SAMPLES, # only used for loupe browser
        # condition = [lane]*3 #CONDITION, ### ->>>>
    run:
        with open(output[0],'w') as out:
            # only add 'batch' when you need Chemistry Batch Correction
            # batch is optional, see docs if you have bath info
            # https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate
            out.write("sample_id,molecule_h5,run\n")
            for s, b in zip(params.samples, params.ids):
                out.write(f"{s},{s}_cnt/outs/molecule_info.h5,{b}\n")


rule cellranger_aggr:
    """
    see the docs for cellranger aggr
    https://kb.10xgenomics.com/hc/en-us/articles/115004217543-How-does-cellranger-aggr-normalize-for-sequencing-depth-among-multiple-libraries-

    when normalized=mapped, 
    the default behavior of cell rangers will subsamples reads from higher-depth libraries 
    until all libraries have an equal number of confidently mapped reads per cell. 
    but this does have the potential for discarding a lot of data

    The seurat authors do not suggest: normalized=mapped, since we'll do cell level normalization (LogNormalize) later
    although we could do normalized=mapped first then run seruat pipeline.
    see here: https://github.com/satijalab/seurat/issues/672
    """
    input:
        mh5 = expand("{sample}_cnt/outs/molecule_info.h5", sample=SAMPLES),
        aggr = "%s.aggr.csv"%PROJECT
    output:
        barcode = "%s_aggr/outs/count/filtered_feature_bc_matrix/barcodes.tsv.gz"%PROJECT,
        features = "%s_aggr/outs/count/filtered_feature_bc_matrix/features.tsv.gz"%PROJECT,
        mtx = "%s_aggr/outs/count/filtered_feature_bc_matrix/matrix.mtx.gz"%PROJECT,
        mh5 =  "%s_aggr/outs/count/filtered_feature_bc_matrix.h5"%PROJECT,
        agg ="%s_aggr/outs/aggregation.csv"%PROJECT,
    params:
        project = PROJECT,
    log: "logs/%s.aggr.log"%PROJECT 
    threads: 32
    run:
        # rm outputdir to solve "pipestance directory" error from cellranger
        shell("rm -rf {params.project}_aggr")
        shell("cellranger aggr --id={params.project}_aggr --nosecondary " 
              "--csv={input.aggr} " 
              "--localcores={threads} "
              "--normalize=none > {log}") # if =mapped, then normalizes for effective sequencing depth by subsampling the reads.
