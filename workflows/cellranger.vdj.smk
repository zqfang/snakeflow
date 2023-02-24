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
TRANSCRIPTOM10X = "/home/groups/bpulend/zqfang/refdata_mm10_YFP_CRE"
VDJ = 

# automatic search sample names using regex
pat = re.compile(r"(.+)_(S\d{1})_(L\d{3})_(\w{1}\d{1})_\d{3}.fastq.gz")
# name, capture, lane, reads = pat.search(f).groups()
files = glob.glob(os.path.join(FASTQ_DIR, "*fastq.gz"))
names = [pat.search(f).groups() for f in files]
# remove duplicate entries
samples = sorted(list(set([f[0].split("/")[-1] for f in names])))

PROJECT = "iPSC-NK"
SAMPLES = ["CiNK", "FiNK", "PBNK", "PBNKIL15"] #samples 
### get fastqs for Gene Expression
FASTQ_GE = {s:[os.path.join(FASTQ_DIR, f) for f in os.listdir(FASTQ_DIR) if f.startswith(s)] for s in SAMPLES}

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
rule prepare_config:
    input: 
        ref_tx=TRANSCRIPTOM10X,
        ref_vdj=VDJ,
    output: expand("{sample}.multi.config.csv", sample=SAMPLES)
    params:
        ref_tx=TRANSCRIPTOM10X,
        ref_vdj=VDJ,   
        samp_gex = SAMPLES,
        samp_vdj = [s + "-BCR" for s in SAMPLES],
        fastq = FASTQ_DIR
    run:
        for sg, sb in zip(params.samp_gex, params.vdj):
            config="""[gene-expression]
                    reference,{ref},
                    [vdj]
                    reference,{vdj},
                    [libraries]
                    fastq_id,fastqs,feature_types
                    {s1},{fqdir},gene expression
                    {s2},{fqdir},vdj-b
                    """.format(ref=params.ref_tx, vdj=params.ref_vdj, 
                               fqdir=params.fastq, s1=sg, s2=sb)
            with open(output[0], "w") as out:
                out.write(config)
    
rule cellranger_multi:
    input:
        fastq_ge = lambda wildcards: FASTQ_GE[wildcards.sample],
        transcriptom = os.path.join(TRANSCRIPTOM10X, "fasta/genome.fa"),
        vdj = os.path.join(VDJ, "fasta/regions.fa"),
        config = "{sample}.multi.config.csv"
    output:
        h5 = "{sample}/outs/per_sample_outs/{sample}/sample_filtered_feature_bc_matrix.h5",
    threads: 12
    log: "logs/{sample}.count.log"
    version: "0.1"
    params:
        extra=" --expect-cells=6000 --nosecondary", # loading cells 9000
        fastqs=FASTQ_DIR,
        transcriptom=TRANSCRIPTOM10X,
    run:
        # rm outputdir to solve "pipestance directory" error from cellranger
        shell("rm -rf {wildcards.sample}") # --localmem=48
        shell("cellranger multi --id={wildcards.sample} --csv={input.config} --localcores={threads} ")  

rule aggregate_list:
    input: expand("{sample}/outs/per_sample_outs/{sample}/sample_filtered_feature_bc_matrix.h5", sample=SAMPLES)
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
            out.write("sample_id,sample_outs,donor,origin\n")
            for s, b in zip(params.samples, params.ids):
                out.write(f"{s},{s}/outs/per_sample_outs/{s},{s},{s}\n")


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
