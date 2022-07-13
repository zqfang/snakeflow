# snakemake pipeline for 10x genomics single cell CITE-seq
# @author: Zhuoqing Fang
# @email: maxzqfang@stanford.edu
# @date: May 7th, 2020
import glob, re, os

WORKSPACE = "/data/bases/fangzq/20220412_COVID_Prabhu"
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
HASHTAG = os.path.join(WORKSPACE, "TotalseqA.barcodes.csv")
FASTQ_DIR = os.path.join(WORKSPACE, "raw_fastq")

# automatic search sample names using regex
pat = re.compile(r"(.+)_(S\d+)_(L\d+)_(\w\d)_001.fastq.gz")
# name, capture, lane, reads = pat.search(f).groups()
files = glob.glob(os.path.join(FASTQ_DIR, "*fastq.gz"))
names = [pat.search(f).groups() for f in files]
# remove duplicate entries
samples = sorted(list(set([f[0].split("/")[-1] for f in names])))


PROJECT = "AsymptomaticCOVID"
SAMPLES = [s[:-3] for s in samples]
SAMPLES_HASH=samples

### get fastqs for Gene Expression and Antibody Capture
FASTQ_GE = {s:[os.path.join(FASTQ_DIR, s, f) for f in os.listdir(os.path.join(FASTQ_DIR, s)) if pat.search(f)] for s in SAMPLES}
FASTQ_FB = {s[:-3]: [os.path.join(FASTQ_DIR, f) for f in os.listdir(FASTQ_DIR) if pat.search(f) and f.find("s") > 0 ]   for s in SAMPLES_HASH}


# cellranger
TRANSCRIPTOM10X = "/home/fangzq/genome/10XRefdata/refdata-gex-GRCh38-2020-A" 

# # vireo
# SNV = "/home/fangzq/genome"
# KNOWN_SNV = os.path.join(SNV, "genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf")
# DONOR_NUM_PER_SAMPLE = 3 # per 10x run samples/patients ...


#### Target output files #####
MTX =  expand("{sample}_cnt/outs/filtered_feature_bc_matrix/matrix.mtx.gz",sample=SAMPLES)
BARCODE = expand("{sample}_cnt/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", sample=SAMPLES)
FEATURES = expand("{sample}_cnt/outs/filtered_feature_bc_matrix/features.tsv.gz",sample=SAMPLES)
AGGR="{project}_aggr/outs/count/filtered_feature_bc_matrix.h5".format(project = PROJECT)
BAM = expand("{sample}_cnt/outs/possorted_genome_bam.bam", sample=SAMPLES)
VIREO = expand("{sample}_vireo/results.demux/donor_ids.tsv", sample=SAMPLES)

#### Rules ########
rule target: 
    input: MTX, BARCODE, FEATURES, BAM, AGGR, #VIREO


# rule prepare_count:
#     output: 
#         expand("{sample}.libraries.csv", sample=SAMPLES)
#     params:
#         samples=SAMPLES,
#         fastqs=FASTQ_DIR,
#         capture="HASH", # tag name represent which file is hashtag fastq
#         files=files,
#     run:
#         for sample in params.samples:
#             hashtag = params.files[sample]
#             out = sample+".libraries.csv"
#             with open(out, 'w') as w:
#                 w.write("fastqs,sample,library_type\n")
#                 for k, v in hashtag.items():
#                     if k == params.capture:
#                         w.write("%s,%s,Antibody Capture\n"%(params.fastqs, sample))
#                     else:
#                         w.write("%s,%s,Gene Expression\n"%(params.fastqs, sample))



# snakemake creates the folder if it does not exisit. but cellranger fails if the folder is 
# not created by itself. 

# trick 1: rm {sample}_cnt/_lock file and re-run pipeline, it will be OK.
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/troubleshooting
# trick 2: rm output folder before cellranger run
rule citeseq_count:
    input:
        fastq_ge = lambda wildcards: FASTQ_GE[wildcards.sample],
        fastq_fb = lambda wildcards: FASTQ_FB[wildcards.sample],
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
    threads: 16
    log: "logs/{sample}.count.log"
    version: "0.1"
    params:
        extra=" --expect-cells=20000 --nosecondary", # --expect-cells=12000
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
        ids = SAMPLES, # only used for loupe browser
        # condition = [lane]*3 #CONDITION, ### ->>>>
    run:
        #run = [s[:4] for s in SAMPLES]
        with open(output[0],'w') as out:
            # only add 'batch' when you need Chemistry Batch Correction
            # batch is optional, see docs if you have bath info
            # https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate
            out.write("library_id,molecule_h5,run\n")
            for s, b in zip(params.samples, params.ids):
                out.write(f"{s},{s}_cnt/outs/molecule_info.h5,{b}\n")




rule citeseq_aggr:
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
        shell("cellranger aggr --id={params.project}_aggr --nosecondary " +\ 
              "--csv={input.aggr} " +\
              "--localcores={threads} "+\
              "--normalize=none > {log}") # if =mapped, then normalizes for effective sequencing depth by subsampling the reads.

# demuxlet --sam ${sample}citeseq/outs/possorted_genome_bam.bam  \
#          --vcf /home/fangzq/projects/merged.sorted.vcf \
#          --out ${sample}.demuxlet \
#          --min-mac 1 --field GT --alpha 0  \
#          --group-list /home/fangzq/projects/SFGF-Peltz-YG-13709/${sample}citeseq/outs/filtered_feature_bc_matrix/barcodes.tsv.gz

# call single cell SNP

# a single BAM/SAM file, e.g., from cellranger, 
# a list of cell barcodes, a VCF file for common SNPs. 
# This mode is recommended comparing to mode 2, if a list of common SNP is known, e.g., human (see Candidate SNPs below)
# rule cell_snp:
#     input:
#         bam="{sample}_cnt/outs/possorted_genome_bam.bam",
#         barcode="{sample}_cnt/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
#         vcf=KNOWN_SNV,
#     output:
#         "{sample}_vireo/CELL_DATA/cellSNP.cells.vcf.gz",
#         "{sample}_vireo/CELL_DATA/cellSNP.base.vcf.gz"
#     threads: 12
#     log: "logs/{sample}.cellsnp.log"
#     shell:
#         "cellSNP -s {input.bam} -b {input.barcode} -R {input.vcf} "
#         "-p {threads} --minMAF 0.1 --minCOUNT 20 "
#         "-O {wildcards.sample}_vireo/CELL_DATA > {log}"

# # # deconvolution without known genotype
# # You have to know how many individals, e.g. n_donor=3
# rule vireo_demux:
#     input: 
#         "{sample}_vireo/CELL_DATA/cellSNP.cells.vcf.gz",
#         "{sample}_vireo/CELL_DATA/cellSNP.base.vcf.gz",
#     output: 
#         "{sample}_vireo/results.demux/donor_ids.tsv",
#         "{sample}_vireo/results.demux/GT_donors.vireo.vcf.gz",
#     params:
#         donor_num = DONOR_NUM_PER_SAMPLE,
#         CELL_DATA = "{sample}_vireo/CELL_DATA"
#         # with known genotype, add below argument
#         # known_gt = "-d {DONOR_GT_FILE} -t GT" # DONOR_GT_FILE(vcf) obtained from input samples
#     log: "logs/{sample}.vireo.log"
#     shell:
#         "vireo -c {params.CELL_DATA} -N {params.donor_num} -o {wildcards.sample}_vireo/results.demux > {log}"

# # Mode 3: pileup a list of SNPs for one or multiple BAM/SAM files
# # rule cell_SNP3:
# #     input:
# #         bam=expand("{sample}_cnt/outs/possorted_genome_bam.bam", sample=SAMPLES)
# #         barcode=expand("{sample}_cnt/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", sample=SAMPLES)
# #         vcf=
# #     output:
# #     params:
# #         extra="--UMItag None" # only include this if smart-seq or bulk rna-seq
# #         bams=",".join(expand("{sample}_cnt/outs/possorted_genome_bam.bam", sample=SAMPLES)),
# #         ids=",".join(SAMPLES)
# #     log:
# #     shell:
# #         "cellSNP -s {params.bams} -I {params.ids} -o $OUT_FILE -R {input.vcf} -p 20"

# # https://github.com/10XGenomics/vartrix
# rule vatrix:
#     input:
#         bam="{sample}_cnt/outs/possorted_genome_bam.bam",
#         barcode="{sample}_cnt/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
#         vcf=KNOWN_SNV,
#         genome=os.path.join(TRANSCRIPTOM10X, "fasta/genome.fa")
#     output:
#         "{sample}.vatrix/matrix.mtx"
#     threads: 12
#     shell:
#         "vartrix -v {input.vcf} -b {input.bam} -f {input.genome} -c {input.barcode} "
#         "--threads {threads} -o {wildcards.sample}.vatrix"

# rule vcf2SNVloci:
#     input: KNOWN_SNV
#     output: "SNV.loci.txt"
#     shell:
#         "awk '{{print $1,$2}}' {input} | sed -i 's/\s/:/g' > {output}"
