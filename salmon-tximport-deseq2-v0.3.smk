from os.path import join, isfile
from itertools import combinations

include: "rules/common.smk"

# configfile: 'config.yml'
workdir: config['workdir']



################### globals #############################################

# Full path to an uncompressed FASTA file with all chromosome sequences.
CDNA = config['cdna']

# Full path to a folder that holds all of your FASTQ files.
FASTQ_DIR = config['fastq_dir']
READ_LEN = config['read_length']
PAIRED = config['paired']
# Full path to a Genome.
GENOME = config['genome']
#CDNA =           join(GENOME,"gencode.v25.transcripts.fa")
# genome sequence
FASTA_REF =      config['dna']
# index_dir
SALMON_INDEX_DIR=config['salmon_index']
# index basename
INDEX_PREFIX = 'hg38'
# gtf
GTF_FILE =       config['gtf']
GTF_Genes =      GTF_FILE.rstrip(".gtf")+".extracted.genes.annotation.txt"
GTF_Trans =      GTF_FILE.rstrip(".gtf")+".extracted.transx2gene.txt"
############ Samples ##################
# A Snakemake regular expression matching the forward mate FASTQ files.
# the part in curly brackets {} will be saved, so the variable SAMPLES
# is a list of strings #['Sample1','Sample2'].

#notice that SAMPLES, has a trailing comma.
#you must include this trailing comma, or else the code won’t work correctly.

#SAMPLES, = glob_wildcards(join(FASTQ_DIR, '{sample, SRR[^/]+}_R1.fastq.gz'))
SAMPLES,SAMPLES_ALIAS,GROUP,TIME = parse_samples(config['sample_meta'])
uGroup=unique(GROUP)

# Patterns for the 1st mate and the 2nd mate using the 'sample' wildcard.
PATTERN_R1 = config['read_pattern']['r1']
PATTERN_R2 = config['read_pattern']['r2']

# dirs
DIRS = ['qc','mapped','alternative_splicing', 'gene_expression',
        'differential_expression','logs','temp']
# go domain
GO_DOMAIN = config['enrichr_library']
# cutoff
LOG2FC = config['log2fc']
FDR = config['fdr']
########### Target output files #################
SALMON_INDEX = expand(SALMON_INDEX_DIR+"/{prefix}.bin", prefix=['pos','mphf'])
SALMON_QUANT_Trans = expand("salmon/{sample}/quant.sf", sample=SAMPLES)
SALMON_QUANT_Genes = expand("salmon/{sample}/quant.genes.sf", sample=SAMPLES)

# RAW_COUNTS ="quant/sample.raw.counts.txt"
SAMPLE_TPM_ANNO = "quant/gene_expression.TPM.annotated.csv"
SAMPLE_TXTPM_ANNO ="quant/transcripts_expression.TPM.annotated.csv"

SALMON_OBJ= "quant/txi.salmon.RData"
DESEQ_DDS = "quant/deseq2.dds.RData"
DESEQ_NTD = "quant/deseq2.ntd.Rdata"
DESEQ_RES = ["differential_expression/diff_{t}_vs_{c}/diff_{t}_vs_{c}_results.txt".format(t=j, c=i)
             for i, j in combinations(uGroup, 2)]
DESEQ_ANNO = [res.replace(".txt", ".annotated.xls") for res in DESEQ_RES]
DESEQ_HEATMAP = ["differential_expression/diff_{t}_vs_{c}/diff_{t}_vs_{c}_all.degs.pdf".format(t=j, c=i)
               for i, j in combinations(uGroup, 2)]

GSEA_RES=["differential_expression/GO/GSEA_{treat}_vs_{ctrl}/%s/gseapy.gsea.gene_set.report.csv"%domain for domain in GO_DOMAIN]
GSEA_FINAL=["differential_expression/GO/GSEA_%s_vs_%s/KEGG_2016/gseapy.gsea.gene_set.report.csv"%(j, i) for i, j in combinations(uGroup, 2)]
#Enrichr = ["GO/Enrichr_{treat}_vs_{ctrl}/{domain}_{types}/{domain}.{type}.enrichr.reports.txt",type=["all","up","down"]
#GSEA_OUT = [ GSEA_RES.format(treat=uGroup[i], ctrl=uGroup[i-1]) for i in range(1, len(uGroup))]
################## Rules #######################################


rule target:
    input: DESEQ_DDS,  #DESEQ_ANNO,
           SAMPLE_TPM_ANNO, SAMPLE_TXTPM_ANNO,
           #DESEQ_RES, DESEQ_HEATMAP, GSEA_FINAL

rule salmon_index:
    input: CDNA
    output: 
        bin = SALMON_INDEX, 
        #fa = join(SALMON_INDEX_DIR,"ref_k31_fixed.fa")
    threads: 8
    params:
        genome_dir=config['genome'],
        cdna=CDNA,
        outdir=SALMON_INDEX_DIR,
        extra=" --gencode "
    shell:
        "salmon index {params.extra} -t {input} -i {params.outdir}"
#"salmon index -i /genome/{params.outdir} -t /genome/{params.cdna} {params.extra}"
########## notes on salmon quant ###################################################################

###   <LIBTYPE>
### A  (automatically infer the library type)
### IU (an unstranded paired-end library where the reads face each other)
### SF (a stranded single-end protocol where the reads come from the forward strand)
### more types visit: http://salmon.readthedocs.io/en/latest/salmon.html#quasi-mapping-based-mode-including-lightweight-alignment
"""
salmon quant

    -l (depends on the lib type, ISR for truseq stranded, equivalent to tophat -fr-firststrand)
    -p (the number of available cores on the instance)
    -g (the gene level counts will be part of output)
    --incompatPrior 0 (we don’t want reads incompatible with the library type)
    --fldMean 250 (for single-ended reads only, kallisto default, can change if there is info about lib prep)
    --fldSD 25 (for single-ended reads only, kallisto default, can change if there is info about lib prep)
    --numBootstraps 100 (maybe good for samples without technical replicates)
    --seqBias (this option not for single-end )
    --gcBias (this option not for single-end )
    --writeUnmappedNames

    Note:
       Choose not to use --useVBOpt based on documentation and this link
       https://groups.google.com/forum/#!topic/sailfish-users/-LBZD4aoJSc
       The behavior will be more like Kallisto and RSEM instead of BitSeq
"""
############ salmon quant start ####################################################################
rule salmon_quant:
    input:
        index=SALMON_INDEX,
        #index_dir=SALMON_INDEX_DIR,
        gtf=GTF_FILE,
        r1=join(FASTQ_DIR, PATTERN_R1),
        r2=join(FASTQ_DIR, PATTERN_R2)
    output:
        sf="salmon/{sample}/quant.sf",
        #outdir=directory("salmon/{sample}") # Directories as outputs
    threads: 8
    params:
        r1=join(FASTQ_DIR, PATTERN_R1),
        r2=join(FASTQ_DIR, PATTERN_R2),
        index_dir=SALMON_INDEX_DIR,
        outdir=join(config['workdir'],"salmon/{sample}"),
        extra_paried=" --incompatPrior 0  --numBootstraps 100 --seqBias --gcBias --writeUnmappedNames",
        #extra_single=" --fldMean 250 --fldSD 25 --incompatPrior 0  --numBootstraps 100 --writeUnmappedNames"
    log: "logs/salmon/{sample}_salmons_quant.log"
    shell:
        "salmon quant -l A -i {params.index_dir} -1 {params.r1} -2 {params.r2} "
        "-p {threads} -o {params.outdir} {params.extra_paried} &> {log}"
rule tximport:
    '''used for kallisto, Salmon, Sailfish, and RSEM. see:
    http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
    ###
    good tutorial to look:
    https://github.com/crazyhottommy/RNA-seq-analysis/blob/master/salmon_kalliso_STAR_compare.md#counts-versus-tpmrpkmfpkm
    '''
    input:
        quant=expand("salmon/{sample}/quant.sf", sample=SAMPLES),
        tx2gene=GTF_Trans
    output:
        tpm=temp("quant/gene_expression.TPM.txt"),
        txtpm=temp("quant/transcripts_expression.TPM.txt"),
        #counts=RAW_COUNTS,
        image=SALMON_OBJ
    params:
        ids =",".join(SAMPLES)
    threads: 1
    script:
        "scripts/runTximport.R"


# rule salmon_quantmerge:
#     input:
#         expand("salmon/{sample}/quant.sf", sample=SAMPLES),
#         expand("salmon/{sample}/quant.genes.sf", sample=SAMPLES)
#     output:
#         tpm=SAMPLE_TPM,
#         txtpm=SAMPLE_TXTPM,
#         counts=RAW_COUNTS
#     params:
#         ids =" ".join(SAMPLES)
#     shell:
#         """salmon quantmerge --quants {params.ids}  --genes -o {output.tpm}
#            salmon quantmerge --quants {params.ids}  --genes -c numreads -o {output.counts}
#            salmon quantmerge --quants {params.ids}  -o {output.txtpm}
#         """

        
rule deseq2:
    input:
        image="quant/txi.salmon.RData",#SALMON_TPM
    output:
        res=DESEQ_RES,
        ddsimage="quant/deseq2.dds.RData", #DESEQ_DDS
        ntdimage="quant/deseq2.ntd.RData", #DESEQ_NTD
    params:
        group=" ".join(GROUP),#used for grouping each sample, to dectect degs.
        time=" ".join(TIME),
        alias=" ".join(SAMPLES_ALIAS)
    threads: 8
    script:
        "scripts/runDESeq2.R"

rule gtf_extract:
    input: GTF_FILE
    output:
        gene_anno=GTF_Genes,
        tx2gene = GTF_Trans
    script:
        "scripts/extractGTF.py"

rule anno_DEGs:
    input: 
        "differential_expression/diff_{treat}_vs_{ctrl}/diff_{treat}_vs_{ctrl}_results.txt"
    output: 
        "differential_expression/diff_{treat}_vs_{ctrl}/diff_{treat}_vs_{ctrl}_results.annotated.xls"
    params:
        gene_anno=GTF_Genes,
        tpm=SAMPLE_TPM,
        alias=SAMPLES_ALIAS,
        samples=SAMPLES,
        group=GROUP,
        treat="{treat}",
        ctrl="{ctrl}",
        log2fc=LOG2FC,
        padj=FDR
    script:
        "scripts/annotateDEGs.py"

rule pheatmap_degs:
    input:
        degstab="differential_expression/diff_{treat}_vs_{ctrl}/diff_{treat}_vs_{ctrl}_results.txt",
        image="quant/deseq2.ntd.RData"
    output:
        "differential_expression/diff_{treat}_vs_{ctrl}/diff_{treat}_vs_{ctrl}_all.degs.pdf",
        "differential_expression/diff_{treat}_vs_{ctrl}/diff_{treat}_vs_{ctrl}_top20genes.pdf"
    params:
        treat="{treat}",
        ctrl="{ctrl}",
        padj=FDR,
        topgene=20,
    script:
        "scripts/pheatmapDEGs.R"

rule anno_samples:
    input:
        GTF_Genes,
        GTF_Trans,
        "quant/gene_expression.TPM.txt",
        "quant/transcripts_expression.TPM.txt",
    output:
        "quant/gene_expression.TPM.annotated.csv",
        "quant/transcripts_expression.TPM.annotated.csv",
    params:
        group=GROUP,
        alias=SAMPLES_ALIAS,
        samples=SAMPLES,
    script:
        "scripts/annotateTPMs.py"


rule GSEA_Enrichr:
    input:
        "differential_expression/diff_{treat}_vs_{ctrl}/diff_{treat}_vs_{ctrl}_results.annotated.xls"
    output:
        GSEA_RES,
        #directory("differential_expression/GO/GSEA_{treat}_vs_{ctrl}",
        directory("differential_expression/GO/Enrichr_{treat}_vs_{ctrl}")
    params:
        treat="{treat}",
        ctrl="{ctrl}",
        log2fc=LOG2FC,
        padj=FDR,
        go=GO_DOMAIN,
    script:
        "scripts/gseaEnrichr.py"
