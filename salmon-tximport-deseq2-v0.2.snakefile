from os.path import join, isfile

configfile: 'config.yml'
#Workding directory
workdir: config['workdir']

def unique(seq):
    """Remove duplicates from a list in Python while preserving order.
    :param seq: a python list object.
    :return: a list without duplicates while preserving order.
    """    
    seen = set()
    seen_add = seen.add

    return [x for x in seq if x not in seen and not seen_add(x)]


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
FASTA_REF =      config['fasta']
# index_dir
SALMON_INDEX_DIR=config['index_dir']
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

if isfile(config['samples']['coldata']):
    SAMPLES=[]
    SAMPLES_ALIAS=[]
    GROUP=[]
    TIME=[]
    with open(config['samples']['coldata']) as f:
        lines = f.readlines()
    for line in lines:
        item = line.rstrip("\n").split(" ")
        SAMPLES.append(item[0])
        SAMPLES_ALIAS.append(item[1])
        GROUP.append(item[2])
        TIME.append(item[3])
else:
    SAMPLES = config['samples']['name'].split()
    SAMPLES_ALIAS = config['samples']['alias'].split()
    GROUP=config['samples']['group'].split()
    TIME=config['samples']['time'].split()

uGroup=unique(GROUP)

# Patterns for the 1st mate and the 2nd mate using the 'sample' wildcard.
#PATTERN_R1 = '{sample}_R1.fastq.gz'
#PATTERN_R2 = '{sample}_R2.fastq.gz'
PATTERN_R1 = config['read_pattern']['r1']
PATTERN_R2 = config['read_pattern']['r2']

# dirs 
DIRS = ['qc','mapped','counts','alternative_splicing', 'gene_expression',
        'differential_expression','logs','temp']
# go domain
GO_DOMAIN = ['GO_Cellular_Component_2015','GO_Molecular_Function_2015',
             'GO_Biological_Process_2015','Human_Phenotype_Ontology',
             'MSigDB_Oncogenic_Signatures','WikiPathways_2016',
             'KEGG_2016']

########### Target output files #################

FASTQC = expand("qc/fastqc/{sample}_{read}_fastqc.{suf}", sample=SAMPLES,read=['R1','R2'],suf=['html','zip'])

SALMON_INDEX = expand(SALMON_INDEX_DIR+"/{prefix}.bin", prefix=['hash','rsd','sa','txpInfo'])
SALMON_QUANT_Trans = expand("salmon/{sample}/quant.sf", sample=SAMPLES)
SALMON_QUANT_Genes = expand("salmon/{sample}/quant.genes.sf", sample=SAMPLES)

RAW_COUNTS ="counts/sample.raw.counts.txt"
SAMPLE_TPM ="gene_expression/gene_expression.TPM.txt"
SAMPLE_TPM_ANNO = "gene_expression/gene_expression.TPM.annotated.csv"
ROBJ_DESeq ="salmon/txi.salmon.RData"   
DESEQ_RES = ["differential_expression/diff_%s_vs_%s_results.txt"%(uGroup[i], uGroup[i-1]) for i in range(1, len(uGroup))]
DESEQ_ANNO = [res.replace(".txt", ".annotated.xls") for res in DESEQ_RES]
GSEA_RES=["GO/GSEA_{treat}_vs_{ctrl}/%s/gseapy.prerank.gene_sets.report.csv"%domain for domain in GO_DOMAIN]
#Enrichr = ["GO/Enrichr_{treat}_vs_{ctrl}/{domain}_{types}/{domain}_{type}_enrichr.reports.txt",type=["all","up","down"]
#GSEA_OUT = [ GSEA_RES.format(treat=uGroup[i], ctrl=uGroup[i-1]) for i in range(1, len(uGroup))]
################## Rules #######################################


rule target:
    input: FASTQC, RAW_COUNTS, ROBJ_DESeq, DESEQ_ANNO, SAMPLE_TPM_ANNO,
           DESEQ_RES, "GO/GSEA_HDE_vs_Ctrl/KEGG_2016/gseapy.prerank.gene_sets.report.csv"


rule fastqc:
    input:
        #lambda wildcards: join(FASTQ_DIR, wildcards.prefix +".fastq.gz")
        join(FASTQ_DIR, "{prefix}.fastq.gz")
    output:
        "qc/fastqc/{prefix}_fastqc.html",
        "qc/fastqc/{prefix}_fastqc.zip",
    params:
        outdir="qc/fastqc"

    shell: "fastqc -o {params.outdir} {input}"

rule salmon_index:
    input: CDNA
    output: SALMON_INDEX
    threads: 8
    params: 
        genome_dir=config['genome'],
        cdna=CDNA.split("/")[-1],
        outdir=SALMON_INDEX_DIR.split("/")[-1],
        extra=" --gencode --type quasi -k 31"
    shell: 
        #"salmon index {params.extra} -t {input} -i {params.outdir}"
        "docker run  -v {params.genome_dir}:/genome 435326e03520  "
        "salmon index -i /genome/{params.outdir} -t /genome/{params.cdna} {params.extra}"
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
        gtf=GTF_FILE,
        r1=join(FASTQ_DIR, PATTERN_R1), 
        r2=join(FASTQ_DIR, PATTERN_R2)
    output: 
        "salmon/{sample}/quant.sf",
    threads: 8
    params:
        r1=join(FASTQ_DIR.split("/")[-1], PATTERN_R1), 
        r2=join(FASTQ_DIR.split("/")[-1], PATTERN_R2),
        workdir=config['workdir'],
        index_dir=SALMON_INDEX_DIR,
        outdir="salmon/{sample}",
        extra_paried=" --incompatPrior 0  --numBootstraps 100 --seqBias --gcBias --writeUnmappedNames",
        #extra_single=" --fldMean 250 --fldSD 25 --incompatPrior 0  --numBootstraps 100 --writeUnmappedNames"
    log: "logs/salmon/{sample}_salmons_quant.log"
    shell:        
        "docker run -v {params.index_dir}:/index  -v {params.workdir}:/data 435326e03520 "
        "salmon quant -i /index -1 /data/{params.r1} -2 /data/{params.r2} "
        "-l A -o /data/{params.outdir} {params.extra_paried}"
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
        tpm=SAMPLE_TPM,
        counts=RAW_COUNTS,
        image=ROBJ_DESeq       
    params:
        ids =",".join(SAMPLES)
    threads: 1
    script:
        "scripts/runTximport.R"


rule deseq2:
    input: 
        image=ROBJ_DESeq,
    output: 
        res=DESEQ_RES
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

rule annotate_diGenes:
    input: 
        "differential_expression/diff_{treat}_vs_{ctrl}_results.txt"
    output: 
        "differential_expression/diff_{treat}_vs_{ctrl}_results.annotated.xls"

    params:
        gene_anno=GTF_Genes,
        tpm=SAMPLE_TPM,
        alias=SAMPLES_ALIAS,
        samples=SAMPLES,
        group=GROUP,
        log2fc=1,
        padj=0.05
    script:
        "scripts/annotateDEGs.py"


rule annotate_samples:
    input: SAMPLE_TPM
    output: SAMPLE_TPM_ANNO
    params:
        gene_anno=GTF_Genes,
        group=GROUP,
        alias=SAMPLES_ALIAS,
        samples=SAMPLES,
    run:
        from pandas import read_table, concat

        anno = read_table(params.gene_anno, index_col='gene_id')
        tpm = read_table(input[0], index_col=0)
        tpm = tpm[params.samples]
        tpm.columns = ["TPM.%s.%s"%(g,a) for a, g, s in zip(params.alias, params.group, params.samples)]

        merge = concat([anno, tpm], axis=1, join='inner')
        merge.index.name ='gene_id'
        merge.to_csv(output[0])


rule GSEA_Enrichr:
    input: 
        "differential_expression/diff_{treat}_vs_{ctrl}_results.annotated.xls"

    output:
        #gsea="GO/GSEA_{treat}_vs_{ctrl}/KEGG_2016/gseapy.prerank.gene_sets.report.csv"
        enrichr=
        #gsea=GSEA_RES 
        GSEA_RES   
    threads: 8
    params:
        log2fc=1,
        padj=0.05,
        go=GO_DOMAIN,
    script:
        "scripts/gseaEnrichr.py"


    