#this is snakemake script used for RNA-seq

from os.path import join, basename, dirname
from snakemake.utils import R
from snakemake.shell import shell
from pandas import read_table, read_excel, concat, ExcelWriter
from numpy import log10
import json


configfile: 'config.yml'
#Workding directory
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
SAMPLES = config['samples']['name'].split()
GROUP=config['samples']['group'].split()
TIME=config['samples']['time'].split()



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
SAMPLE_TPM ="gene_expression/gene_expression.TPM.txt",
SAMPLE_TPM_ANNO = "gene_expression/gene_expression.TPM.annotated.csv"
SAMPLE_DIFF_ANNO = "differential_expression/differential_expression_annotated.xls"
ROBJ_DESeq ="salmon/txi.salmon.RData"   
DESEQ_RES = "differential_expression/deseq2.results.txt"
ENRICHR_BAR = expand("differential_expression/Enrichr_{domain}_{types}/enrichr.reports.{types}.pdf",
                       domain=GO_DOMAIN, types=['all','up','down'])
GSEA_PRERANK = expand("differential_expression/GSEA_prerank_{domain}/gseapy.prerank.reports.csv", domain=GO_DOMAIN)
################## Rules #######################################


rule target:
    input: FASTQC, RAW_COUNTS, SAMPLE_TPM, SAMPLE_TPM_ANNO, ROBJ_DESeq,
           DESEQ_RES, GSEA_PRERANK
           #ENRICHR_BAR


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
        outdir=SALMON_INDEX_DIR，
        extra=" --gencode --type quasi -k 31"
    shell: 
        "salmon index {params.extra} -t {input} -i {params.outdir}"

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
        "salmon/{sample}/quant.genes.sf"
    threads: 8
    params:
        paired=PAIRED,
        index_dir=SALMON_INDEX_DIR,
        outdir="salmon/{sample}",
        extra_paried=" --incompatPrior 0  --numBootstraps 100 --seqBias --gcBias --writeUnmappedNames",
        extra_single=" --fldMean 250 --fldSD 25 --incompatPrior 0  --numBootstraps 100 --writeUnmappedNames"
    log: "logs/salmon/{sample}_salmons_quant.log"
    shell:
        "salmon quant -i {params.index_dir} -l A -1 {input.r1} -2 {input.r2} "
        " -g {input.gtf} -p {threads} -o {params.outdir} {params.extra_paried} &> {log} "
            
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
        res=DESEQ_RES,
    params:
        group=GROUP,#used for grouping each sample, to dectect degs.
    script:
        "scripts/runDESeq2.R"

rule gtf_extract:
    input: GTF_FILE
    output: 
        gene_anno=GTF_Genes,
        tx2gene = GTF_Trans
    run:
        lines_seen = set()
        tx2gene_seen = set()

        out1 = open(output.gene_anno, 'w') 
        out2 = open(output.tx2gene, 'w')
        out1.write("gene_id\tgene_name\tgene_type\tgene_status\n")
        out2.write("tx_id\tgene_id\n")
        for line in open(input[0]):
            if line.startswith('#'): continue
            attr = dict(item.strip().split(' ')  for item in line.split('\t')[8].strip('\n').split(';') if item)
            line_1st = attr['gene_id'].strip('\"')+'\t'+ attr['gene_name'].strip('\"')+'\t'
            line_2nd = attr['gene_type'].strip('\"')+'\t'+attr['gene_status'].strip('\"')+'\n'
            line_out = line_1st + line_2nd

            if line_out not in lines_seen:
                out1.write(line_out)
                lines_seen.add(line_out)
            if 'transcript_id' in attr:
                tx2gene_out = attr['transcript_id'].strip('\"')+'\t'+ attr['gene_id'].strip('\"')+'\n'
                if tx2gene_out not in tx2gene_seen:
                    out2.write(tx2gene_out)
                    tx2gene_seen.add(tx2gene_out)
        out1.close()
        out2.close()


rule annotate_genes:
    input: GTF_Genes, SAMPLE_TPM, DESEQ_RES
    output: 
        sample_anno=SAMPLE_TPM_ANNO,
        diff_anno=SAMPLE_DIFF_ANNO,
    params:
        log2fc=1,
        padj=0.05
    run:
        anno = read_table(input.annotation, index_col='gene_id')

        tpm = read_table(input.tpm, index_col=0)
        tpm.columns = ['TPM.'+col for col in tpm.columns]
        deseq = read_table(input.deseq, index_col=0)
        #merge results
        merge = concat([anno, deseq, tpm], axis=1, join='inner')
        merge.index.name ='gene_id'
        merge.to_csv(output.sample_anno)

        sig_deg = merge[(merge['log2FoldChange'].abs()> params.log2fc) & (merge['padj'] < params.padj) ]
        sig_deg = sig_deg.sort_values('padj',axis=0)
        sig_deg.loc[:,'up_down'] = sig_deg.log2FoldChange.apply(lambda x : 'up' if x > 0 else 'down' )
         

        sig_deg_up = sig_deg[sig_deg['up_down'] == 'up']
        sig_deg_dw = sig_deg[sig_deg['up_down'] == 'down'] 

        writer = ExcelWriter(output.diff_anno)
        sig_deg.to_excel(writer,sheet_name="sig-log2fc%s-padj%s"%(params.log2fc, params.padj))
        sig_deg_up.to_excel(writer,sheet_name="sig-up",)
        sig_deg_dw.to_excel(writer,sheet_name="sig-down",)
        writer.save()

rule GSEA_Enrichr:
    input: 
        diff=SAMPLE_DIFF_ANNO
    output:
        gsea=GSEA_PRERANK,
        enrichr=expand("differential_expression/Enrichr_{domain}_{types}/enrichr.reports.{types}.txt",
                       domain=GO_DOMAIN, types=['all','up','down'])
    params:
        log2fc=1,
        padj=0.05,
        go=GO_DOMAIN,
    run:
        writer = input.diff
        sig_deg = read_excel(writer,sheet_name="sig-fc%s-padj%s"%(params.log2fc, params.padj))
        sig_deg_up= read_excel(writer,sheet_name="sig-up",)
        sig_deg_dw= read_excel(writer,sheet_name="sig-down",)
        
        degs_sig = [deg.gene_name.squeeze().tolist() for deg in[sig_deg, sig_deg_up,sig_deg_dw]]

        sig_deg_gsea = sig_deg[['gene_name','log2FoldChange']]
        sig_deg_gsea_sort = sig_deg_gsea.sort_values('log2FoldChange',ascending=False)
        sig_deg_gsea_sort = sig_deg_gsea_sort.reset_index(drop=True)


        import gseapy as gp
        for domain in params.go:
            prerank = gp.prerank(rnk=sig_deg_gsea_sort, gene_sets=domain, pheno_pos='', pheno_neg='', min_size=15, max_size=500, 
                                 outdir='differential_expression/GSEA_prerank_'+domain)
        for domain in params.go:
            for glist, gl_type in zip(degs_sig, ['all','up','down']):
                enrichr = gp.enrichr(gene_list=glist, gene_sets=domain, 
                                     no_plot=True,  description=gl_type,
                                     outdir='differential_expression/Enrichr_%s_%s'%(domain, gl_type))


rule barplot:
    """Enrichr Results plotting"""
    input: "differential_expression/Enrichr_{domain}/enrichr.reports.{description}.txt"
    output:
        png="differential_expression/Enrichr_{domain}/enrichr.reports.{description}.png",
        pdf="differential_expression/Enrichr_{domain}/enrichr.reports.{description}.pdf",
    run:
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
        from matplotlib.figure import Figure
        d = read_table(input[0])
        d['logAP'] = -log10(d['Adjusted P-value']) 
        d = d.sort_values('logAP', ascending=False)
        dd = d.head(10).sort_values('logAP')
        fig = Figure(figsize=(12,6))
        canvas = FigureCanvas(fig)
        ax = fig.add_subplot(111)
        bar = dd.plot.barh(x='Term', y='logAP', color="salmon", alpha=0.75, edgecolor='none',fontsize=32, ax=ax)
        bar.set_xlabel("-log$_{10}$ Adjust P-value", fontsize=32)
        bar.set_ylabel("")
        bar.set_title(f.split("/")[-1].split(".")[-2],fontsize=32)
        bar.legend(loc=4)
        fig.savefig(output.png, bbox_inches='tight')
        fig.savefig(output.pdf, bbox_inches='tight')




    