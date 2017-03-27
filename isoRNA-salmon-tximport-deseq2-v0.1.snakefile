#this is snakemake script used for RNA-seq

from os.path import join, basename, dirname
from snakemake.utils import R
from snakemake.shell import shell
from pandas import read_table, read_excel, concat, ExcelWriter
import json


configfile: 'config.yml'
#Workding directory
workdir: config['workdir']

################### globals #############################################

# Full path to an uncompressed FASTA file with all chromosome sequences.
CDNA = config['cdna']

# Full path to a folder that holds all of your FASTQ files.
FASTQ_DIR = config['fastq_dir']
READ_LEN = congfig['read_length']
PAIRED = congfig['paired']
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

FASTQC = expand("qc/fastqc/{prefix}_fastqc.{suf}", prefix=SAMPLES,suf=['html','zip'])

SALMON_INDEX = expand(SALMON_INDEX_DIR+"/{prefix}.bin", prefix=['hash','rsd','sa','txpinfo'])
SALMON_QUANT_Trans = expand("salmon/{sample}/quant.sf", sample=SAMPLES)
SALMON_QUANT_Genes = expand("salmon/{sample}/quant.genes.sf", sample=SAMPLES)

RAW_COUNTS ="counts/sample.raw.counts.txt"
SAMPLE_TPM ="gene_expression/gene_expression.TPM.txt",
SAMPLE_TPM_ANNO = "gene_expression/gene_expression.TPM.annotated.csv"
SAMPLE_DIFF_ANNO = "differential_expression/differential_expression_annotated.xls"
ROBJ_DESeq ="salmon/txi.salmon.RData"   
DESEQ_RES = "differential_expression/deseq2.results.txt"
################## Rules #######################################


rule target:
    input: FASTQC, RAW_COUNTS, SAMPLE_TPM, SAMPLE_TPM_ANNO, ROBJ_DESeq,DESEQ_RES


rule fastqc:
    input:
        "fastq/{prefix}.fastq.gz",
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
        outdir=SALMON_INDEX_DIR
    shell: 
        "salmon index -t {input} -i {params.outdir} --type quasi -k 31"

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
    run:
        if paired:
            #paried end
            shell("""
                  salmon quant -i {params.index_dir} -l A -1 {input.r1} -2 {input.r2} -g {input.gtf} -p {threads} -o {params.outdir} {params.extra_paried} &> {log}
                  """)
        else:
            # single end
            shell("""
                 salmon quant -i {params.index_dir}  -l A -r {input.r1} -p {threads} -g {input.gtf} -o {params.outdir} {params.extra_signle} &> {log}
                 """)
            
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
    threads: 8
    run:
        R("""
          library(tximport)
          library(readr)
          
          # hg38, create a tx2gene.txt table
          #library(EnsDb.Hsapiens.v86)
          #edb <- EnsDb.Hsapiens.v86
          #Tx.ensemble <- transcripts(edb,
                                     columns = c("tx_id", "gene_id", "gene_name"),
                                     return.type = "DataFrame")
          #nrow(Tx.ensemble)
          #tx2gene<- Tx.ensemble[,c(1,2)]
          #write.table(tx2gene, "{output.tx2gene}", col.names = F, row.names = F, sep="\t", quote=F)

          tx2gene <- read.table("{input.tx2gene}", header= T, sep="\t", stringsAsFactors = F)
          samples <- unlist(strsplit("{params.ids}",","))

          salmon.files <- file.path('salmon',samples, "quant.sf")
          names(salmon.files) <- samples
          all(file.exists(salmon.files))

          #aggregate transcripts to genes 
          txi.salmon <- tximport(salmon.files, type = "salmon", 
                                 tx2gene = tx2gene, reader = read_tsv,
                                 countsFromAbundance = "lengthScaledTPM")

          salmon.counts<- txi.salmon$counts
          salmon.counts<- as.data.frame(salmon.counts)
          #salmon.counts$gene_name<- rownames(salmon.counts)
          write.table(salmon.counts, "{output.counts}", sep="\t", quote=F)

          salmon.TPM<- txi.salmon$abundance
          salmon.TPM<- as.data.frame(salmon.TPM)
          #salmon.TPM$gene_name<- rownames(salmon.TPM)
          write.table(salmon.TPM, "{output.tpm}", sep="\t", quote=F)
          save.image(file="{output.image}")
          """)

rule deseq2:
    input: 
        image=ROBJ_DESeq,
    output: 
        res=DESEQ_RES,
    params:
        group=GROUP,#used for grouping each sample, to dectect degs.
    run:
        R("""
          library(DESeq2)
          load("{input.image}")

          #assign each sample to differrent group.
          group=unlist(strsplit("{group.group}", " "))
          sampleTable <- data.frame(condition = factor(group))
          rownames(sampleTable) <- colnames(txi.salmon$counts)

          #run DESeq2
          dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~condition)
          dds <- DESeq(dds)
          res <- results(dds,addMLE=TRUE)
          resOrdered <- res[order(res$padj),]
          resOrdered = as.data.frame(resOrdered)
          write.table(resOrdered, file="{output.res}",sep="\t")

          """)

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
            tx2gene_out = attr['transcript_id'].strip('\"')+'\t'+ attr['gene_id'].strip('\"')+'\n'
            if line_out not in lines_seen:
                out1.write(line_out)
                lines_seen.add(line_out)
            if tx2gene_out not in tx2gene_seen:
                out2.write(tx2gene_out)
                tx2gene_seen.add(tx2gene_out)
        out1.close()
        out2.close()


rule annotation:
    input: 
        annotation=GTF_Genes, 
        tpm=SAMPLE_TPM, 
        deseq=DESEQ_RES
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

        sig_deg = merge[(merge['log2FoldChange'].abs()> params.log2fc) & (merge1['padj'] < params.padj) ]
        sig_deg = sig_deg.sort_values('padj',axis=0)
        sig_deg['up_down'] = sig_deg.log2FoldChange.apply(lambda x : 'up' if x > 0 else 'down' )
         

        sig_deg_up = sig_deg[sig_deg['up_down'] == 'up']
        sig_deg_dw = sig_deg[sig_deg['up_down'] == 'down'] 

        writer = ExcelWriter(output.diff_anno)
        sig_deg.to_excel(writer,sheet_name="sig-fc%s-padj%s"%(params.log2fc, params.padj))
        sig_deg_up.to_excel(writer,sheet_name="sig-up",)
        sig_deg_dw.to_excel(writer,sheet_name="sig-down",)
        writer.save()

rule GSEA_prerank:
    input: 
        go=GO_DOMAIN,
        diff=SAMPLE_DIFF_ANNO
    output:
    params:
        log2fc=1,
        padj=0.05
    run:
        writer = input.diff
        sig_deg = read_excel(writer,sheet_name="sig-fc%s-padj%s"%(params.log2fc, params.padj))
        sig_deg_up= read_excel(writer,sheet_name="sig-up",)
        sig_deg_dw= read_excel(writer,sheet_name="sig-down",)

        deg_all = sig_deg.gene_name.squeeze().tolist()
        deg_up = sig_deg_up.gene_name.squeeze().tolist()
        deg_dw = sig_deg_dw.gene_name.squeeze().tolist()


        sig_deg_gsea = sig_deg[['gene_name','log2FoldChange']]
        sig_deg_gsea_sort = sig_deg_gsea.sort_values('log2FoldChange',ascending=False)
        sig_deg_gsea_sort = sig_deg_gsea_sort.reset_index(drop=True)


        import gseapy as gp
        for domain in input.go:
            prerank = gp.prerank(rnk=sig_deg_gsea_sort, gene_sets=domain, pheno_pos='', pheno_neg='', min_size=15, max_size=500, 
                                 outdir='GSEA_prerank_'+domain)
        for domain in input.go:
            for glist, gl_type in zip([deg_all, deg_up, deg_dw],['all','up','down']):
                enrichr = gp.enrichr(glist=glist, gene_sets=domain, no_plot=True, outdir='Enrichr_%s_%s'%(domain, gl_type))





    