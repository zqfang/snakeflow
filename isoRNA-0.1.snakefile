from os.path import join
from snakemake.utils import R
from snakemake.shell import shell

############# Globals ######################################

# Full path to a folder that holds all of your FASTQ files.
FASTQ_DIR = 'fastq'
READ_LEN = 150
# Full path to a Genome.
GENOME = "/Users/bioninja/genome"
# genome sequence
FASTA_REF =      GENOME+"/hg19_ucsc_fasta/hg19.fa"
# index dir
STAR_REFDIR =    GENOME+"/starIndex/"
HISAT2_REFDIR=   GENOME+"/hisat2Indexes/hg19_ucsc"
# index basename
INDEX_PREFIX =   'hg19'
# gtf
GTF_FILE =       GENOME+"/gtf/gencode.v19.annotation.gtf"
GTF_EXTRACT =    GTF_FILE.rstrip(".gtf")+".extracted.txt"

# other
RRNA =           GENOME+"/Annotation/GRCm38_rRNA.list"
RSEQC_ANNO =     GENOME+"/rseqc_ann/hg19.HouseKeepingGenes.bed"

####### Tools dir#################
RMATS_DIR="/picb/external/isotex/program/rMATS.3.2.5/RNASeq-MATS.py"


############ Samples ##################
# A Snakemake regular expression matching the forward mate FASTQ files.
# the part in curly brackets {} will be saved, so the variable SAMPLES
# is a list of strings #['Sample1','Sample2'].

#notice that SAMPLES, has a trailing comma.
#you must include this trailing comma, or else the code won’t work correctly.

SAMPLES, = glob_wildcards(join(FASTQ_DIR, '{sample, [^/]+}_R1.fastq.gz'))

# Patterns for the 1st mate and the 2nd mate using the 'sample' wildcard.
PATTERN_R1 = '{sample}_R1.fastq.gz'
PATTERN_R2 = '{sample}_R2.fastq.gz'

# group assignment for each sample
REPLICATES = [0,1,2]
GROUP  = "ES2 ES2 ES8 ES8".split()
GROUP1 = "2D0 2D0".split()
GROUP2 = "8D0 8D0".split()
#PERRETY_SAMPLE = expand("mapped/{sample}.{replicate}.bam", sample=SAMPLES, replicate=[0, 1])
#SAMPLES,GROUP,IDS,= glob_wildcards(join(FASTQ_DIR, '{sample}_{group}_{id}_R1.fastq.gz'))

# gene ontology library
GO_DOMAINS = ['biological_process','cellular_component','molecular_function']

# dirs 
DIRS = ['qc','mapped','counts','alternative_splicing', 'gene_expression',
        'differential_expression','logs','temp']

########### Target output files #################
MULTIQC = 'qc/multiqc_report.html'
DESEQ2_CNT = "counts/All.raw.counts.for.Deseq2.txt"
STRTIEQ = ['gene_expression/'+f+'.tab' for f in SAMPLES]
STRTIE_COUNTS = "counts/gene_count_matrix.csv"
STRTIE_COMPILE = expand("gene_expression/gene_expression_table_annotated.{suf}.csv", suf=['full','tpm','fpkm'])
BALLGOWN = ["gene_expression/ballgown_transcripts_expression_table.csv",
            "gene_expression/ballgown_gene_expression_table.csv"]
################## Rules #######################################

rule target:
    input: STRTIEQ, DESEQ2_CNT, STRTIE_COMPILE, STRTIE_COUNTS, BALLGOWN,  MULTIQC


rule fastqc:
    input:
        "fastq/{prefix}.fastq.gz",
    output:
        "qc/fastqc/{prefix}_fastqc.html",
        "qc/fastqc/{prefix}_fastqc.zip",

    shell: "fastqc -o qc/fastqc {input}"

rule hisat2_index:
    input:
        fasta = FASTA_REF,
        gtf = GTF_FILE,
    output: expand(join(HISAT2_REFDIR,INDEX_PREFIX)+".{ids}.ht2",ids=range(1,9))
    params:
        basename=join(HISAT2_REFDIR, INDEX_PREFIX)
    log: "logs/hisat2/hisat2.index.build.log"
    threads: 12
    shell: "hisat2-build -f {input.fasta}  -p {threads}  {params.basename} &> {log}"


rule hisat2_extract_splicesites:
    input: GTF_FILE
    output:
        splice = join(HISAT2_REFDIR, 'splicesites.txt'),
        exon =   join(HISAT2_REFDIR, 'exon.txt')
    threads: 12
    shell:
        """
        hisat2_extract_splice_sites.py {input} > {output.splice}
        hisat2_extract_exons.py {input} > {output.exon}
        """

rule hisat2_align:
    input:
        index=expand(join(HISAT2_REFDIR,INDEX_PREFIX)+".{ids}.ht2", ids=range(1,9)),
        site = join(HISAT2_REFDIR, "splicesites.txt"),
        r1 = join(FASTQ_DIR, PATTERN_R1),
        r2 = join(FASTQ_DIR, PATTERN_R2)
    output:
        temp('mapped/{sample}.bam')
    log:
        "logs/hisat2/{sample}.align.log"
    threads: 12
    params:
        ref = join(HISAT2_REFDIR, INDEX_PREFIX),
        extra="--min-intronlen 1000 --dta -t"
    shell:
        "(hisat2 {params.extra} --threads {threads} -x {params.ref}"
        " -1 {input.r1} -2 {input.r2}  --known-splicesite-infile {input.site}"
        " | samtools view -Sb -@ {threads}  -o {output} - ) 2> {log}"

rule bam_sort:
    input: "mapped/{sample}.bam"
    output: protected("mapped/{sample}.sorted.bam")
    threads: 12
    shell: "samtools sort -T mapped/{wildcards.sample} -@ {threads} -O bam  {input} > {output}"

rule bam_index:
    input: "mapped/{sample}.sorted.bam"
    output: "mapped/{sample}.sorted.bam.bai"
    shell: "samtools index {input}"

rule bam_stats:
    input:
        bam="mapped/{sample}.sorted.bam",
        bai="mapped/{sample}.sorted.bam.bai"
    output:
        "logs/bamstats/{sample}.bam.stats.txt"
    shell:
        """
		samtools idxstats {input.bam} > {output}
		rm mapped/{wildcards.sample}.bam
		"""

rule geneBody_coverage:
    input:
        bam="mapped/{sample}.sorted.bam",
	bai="mapped/{sample}.sorted.bam.bai",
        anno=RSEQC_ANNO
    output:
        "qc/rseqc/{sample}.geneBodyCoverage.r",
        "qc/rseqc/{sample}.geneBodyCoverage.txt"
    log:"logs/rseqc/{sample}.geneBodyCoverage.log"
    run:
        shell("geneBody_coverage.py -r {input.anno} -i {input.bam}  -o qc/rseqc/{wildcards.sample} &> {log}")

#quantification
rule stringtie:
    input:
        gtf= GTF_FILE,
        bam="mapped/{sample}.sorted.bam",
	bai="mapped/{sample}.sorted.bam.bai"
    output:
        anno="gene_expression/{sample}/{sample}.gtf",
        tab="gene_expression/{sample}.tab"
    threads: 12
    params:
        extra="-e -B"
    shell:
        "stringtie {params.extra} -G {input.gtf} -p {threads} -A {output.tab} -o {output.anno}  {input.bam}"


rule ballgown:
    input: expand("gene_expression/{sample}/{sample}.gtf", sample=SAMPLES)
    output:
        transx="gene_expression/ballgown_transcripts_expression_table.csv",
        genex="gene_expression/ballgown_gene_expression_table.csv",
    params:
        ids =",".join(expand("gene_expression/{sample}",sample=SAMPLES)),
    run:
        R("""
          library("ballgown")
          samples = unlist(strsplit("{params.ids}",","))
          bg = ballgown(samples = samples, meas='all')
          whole_tx_table = texpr(bg, 'all')
          gene_expression = gexpr(bg)
          write.csv(whole_tx_table,file="{output.transx}")
          write.csv(gene_expression,file="{output.genex}")
          """)

##### Read Count #####
rule stringtie_counts:
    input:
        source="scripts/preDEseq.py",
        gtf=expand("gene_expression/{sample}/{sample}.gtf", sample=SAMPLES)
    output: "counts/gene_count_matrix.csv"
    params:
        extra="-l %s"%READ_LEN
    shell: "python {input.source} -i gene_expression {params.extra} -g {output}"

rule htseq:
    input: 
        bam="mapped/{sample}.sorted.bam", 
	bai="mapped/{sample}.sorted.bam.bai",
	gtf=GTF_FILE
    output: "counts/{sample}.htseq.tsv"
    log: "logs/htseq/{sample}.htseq-count.log"
    threads: 1
    shell: "htseq-count -r pos -s no -f bam {input.bam} {input.gtf}  > {output} 2> {log}"

rule complie_htseq:
    input: cnt=expand("counts/{sample}.htseq.tsv", sample=SAMPLES)
    output: deseq="counts/All.raw.counts.for.Deseq2.txt"
    run:
        from pandas import read_table, concat
        count_df=[]
        for count in input.cnt:
            cnt = read_table(count, index_col=0, header=None)
            cnt.columns = [count.split("/")[-1].rstrip(".htseq.tsv")]
            cnt = cnt.sort_index()
            count_df.append(cnt)
        merge_cnt = concat(count_df, axis=1)
        merge_cnt.to_csv(output.deseq)

####complied stringtie output####
rule extract_gtf:
    input: GTF_FILE
    output: GTF_EXTRACT
    run:
        lines_seen = set()
        with open(output[0], 'w') as out:
            out.write("gene_id\tgene_name\tgene_type\tgene_status\n")
            for line in open(input[0]):
                if line.startswith('#'): continue
                attr = dict(item.strip().split(' ')  for item in line.split('\t')[8].strip('\n').split(';') if item)
                line_1st = attr['gene_id'].strip('\"')+'\t'+ attr['gene_name'].strip('\"')+'\t'
                line_2nd = attr['gene_type'].strip('\"')+'\t'+attr['gene_status'].strip('\"')+'\n'
                line_out = line_1st + line_2nd
                if line_out not in lines_seen:
                    out.write(line_out)
                    lines_seen.add(line_out)

rule compile_stringtie:
    """
    Compile stringtie output for all samples.
    """
    input:
        annotation = GTF_EXTRACT,
        filelist = expand("gene_expression/{sample}.tab", sample=SAMPLES)
    output:
        full = "gene_expression/gene_expression_table_annotated.full.csv",
        tpm  = "gene_expression/gene_expression_table_annotated.tpm.csv",
        fpkm = "gene_expression/gene_expression_table_annotated.fpkm.csv",
    run:
        from pandas import read_table, concat
        anno = read_table(input.annotation)
        frames = []
        for f in input.filelist:
            name = f.split("/")[-1].strip(".tab")
            frame = read_table(f, index_col="Gene ID")
            frame = frame.sort_index()
            frame.rename(columns={'Coverage': 'Coverage.'+name,
                                  'FPKM': 'FPKM.'+name, 'TPM':'TPM.'+name},
                         inplace=True)
            frames.append(frame)

        #when concat dataframe, differrent dataframes will ordered by their index by default
        result = concat(frames, axis=1)
        df = result.loc[:,~result.columns.duplicated()]
        df_merge = anno.merge(df, left_on='gene_id',right_index=True, how='right')
        df_merge.drop('Gene Name', axis=1, inplace=True)
        df_merge.to_csv(output.full, index=False)
        col_tpm = ['gene_id','gene_name'] + [item for item in df_merge.columns if item.startswith("TPM.")]
        col_fpkm = ['gene_id','gene_name'] + [item for item in df_merge.columns if item.startswith("FPKM.")]
        df_merge[col_tpm].to_csv(output.tpm, index=False)
        df_merge[col_fpkm].to_csv(output.fpkm, index=False)

# Alternative splicing
rule rmats:
    input:
        bam=expand("mapped/{sample}.sorted.bam", sample=SAMPLES),
        gtf=GTF_FILE
    output:
        expand("alternative_splicing/{sample1}_vs_{sample2}/MATS_output/{type}",
                sample1='GROUP1', sample2='GROUP2', type=['A3SS','A5SS','MXE','RI','SE'])
    params:
        b1=",".join(expand("mapped/{sample}.sorted.bam", sample=GROUP1)),
        b2=",".join(expand("mapped/{sample}.sorted.bam", sample=GROUP2)),
        rmats=RMATS_DIR,
        prefix="alternative_splicing/{sample1}_vs_{sample2}".format(sample1='GROUP1', sample2='GROUP2'),
        extra="-t paried -len %s -a 1 -c 0.0001 -analysis U -novelSS 0"%READ_LEN
    shell:
        "python {params.rmats} -b1 {params.b1} -b2 {params.b2} "
        "-gtf {input.gtf}  -o {params.prefix} {params.extra}"


##### Aggreate QC ##########
rule multiqc:
    input:
        expand("qc/fastqc/{sample}_{read}_fastqc.{suf}", sample=SAMPLES, read=['R1','R2'], suf=['zip','html']),
        expand("qc/rseqc/{sample}.geneBodyCoverage.{suf}", sample=SAMPLES, suf=['r','txt']),
	expand("logs/bamstats/{sample}.bam.stats.txt",sample=SAMPLES),
	expand("counts/{sample}.htseq.tsv",sample=SAMPLES)
    output: html='qc/multiqc_report.html'
    params:
        analysis_dir='.',
	outdir="qc",
        extra="--config multiqc_config.yaml"
    shell: "multiqc --quiet --force -p -o {params.outdir} {params.analysis_dir}"
