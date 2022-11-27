import os
from os.path import join, isfile
from itertools import combinations
from snakemake.shell import shell
import pandas as pd


include: "rules/common.smk"

#configfile: 'config.yml'
workdir: config['workdir']

################### globals #############################################

# Full path to an uncompressed FASTA file with all chromosome sequences.
CDNA = config['cdna']

# Full path to a folder that holds all of your FASTQ files.
FASTQ_DIR = config['fastq_dir']
READ_LEN = config['read_length']
PAIRED = 'paired' if config['paired'] else 'single'
# Full path to a Genome.
GENOME = config['genome']
#CDNA =           join(GENOME,"gencode.v25.transcripts.fa")
# genome sequence
FASTA_REF =      config['dna']
# index_dir
STAR_REFDIR =    config['star_index']
# index basename
INDEX_PREFIX = 'hg38'
# gtf
GTF_FILE =       config['gtf']
# transcriptiome cdna
CDNA = config['cdna'] 

#rseqc_annotation
RSEQC_ANNO  = config['rseqc']
############ Samples ##################
# A Snakemake regular expression matching the forward mate FASTQ files.
# the part in curly brackets {} will be saved, so the variable SAMPLES
# is a list of strings #['Sample1','Sample2'].

#notice that SAMPLES, has a trailing comma.
#you must include this trailing comma, or else the code won’t work correctly.
#SAMPLES, = glob_wildcards(os.path.join(FASTQ_DIR, '{sample, SRR[^/]+}_R1.fastq.gz'))

SAMPLES, SAMPLES_ALIAS,GROUP,TIME = parse_samples(config['sample_meta'])
# #rMATS
# uGroup=unique(GROUP)
uGroup = SAMPLES
RMATS_DICT = [[] for i in range(len(uGroup))]
# for i,g in enumerate(GROUP):
#     for j, u in enumerate(uGroup):
#         if g == u:
#             RMATS_DICT[j].append(SAMPLES[i])

PATTERN_U = config['read_pattern']['u']
PATTERN_R1 = config['read_pattern']['r1']
PATTERN_R2 = config['read_pattern']['r2']
# FASTQ_GE = {s:[os.path.join(FASTQ_DIR, f) 
#                for f in os.listdir(FASTQ_DIR) 
#                if (f.startswith(s) and (f.endswith(("fq", "fq.gz", "fastq", "fastq.gz")))) 
#                ] 
#                for s in SAMPLES}
FASTQ_GE = get_sample_fastqs(SAMPLES, FASTQ_DIR)

# dirs
DIRS = ['qc','BAMs','AS_rMATS', 'quant', 'degs', 'logs','temp']


########### Target output files #################

RMATS_TURBO =expand("AS_rMATS/{t}.MATS.JCEC.txt", t=['A3SS','A5SS','MXE','RI','SE'])
BIGWIG = expand("igv/{sample}.sorted.bw", sample=SAMPLES)
COUNTS = "quant/count_matrix.txt"
QUANT = ["quant/gene_expression.TPM.annotated.csv",
        "quant/transcript_expression.TPM.annotated.csv"]

MULTIQC = 'multiqc/multiqc_report.html'
################## Rules #######################################

rule target:
    input: COUNTS, QUANT, RMATS_TURBO, MULTIQC#BIGWIG, # RMATS_TURBO

# rule fastqc:
#     input:  
#         fastqs = lambda wildcards: FASTQ_GE[wildcards.sample]['fastq']
#     output: 
#         "qc/fastqc/{sample}_1_fastqc.zip", 
#         "qc/fastqc/{sample}_2_fastqc.zip"
#     log:    "logs/fastqc/{sample}_fastqc"
#     threads: 2
#     params : jobname = "{sample}"
#     message: "fastqc {input}: {threads}"
#     shell:
#         # fastqc works fine on .gz file as well:  
#         """
#         fastqc --threads {threads} -o qc/fastqc -f fastq --noextract {input.fastqs} 2> {log}
#         """

rule star_index:
    input:
        fasta = FASTA_REF,
        gtf = GTF_FILE,
    output:
        os.path.join(STAR_REFDIR, "SAindex"),
        os.path.join(STAR_REFDIR, "geneInfo.tab"),
        os.path.join(STAR_REFDIR, "sjdbList.out.tab"),
        
    params:
        outdir = STAR_REFDIR,
        read_length = int(READ_LEN) -1 
    threads: 16
    log: "logs/star/star.index.log"
    message: "every time you have different read length, you need to re-build index for the read length you've got"
    shell:
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir {params.outdir} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang {params.read_length} "
        "2> {log}"

rule star_align:
    input:
        index=[os.path.join(STAR_REFDIR, "SAindex"),
               os.path.join(STAR_REFDIR, "geneInfo.tab"),
               os.path.join(STAR_REFDIR, "sjdbList.out.tab")],
        fastqs = lambda wildcards: FASTQ_GE[wildcards.sample]['fastq'],
    output:
        'BAMs/{sample}.Aligned.sortedByCoord.out.bam',
        'BAMs/{sample}.ReadsPerGene.out.tab',
        "BAMs/{sample}.Aligned.toTranscriptome.out.bam"
    log:
        "logs/star/{sample}.align.log"
    threads: 12
    params:
        ref = STAR_REFDIR
    shell:
        "STAR --genomeDir {params.ref}  "
        "--readFilesIn   {input.fastqs} "
        "--readFilesCommand zcat --runThreadN {threads} " 
        "--outFileNamePrefix BAMs/{wildcards.sample}.  " 
        "--outSAMtype BAM SortedByCoordinate "
        "--outSAMstrandField intronMotif " 
        "--outSAMunmapped Within "
        "--alignEndsType EndToEnd "
        "--quantTranscriptomeBan Singleend " # salmon, if RESM, use IndelSoftclipSingleend
        "--quantMode TranscriptomeSAM GeneCounts "
        #"--genomeLoad LoadAndRemove "
        "2> {log}"

# rule bam_sort:
#     input: "BAMs/{sample}.bam"
#     output: protected("BAMs/{sample}.sorted.bam")
#     threads: 8
#     shell: "samtools sort -@ {threads} {input} > {output}"

rule bam_index:
    input: 'BAMs/{sample}.Aligned.sortedByCoord.out.bam'
    output: 'BAMs/{sample}.Aligned.sortedByCoord.out.bam.bai'
    shell: "samtools index {input}"

rule bam2bw:
    input:
        bam='BAMs/{sample}.Aligned.sortedByCoord.out.bam',
        bai='BAMs/{sample}.Aligned.sortedByCoord.out.bam.bai'
    output:
        "igv/{sample}.sorted.bw"
    threads: 8
    shell:
        "bamCoverage --normalizeUsing RPKM -p {threads} -b {input.bam} -o {output}"


rule generate_transcriptome:
    input: 
        fasta=FASTA_REF,
        gtf = GTF_FILE,
    output: CDNA
    shell:
        # https://github.com/gpertea/gffread 
        # or download from gencode if you genome and gtf also from gencode
        "gffread -w {output} -g {input.fasta} {input.gtf} "


rule salmon_quant:
    input: 
        bam="BAMs/{sample}.Aligned.toTranscriptome.out.bam",
        cdna=CDNA,
        gtf=GTF_FILE
    output: 
        "salmon_star/{sample}/quant.sf",
        "salmon_star/{sample}/quant.genes.sf"
    threads: 8
    log: "logs/salmon/salmon.quant.{sample}.log"
    shell:
        "salmon quant --threads {threads} "
        "--targets {input.cdna} "
        "--gencode "
        "--geneMap {input.gtf} "
        "--libType A " # -l (depends on the lib type, ISR for truseq stranded, equivalent to tophat -fr-firststrand)
        "--output salmon_star/{wildcards.sample} "
        "--alignments {input.bam} "
        "--seqBias --gcBias "
        "2> {log}"


rule count_matrix:
    input: expand('BAMs/{sample}.ReadsPerGene.out.tab', sample=SAMPLES)
    output: "quant/count_matrix.txt"
    params: 
        samples = SAMPLES
    run:
        header = "gene_id\t" + "\t".join(params.samples) + "\n"
        # retrieve the 2th column of each "ReadsPerGene.out.tab" file + the first column that contains the gene IDs
        shell("""paste {input} | grep -v "_" | awk '{{printf "%s", $1}} {{for (i=2;i<=NF;i+=4) printf "\\t%s", $i; printf "\\n" }}' > {output}
              """)
        # insert header in the first line
        shell("sed  -i '1i %s' {output} "%header)


rule tpm_matrix:
    input: expand("salmon_star/{sample}/quant.genes.sf", sample=SAMPLES)
    output: temp("quant/gene_expression.TPM.txt")
    params: 
        samples = SAMPLES
    run:
        header = "gene_id\t" + "\t".join(params.samples) + "\n"
        # retrieve the 4th column of each "quant.genes.sf" file + the first column that contains the gene IDs
        shell("""paste {input} | grep -v "^Name" | awk '{{printf "%s", $1}} {{for (i=4;i<=NF;i+=5) printf "\\t%s", $i; printf "\\n" }}' > {output}
              """)
        # insert header in the first line
        shell("sed  -i '1i %s' {output} "%header)

rule tpm_matrix_tx:
    input: expand("salmon_star/{sample}/quant.sf", sample=SAMPLES)
    output: temp("quant/transcript_expression.TPM.txt")
    params: 
        samples = SAMPLES
    run:
        header = "tx_id\t" + "\t".join(params.samples) + "\n"
        # retrieve the 4th column of each "quant.genes.sf" file + the first column that contains the gene IDs
        shell("""paste {input} | grep -v "^Name" | awk '{{printf "%s", $1}} {{for (i=4;i<=NF;i+=5) printf "\\t%s", $i; printf "\\n" }}' > {output}
              """)
        # insert header in the first line
        shell("sed  -i '1i %s' {output} "%header)


rule gtf_extract:
    input: GTF_FILE
    output:
        gene_anno="gene_anno.txt",
        tx2gene = "tx_anno.txt"
    script:
        "scripts/extractGTF.py"

rule anno_samples:
    input:
        "gene_anno.txt",
        "tx_anno.txt",
        "quant/gene_expression.TPM.txt",
        "quant/transcript_expression.TPM.txt",
    output:
        "quant/gene_expression.TPM.annotated.csv",
        "quant/transcript_expression.TPM.annotated.csv",
    params:
        samples=SAMPLES,
    script:
        "scripts/annotateTPMs.py"


rule rMATS_pre:
    """prepared bam and gtf files for rmats docker image"""
    input:
        bam=expand("BAMs/{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        #gtf=GTF_FILE
    output:
        #groups= ["temp/rmats/%s_vs_%s.rmats.txt"%(j, i) for i, j in combinations(uGroup, 2)],
        groups = "AS_rMATS/rmats.samples.all.txt"
    run:
        with open(output[0],'w') as out:
            ss = ",".join(input.bam)
            out.write(ss+"\n")
        # for u, g in zip(params.ugroup, params.ugsamples):
        #     out = open("temp/rmats/b_%s.txt"%u, 'w')
        #     temp = ["BAMs/%s.sorted.bam"%sample for sample in g]
        #     line=",".join(temp)
        #     out.write(line)
        #     out.close()
        # for i, j in combinations(params.ugroup, 2):
        #     outname = "temp/rmats/%s_vs_%s.rmats.txt"%(j,i)
        #     out2 = open(outname,'w')
        #     out2.write("temp/rmats/b_%s.txt\n"%j)
        #     out2.write("temp/rmats/b_%s.txt\n"%i)
        #     out.close()
        #shell("cp {input.gtf} {output.gtf_tmp}")

rule rMATS_turbo:
    input:
        bam=expand("BAMs/{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        bai=expand("BAMs/{sample}.Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES),
        gtf = GTF_FILE, 
        groups = "AS_rMATS/rmats.samples.all.txt"
    output:
        "AS_rMATS/SE.MATS.JCEC.txt",
        "AS_rMATS/A3SS.MATS.JCEC.txt",
        "AS_rMATS/A5SS.MATS.JCEC.txt",
        "AS_rMATS/RI.MATS.JCEC.txt",
        "AS_rMATS/MXE.MATS.JCEC.txt"
    threads: 16
    log: "logs/rMATS-turbo/rMATS.turbo.log",
    params:
        prefix="AS_rMATS",
        extra=" -t %s --readLength %s --anchorLength 1 "%(PAIRED, READ_LEN),
        wkdir= config['workdir'],
        #gtf = join("temp", GTF_FILE.split("/")[-1])
    message: "quantify alternative splicing events in all samples using rMATS. Just like gene expression TPM values."
    shell:
        # conda install rmats
        "rmats.py --b1 {input.groups} "
        "--gtf {input.gtf} --od {params.prefix} "
        "--nthread {threads} --statoff {params.extra} "
        "--tmp {params.prefix}/splicing_graph &> {log}"


# rule rMATS_anno:
#     input:
#         "differential_expression/diff_{treat}_vs_{ctrl}/diff_{treat}_vs_{ctrl}_results.annotated.xls",
#         "alternative_splicing/rMATS.{treat}_vs_{ctrl}/SE.MATS.JCEC.txt",
#         "alternative_splicing/rMATS.{treat}_vs_{ctrl}/A3SS.MATS.JCEC.txt",
#         "alternative_splicing/rMATS.{treat}_vs_{ctrl}/A5SS.MATS.JCEC.txt",
#         "alternative_splicing/rMATS.{treat}_vs_{ctrl}/RI.MATS.JCEC.txt",
#         "alternative_splicing/rMATS.{treat}_vs_{ctrl}/MXE.MATS.JCEC.txt"
#     output:
#         "alternative_splicing/rMATS.{treat}_vs_{ctrl}_sig/SE.MATS.JCEC.sig.txt",
#         "alternative_splicing/rMATS.{treat}_vs_{ctrl}_sig/A3SS.MATS.JCEC.sig.txt",
#         "alternative_splicing/rMATS.{treat}_vs_{ctrl}_sig/A5SS.MATS.JCEC.sig.txt",
#         "alternative_splicing/rMATS.{treat}_vs_{ctrl}_sig/RI.MATS.JCEC.sig.txt",
#         "alternative_splicing/rMATS.{treat}_vs_{ctrl}_sig/MXE.MATS.JCEC.sig.txt",
#         "alternative_splicing/rMATS.{treat}_vs_{ctrl}_sig/Skip_Exons/SE.MATS.JCEC.sig.annotated.csv",
#     params:
#         indir="alternative_splicing/rMATS.{treat}_vs_{ctrl}",
#         outdir="alternative_splicing/rMATS.{treat}_vs_{ctrl}_sig",
#         go=config['enrichr_library'],
#         rbps=config['rbps']
#     script:
#         "scripts/annotateRMATS.py"

rule bam_stats:
    input:
        bam='BAMs/{sample}.Aligned.sortedByCoord.out.bam',
        bai='BAMs/{sample}.Aligned.sortedByCoord.out.bam.bai'
    output: "qc/rseqc/{sample}.bamstats.txt"
    shell: "bam_stat.py -i {input.bam} > {output}"

rule geneBody_coverage:
    input:
        bam='BAMs/{sample}.Aligned.sortedByCoord.out.bam',
        bai='BAMs/{sample}.Aligned.sortedByCoord.out.bam.bai',
        anno=RSEQC_ANNO['housekeep']
    output:
        "qc/rseqc/{sample}.geneBodyCoverage.r",
        "qc/rseqc/{sample}.geneBodyCoverage.txt"
    log:"logs/rseqc/{sample}.geneBodyCoverage.log"
    shell:
        "geneBody_coverage.py -r {input.anno} -i {input.bam}  -o qc/rseqc/{wildcards.sample} &> {log}"

rule read_distribution:
    input:
        bam='BAMs/{sample}.Aligned.sortedByCoord.out.bam',
        bai='BAMs/{sample}.Aligned.sortedByCoord.out.bam.bai',
        bed=RSEQC_ANNO['refseq']
    output:
        "qc/rseqc/{sample}.readDistribution.txt"
    shell:
        "read_distribution.py -i {input.bam} -r {input.bed} > {output}"

rule collopase_annotation:
    input: GTF_FILE
    output: "temp/gencode.collopased.annotation.gtf"
    shell:
        ## note: works for gencode gtf
        "python3 collapse_annotation.py {input} {output}"

rule rnaseqc:
    input:
        bam='BAMs/{sample}.Aligned.sortedByCoord.out.bam',
        bai='BAMs/{sample}.Aligned.sortedByCoord.out.bam.bai',
        gtf="temp/gencode.collopased.annotation.gtf"
    output:
        "qc/rnaseqc/{sample}.coverage.tsv",
        "qc/rnaseqc/{sample}.metrics.tsv"
    params:
        outdir="qc/rnaseqc/{sample}"
    shell:
        # conda install rna-seqc
        "rnaseqc {input.gtf} {input.bam} --coverage --sample {wildcards.sample} {params.outdir}"

rule multiqc:
    """Aggreate QC """
    input:
        expand("qc/rseqc/{sample}.geneBodyCoverage.txt", sample=SAMPLES),
        expand("qc/rseqc/{sample}.readDistribution.txt", sample=SAMPLES),
        expand("salmon_star/{sample}/quant.genes.sf", sample=SAMPLES),
        #expand("qc/fastqc/{sample}_{r}_fastqc.zip", sample=SAMPLES, r=[1,2]),
    output: html='multiqc/multiqc_report.html'
    params:
        analysis_dir=["qc", "BAMs", "salmon_star"],
        extra="--config multiqc_config.yaml"
    shell: "multiqc --interactive --quiet --force -p -o multiqc {params.analysis_dir}"