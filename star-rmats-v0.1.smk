from os.path import join, isfile
from itertools import combinations
from snakemake.shell import shell

include: "rules/common.smk"

#configfile: 'config.yml'
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
STAR_REFDIR =config['star_index']
# index basename
INDEX_PREFIX = 'hg38'
# gtf
GTF_FILE =       config['gtf']

############ Samples ##################
# A Snakemake regular expression matching the forward mate FASTQ files.
# the part in curly brackets {} will be saved, so the variable SAMPLES
# is a list of strings #['Sample1','Sample2'].

#notice that SAMPLES, has a trailing comma.
#you must include this trailing comma, or else the code won’t work correctly.

#SAMPLES, = glob_wildcards(join(FASTQ_DIR, '{sample, SRR[^/]+}_R1.fastq.gz'))

SAMPLES,SAMPLES_ALIAS,GROUP,TIME = parse_samples(config['sample_meta'])
# #rMATS
# uGroup=unique(GROUP)
uGroup = SAMPLES
RMATS_DICT = [[] for i in range(len(uGroup))]
# for i,g in enumerate(GROUP):
#     for j, u in enumerate(uGroup):
#         if g == u:
#             RMATS_DICT[j].append(SAMPLES[i])

# Patterns for the 1st mate and the 2nd mate using the 'sample' wildcard.
#PATTERN_R1 = '{sample}_R1.fastq.gz'
#PATTERN_R2 = '{sample}_R2.fastq.gz'
PATTERN_U = config['read_pattern']['u']
PATTERN_R1 = config['read_pattern']['r1']
PATTERN_R2 = config['read_pattern']['r2']


FASTQ_GE = {s:[os.path.join(FASTQ_DIR, f) 
               for f in os.listdir(FASTQ_DIR) 
               if (f.startswith(s) and (f.endswith(".fastq") or f.endswith(".gz") or f.endswith(".fq"))) 
               ] 
               for s in SAMPLES}
#SAMPLES, GROUP, IDS,= glob_wildcards(join(FASTQ_DIR, '{sample}_{group}_{id}_R1.fastq.gz'))
# dirs
DIRS = ['qc','BAMs','alternative_splicing', 'gene_expression',
        'differential_expression','logs','temp']


########### Target output files #################
RMATS_TEMP=["alternative_splicing/rMATS.%s_vs_%s_sig/{type}.MATS.JCEC.sig.txt"%(j, i) for i, j in combinations(uGroup, 2)]
RMATS_TURBO =[temp.format(type=t) for temp in RMATS_TEMP for t in ['A3SS','A5SS','MXE','RI','SE']]
BIGWIG = expand("igv/{sample}.sorted.bw", sample=SAMPLES)
################## Rules #######################################

rule target:
    input: "counts/count_matrix.txt", #BIGWIG, # RMATS_TURBO

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
        read_length = READ_LEN
    threads: 12
    shell:
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir {params.outdir} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang {params.read_length} "

rule star_align:
    input:
        index=[os.path.join(STAR_REFDIR, "SAindex"),
               os.path.join(STAR_REFDIR, "geneInfo.tab"),
               os.path.join(STAR_REFDIR, "sjdbList.out.tab")],
        fastqs = lambda wildcards: FASTQ_GE[wildcards.sample],
    output:
        'BAMs/{sample}.Aligned.sortedByCoord.out.bam',
        'BAMs/{sample}.ReadsPerGene.out.tab'
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
        "--alignEndsType EndToEnd "
        "--quantMode GeneCounts "
        #"--genomeLoad LoadAndRemove "
        "2> {log}"

# STAR --genomeDir {params.ref}  \
#      --readFilesIn   {input.fastqs} \
#      --readFilesCommand zcat \ 
#      --runThreadN {threads} \
#      --outFileNamePrefix BAMs/{wildcards.sample}  \
#      --outSAMtype BAM SortedByCoordinate \
#      --outSAMstrandField intronMotif \
#      --alignEndsType EndToEnd \
#      --outSAMunmapped Within \
#      --quantMode GeneCounts \
#      --genomeLoad LoadAndRemove

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

rule count_matrix:
    input: expand('BAMs/{sample}.ReadsPerGene.out.tab', sample=SAMPLES)
    output: "counts/count_matrix.txt"
    params: 
        samples = SAMPLES
    run:
        header = "gene_id\t" + "\t".join(params.samples) + "\n"
        # retrieve the 2th column of each "ReadsPerGene.out.tab" file + the first column that contains the gene IDs
        shell("""paste {input} | grep -v "_" | awk '{{printf "%s\\t", $1}} {{for (i=2;i<=NF;i+=4) printf "%s\\t", $i; printf "\\n" }}' > {output}
              """)
        # insert header in the first line
        shell("sed  -i '1i %s' {output} "%header)



# rule rMATS_pre:
#     """prepared bam and gtf files for rmats docker image"""
#     input:
#         bam=expand("BAMs/{sample}.sorted.bam", sample=SAMPLES),
#         bai=expand("BAMs/{sample}.sorted.bam.bai", sample=SAMPLES),
#         gtf=GTF_FILE
#     output:
#         groups= ["temp/rmats/%s_vs_%s.rmats.txt"%(j, i) for i, j in combinations(uGroup, 2)],
#         gtf_tmp = join("temp", GTF_FILE.split("/")[-1])
#     params:
#         ugsamples=RMATS_DICT,
#         ugroup=uGroup,
#     run:
#         for u, g in zip(params.ugroup, params.ugsamples):
#             out = open("temp/rmats/b_%s.txt"%u, 'w')
#             temp = ["/data/BAMs/%s.sorted.bam"%sample for sample in g]
#             line=",".join(temp)
#             out.write(line)
#             out.close()
#         for i, j in combinations(params.ugroup, 2):
#             outname = "temp/rmats/%s_vs_%s.rmats.txt"%(j,i)
#             out2 = open(outname,'w')
#             out2.write("temp/rmats/b_%s.txt\n"%j)
#             out2.write("temp/rmats/b_%s.txt\n"%i)
#             out.close()
#         shell("cp {input.gtf} {output.gtf_tmp}")

# rule rMATS_turbo:
#     input:
#         bam=expand("BAMs/{sample}.sorted.bam", sample=SAMPLES),
#         bai=expand("BAMs/{sample}.sorted.bam.bai", sample=SAMPLES),
#         gtf = join("temp", GTF_FILE.split("/")[-1]),
#         rmats="temp/rmats/{treat}_vs_{ctrl}.rmats.txt"
#     output:
#         "alternative_splicing/rMATS.{treat}_vs_{ctrl}/SE.MATS.JCEC.txt",
#         "alternative_splicing/rMATS.{treat}_vs_{ctrl}/A3SS.MATS.JCEC.txt",
#         "alternative_splicing/rMATS.{treat}_vs_{ctrl}/A5SS.MATS.JCEC.txt",
#         "alternative_splicing/rMATS.{treat}_vs_{ctrl}/RI.MATS.JCEC.txt",
#         "alternative_splicing/rMATS.{treat}_vs_{ctrl}/MXE.MATS.JCEC.txt"
#     threads: 8
#     log: "logs/rMATS-turbo/{treat}_vs_{ctrl}.rMATS.turbo.log",
#     params:
#         prefix="alternative_splicing/rMATS.{treat}_vs_{ctrl}",
#         extra=" -t %s --readLength %s --anchorLength 1 "%(PAIRED, READ_LEN),
#         wkdir= config['workdir'],
#         gtf = join("temp", GTF_FILE.split("/")[-1])
#     run:
#         # blacklist to skip
#         if isfile("temp/blacklist.txt"):
#             with open("temp/blacklist.txt") as black:
#                 blacklist = [ bla.strip().split("/")[-1] for bla in black]
#             # groups you want to skip
#             bk = "diff_%s_vs_%s_results.annotated.xls"%(wildcards.treat, wildcards.ctrl)
#             if bk in blacklist:
#                 for ast in output:
#                     shell("touch %s"%ast)
#                 return

#         shell("""docker run -v {params.wkdir}:/data rmats:turbo01 \
#                  --b1 /data/temp/rmats/b_{wildcards.treat}.txt --b2 /data/temp/rmats/b_{wildcards.ctrl}.txt \
#                  --gtf /data/{params.gtf} --od /data/{params.prefix} \
#                  --nthread {threads} --tstat {threads} {params.extra} &> {log}""")


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
