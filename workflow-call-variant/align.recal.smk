rule bwa_index:
    input: GENOME,
    output: ""
    params: 
        prefix="GRCm38"
    log: "logs/bwa.index.log"
    shell:
        "bwa index -a bwtsw -p {params.prefix} {input} 2> {log}"

rule bwa_mem:
    input:
        r1="fastq/{sample}_R1.fastq.gz",
        r2="fastq/{sample}_R2.fastq.gz",,
        index="GRCm38"
    ouput: temp(".sam")
    threads: 8
    log: "logs/bwa.mem.log"
    params:
        # PL has to be one of ILLUMINA，SLX，SOLEXA，SOLID，454，LS454，COMPLETE，
        # PACBIO，IONTORRENT，CAPILLARY，HELICOS, or UNKNOWN
        # ID：输入reads集的ID号； LB： reads集的文库名； SM：样本名称； PL：测序平台
        RG="@RG\tID:SAMPLE_ID\tPL:illumina\tSM:{sample}"
    shell:
        "bwa mem -t {threads} -R '{param.RG}' {input.index} "
        "{input.r1} {input.r2} 1> {output} 2> {log}"
        
rule sam_sort:
    input: "BAM/{sample}.sam"
    output: temp("BAM/{sample}.bam")
    params:
        tmpdir=TMPDIR,
        java_ops= "-Xmx32G -Djava.io.tmpdir=%s"%TMPDIR
    shell:
        "gatk --java-options '{params.java_ops}' SortSam -I {input} -O {output} "
        "-SO coordinate --CREATE_INDEX true"

rule mark_dups:
    input:  "BAM/{sample}.bam"
    output: 
        bam = temp("{sample}.marked.bam")
        metric = temp("{sample}.marked.bam.metrics")
    params:
        extra = "--CREATE_INDEX" # if skip fixMateInformation, you need this
    shell:
        "gatk MarkDuplicates -I {input} -O {output.bam} -M {output.metrics}"

rule fix_mate:
    input: "BAM/sample}.marked.bam"
    output: temp("BAM/{sample}.marked.fixed.bam")
    shell:
        "gatk FixMateInformation -I {input} -O {output.bam} "
        "-SO coordinate --CREATE_INDEX true"

rule genome_dict:
    input: GENOME
    output: GENOME + ".dict"
    shell:
        "gatk CreateSequenceDictionary -R {input} -O {output}"

## dbSNP, dbINDEL need to be indexed
# rule vcf_index:
#     input: dbINDEL, dbSNP
#     output: dbINDEL +".idx", dbSNP +".idx"
#     shell:
#         """gatk IndexFeatureFile -F {input[0]}
#            gatk IndexFeatureFile -F {input[1]}
#         """

rule baseRecalibrate:
    input: 
        genome=GENOME,
        bam = "BAM/{sample}.marked.fixed.bam",
        bai = "BAM/{sample}.marked.fixed.bai",
    output: 
        "BAM/{sample}.marked.fixed.table"
    params:
        dbsnp = dbSNP,
        dbindel = dbINDEL,
    shell:
        "gatk BaseRecalibrator -R {input.genome} -I {input.bam} -O {output} "
        "--known-sites {params.dbsnp} --known-sites {params.dbindel}"

rule applyBQSR:
    input: 
        genome=GENOME
        bam = "BAM/{sample}.marked.fixed.bam",
        bai = "BAM/{sample}.marked.fixed.bai",
        table = "BAM/{sample}.marked.fixed.table",
    output: 
        protected("BAM/{sample}.marked.fixed.BQSR.bam")
    params:
        dbsnp=dbSNP,
        dbindel=dbINDEL
    shell:
        "gatk ApplyBQSR -R {input.genome} -I {input.bam} "
        "--bqsr-recal-file {input.table} -O {output}"
