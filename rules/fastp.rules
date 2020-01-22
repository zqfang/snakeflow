rule trim_galore:
    input: 
        r1=join(FASTQ_RAW,"{sample}_R1.fastq.gz"),
        r2=join(FASTQ_RAW,"{sample}_R2.fastq.gz"),
    output:
        r1=join(FASTQ_CLEAN, "{sample}_R1_val_1.fastq.gz"),
        r2=join(FASTQ_CLEAN, "{sample}_R2_val_2.fastq.gz"),

    log: "logs/trim-galore/${sample}.trimgalore.log"
    params:
        outdir="fastq_clean",
        fastqc_outdir="qc/fastqc_clean"
    shell:
        "trim_galore -q 20 --phred33  --stringency 1 -e 0.1  --length 50 --paired " 
        "-o {params.outdir}  {input.r1} {input.r2} > {log}"

rule trimmonatic:
    input:
        r1=join(FASTQ_RAW, "{sample}_R1.fastq.gz"),
        r2=join(FASTQ_RAW, "{sample}_R2.fastq.gz"),
        adap=config['adaptors']['illumina']
    output:
        r1=join(FASTQ_CLEAN, "{sample}_R1_val_1.fastq.gz"),
        r2=join(FASTQ_CLEAN, "{sample}_R2_val_2.fastq.gz"),
        r1up=join(FASTQ_CLEAN, "{sample}_R1_val_upparied_1.fastq.gz"),
        r2up=join(FASTQ_CLEAN, "{sample}_R1_val_upparied_2.fastq.gz)"
    log: "logs/trimmonatic/${sample}.trimmonatic.log"
    thread: 8
    shell:
        "trimmomatic PE -threads {thread} -phred33 {output.r1} {ouput.upr1} {output.r2} {ouput.upr2} "
        "LLUMINACLIP:{input.adap}:2:30:10 LEADING:10 TRAILING:20 SLIDINGWINDOW:4:25 MINLEN:50 "
        "| tee {log} "

rule fastp:
    input:
        r1=join(FASTQ_RAW, "{sample}_R1.fastq.gz"),
        r2=join(FASTQ_RAW, "{sample}_R2.fastq.gz"),   
    output:
        r1=join(FASTQ_CLEAN, "{sample}_R1.fastq.gz"),
        r2=join(FASTQ_CLEAN, "{sample}_R2.fastq.gz"),
    tread: 8
    params:
    log: "logs/fastp/${sample}.fastp.log"
    shell:
        "fastp -q 20 --thread {thread} -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} --html {wildcards.sample}.html &> {log}" 
