import os
workdir: "/u04/sv_calling/bams"

with open("/u04/sv_calling/fastq.files.txt", 'r') as fastq:
    fastqs_r1 = fastq.read().strip().split()
with open("/u04/sv_calling/fq.files.txt", 'r') as fastq:
    fq_r1 = fastq.read().strip().split()
#fastqs_r2 = fastqs_r1.replace("_1.fastq.gz", "_2.fastq.gz")

samples = []
ids = []

for fq in fastqs_r1:
    tempx = fq.split("/")        
    if tempx[-3] == 'BTBR': # BTBR
        samples.append(tempx[-3])
    else:
        samples.append(tempx[-2])
    ids.append(tempx[-1].replace("_1.fastq.gz", ""))

for fq in fq_r1:
    tempx = fq.split("/")        
    samples.append(tempx[-2])
    ids.append(tempx[-1].replace("_1.fq.gz", ""))


OUTPUT_BAM = ["/u04/sv_calling/bams/%s/%s.bam"%(s, sid) for s, sid in zip(samples, ids)]

GENOME =  "/home/mzheng/NGS/ref_B38/GRCm38_68.fa"


rule target:
    input: OUTPUT_BAM

rule align:
    input:
        r1="/home/mzheng/rawData_NGS/{sample}/{ID}_1.fastq.gz",
        r2="/home/mzheng/rawData_NGS/{sample}/{ID}_2.fastq.gz",
        genome=GENOME, # bwa indexed
    output: 
        bam="/u04/sv_calling/bams/{sample}/{ID}.bam",
        discordants="/u04/sv_calling/bams/{sample}/{ID}.discordants.bam",
        splitters="/u04/sv_calling/bams/{sample}/{ID}.splitters.bam",
    params:
        outprefix="/u04/sv_calling/bams/{sample}/{ID}",
        # PL has to be one of ILLUMINA，SLX，SOLEXA，SOLID，454，LS454，COMPLETE，
        # PACBIO，IONTORRENT，CAPILLARY，HELICOS, or UNKNOWN
        # ID：输入reads集的ID号； LB： reads集的文库名； SM：样本名称； PL：测序平台
        RG="@RG\tID:{sample}.{ID}\tPL:illumina\tSM:{sample}\tLB:{ID}"
    run:
        shell("/home/fangzq/github/speedseq/bin/speedseq align -o {params.outprefix} -R '{params.RG}' {input.genome} {input.r1} {input.r2}")
        shell("echo -e {params.RG} >> /u04/sv_calling/bams/{sample}/all.LBIDs.txt")

rule align2:
    input:
        r1="/home/mzheng/rawData_NGS/Sanger/{sample}/{ID}_1.fastq.gz",
        r2="/home/mzheng/rawData_NGS/Sanger/{sample}/{ID}_2.fastq.gz",
        genome=GENOME, # bwa indexed
    output: 
        bam="/u04/sv_calling/bams/{sample}/{ID}.bam",
        discordants="/u04/sv_calling/bams/{sample}/{ID}.discordants.bam",
        splitters="/u04/sv_calling/bams/{sample}/{ID}.splitters.bam",
    params:
        outprefix="/u04/sv_calling/bams/{sample}/{ID}",
        # PL has to be one of ILLUMINA，SLX，SOLEXA，SOLID，454，LS454，COMPLETE，
        # PACBIO，IONTORRENT，CAPILLARY，HELICOS, or UNKNOWN
        # ID：输入reads集的ID号； LB： reads集的文库名； SM：样本名称； PL：测序平台
        RG="@RG\tID:{sample}.{ID}\tPL:illumina\tSM:{sample}\tLB:{ID}"
    run:
        shell("/home/fangzq/github/speedseq/bin/speedseq align -o {params.outprefix} -R '{params.RG}' {input.genome} {input.r1} {input.r2}")
        shell("echo -e {params.RG} >> /u04/sv_calling/bams/{sample}/all.LBIDs.txt")


rule align22:
    input:
        r1="/home/mzheng/rawData_NGS/{sample}/{ID}_1.fq.gz",
        r2="/home/mzheng/rawData_NGS/{sample}/{ID}_2.fq.gz",
        genome=GENOME, # bwa indexed
    output: 
        bam="/u04/sv_calling/bams/{sample}/{ID}.bam",
        discordants="/u04/sv_calling/bams/{sample}/{ID}.discordants.bam",
        splitters="/u04/sv_calling/bams/{sample}/{ID}.splitters.bam",
    params:
        outprefix="/u04/sv_calling/bams/{sample}/{ID}",
        # PL has to be one of ILLUMINA，SLX，SOLEXA，SOLID，454，LS454，COMPLETE，
        # PACBIO，IONTORRENT，CAPILLARY，HELICOS, or UNKNOWN
        # ID：输入reads集的ID号； LB： reads集的文库名； SM：样本名称； PL：测序平台
        RG="@RG\tID:{sample}.{ID}\tPL:illumina\tSM:{sample}\tLB:{ID}"
    run:
        shell("/home/fangzq/github/speedseq/bin/speedseq align -o {params.outprefix} -R '{params.RG}' {input.genome} {input.r1} {input.r2}")
        shell("echo -e {params.RG} >> /u04/sv_calling/bams/{sample}/all.LBIDs.txt")


rule align3:
    input:
        r1="/home/mzheng/rawData_NGS2/Sanger/{sample}/{ID}_1.fastq.gz",
        r2="/home/mzheng/rawData_NGS2/Sanger/{sample}/{ID}_2.fastq.gz",
        genome=GENOME, # bwa indexed
    output: 
        bam="/u04/sv_calling/bams/{sample}/{ID}.bam",
        discordants="/u04/sv_calling/bams/{sample}/{ID}.discordants.bam",
        splitters="/u04/sv_calling/bams/{sample}/{ID}.splitters.bam",
    params:
        outprefix="/u04/sv_calling/bams/{sample}/{ID}",
        # PL has to be one of ILLUMINA，SLX，SOLEXA，SOLID，454，LS454，COMPLETE，
        # PACBIO，IONTORRENT，CAPILLARY，HELICOS, or UNKNOWN
        # ID：输入reads集的ID号； LB： reads集的文库名； SM：样本名称； PL：测序平台
        RG="@RG\tID:{sample}.{ID}\tPL:illumina\tSM:{sample}\tLB:{ID}"
    run:
        shell("/home/fangzq/github/speedseq/bin/speedseq align -o {params.outprefix} -R '{params.RG}' {input.genome} {input.r1} {input.r2}")
        shell("echo -e {params.RG} >> /u04/sv_calling/bams/{sample}/all.LBIDs.txt")

#BTBR
rule align4:
    input:
        r1="/home/mzheng/rawData_NGS2/Sanger/BTBR/shuffled_version_201604/{ID}_1.fastq.gz",
        r2="/home/mzheng/rawData_NGS2/Sanger/BTBR/shuffled_version_201604/{ID}_2.fastq.gz",
        genome=GENOME, # bwa indexed
    output: 
        bam="/u04/sv_calling/bams/BTBR/{ID}.bam",
        discordants="/u04/sv_calling/bams/BTBR/{ID}.discordants.bam",
        splitters="/u04/sv_calling/bams/BTBR/{ID}.splitters.bam",
    params:
        outprefix="BTBR/{ID}",
        # PL has to be one of ILLUMINA，SLX，SOLEXA，SOLID，454，LS454，COMPLETE，
        # PACBIO，IONTORRENT，CAPILLARY，HELICOS, or UNKNOWN
        # ID：输入reads集的ID号； LB： reads集的文库名； SM：样本名称； PL：测序平台
        RG="@RG\tID:BTBR.{ID}\tPL:illumina\tSM:BTBR\tLB:{ID}"
    run:
        shell("/home/fangzq/github/speedseq/bin/speedseq align -o {params.outprefix} -R '{params.RG}' {input.genome} {input.r1} {input.r2}")
        shell("echo -e {params.RG} >> /u04/sv_calling/bams/BTBR/all.LBIDs.txt")