import glob, os

workdir: "/data/bases/fangzq/strains"
GENOME = "/home/fangzq/genome/mouse/GRCm38_68.fa"
BAMs = glob.glob("./*/output.GATKrealigned.Recal.bam")
STRAINS = sorted([ b.split("/")[1] for b in BAMs])
COVERAGE = expand("{sample}/{sample}.avg.coverage.txt", sample=STRAINS) 

rule target:
    input: "samples.depth.txt"


rule genome_size:
    """
    use this to get genome size instead of awk's NR
    """
    input:
        bam="{sample}/output.GATKrealigned.Recal.bam",
        bami="{sample}/output.GATKrealigned.Recal.bam.bai"
    output:
         temp("{sample}/{sample}.chrom.size")
    shell:
        "samtools view -H {input.bam} | grep -P '^@SQ' | "
        "cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}' > {output}"

        
# rule Coverage_single:
#     input:
#         genome=GENOME,
#         bam="{sample}/output.GATKrealigned.Recal.bam",
#         bami="{sample}/output.GATKrealigned.Recal.bam.bai"
#     output:
#          "{sample}/{sample}.avg.coverage.txt"
#     shell:
#         ## positions with zero coverage omitted if not -a
#         "samtools depth -a --reference {input.genome} {input.bam} | "
#         "awk '{{sum+=$3; if($3>0) total+=1}} END "
#         "{{print \"Depth (of covered bases) =\"sum/total; "
#         "print \"Depth (of genome size) =\"sum/NR; "
#         "print \"Breath = \"(total/NR)*100}}' "
#         " > {output}"

rule Coverage:
    input: 
        genome=GENOME,
        bam=expand("{sample}/output.GATKrealigned.Recal.bam", sample=STRAINS),
        bami=expand("{sample}/output.GATKrealigned.Recal.bam.bai", sample=STRAINS)
    output: "samples.depth.txt"
    params:
        samples=",".join(STRAINS),
    shell:
        "samtools depth -a --reference {input.genome} {input.bam}  | "
        " awk -v sample={params.samples} '{{ for (i=3;i<=NF;i++) {{ sum[i]+=$i; if($i>0) total[i]+=1;}} }} "
        " END {{ split(sample,arr,\",\"); print \"Sample\tDepthOfCovered\tDepthOfGenome\tBreath\"; " 
        " for (i in sum) print arr[i], sum[i]/total[i], sum[i]/NR, total[i]/NR*100 }}' "
        " > {output} "
        