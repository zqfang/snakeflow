### pacbio sequencing
import glob, os

workdir: "/data/bases/fangzq/20200815_SV/PacBio_project"

GENOME = "/home/fangzq/genome/mouse/GRCm38_68.fa"
FASTQs = glob.glob("RAW-DATA/*fastq.gz")
SAMPLES = sorted([b.split("/")[-1].split(".")[0] for b in glob.glob("BAMs/*sorted.MD.bam")])

### Sniffiles parameters #####
MIN_HOMO_AF = 0.7
MIN_HET_AF = 0.3
MIN_SUPPORT_READ = 8
MIN_SVLEN = 50


### Output ###################
SVs = expand("VCFs/{sample}.raw.vcf", sample=SAMPLES)
MERGED_CALLSET = "VCFs/merged_gt_SURVIVOR_1kbpdist.vcf"
MERGED_FILTER = "VCFs/merged_gt_SURVIVOR_1kbpdist.filter.vcf"
MERGED_VEP =  "VEPs/merged_gt_SURVIVOR_1kbpdist.pass.txt.gz"
VEP_OUT = expand("VEPs/{sample}.pass.vep.txt.gz", sample=SAMPLES)
BAM_STATS = expand("BAMs/{sample}.bamstat.txt", sample=SAMPLES)
COVERAGE = "BAMs/samples.coverage.txt"

rule target:
    input: SVs, MERGED_CALLSET, VEP_OUT, MERGED_VEP, BAM_STATS, COVERAGE

rule ngmlr:
    input: 
        ref= GENOME,
        fastq="{sample}.fastq.gz",
    output: "BAMs/{sample}.ngmlr.bam"
    threads: 12
    log: "logs/{sample}.ngmlr.log"
    params:
        extra = "-x ont" # nanopore 
    shell:
        "ngmlr -t {threads} -r {input.ref} -q {input.fastq} | "
        "samtools view -Sbh -@ {threads}  -o {output} - ) 2> {log}"


rule minimap2:
    input: 
        ref= GENOME,
        fastq="{sample}.fastq.gz"
    output: "BAMs/{sample}.bam"
    threads: 12
    log: "logs/{sample}.minimap2.log"
    shell:
        "(minimap2 -c -a --MD -x map-pb -t {threads} {input.ref} {input.fastq} |"
        "samtools view -Sbh -@ {threads}  -o {output} - ) 2> {log}"


rule bam_sort:
    input: "BAMs/{sample}.bam"
    output: protected("BAMs/{sample}.sorted.MD.bam")
    threads: 8
    shell: "samtools sort -@ {threads} {input} > {output}"

rule bam_index:
    input: "BAMs/{sample}.sorted.MD.bam"
    output: "BAMs/{sample}.sorted.MD.bam.bai"
    shell: "samtools index {input}"

# give mosdepth a try, much faster than samtools depth
# rule Coverage:
#     input:
#         genome=GENOME,
#         bam="BAMs/{sample}.sorted.MD.bam",
#         bami="BAMs/{sample}.sorted.MD.bam.bai"
#     output:
#          "BAMs/{sample}.avg.coverage.txt"
#     shell:
#         ## positions with zero coverage omitted if not -a
#         "samtools depth -a --reference {input.genome} {input.bam} | "
#         "awk '{{sum+=$3; if($3>0) total+=1}} END "
#         "{{print \"Depth (of covered bases) =\"sum/total; "
#         "print \"Depth (of genome size) =\"sum/NR; "
#         "print \"Breath = \"(total/NR)*100}}' "
#         #"awk '{{sum+=$3}} END {{ print \"Average Depth = \"sum/NR}}' "
#         " > {output}"

rule Coverage:
    input: 
        genome=GENOME,
        bam=expand("BAMs/{sample}.sorted.MD.bam", sample=SAMPLES),
        bami=expand("BAMs/{sample}.sorted.MD.bam.bai", sample=SAMPLES),
    output: COVERAGE
    params:
        samples=",".join(SAMPLES),
    shell:
        "samtools depth -a --reference {input.genome} {input.bam}  | "
        " awk -v sample={params.samples} '{{ for (i=3;i<=NF;i++) {{ sum[i-2]+=$i; if($i>0) total[i-2]+=1;}} }} "
        " END {{ split(sample,arr,\",\"); print \"Sample DepthOfCovered DepthOfGenome Breath\"; " 
        " for (i in sum) print arr[i], sum[i]/total[i], sum[i]/NR, total[i]/NR*100 }}' "
        " > {output} "

rule bam_stats:
    input: 
        genome=GENOME,
        bam="BAMs/{sample}.sorted.MD.bam",
        bami="BAMs/{sample}.sorted.MD.bam.bai"
    output: 
        "BAMs/{sample}.bamstat.txt"
    threads: 8
    shell:
        "samtools stats --reference {input.genome} -i 20000 "
        "--threads {threads} {input.bam} | grep '^SN' > {output}"


rule sniffles_single:
    input: "BAMs/{sample}.sorted.MD.bam"
    output: protected("VCFs/{sample}.raw.vcf")
    threads: 8
    log: "logs/{sample}.sniffles.single.log"
    params:
        svlen = MIN_SVLEN,
        sr = MIN_SUPPORT_READ,
        homo_af = MIN_HOMO_AF, ## if DV / (DR + DV) > homo_af, assign 1/1
        het_af = MIN_HET_AF, ## if DV / (DR + DV) > homo_af, assign 0/1, else 0/0
        binary = "/home/fangzq/github/Sniffles/bin/sniffles-core-1.0.12"
    shell:
        "{params.binary}/sniffles --genotype -m {input} -v {output} -t {threads} "
        "--min_length {params.svlen} --min_support {params.sr} "
        "--min_homo_af {params.homo_af} --min_het_af {params.het_af} "
        "--cluster 1>{log}"


rule list_raw:
    input: expand("VCFs/{sample}.raw.vcf", sample=SAMPLES),
    output: temp("VCFs/all.raw.txt")
    params:
        samples = expand("VCFs/{sample}.raw.vcf", sample=SAMPLES)
    run:
        with open(output[0], 'w') as out:
            for sample in input:
                out.write(sample+"\n")

rule SURVIOR_merge:
    input:
        vcf=expand("VCFs/{sample}.raw.vcf", sample=SAMPLES),
        txt="VCFs/all.raw.txt"
    output: temp("VCFs/merged.temp.vcf")
    shell:
        "SURVIVOR merge {input.txt} 1000 1 1 -1 -1 -1 {output}"

rule sniffles_joint:
    input: 
        bam = "BAMs/{sample}.sorted.MD.bam",
        ivcf = "VCFs/merged.temp.vcf"
    output: "VCFs/{sample}.joint.vcf"
    threads: 8
    log: "logs/{sample}.sniffles.joint.log"
    params:
        svlen = MIN_SVLEN,
        sr = MIN_SUPPORT_READ,
        homo_af = MIN_HOMO_AF, ## if DV / (DR + DV) > homo_af, assign 1/1
        het_af = MIN_HET_AF, ## if DV / (DR + DV) > homo_af, assign 0/1, else 0/0
        binary = "/home/fangzq/github/Sniffles/bin/sniffles-core-1.0.12"
    shell:
        "{params.binary}/sniffles --genotype -m {input.bam} -v {output} -t {threads} "
        "--min_length {params.svlen} --min_support {params.sr} "
        "--min_homo_af {params.homo_af} --min_het_af {params.het_af} "
        "--cluster --Ivcf {input.ivcf} 1>{log}"

rule list_joint:
    input: 
        expand("VCFs/{sample}.joint.vcf", sample=SAMPLES),
    output: 
        temp("VCFs/all.joint.txt"),
        "sample.order.txt",
    run:
        with open(output[0], 'w') as out:
            for sample in input:
                out.write(sample+"\n")
        with open(output[1], 'w') as s:
            s.write("\n".join(SAMPLES) +"\n")

rule SURVIOR_merged_callset:
    input:
        vcf=expand("VCFs/{sample}.joint.vcf", sample=SAMPLES),
        txt="VCFs/all.joint.txt"
    output: temp("VCFs/merged_gt_temp.vcf")
    shell:
        "SURVIVOR merge {input.txt} 1000 -1 1 -1 -1 -1 {output}"


# rule bcftools_filter:
#     input: MERGED_CALLSET
#     output:
#     params: 
#     shell:
#         "bcftools view -f PASS,.  "
#         "-i 'GT==\"AA\" && (INFO/SVLEN>50 || INFO/SVLEN<-50) && INFO/SVTYPE!=\"BND\"' "
#         "-o {output} {input} "

rule bcftools_reheader:
    """
    handy tool to rename sample names, see docs 
    here: http://samtools.github.io/bcftools/bcftools-man.html#reheader
    """
    input: 
        vcf="VCFs/merged_gt_temp.vcf",
        samples="sample.order.txt"
    output:
        MERGED_CALLSET
    shell:
        "bcftools reheader -s {input.samples} -o {output} {input.vcf}"

rule ensemble_vep:
    input: "VCFs/{sample}.raw.vcf",
    output: "VEPs/{sample}.pass.vep.txt.gz"
    threads: 8
    params:
        svlen = 50,
        genome_build="GRCm38",
        species="mus_musculus",
        BIN= "/home/fangzq/github/ensembl-vep"
    shell:
        "bcftools view -f PASS,.  "
        "-i 'GT==\"AA\" && INFO/SVTYPE!=\"BND\"' {input} | "
        "{params.BIN}/vep -a {params.genome_build} --species {params.species} "
        "--fork {threads} --offline --uniprot --cache --format vcf --force_overwrite -overlaps "
        "--plugin TSSDistance --domains --plugin StructuralVariantOverlap "
        "--plugin phenotypes --symbol --gencode_basic --gene_phenotype MGI"
        "--nearest gene --regulatory --distance 5000 --individual all "
        "--no_check_variants_order --max_sv_size 100000 " # --check_svs " 
        "--compress_output gzip -o {output} --tab"

rule merged_callset_filtering:
    input: MERGED_CALLSET
    output: MERGED_FILTER
    run:
        out = open(output[0], 'w')
        header = []
        ref_idx = -1
        ref = "C57BL6J"
        with open(input[0], 'r') as callset:
            for line in callset:
                if line.startswith("#"):
                    out.write(line)
                    if line.startswith("#CHROM"):
                        header += line.strip("#\n").split("\t")
                        ref_idx = header.index(ref)
                    continue
                record = line.strip().split("\t")
                if record[6] != "PASS":
                    continue
                if record[0] in ['Y','MT']:
                    continue
                aa = 0
                rr = 0
                ra = 0
                for gt in record[9:]:
                    if gt.startswith("1/1"):
                        aa +=1
                    if gt.startswith("0/0") or  gt.startswith("./."):
                        rr +=1   
                    if gt.startswith("0/1"):
                        ra += 1 
                if rr == len(record[9:]): # remove entries with if all samples are  0/0 (./.)
                    continue
                if aa < 1:  # at least one 1/1
                    continue  
                sample_size = len(record[9:])
                if ra == sample_size or aa == sample_size or rr == sample_size :
                # if all samples are same genotype to reference genome (C57BL6J) , it's atifical SVs
                    continue

                # if the genotype of reference genome (C57BL6J) is 1/1 , it's atifical SVs 
                if (ref_idx != -1) and record[ref_idx].startswith("1/1"):
                    continue

                if record[7].find("SVTYPE=TRA") >= 0:
                    continue
                if record[7].find("SVTYPE=BND") >= 0:
                    continue                
                out.write(line)
        out.close()
    
rule ensemble_vep_joint:
    input: MERGED_FILTER,
    output: MERGED_VEP
    threads: 8
    params:
        svlen = 50,
        genome_build="GRCm38",
        species="mus_musculus",
        BIN= "/home/fangzq/github/ensembl-vep"
    shell:
        "{params.BIN}/vep -a {params.genome_build} --species {params.species} "
        "--fork {threads} --offline --uniprot --cache --format vcf --force_overwrite -overlaps "
        "--plugin TSSDistance --domains --plugin StructuralVariantOverlap "
        "--plugin phenotypes --symbol --gencode_basic --gene_phenotype MGI"
        "--nearest gene --regulatory --distance 5000 --individual all "
        "--no_check_variants_order --max_sv_size 500000 " # --check_svs " 
        "--compress_output gzip -o {output} --tab -i {input}"




# bcftools view -f PASS -i 'INFO/SVTYPE!="BND"' strains39.all.vcf.gz | \
# /home/fangzq/github/ensembl-vep/vep -a GRCm38 --species mus_musculus \
# --compress_output gzip -o {output} --tab \
# --fork 4 --offline --uniprot --cache --format vcf --force_overwrite -overlaps \
# --plugin TSSDistance --domains --plugin StructuralVariantOverlap \
# --plugin phenotypes --symbol --gencode_basic --gene_phenotype MGI \
# --nearest gene --regulatory --distance 5000 --individual all \
# --no_check_variants_order --max_sv_size 100000
