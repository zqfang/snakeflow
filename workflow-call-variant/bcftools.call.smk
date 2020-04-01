

# samtools-bcftools-calling
rule faidx: 
    input: GENOME
    output: GENOME + ".fai"
    shell:
        "samtools faidx input.fa"

rule getChromSize:
    input: rule.faidx.output
    output: temp(expand("chr{i}.txt",i=CHROMOSOME))
    params:
        #strng='{print "chr"$1":1-"$2}'
        strng='{print $1":1-"$2}'
    shell:
        "awk -F'\t' '{params.string}' {input} | while read line;"
        "do echo -e $line > region.$line done"

rule bcftoolsCalling:
    input:
        genome = GENOME,
        bam = expand("BAM/{sample}.marked.fixed.BQSR.bam", sample=STRAINS),
        chroms_region="region.{region}", # chrom region file
    output: "VCFs/combined.{region}.raw.vcf"
    params:
        chrom = lambda wildcards, output: output[0].split(":")[0] # <- "Y"
        region="{region}", # <- "Y:1-91744698"
        bam = " ".join(expand("BAM/{sample}.marked.fixed.BQSR.bam", sample=STRAINS))
    threads: 8
    shell:
        "bcftools mpileup -t DP,AD,ADF,ADR,SP,INFO/AD --threads {threads} " 
        "-E -Q 0 -p -m3 -F0.25 -d500 -r {params.region} "
        "-Ou -f {input.genome} {params.bam} | "
        "bcftools call --threads {threads} -mv -f GQ,GP -Ov  > {output}"

rule tabix:
    input: "VCFs/combined.{region}.raw.vcf"
    output: "VCFs/combined.{region}.raw.vcf.tbi"
    shell:
        "tabix -p vcf {input}"

rule bcftools_stats:
    input: 
        genome=GENOME,
        vcf="VCFs/combined.{region}.raw.vcf"
    output: "VCFs/combined.{region}.raw.vcf.stats"
    shell:
        "bcftools stats -F {input.genome} -s - {input.vcf} > {output}"

rule bcftools_plot:
    input: "VCFs/combined.{region}.raw.vcf.stats"
    output: "VCFs/combined.{region}.raw.vcf.stats.pdf"
    shell:
        "plot-vcfstats -p plots {input}"

rule bcfcall_filtering:
    input: 
        vcf="VCFs/combined.{region}.raw.vcf",
        vcfi="VCFs/combined.{region}.raw.vcf.tbi"
    output: "VCFs/combined.{region}.filter.vcf.gz"
    shell: 
        "bcftools filter -Oz -o {ouput} -s LOWQUAL -i'%QUAL>10' {input}"