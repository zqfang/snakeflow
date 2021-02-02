
import os
############### Globals ########################
# configfile: "config.yaml"
workdir: config['BCFTOOLS']['WORKSPACE']

GENOME = config['GENOME']
dbSNP = config['dbSNP']
BAM_DIR = config['BAM_DIR']
STRAINS = sorted(config['STRAINS'])

#CHROMSOME = [ str(c) for c in range(1,20)] + ["X", "Y", "MT"]
CHROMSOME = ['1'] + [ str(c) for c in range(10,20)] + [ str(c) for c in range(2,10)]+ ["X", "Y"]

# OUTPUT
VCF_PASS = expand("VCFs/combined.chr{i}.hardfilter.vcf.gz", i=CHROMSOME)
VEP_ANNO = expand("VEP/combined.chr{i}.hardfilter.vep.txt.gz", i=CHROMSOME)
VCF_RAW = expand("VCFs/combined.chr{i}.raw.vcf.gz", i=CHROMSOME)
VCF_STATS = expand("VCFs/combined.chr{i}.raw.vcf.gz.stats", i=CHROMSOME)

SNPDB = expand("SNPs/chr{i}.txt", i=CHROMSOME)

###########################################################################################
rule target:
    input: SNPDB,#VCF_RAW, VCF_PASS, SNPDB, #VEP_ANNO, #VCF_STATS

# samtools-bcftools-calling
rule faidx: 
    input: GENOME
    output: GENOME + ".fai"
    shell:
        "samtools faidx input.fa"

rule getChromSize:
    input: GENOME + ".fai"
    output: temp(expand("chr{i}.tmp",i=CHROMSOME))
    params:
        cmd='{print $1":1-"$2}',
        chrom=CHROMSOME
    #shell:
        # snakemake use bash -c mode 
        # Escape certain characters, such as \t by \\t, $ by \$, and { by {{.
        # Use triple quotation marks to surround the command line call.
        # """awk -F"\\t" '{{print $1":1-"$2}}' {input} | while read line;
        #    do echo -e $line > region.$line done
        # """
        # this works too
        # "awk -F'\\t' '{params.cmd}' {input} | while read line; "
        # "do echo -e $line > region.$line done"
    run:
        with open(input[0]) as cords:
            for c in cords:
                c1 = c.strip().split("\t")
                cmd = f"echo -e '{c1[0]}:1-{c1[1]}' > chr{c1[0]}.tmp"
                shell(cmd)
            
rule bcftools_call:
    input:
        genome = GENOME,
        bam = expand(os.path.join(BAM_DIR, "{sample}/output.GATKrealigned.Recal.bam"), sample=STRAINS),
        #bam = expand("BAM/{sample}.marked.fixed.BQSR.bam", sample=STRAINS),
        chroms_region="chr{i}.tmp"
    output: protected('VCFs/combined.chr{i}.raw.vcf.gz')
    params:
        #bam = " ".join(expand("BAM/{sample}.marked.fixed.BQSR.bam", sample=STRAINS))
        bam=" ".join(expand(os.path.join(BAM_DIR, "{sample}/output.GATKrealigned.Recal.bam"), sample=STRAINS))
    run:
        with open(input.chroms_region) as reg:
            region = reg.read().strip()
        ## cmds same to mouse genome project, but different here: "bcftools call -mv -f GQ,GP -p 0.99"
        ## see here: https://doi.org/10.1186/s13059-016-1024-y
        ## we used the updated annotations to replace the deprecated ones
        cmd = "bcftools mpileup "+\
              "-a DP,AD,ADF,ADR,SP,INFO/AD "+\  
              "-E -F0.25 -Q0 -p -m3 -d500 -r %s "%region +\
              "-Ou -f {input.genome} {params.bam} | " +\  
              "bcftools call -mv -f GQ,GP -Oz  > {output}"
        shell(cmd)

# rule bcftools_reheader:
#     """
#     handy tool to rename sample names, see docs 
#     here: http://samtools.github.io/bcftools/bcftools-man.html#reheader
#     """
#     input: 
#         vcf="VCFs/combined.chr{i}.raw.vcf.gz",
#         samples="sample.new.name.txt"
#     output:
#         "VCFs/combined.chr{i}.raw.vcf.gz"
#     shell:
#         "bcftools reheader -s -o {output} {input}"


rule bcftools_norm:
    """
    Left-align and normalize indels
    """
    input:  "VCFs/combined.chr{i}.raw.vcf.gz"
    output: temp("VCFs/combined.chr{i}.normed.vcf.gz")
    shell:
        "bcftools norm -d none -s -m+indels -Oz -o {output} {input} "


rule tabix:
    input: "VCFs/combined.{chr}.normed.vcf.gz"
    output: 
        temp("VCFs/combined.{chr}.normed.vcf.gz.tbi")
    shell:
        "tabix -p vcf {input} "
        
rule bcftools_stats:
    input: 
        genome=GENOME,
        vcf="VCFs/combined.{chr}.normed.vcf.gz",
        vcfi="VCFs/combined.{chr}.normed.vcf.gz.tbi"
    output: 
        "VCFs/combined.{chr}.normed.vcf.gz.stats"
    shell:
        "bcftools stats -F {input.genome} -s - {input.vcf} > {output}"

rule bcftools_plot:
    input: "VCFs/combined.{chr}.normed.vcf.gz.stats"
    output: "figures/combined.{chr}.normed.vcf.gz.stats.pdf"
    shell:
        "plot-vcfstats -p figures -T {wildcards.chr} {input}"

rule bcftools_filter:
    input: 
        vcf="VCFs/combined.{chr}.normed.vcf.gz",
        vcfi="VCFs/combined.{chr}.normed.vcf.gz.tbi"
    output: 
        protected("VCFs/combined.{chr}.hardfilter.vcf.gz")
    params:
        snp="-i 'TYPE=\"snp\" && %QUAL>20' ",  # if only snp need
        expression = "%QUAL<20 || INFO/MQ < 20 || "+\
                     "INFO/SB < 0.0001 || INFO/RPB < 0.0001 || "+\
                     "INFO/BQB < 0 || INFO/VDB < 0 || " +\
                     "%MAX(FORMAT/AD)<=3 || %MAX(FORMAT/AD)/%MAX(FORMAT/DP)<=0.3 " 
        ## MIN(DP)>5 ? both FORMAT and INFO have DP
        ## DP:  sum of the DP value over all samples 
    shell: 
        # select snp and indels
        "bcftools filter -Oz -o {output} -s LowQual -g3 -G10 " 
        "-e '%QUAL<20 || %MAX(FORMAT/AD)<=3 || %MAX(FORMAT/AD)/%MAX(FORMAT/DP)<=0.3' " 
        "{input.vcf}"


rule strainOrder:
    output: "strain.order.snpdb.txt"
    run:
        with open(output[0], 'w') as s:
            s.write("\n".join(STRAINS) +"\n")
        
rule snp2NIEHS:
    input:  
        strain = "strain.order.snpdb.txt",
        vcf = "VCFs/combined.chr{i}.hardfilter.vcf.gz",
    output: 
        protected("SNPs/chr{i}.txt")
    params:
        qual = config['BCFTOOLS']['qual'], 
        het = config['BCFTOOLS']['phred_likelihood_diff'],
        ad = config['BCFTOOLS']['allele_depth'],
        ratio = config['BCFTOOLS']['allele_mindepth_ratio'],
        mq = config['BCFTOOLS']['mapping_quality'],
        sb = config['BCFTOOLS']['strand_bias_pvalue'], 
        BIN = config['HBCGM']['BIN']# path to haplomap binary
    log: "logs/combined.chr{i}.snp2niehs.log"
    shell:
        "bcftools view -f .,PASS -v snps {input.vcf} | "
        "{params.BIN}/haplomap niehs -o {output} -a {params.ad} -r {params.ratio} "
        "-q {params.qual} -p {params.het} -m {params.mq} -b {params.sb} "
        "-s {input.strain} -v > {log}"


rule variantEeffectPrediction:
    """emsemble-vep"""
    input: 
        vcf="VCFs/combined.{chrom}.hardfilter.vcf.gz",
        reference=GENOME,
    output: "VEP/combined.{chrom}.hardfilter.pass.vep.txt.gz"
    params:
        #genome_build = " -a GRCm38 --species mus_musculus ",
        genome_build = config['VEP']['GENOME_BUILD'],
        VEPBIN = config['VEP']['BIN'],
        extra=" --dir_cache "  + config['VEP']['CACHE_DIR']
    threads: 1
    shell:
        ## emsemble-vep
        # https://github.com/Ensembl/ensembl-vep
        "bcftools view -f .,PASS {input.vcf} | "
        "{params.VEPBIN}/vep --fasta {input.reference} {params.genome_build} "
        "--format vcf --fork {threads} --hgvs --force_overwrite "
        "--uniprot --domains --symbol --regulatory --distance 1000 --biotype "
        "--gene_phenotype MGI --check_existing  --pubmed --numbers "
        "--offline --cache --variant_class "
        "--gencode_basic --no_intergenic --individual all "
        "-o {output} --tab --compress_output gzip"