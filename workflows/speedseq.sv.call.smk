
import os, glob
############### Globals ########################
workdir: "/data/bases/fangzq/20200815_SV"
SAMPLES = ['129P2', '129S1', '129S5', 'AKR', 'A_J', 'B10', 
        'BPL', 'BPN', 'BTBR', 'BUB', 'B_C', 'C3H', 'C57BL10J',
        'C57BL6NJ', 'C57BRcd', 'C57LJ', 'C58', 'CBA', 'CEJ', 
        'DBA', 'DBA1J', 'FVB', 'ILNJ', 'KK', 'LGJ', 'LPJ', 
        'MAMy', 'MRL','NOD', 'NON', 'NOR', 'NUJ', 'NZB', 'NZO', 'NZW', 
        'PJ', 'PLJ', 'RFJ', 'RHJ', 'RIIIS', 'SEA', 'SJL', 'SMJ', 'ST', 'SWR', 'TALLYHO', 'RBF'] + \
        ['CAST', 'MOLF', 'PWD','PWK', 'SPRET', 'WSB']  # <- wild derived

BAM_DIR= "/data/bases/fangzq/strains"
GENOME = "/home/fangzq/genome/mouse/GRCm38_68.fa"
SPEESEQ = "/home/fangzq/github/speedseq/bin"
EXCLUDE = "/home/fangzq/github/snakeflow/data/mouse.mm10.excl.bed"

### use svtools env ####
##  snakemake --use-conda --conda-prefix /home/fangzq/miniconda/envs/sv


################################ rules ###################################

rule target:
    input: 
        #expand("bams/{sample}.realign.bam", sample=SAMPLES),
        "svtools/merged.sv.pruned.vcf.gz",
        #"svtools/output.largesample.vcf.gz",
        "svtools/merged.sv.pruned.annot.txt"

rule realign:
    input:
        bam = os.path.join(BAM_DIR, "{sample}/output.GATKrealigned.Recal.bam"),
        genome = GENOME
    output:
        bam = protected("bams/{sample}.realign.bam"),
        discordant_bam= protected("bams/{sample}.realign.discordants.bam"),
        splitter_bam = protected("bams/{sample}.realign.splitters.bam"), 
    params:
        outprefix = "bams/{sample}.realign",
        BIN=SPEESEQ,
    threads: 8
    shell:
        "{params.BIN}/speedseq realign -o {params.outprefix} "
        "-t {threads} -M 16  {input.genome} {input.bam}" 

## multisample sv calling
# rule callStructVariants:
#     input:
#         bam = expand("bams/{sample}.realign.bam", sample=SAMPLES),
#         discordant_bam=expand("bams/{sample}.realign.discordants.bam",sample=SAMPLES),
#         splitter_bam =expand("bams/{sample}.realign.splitters.bam", sample=SAMPLES),
#         genome=GENOME,
#     output:
#         "strains.analysis.sv.vcf.gz",
#         "strains.analysis.sv.vcf.gz.tbi"
#     params:
#         outprefix="strains.analysis",
#         BIN=SPEESEQ,
#     shell:
#         "{params.BIN}/speedseq sv -o {params.outprefix} -B {input.bam} "
#         "-D {input.discordant_bam} -S {input.splitter_bam}" 
#         "-v -d -P -g -k"  

## single sample sv calling
rule callStructVariantsSingle:
    input:
        bam = "bams/{sample}.realign.bam",
        discordant_bam="bams/{sample}.realign.discordants.bam",
        splitter_bam ="bams/{sample}.realign.splitters.bam", 
        exclude=EXCLUDE, 
        genome=GENOME,
    output:
        protected("lumpy/{sample}.sv.vcf.gz"),
        protected("lumpy/{sample}.sv.vcf.gz.tbi"),
        "lumpy/{sample}.sv.{sample}.realign.bam.readdepth.bed",
        "SVTEMP_{sample}/cnvnator-temp/{sample}.realign.bam.hist.root",
    params:
        outprefix="lumpy/{sample}",
        BIN=SPEESEQ,
        tempdir="SVTEMP_{sample}",
    log: "logs/svtools/sv/{sample}.log"
    shell:
        # regions to exclude for mouse, yeast, drosophila: 
        # https://github.com/dellytools/delly/tree/master/excludeTemplates
        #  excluded regions of the reference genome with consistently 
        #  high sequencing depth over multiple individuals, 
        #  since high depth is indicative of artifacts in the reference assembly
        # We used BEDTools v2.17.0 to generate a BED graph of aggregate per-base coverage 
        # from all 17 individuals [23]. 
        # The mode and standard deviation of the aggregate depth were calculated 
        # separately for the autosomes and sex chromosomes. 
        # Any regions with depth exceeding 2 * mode + 3 standard deviations were excluded from our analyses. 
        # (We chose to double the mode to allow inclusion of duplicated copy number variant regions.)    
        # abnormally high coverage, the mitochondrial genome, the decoy genome and the genome of the Epstein-Barr virus (EBV)
        "{params.BIN}/speedseq sv -o {params.outprefix} -R {input.genome} "
        "-B {input.bam} -D {input.discordant_bam} -S {input.splitter_bam} "
        "-T {params.tempdir} -v -d -P -g -k 2> {log}"   # -x {input.exclude}


### use svtools to create a callset
rule svtools_lsort:
    # takes a space separated list of all of the LUMPY VCF files
    input: expand("lumpy/{sample}.sv.vcf.gz", sample=SAMPLES)
    output: temp("svtools/sorted.vcf.gz")
    # conda: "envs/svtools.yaml"
    shell:
        "svtools lsort {input} | bgzip -c > {output}"

rule svtools_lmerge:
    # merge variant calls likely representing the same variant in the sorted VCF
    input: "svtools/sorted.vcf.gz"
    output: "svtools/merged.vcf.gz"
    #conda: "envs/svtools.yaml"
    shell:
        "zcat {input} | svtools lmerge -i /dev/stdin -f 20 | "
        "bgzip -c > {output}"

rule svtools_genotype:
    # genotype all samples for all variants present in the merged set
    input: 
        vcf="svtools/merged.vcf.gz",
        bam="bams/{sample}.realign.bam",
    output: "gt/{sample}.vcf"
    # conda:  "envs/svtools.yaml"
    params: BIN=SPEESEQ
    shell:    
        "zcat {input.vcf} | {params.BIN}/vawk --header '{{ $6=\".\"; print }}'  "
        "| svtools genotype -B {input.bam} -l {input.bam}.json  "
        "| sed 's/PR...=[0-9\.e,-]*\(;\)\{{0,1\}}\(\\t\)\{{0,1\}}/\2/g' - "
        "> {output}"

rule create_coordinates:
    input: "svtools/merged.vcf.gz",
    output: 
        vcf=temp("svtools/merged.vcf"),
        coord=temp("coordinates"),
    # conda: "envs/svtools.yaml"
    shell: 
        """
        zcat {input} > {output.vcf}
        create_coordinates -i {output.vcf} -o {output.coord}
        """

rule svtools_copynumber:
    # https://github.com/hall-lab/svtools/blob/master/Tutorial.md 
    # need DNVnator, details see the link
    input:  
        gt="gt/{sample}.vcf",
        coord="coordinates" ,
        cnvnator="SVTEMP_{sample}/cnvnator-temp/{sample}.realign.bam.hist.root" 
    output: 
        "cn/{sample}.cn.vcf"
    # conda: "envs/svtools.yaml"
    params: BIN=SPEESEQ
    shell:
        ## annotate variants with copynumber from cnvnator
        ## CNVnator is run as part of speedseq sv
        "svtools copynumber --cnvnator {params.BIN}/cnvnator -s {wildcards.sample} "
        "-w 100 -r {input.cnvnator} -c {input.coord} -i {input.gt} "
        "> {output}"

rule copynumber_list:
    input: 
        vcf = expand("cn/{sample}.cn.vcf", sample=SAMPLES)
    output: temp("cn/cn.list")
    # conda: "envs/svtools.yaml"
    run:
        with open(output[0], 'w') as cn:
            cn.write("\n".join(input.vcf))
        

rule svtools_vcfpaste:
# re-assemble a cohort-level VCF file
    input: 
        vcf="svtools/merged.vcf",
        cnlist="cn/cn.list"
    # conda: "envs/svtools.yaml"
    output:
        temp("svtools/merged.sv.gt.cn.vcf.gz")
    shell:
        "svtools vcfpaste -m {input.vcf} -f {input.cnlist} -q "
        "| bgzip -c > {output}"

rule svtools_prune:
# filter out additional variants deemed to be identical
    input: "svtools/merged.sv.gt.cn.vcf.gz"
    output: protected("svtools/merged.sv.pruned.vcf.gz")
    # conda: "envs/svtools.yaml"
    shell:
        "zcat {input} "
        "| svtools afreq "
        "| svtools vcftobedpe "
        "| svtools bedpesort "
        "| svtools prune -s -d 100 -e 'AF' "
        "| svtools bedpetovcf "
        "| bgzip -c > {output}"

# svtools classify:
#  refine genotypes and SV types
# https://github.com/hall-lab/svtools/blob/master/Tutorial.md

        
# All svtools classify commands require a BED file of repeats for classifying Mobile Element Insertions (MEI). 
# This can be created from the UCSC genome browser.
# example for hg19

# curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz \
# | gzip -cdfq \
# | awk '{ gsub("^chr", "", $6); if ($3<200) print $6,$7,$8,$12"|"$13"|"$11,$3,$10 }' OFS="\t" \
# | sort -k1,1V -k2,2n -k3,3n \
# | awk '$4~"LINE" || $4~"SINE" || $4~"SVA"' \
# | bgzip -c > repeatMasker.recent.lt200millidiv.LINE_SINE_SVA.b37.sorted.bed.gz

# hg38
# curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz \
# | gzip -cdfq \
# | awk '{ if ($3<200) print $6,$7,$8,$12"|"$13"|"$11,$3,$10 }' OFS="\t" \
# | sort -k1,1V -k2,2n -k3,3n \
# | awk '$4~"LINE" || $4~"SINE" || $4~"SVA"' \
# | bgzip -c > repeatMasker.recent.lt200millidiv.LINE_SINE_SVA.GRCh38.sorted.bed.gz


# mm10
# curl -s http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz \
rule repeat_elements:
    output: "repeatMasker.recent.lt200millidiv.LINE_SINE_SVA.mm10.sorted.bed.gz" 
    params:
        genome_build = "mm10", # hg38, hg19. mm10
    shell:
        "curl -s http://hgdownload.cse.ucsc.edu/goldenPath/{params.genome_build}/database/rmsk.txt.gz "
        "| gzip -cdfq "
        "| awk '{{ gsub(\"^chr\", \"\", $6); if ($3<200) print $6,$7,$8,$12\"|\"$13\"|\"$11,$3,$10 }}' OFS=\"\\t\" "
        "| sort -k1,1V -k2,2n -k3,3n "
        "| awk '$4~\"LINE\" || $4~\"SINE\" || $4~\"SVA\"' "
        "| bgzip -c > {output}"



#svtools requires a tab-delimited file to specify the number of X chromosome copies
rule chrXnumbers:
    output: "strain.chrXnum.txt"
    params:
        strain=SAMPLES,
        #cs=SAMPLES_cn,
    run:
        # for inbred mice sequenced by us, all are females, so 2 copies.
        with open(output[0],'w') as out:
            for s in params.strain:
                out.write(f"{s}\t2\n")
            
# a file of high-quality, simple deletions and duplications
# https://github.com/hall-lab/svtools/blob/master/Classifier.md
# # see the links to create custom training data
# rule trainingdata:
#     output: "training_vars.bedpe.gz"
#     shell:
#         "curl -O https://raw.githubusercontent.com/ernfrid/svtools/classifier_docs/resources/training_vars.bedpe.gz"


# Creating a set of high-quality, simple deletions and duplications for classifier training
# A VCF from a large cohort (>30 samples) should be used.
# Find the mean per-site copy number and overall percentiles for deletions and duplications
rule find_cn:
    input: "svtools/merged.sv.pruned.vcf.gz"
    output:
        mean = "per_site.means.txt",
        percenttile = "overall_percentiles.txt"
    params:
        BIN = SPEESEQ,
    shell:
        "zcat {input} | {params.BIN}/vawk "
        "--header '{{if((I$SVTYPE==\"DEL\" || I$SVTYPE==\"DUP\" || I$SVTYPE==\"MEI\") && I$AF>0.01 && $1!=\"X\" && $1!=\"Y\") print $0}}'  "
        "| perl ../scripts/mean_cn.pl 1> {output.mean} 2> {output.percentile}"


rule extract_highquality_dups_dels:
    input: 
        mean = "per_site.means.txt",
        percenttile = "overall_percentiles.txt",        
    output: "training_vars.bedpe.gz"
    params:
    shell:
        "cat {input.mean} | cut -f -8 "
        "| zjoin -a stdin -b <(cat {input.percentile} | cut -f -8 ) -1 2 -2 1 "
        "| awk '{{if($5>$11 && $5<$12 && $8>$15 && $8<$16) print $0}}' "
        "| cut -f 1 | zjoin -a <(zcat sv.vcf.gz) -b stdin -1 3 -2 1 "
        "| svtools vcftobedpe | bgzip -c > {output}"

rule training_variants:
    input: 
        vcf="svtools/merged.sv.pruned.vcf.gz",
        vars="training_vars.bedpe.gz"
    output: "training.vars.vcf.gz"
    # conda: "envs/svtools.yaml"
    shell:
        "zcat {input.vcf} "
        "| svtools vcftobedpe  "
        "| svtools varlookup -a stdin -b {input.vars} -c HQ -d 50 "
        "| svtools bedpetovcf "
        "| svtools vcfsort "
        "| vawk --header '{{ if(I$HQ_AF>0) print $0 }}' "
        "| bgzip -c > {output}"


rule svtools_classify_naivebayes:
# https://github.com/hall-lab/svtools/blob/master/Classifier.md
    input:
        vcf="svtools/merged.sv.pruned.vcf.gz",
        repeatmasker="repeatMasker.recent.lt200millidiv.LINE_SINE_SVA.mm10.sorted.bed.gz",
        chrXcopy="strain.chrXnum.txt",
        train="training.vars.vcf.gz"
    output: "output.naivebayes.vcf.gz"
    # conda: "envs/svtools.yaml"
    shell:
        "zcat {input.vcf} |  svtools classify "
        "-g {input.chrXcopy} -a {input.repeatmasker} "
        "-m naive_bayes -t {input.train} | bgzip -c > {output}"

rule svtools_classify_largesample: 
# large sample mode is recommended when samples > 30
    input:
        vcf="svtools/merged.sv.pruned.vcf.gz",
        repeatmasker="repeatMasker.recent.lt200millidiv.LINE_SINE_SVA.mm10.sorted.bed.gz",
        chrXcopy="strain.chrXnum.txt",
    output: "svtools/output.largesample.vcf.gz"
    # conda: "envs/svtools.yaml"
    shell:
        "zcat {input.vcf} |  svtools classify "
        "-g {input.chrXcopy} -a {input.repeatmasker} "
        "-m large_sample | bgzip -c > {output}"

rule ensemble_vep:
    input: "svtools/merged.sv.pruned.vcf.gz",
    output: "svtools/merged.sv.pruned.annot.txt"
    threads: 8
    params:
        genome_build="GRCm38",
        species="mus_musculus",
        BIN= "/home/fangzq/github/ensembl-vep"
    shell:
        "{params.BIN}/vep -i {input} -o {output} --tab -a {params.genome_build} --species {params.species} "
        "--fork {threads} --offline --uniprot --cache --format vcf --force_overwrite -overlaps "
        "--plugin TSSDistance --domains --plugin StructuralVariantOverlap "
        "--plugin phenotypes --plugin miRNA --symbol "
        "--nearest gene --regulatory --distance 5000 "
        "--no_check_variants_order --dont_skip –max_sv_size 50000" # --check_svs 

rule annotSV:
    input: "svtools/merged.sv.pruned.vcf.gz",
    output: "svtools/merged.sv.pruned.annotSV.tsv"
    threads: 8
    params:
        genome_build="mm10",
        species="mus_musculus",
        BIN= "/home/fangzq/github/AnnotSV/bin"
    shell:
        "{params.BIN}/AnnotSV -SvinputFile {input} -genomeBuild {params.genome_build} "
        "-outputFile {output} -overwrite -promoterSize 2000"