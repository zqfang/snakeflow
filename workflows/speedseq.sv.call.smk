
import os, glob
############### Globals ########################
workdir: "/data/bases/fangzq/20200815_SV"
# SAMPLES = ['129P2', '129S1', '129S5', 'AKR', 'A_J', 'B10', 
#         'BPL', 'BPN', 'BTBR', 'BUB', 'B_C', 'C3H', 'C57BL10J',
#         'C57BL6NJ', 'C57BRcd', 'C57LJ', 'C58', 'CBA', 'CEJ', 
#         'DBA', 'DBA1J', 'FVB', 'ILNJ', 'KK', 'LGJ', 'LPJ', 
#         'MAMy', 'MRL','NOD', 'NON', 'NOR', 'NUJ', 'NZB', 'NZO', 'NZW', 
#         'PJ', 'PLJ', 'RFJ', 'RHJ', 'RIIIS', 'SEA', 'SJL', 'SMJ', 'ST', 'SWR', 'TALLYHO', 'RBF'] + \
#         ['CAST', 'MOLF', 'PWD','PWK', 'SPRET', 'WSB']  # <- wild derived except MRL

SAMPLES = ['129P2', '129S1', '129S5', 'A_J', 'B10', 
        'BPL', 'BPN', 'BTBR', 'BUB', 'B_C', 'C3H', 'C57BL10J',
        'C57BL6NJ', 'C57BRcd', 'C57LJ', 'C58', 'CBA', 'CEJ', 
        'DBA', 'ILNJ', 'KK', 'LGJ', 'LPJ', 
        'MAMy', 'MRL','NOD', 'NON', 'NUJ', 'NZB', 'NZO', 'NZW', 
        'PJ', 'PLJ', 'RFJ', 'RHJ', 'RIIIS', 'SEA', 'SJL', 'SMJ', 'ST', 'SWR'] + \
        ['CAST', 'MOLF', 'SPRET', 'WSB']  # <- wild derived except MRL
BAM_DIR= "/data/bases/fangzq/strains"
GENOME = "/home/fangzq/genome/mouse/GRCm38_68.fa"
SPEESEQ = "/home/fangzq/github/speedseq/bin/"
################################ rules ###################################

rule target:
    input: 
        expand("bams/{sample}.realign.bam", sample=SAMPLES),
        "strains.analysis.sv.vcf.gz",
        "svtools/merged.sv.pruned.vcf.gz"

rule realign:
    input:
        bam = os.path.join(BAM_DIR, "{sample}/output.GATKrealigned.Recal.bam"),
        genome = GENOME
    output:
        bam = "bams/{sample}.realign.bam",
        discordant_bam="bams/{sample}.realign.discordants.bam",
        splitter_bam ="bams/{sample}.realign.splitters.bam", 
    params:
        outprefix = "bams/{sample}.realign",
        BIN=SPEESEQ,
    threads: 8
    shell:
        "{params.BIN}/speedseq realign -o {params.outprefix} "
        "-t {threads} -M 16  {input.genome} {input.bam}" 

## multisample sv calling
rule callStructVariants:
    input:
        bam = expand("bams/{sample}.realign.bam", sample=SAMPLES),
        discordant_bam=expand("bams/{sample}.realign.discordants.bam",sample=SAMPLES),
        splitter_bam =expand("bams/{sample}.realign.splitters.bam", sample=SAMPLES),
        genome=GENOME,
    output:
        "strains.analysis.sv.vcf.gz",
        "strains.analysis.sv.vcf.gz.tbi"
    params:
        #bed2exclude="-x {input.annotation}",
        bams=",".join(expand("bams/{sample}.bam", sample=SAMPLES)),
        discordants=",".join(expand("bams/{sample}.discordants.bam",sample=SAMPLES)),
        splitters=",".join(expand("bams/{sample}.splitters.bam", sample=SAMPLES)),
        outprefix="strains.analysis",
        BIN=SPEESEQ,
    shell:
        "{params.BIN}speedseq sv -o {params.outprefix} -B {params.bam} "
        "-D {params.discordant_bam} -S {params.splitter_bam}" 

## single sample sv calling
rule callStructVariantsSingle:
    input:
        bam = "bams/{sample}.realign.bam",
        discordant_bam="bams/{sample}.realign.discordants.bam",
        splitter_bam ="bams/{sample}.realign.splitters.bam", 
        genome=GENOME,
    output:
        "sv/{sample}.sv.vcf.gz",
        "sv/{sample}.sv.vcf.gz.tbi"
    params:
        #bed2exclude="-x {input.annotation}",
        outprefix="sv/{sample}.sv",
        BIN=SPEESEQ
    shell:
        "{params.BIN}/speedseq sv -o {params.outprefix} -B {input.bam} "
        "-D {input.discordant_bam} -S {input.splitter_bam} -d -g" 


### use svtools to create a callset
rule svtools_lsort:
    # takes a space separated list of all of the LUMPY VCF files
    input:
        vcf = expand("bams/{sample}.realign.bam", sample=SAMPLES),
    output:
        "svtools/sorted.vcf.gz"
    conda: 
    shell:
        "svtools lsort {input} | bgzip -c > {output}"

rule svtools_lmerge:
    # merge variant calls likely representing the same variant in the sorted VCF
    input: "svtools/sorted.vcf.gz"
    output: "svtools/merged.vcf.gz"
    shell:
        "zcat {input} | svtools lmerge -i /dev/stdin -f 20 | "
        "bigzip -c > {output}"

rule svtools_genotype:
    # genotype all samples for all variants present in the merged set
    input: 
        vcf="svtools/merged.vcf.gz",
        bam="bams/{sample}.realign.bam",
    output:
        "gt/{sample}.vcf"
    shell:
        
        "zcat {input.vcf} | vawk --header '{{ $6=\".\"; print }}'  "
        "| svtools genotype -B {input.bam} -l {input.bam}.json  "
        "| sed 's/PR...=[0-9\.e,-]*\(;\)\{{0,1\}}\(\\t\)\{{0,1\}}/\2/g' - "
        "> {output}"

rule create_coordinates:
    input: "svtools/merged.vcf.gz",
    output: 
        vcf="svtools/merged.vcf",
        coord="coordinates"
    shell: 
        """
        zcat {input} > {output.vcf}
        create_coordinates -i svtools/merged.vcf -o {output.coord}
        """


rule svtools_copynumber:
    # https://github.com/hall-lab/svtools/blob/master/Tutorial.md 
    # need DNVnator, details see the link
    input:  
        gt="gt/{sample}.vcf",
        coord="coordinates" ,
        cnvnator="/temp/cnvnator-temp/{sample}.bam.hist.root"
    output: 
        "cn/{sample}.vcf"
    params:
        BIN=SPEESEQ
    shell:
        ## annotate variants with copynumber from cnvnator
        ## CNVnator is run as part of speedseq sv
        "svtools copynumber --cnvnator {params.BIN}/cnvnator -s {wildcards.sample} "
        "-w 100 -r {input.cnvnator} \
        "-c {input.coord} -i {input.gt} \
        "> {output}"

rule generate_cnlist:
    input: expand("cn/{sample}.vcf", sample=SAMPLES)
    output: "cn/cn.list"
    shell:
        "ls -1 cn/*vcf > {output}"
        

rule svtools_vcfpaste:
# re-assemble a cohort-level VCF file
    input: 
        vcf="svtools/merged.vcf",
        cnlist="cn/cn.list"
    output:
        "svtools/merged.sv.gt.cn.vcf.gz"
    shell:
        "svtools vcfpaste -m {input.vcf} -f {input.cnlist} -q "
        "| bgzip -c > {output}"

rule svtools_prune:
# filter out additional variants deemed to be identical
    input: "svtools/merged.sv.gt.cn.vcf.gz"
    output: "svtools/merged.sv.pruned.vcf.gz"
    conda: "/home/fangzq/miniconda/envs/sv"
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


#svtools requires a tab-delimited file to specify the number of X chromosome copies
rule chrXnumbers:
    output: "ceph.sex.txt"
    params:
        strain=SAMPLES,
        cs=SAMPLES_cn,
    run:
        with open(output[0],'w') as out:
            for s, c in zip (params.strain, params.cs):
                out.write(f"{s}\t{c}\n")
            
# a file of high-quality, simple deletions and duplications
# https://github.com/hall-lab/svtools/blob/master/Classifier.md
# see the links to create custom training data
rule trainingdata:
    output: "training_vars.bedpe.gz"
    shell:
        "curl -O https://raw.githubusercontent.com/ernfrid/svtools/classifier_docs/resources/training_vars.bedpe.gz"

rule training_variants:
    input: 
        vcf="svtools/merged.sv.pruned.vcf.gz",
        vars="training_vars.bedpe.gz"
    output: "training.vars.vcf.gz"
    shell:
        "zcat {input.vcf} "
        "| svtools vcftobedpe  "
        "| svtools varlookup -a stdin -b {input.vars} -c HQ -d 50 "
        "| svtools bedpetovcf "
        "| svtools vcfsort "
        "| vawk --header '{{ if(I$HQ_AF>0) print $0 }}' "
        "| bgzip -c > {output}"


rule svtools_classify:
    input:
        vcf="svtools/merged.sv.pruned.vcf.gz",
        repeatmasker="repeatMasker.recent.lt200millidiv.LINE_SINE_SVA.b37.sorted.bed.gz",
        chrXnum="ceph.sex.txt",
        train="training.vars.vcf.gz"
    output:
        "output.nb.vcf.gz"
    shell:
        "zcat {input.vcf} |  svtools classify "
        "-g {input.chrXnum} -a {input.repeatmasker} "
        "-m naive_bayes -t {input.train} | bgzip -c > {output}"