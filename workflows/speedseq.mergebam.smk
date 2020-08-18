import os, glob
############### Globals ########################

#configfile: "config.yaml"
workdir: config['SV']['WORKSPACE']

GENOME = config['GENOME']
SAMPLES = ['129P2', '129S1', '129S5', 'AKR', 'A_J', 'B10', 
        'BPL', 'BPN', 'BTBR', 'BUB', 'B_C', 'C3H', 'C57BL10J',
        'C57BL6NJ', 'C57BRcd', 'C57LJ', 'C58', 'CBA', 'CEJ', 
        'DBA', 'DBA1J', 'FVB', 'ILNJ', 'KK', 'LGJ', 'LPJ', 
        'MAMy', 'MRL','NOD', 'NON', 'NOR', 'NUJ', 'NZB', 'NZO', 'NZW', 
        'PJ', 'PLJ', 'RFJ', 'RHJ', 'RIIIS', 'SEA', 'SJL', 'SMJ', 'ST', 'SWR', 'TALLYHO', 'RBF'] + \
         ['CAST', 'MOLF', 'PWD','PWK', 'SPRET', 'WSB']  # <- wild derived except MRL
CHROMSOME = ['1'] + [ str(c) for c in range(10,20)] + [ str(c) for c in range(2,10)]+ ["MT", "X", "Y"]
# OUTPUT

# output files
SVDB = expand("SVs/chr{i}.txt", i=CHROMSOME)

################################ rules ###################################

rule target:
    input: SVDB

rule align:
    input:
        r1="{sample}/{ID}_1.fastq.gz",
        r2="{sample}/{ID}_2.fastq.gz",
        genome=GENOME, # bwa indexed
    output: 
        bam="bams/{sample}/{ID}.bam",
        discordants="bams/{sample}/{ID}.discordants.bam",
        splitters="bams/{sample}/{ID}.splitters.bam",
        IDList = "bams/{sample}/{ID}.IDList.txt"
    params:
        outprefix="bams/{sample}",
        # PL has to be one of ILLUMINA，SLX，SOLEXA，SOLID，454，LS454，COMPLETE，
        # PACBIO，IONTORRENT，CAPILLARY，HELICOS, or UNKNOWN
        # ID：输入reads集的ID号； LB： reads集的文库名； SM：样本名称； PL：测序平台
        RG="@RG\tID:{sample}.{ID}\tPL:illumina\tSM:{sample}\tLB:{ID}"
    run:
        shell("speedseq align -o {params.outprefix} -R {params.RG} {genome} {input.r1} {input.r2}")
        shell("echo -e {ID} > {output.IDList}")

rule collectIDList:
    output:
        "bams/{sample}/allIDLists.txt"
    shell:
        "cat bams/{wildcards.sample}/*.IDList.txt > {output}"

# ## merge samples
rule mergeBam:
    input: "bams/{sample}/allIDLists.txt" # expand("{{sample}}/{ID}.bam", ID=SAMPLE_LB_DICT[wilcard.sample]])
    output: 
        b = "bams/{sample}.merged.bam",
        d = "bams/{sample}.merged.discordants.bam",
        s = "bams/{sample}.merged.splitters.bam",
    run:
        with open(input[0], 'r') as libs:
            ids = libs.read().strip().split()
        bams = " ".join(["bams/{sample}/%s.bam"%i for i in ids])
        bamds = " ".join(["bams/{sample}/%s.discordants.bam"%i for i in ids])
        bamss = " ".join(["bams/{sample}/%s.splitters.bam"%i for i in ids])
        shell("sambamba merge {output.b} %s"%bams)
        shell("sambamba index {output.b}")
        shell("sambamba merge {output.d} %s"%bamds)
        shell("sambamba index {output.d}")
        shell("sambamba merge {output.s} %s"%bamss)
        shell("sambamba index {output.s}")






