
import os, glob
import pysam

## see the tutorial here: 
## https://github.com/vgteam/vg/wiki/Working-with-a-whole-genome-variation-graph
## https://gtpb.github.io/CPANG18/pages/toy_examples
## https://github.com/vgteam/sv-genotyping-paper

## Commands
workdir: "/data/bases/fangzq/20200815_SV/PacBio_project"
THREADS = 48
GENOME = "/home/fangzq/genome/mouse/GRCm38_68.fa"
VCF_5_STRAIN = "/data/bases/fangzq/20200815_SV/PacBio_project/VCFs/merged_gt_SURVIVOR_1kbpdist.filter.vcf"
VCF_DELETION = "merged_deletion.vcf"
SAMPLES =  [ "BTBR"] #, "BALB","SJL","129S1","AJ", "BTBR"] # ["SJL","129S1","BALB"] #"C57BL6J",
#CHROMSOME = [ str(c) for c in range(1,20)] + ["X", "Y", "MT"]
CHROMOSOME = ['1'] + [ str(c) for c in range(10,20)] + [ str(c) for c in range(2,10)]+ ["X"]
VG = expand("vg/merged.deletion.{g}", sample=SAMPLES, g=['vg','xg','gcsa'])
GT = expand("vg/{sample}.genotypes.vcf", sample=SAMPLES)


######################################################################################################3
rule target:
    input: VG, GT

rule vcf_prepared:
    """
    To run vg graph, the vcf ref, alt alleles should be the same to the reference, 
    not the assembly sequence from the orignial VCF output.
    these vcf will be used to construct 
    a graph containing a reference genome and variation to that reference (i.e. a FASTA and a VCF).
    """
    input:
        genome = GENOME,
        vcf = VCF_5_STRAIN,
    output:
        vcf = VCF_DELETION,
    run:
        vcf = input.vcf
        genome = pysam.FastaFile(input.genome)
        outname = output.vcf
        outfile = open(output.vcf,'w')
        result = []
        with open(vcf) as v:
            for rec in v:
                if rec.startswith("#"): 
                    outfile.write(rec)
                    continue
                line = rec.split("\t")
                chrom = line[0]
                start = line[1]
                ref = line[3]
                alt = line[4]
                info = dict(item.strip().split("=")  for item in line[7].strip('\n').split(';') if item.find("=") >= 0)
                end = info['END']
                if chrom in ["Y", "MT"]: continue 
                #if float(info['AF']) < 0.8: continue
                if 50 < int(info['SVLEN']) < 10000: continue
                if info['SVTYPE'] == 'DEL':
                    ref = genome.fetch(region=f"{chrom}:{start}-{end}")
                    alt = genome.fetch(region=f"{chrom}:{start}-{start}") 
                    line[4] = alt
                    line[3] = ref
                    out = "\t".join(line)
                    outfile.write(out)
                    result.append((chrom, start, end, ref, alt))
        outfile.close()


rule vcf_index:
    input: VCF_DELETION
    output: 
        vcf = VCF_DELETION + ".gz",
        vcfi = VCF_DELETION + ".gz.tbi"
    shell:
        """
        bcftools sort -Oz {input} > {output.vcf}
        tabix -p vcf {output.vcf}
       """

rule vg_chrom:
    output:  temp(expand("vgraph.chr{i}.tmp", i=CHROMOSOME))
    run:
        for v in output:
            os.system("touch %s"%v)

##################### variation graph index map ###########################################
rule vg_construct:
    input:
        genome = GENOME,
        tmp = "vgraph.chr{i}.tmp",
        vcf = VCF_DELETION + ".gz",
        vcfi =  VCF_DELETION + ".gz.tbi",
    output:  
        "vg/merged.deletion.chr{i}.vg"
    threads: THREADS
    shell:
        # or construct the graph using toil-vg
        "vg construct --threads {threads} --reference {input.genome} --vcf {input.vcf} "
        "-R {wlidcards.i} -C " # split into chromosomes
        "--alt-paths --handle-sv > {output} "


rule vg_index:
    input: "vg/merged.deletion.chr{i}.vg"
    output: 
        xg = "vg/merged.deletion.chr{i}.xg",
        gcsa = "vg/merged.deletion.chr{i}.gcsa"
    params:
         kmer_size = 16
    threads: THREADS
    run:
        # store the graph in the xg/gcsa index pair
        ## WARNING: the ./tmp need 2T disk space 
        shell("mkdir -p ./tmp")
        shell("vg index --temp-dir ./tmp --threads {threads} "
              "-x {output.xg} -g {output.gcsa} -k {params.kmer_size} {input}")

rule split_bams:
    input: "../bams/{sample}.realign.bam"
    output: "tmp_bams/{sample}.chr{i}.bam"
    shell:
        "samtools view -h -O BAM {input} {wildcards.i} > {output}"


rule vg_map:
    input: 
        bam = "tmp_bams/{sample}.chr{i}.bam",
        xg = "vg/merged.deletion.chr{i}.xg",
        gcsa = "vg/merged.deletion.chr{i}.gcsa",
    output: "vg/{sample}.chr{i}.mapped.gam"
    threads: 24
    shell:
        ### Note: use fastq.gz if you want => -f 
        # -i, --interleaved  fastq or GAM is interleaved paired-ended
        "vg map -i -b {input.bam} -x {input.xg} -g {input.gcsa} -t {threads}  > {output}.mapped.gam" # -S 0 -u 1 -m 1
        # note -m align mode short, or long


rule vg_sort:
    input: "vg/{sample}.chr{i}.mapped.gam"
    output:
        gai = "vg/{sample}.chr{i}.sorted.gam.gai"
        gam = "vg/{sample}.chr{i}.sorted.gam"
    threads: 24
    shell:
        "vg gamsort -t {threads} -i BTBR_chr1_sort.gam.gai BTBR_chr1_default.gam > BTBR_chr1_sorted.gam"


### NOTE: to view the gam or graph,
# vg chunk required a gam.gai file, so vg gamsort first (very slow)
# vg chunk -x graph.xg -a reads.sorted.gam -g -T -b ./subgraph -c 20 -p ${CHROM}:${START}-${END} > chunk.vg
# vg view -dp ${SOMETHING}.vg -A ${SOMETHING}.gam | dot -Tsvg -o chunk.svg

################ SV genotyping ########################################################
rule vg_pack:
    input: 
        vg = "vg/merged.deletion.chr{i}.vg",
        gam = "vg/{sample}.chr{i}.mapped.gam",
        gai= "vg/{sample}.chr{i}.mapped.gam.gai"
    output: "vg/{sample}.chr{i}.pack"
    threads: THREADS
    shell:
        "vg pack --threads {threads} -x {input.vg} -g {input.gam} -o {output} -Q 5"

rule vg_snarls:
    input: "vg/merged.deletion.chr{i}.vg"
    output: "vg/merged.deletion.chr{i}.snarls"
    threads: THREADS
    shell:
        "vg snarls --threads {threads} {input} > {output}"


rule split_vcf:
    input:
        tmp = "vgraph.chr{i}.tmp",
        vcf = VCF_DELETION + ".gz",
        vcfi = VCF_DELETION + ".gz.tbi"
    output:
        vcf = VCF_DELETION + ".chr{i}.vcf.gz",
        vcfi = VCF_DELETION + ".chr{i}.vcf.gz.tbi"
    run:
        shell("bcftools view {input.vcf} {whildcards.i} | bcftools sort -Oz > {output.vcf}")
        shell("tabix -f -p vcf {output.vcf}")


rule vg_call:
    input:
        vg = "vg/merged.deletion.chr{i}.vg",
        snarls = "vg/merged.deletion.chr{i}.snarls",
        pack = "vg/{sample}.chr{i}.pack",
        vcf = VCF_DELETION + ".chr{i}.vcf.gz",
        vcfi = VCF_DELETION + ".chr{i}.vcf.gz.tbi"
    output:
        gt = "vg/{sample}.chr{i}.genotypes.vcf"
    threads: THREADS
    shell:
        # If an insertion fasta file was needed to construct with -I, it must be passed to call with -i.
        "vg call --threads {threads} --pack {input.pack} --snarls {input.snarls} "
        "--vcf {input.vcf} --sample {wildcards.sample} {input.vg} > {output.gt} " 
        # --vcf must have been used to construct input graph with -a)