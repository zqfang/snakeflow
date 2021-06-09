
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

# rule vg_chrom:
#     output:  temp(expand("vgraph.chr{i}.tmp", i=CHROMOSOME))
#     run:
#         for v in output:
#             os.system("touch %s"%v)

##################### variation graph index map ###########################################
rule vg_construct:
    input:
        genome = GENOME,
        # tmp = "vgraph.chr{i}.tmp",
        vcf = VCF_DELETION + ".gz",
        vcfi =  VCF_DELETION + ".gz.tbi",
    output:  
        #"vg/chr{i}.vg"
        "vg/merged.deletion.vg"
    threads: THREADS
    shell:
        # or construct the graph using toil-vg
        "vg construct --threads {threads} --reference {input.genome} --vcf {input.vcf} "
        # " -R {wlidcards.i} -C "
        "--node-max 32 --alt-paths --flat-alts --handle-sv "
        "| vg mod --until-normal 10 - "
        "| vg mod --chop 32 - "
        "| vg ids --sort - > {output} "


rule vg_prune:
    input:  "vg/merged.deletion.vg",
    output:  "vg/merged.deletion.pruned.vg"
    threads: THREADS
    shell:
        "vg prune --threads {threads} -r {input} > {output}"

rule vg_index:
    input: "vg/merged.deletion.pruned.vg"
    output: 
        xg = "vg/merged.deletion.xg",
        gcsa = "vg/merged.deletion.gcsa"
    params:
         kmer_size = 16
    threads: THREADS
    run:
        # store the graph in the xg/gcsa index pair
        ## WARNING: the ./tmp need 2T space 
        shell("mkdir -p ./tmp")
        shell("vg index --temp-dir ./tmp --threads {threads} --xg-alts "
              "-x {output.xg} -g {output.gcsa} -k {params.kmer_size} {input}")


rule samtools_namesort:
    input: 
        bam = "../bams/{sample}.realign.bam",
        genome = GENOME,
    output: "{sample}.rnst.cram"
    threads: THREADS
    shell:
        # if got "Too many open files" error, please increase -m 
        "samtools sort -m 4G -@ {threads} -n {input.bam} -T tmp.{wildcards.sample} "
        "| samtools view -C --reference {input.genome} -o {output} -"

rule vg_map:
    input: 
        cram="{sample}.rnst.cram",
        xg = "vg/merged.deletion.xg",
        gcsa = "vg/merged.deletion.gcsa",
    output: "vg/{sample}.mapped.gam"
    threads: 24
    shell:
        ### Note: use fastq.gz if you want => -f 
        # # convert the alignments to interleaved fastq and map
        "samtools bam2fq {input.cram} | "
        "vg map -if - -x {input.xg} -g {input.gcsa} -t {threads}  > {output}.mapped.gam" # -S 0 -u 1 -m 1
        # note -m align mode short, or long

# rule count_mapped_reads:
#     input:
#         "vg/{sample}.mapped.gam"
#     output:
#         "vg/stats/{sample}.tsv"
#     run:
#         shell("vg view -a {input} | jq -rc '[.name, if .mapping_quality == null then 0 else .mapping_quality end ] | @tsv' | awk '$2>0' | wc -l > {output}")
#         shell("vg view -a {input} | jq -rc '[.name, if .mapping_quality == null then 0 else .mapping_quality end ] | @tsv' | awk '$2>=10' | wc -l >> {output}")
#         shell("vg view -a {input} | jq -rc '[.name, if .mapping_quality == null then 0 else .mapping_quality end ] | @tsv' | awk '$2>=20' | wc -l >> {output}")
#         shell("vg view -a {input} | jq -rc '[.name, if .mapping_quality == null then 0 else .mapping_quality end ] | @tsv' | awk '$2>=30' | wc -l >> {output}")
#         shell("vg view -a {input} | jq -rc '[.name, if .mapping_quality == null then 0 else .mapping_quality end ] | @tsv' | awk '$2>=40' | wc -l >> {output}")
#         shell("vg view -a {input} | jq -rc '[.name, if .mapping_quality == null then 0 else .mapping_quality end ] | @tsv' | awk '$2>=50' | wc -l >> {output}")
#         shell("vg view -a {input} | jq -rc '[.name, if .mapping_quality == null then 0 else .mapping_quality end ] | @tsv' | awk '$2>=60' | wc -l >> {output}")
#         shell("vg view -a {input} | jq -rc '[.name, if .identity == null then 0 else .identity end ] | @tsv' | awk '$2>=1' | wc -l >> {output}")
#         shell("vg view -a {input} | jq -rc '[.name, if .identity == null then 0 else .identity end ] | @tsv' | awk '$2>=0.9' | wc -l >> {output}")
#         shell("vg view -a {input} | jq -rc '[.name, if .identity == null then 0 else .identity end ] | @tsv' | awk '$2>=0.5' | wc -l >> {output}")
#         shell("vg view -a {input} | wc -l >> {output}")

# rule cat_counts:
#     input:
#         stats=expand("mappings/stats/{sample}.tsv", sample=SAMPLES),
#         samples="../../illumina_reads/samples.txt"
#     output:
#         "mappings/stats/all.tsv"
#     run:
#         for f in input.stats:
#             sample = f.split("/")[2].split(".")[0]
#             shell("cat {f} | tr '\\n' '\\t' >> {output}")
#             shell("echo -e -n \"construct\\t\" >> {output}")
#             shell("grep {sample} {input.samples} | awk '{{ print $3}}' >> {output}")

rule sort_gam:
    input:
        "vg/{sample}.mapped.gam"
    output:
        gam="vg/{sample}.mapped.sorted.gam",
        gai="vg/{sample}.mapped.sorted.gam.gai"
    threads: THREADS
    shell:
        "vg gamsort -t {threads} -i {output.gai} {input} > {output.gam}"


################ SV genotyping ########################################################
rule vg_pack:
    input: 
        vg = "vg/merged.deletion.pruned.vg",
        gam = "vg/{sample}.mapped.sorted.gam",
        gai= "vg/{sample}.mapped.sorted.gam.gai"
    output: "vg/{sample}.pack"
    threads: THREADS
    shell:
        "vg pack --threads {threads} -x {input.vg} -g {input.gam} -o {output} -Q 5"

rule vg_snarls:
    input: "vg/merged.deletion.pruned.vg"
    output: "vg/merged.deletion.pruned.snarls"
    threads: THREADS
    shell:
        "vg snarls --threads {threads} {input} > {output}"

rule vg_call:
    input:
        vg = "vg/merged.deletion.pruned.vg",
        snarls = "vg/merged.deletion.pruned.snarls",
        pack = "vg/{sample}.pack",
        vcf = "{sample}.deletion.vcf.gz",
        vcfi = "{sample}.deletion.vcf.gz.tbi",
    output:
        gt = "vg/{sample}.genotypes.vcf"
    threads: THREADS
    shell:
        # If an insertion fasta file was needed to construct with -I, it must be passed to call with -i.
        "vg call --threads {threads} --pack {input.pack} --snarls {input.snarls} "
        "--vcf {input.vcf} --sample {wildcards.sample} > {output.gt} " # --vcf must have been used to construct input graph with -a)