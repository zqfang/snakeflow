import os, glob
############################# Required ###################################
# set output directory 
#configfile: "config.yaml"
workdir: config['HBCGM']['WORKSPACE']
HBCGM_BIN = config['HBCGM']['BIN']
# MPD trait ids
if 'TRAIT_IDS' in config:
    TRAIT_IDS = config['TRAIT_IDS'] # input by --config TRAIT_IDS="ids.txt"
else:
    TRAIT_IDS = config['HBCGM']['TRAIT_IDS']

# ghmap input
TRAIT_DATA =  config['HBCGM']['TRAIT_DATA']
GENETIC_REL = config['HBCGM']['GENETIC_REL']

# eblock input
STRAIN_ANNO = config['HBCGM']['STRAIN_ANNO']
SNPDB = config['HBCGM']['SNPS_DIR']
ANNOVAR = config['HBCGM']['ANNOVAR'] 
KNOWNGENE_META = config['HBCGM']['KNOWNGENE_META']
KNOWNGENE = config['HBCGM']['KNOWNGENE']
GENE_EXPRS = config['HBCGM']['GENE_EXPRS']
# open chromatin regions input
ATAC_PEAKS = glob.glob(config['HBCGM']['ATAC_PEAKS'])

############################################################################

## trait ids
with open(TRAIT_IDS, 'r') as t:
    IDS_ = t.read().strip().split()
    assert len(IDS_) >= 1
## skip run if already done
pat = config['HBCGM']['WORKSPACE'] + "MPD_{ids}{sex}/chrX.results.txt"
IDS = []
for mnum in IDS_:
    found = []
    for suf in ["","-f","-m"]:
        found.append(os.path.exists(pat.format(ids = mnum, sex=suf)))
    if not any(found):
        IDS.append(mnum)


CHROMOSOMES = [str(i) for i in range (1, 20)] + ['X'] # NO 'Y'
# output files
# SNPDB = expand("SNPs/chr{i}.txt", i=CHROMOSOMES)
HBCGM =  expand("MPD_{ids}/chr{i}.results.txt", ids=IDS, i=CHROMOSOMES)
HBLOCKS = expand("MPD_{ids}/chr{i}.hblocks.txt", ids=IDS, i=CHROMOSOMES)
HBCGM_NONCODING = expand("MPD_{ids}/chr{i}.open_region.bed", ids=IDS, i=CHROMOSOMES)

# rules that not work in a new node
localrules: target, traits, strain2trait  


rule target:
    input: HBCGM, #HBCGM_NONCODING


# rule pheno:
#     input: TRAIT_DATA
#     ouput: os.path.join(OUTPUT_DIR, "mpd.ids.txt")
#     shell: 
#         "cut -d, -f1 {input} | uniq | sed '1d' > {output.txt}"
rule traits: 
    output: temp(expand("MPD_{ids}/strain.{ids}.temp", ids=IDS))
    run:
        for out in output:
            shell("touch %s"%out)

rule strain2trait:
    input: 
        strain = STRAIN_ANNO,
        ids = "MPD_{ids}/strain.{ids}.temp"
    output: 
        "MPD_{ids}/strain.{ids}.txt",
        "MPD_{ids}/trait.{ids}.txt",
    params:
        trait = TRAIT_DATA,
        outdir = config['HBCGM']['WORKSPACE'],
        traitid = "{ids}",
        rawdata = config['HBCGM']['USE_RAWDATA']
    script:
        "../scripts/strain2traits.py"

rule annotateSNPs:
    input:
        strains = "MPD_{ids}/strain.{ids}.txt",
        snps = os.path.join(SNPDB, "chr{i}.txt"), 
        annodb = os.path.join(ANNOVAR, "chr{i}.AA_by_strains.pkl"),
        kgxref = KNOWNGENE_META, 
        knowngene= KNOWNGENE,
    output:
        hgnc = "MPD_{ids}/chr{i}.genename.txt",
        ensemble = "MPD_{ids}/chr{i}.geneid.txt",
    script:
        "../scripts/annotateSNPs.py"

# find haplotypes
rule eblocks:
    input: 
        snps = os.path.join(SNPDB, "chr{i}.txt"),
        gene_anno = "MPD_{ids}/chr{i}.genename.txt",
        strains = "MPD_{ids}/strain.{ids}.txt",
    output: 
        hb = protected("MPD_{ids}/chr{i}.hblocks.txt"),
        snphb = temp("MPD_{ids}/chr{i}.snp.hblocks.txt")
    params:
        bin = HBCGM_BIN,
    log: "logs/MPD_{ids}.chr{i}.eblocks.log"
    shell:
        "{params.bin}/haplomap eblocks -a {input.snps} -g {input.gene_anno} "
        "-s {input.strains} -p {output.snphb} "
        "-o {output.hb} -v > {log}"

# statistical testing with trait data       
rule ghmap:
    input: 
        hb = "MPD_{ids}/chr{i}.hblocks.txt",
        trait = "MPD_{ids}/trait.{ids}.txt",
        gene_exprs = GENE_EXPRS,
        rel = GENETIC_REL,
        relid = GENETIC_REL+".id"
    output: "MPD_{ids}/chr{i}.results.txt"
    params:
        bin = HBCGM_BIN,
        cat = "MPD_{ids}/trait.{ids}.categorical"
    log: "logs/MPD_{ids}.chr{i}.ghmap.log"
    run:
        categorical = "-c" if os.path.exists(params.cat) else ''
        cats = "catogorical" if os.path.exists(params.cat) else ''
        cmd = "{params.bin}/haplomap ghmap %s "%categorical +\
              "-e {input.gene_exprs} -r {input.rel} " +\
              "-p {input.trait} -b {input.hb} -o {output} " +\
              "-n MPD_{wildcards.ids}_%s -v > {log}"%cats
        shell(cmd)

rule open_chrom:
    input: 
        ghmap="MPD_{ids}/chr{i}.results.txt",
        peaks=ATAC_PEAKS,
    output:
        "MPD_{ids}/chr{i}.open_region.bed",
    params:
        peaks=" ".join(ATAC_PEAKS)
    shell:
        "sed '1,3d' {input.ghmap} | cut -f6-8 | "
        # add 'chr' and convert hblocks to 0-based coordinate #BEGIN {{FS = \"\\t\";OFS = \"\\t\" }};
        "awk -F'\\t' -v OFS='\\t' -v s=1 '{{print \"chr\"$1, $2-s, $3}}' | "
        "sort -k1,1 -k2,2n | "
        "bedtools intersect -sorted -a stdin -b {params.peaks} -wa -wb > {output}"

# about sort -k field1[,field2]
# To sort on the first field and then on the second: sort -k1,1 -k2,2
# The arguments field1 and field2 have the form m.n (m,n > 0) 
# and can be followed by one or more of the modifiers b, d, f, i, n, g, M and r,
# which correspond to the sort options. -n: numberic sort


