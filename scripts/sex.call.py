import os, gzip, sys 
import pandas as pd 

def GenotypeCountsOfChrX(vcf, outfile):
    if vcf.endswith("gz"):
        _vcf = gzip.open(vcf, 'rt') 
    else:
        _vcf = open(vcf, 'r')
    
    out = open(outfile, "w")
    strains_all = []
    genotype_all = []
    stats_dict = {}
    for line in _vcf:
        # skip description lines
        if line.startswith("##"): continue
        # write header 
        if line.startswith("#CHROM"): 
            strain_name = line.split("\tFORMAT\t")[-1]
            strains = strain_name.strip().split("\t")
            strains_all += strains
            out.write(strain_name)
            continue   

        newline = line.strip().split("\t")
        if newline[6] != "PASS":
            continue
        newline2 = newline[9:]
        gt = []
        for g in newline2: gt.append(g.split(":")[0])
        out.write("\t".join(gt)+"\n")
        genotype_all.append(gt)

    df = pd.DataFrame(genotype_all)
    df.value_counts()
    out.close()
    _vcf.close()

GenotypeCountsOfChrX("VCFs/combined.chrX.hardfilter.pass.vcf.gz", "chrX.stat")