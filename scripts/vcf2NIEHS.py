import time, sys, os, glob

def vcf2niehs(invcf, outdir, sstrain, chromosome, qual_samtools=50, heterzygote_cutoff = 20):
    
    ## strains
    with open(sstrain, 'r') as ss:
        strains = ss.read().strip().split()
    os.makedirs(outdir, exist_ok=True)

    ##### filtering parameters
    qualCutoffForSamtools = qual_samtools
    cutoffForPLscoreForNonHomozygousAlt = heterzygote_cutoff

    numOfAlt, numGoodAlt, theGoodAlt = 0, 0, 0
    totalVariant, totalInterested, numNonPass, numINDEL, numLowQual, numMultiAlt, numGoodSNP = 0, 0, 0, 0, 0, 0, 0
    nucleotide = {"A": 1, "C":1, "G":1, "T":1} 
    
    # prepared output
    ## chromosome = list(range(10,20)) + list(range(1, 10)) + [ "X", "Y", "MT"]
    chromosome = [chromosome]
    outputdict = { str(k) : open(os.path.join(outdir, f"chr{k}.txt" ), 'w') for k in chromosome }
    # add header
    for chrom, output in outputdict.items():
        output.write("LOCAL_IDENTIFIER\tSS_ID\tCHROMOSOME\tACCESSION_NUM\tPOSITION\tSTRAND\tALLELES\t")
        output.write("\t".join(strains)+ "\n")

    # parse VCF
    # CHROM	POS	ID	REF	ALT	QUAL FILTER	INFO FORMAT	AKR
    lineCount = 0
    for line in open(invcf, 'r'):
        # skip description lines
        if line.startswith("#"): continue
        newline = line.strip().split("\t")

        lineCount +=1
        totalVariant +=1

        if lineCount % 100000 == 0:
            sys.stdout.write(f"line: {lineCount}\n")
            
        ## check vcf format
        if len(newline) < 9: 
            sys.exit("Error! format not recognized")
        # skip chr if not interested 
        if newline[0] not in chromosome: continue
        totalInterested +=1
        
        # FIXME: gatk not pass string is ?
        # gatk -> [PASS, LowQual, '.' ], bcftools -> ['.']
        # gatk -> not filtering has been done when it's '.'
        if newline[6] not in ["PASS", "."]:
            numNonPass +=1
            continue
        ### first, get the ref/alternative alleles
        ref = newline[3]
        # skip indels, bcftool col7 is/not INDEL
        # FIXME: len(ref) > 1 is bad ref or indel ?
        if len(ref) > 1 or newline[7].startswith("INDEL"): 
            numINDEL +=1
            continue
        
        # QUAL: The Phred-scaled probability that a REF/ALT polymorphism exists.
        # QUAL: -10 * log (1-p)
        # it's not a very useful property for evaluating the quality of a variant call 
        # FIXME: the cutoff for GATK ???
        if float(newline[5]) < qualCutoffForSamtools:
            numLowQual +=1
            continue

        ###find the entries for GT & PL
        # GATK
        if newline[8] == 'GT:AD:DP:GQ:PL':
            GTind, PLind = 0, 4
        elif newline[8] == 'GT:AD:DP:GQ:PGT:PID:PL':
            GTind, PLind = 0, 6
        # samtools
        elif newline[8] == "GT:PL:DP:DV:SP:DP4:DPR:GP:GQ":
            GTind, PLind = 0, 1
        else:
            IDS = newline[8].split(":") 
            GTind, PLind = -1, -1
            for i in IDS:
                if i == 'GT': GTind = i
                if i == "PL": PLind = i

            if (GTind == -1) | (PLind == -1):
                sys.exit("Error! NOT PL or GT!") 
        
        chrom = newline[0]
        pos = newline[1]  
        alts = newline[4].split(",")
        numOfAlt = len(alts)
        # FIXME: 
        # if numOfAlt > 1: sys.exit("Bad Alt")

        hasAlt= [0] * len(alts)
        alleles = [-1] * len(strains)
        # CHROM	POS	ID	REF	ALT	QUAL FILTER	INFO FORMAT	AKR
        for s, strain in enumerate(strains):
            ## identify the allele for each strai
            strainFormats = newline[9+s].split(":")

            # if NOT found: GTind = PLind = -1
            if (not strainFormats[GTind]) or (not strainFormats[PLind]):
                alleles[s] = "NN"
                continue
            if strainFormats[GTind] == "./.":
                alleles[s] = "NN"
                continue

            GTs = strainFormats[GTind].split("/")
            if len(GTs) != 2:
                sys.exit(f"Improper GT: {strainFormats[GTind]} in {line}\n")
            if GTs[0] != GTs[1]:
                alleles[s] = "NN" # heterozyte -> NN
                continue
            ## now this is supposed to be a homogyzous call. Check whether the call is of good quality
            
            # TODO: what index do here
            GTs = [int(s) for s in GTs]
            index = (GTs[0]+ 1) * (GTs[0] + 2) // 2 - 1 # 0 or 2
            PLs = strainFormats[PLind].split(",")
            # PLs -> ['0/0', '0/1', '1/1']
            PLs = [float(s) for s in PLs]
            if len(PLs) != ((numOfAlt + 1) * (numOfAlt + 2) // 2):
                sys.exit(f"Improper PL found: {strainFormats[PLind]}\n{line}")
            
            minScore = 10000000000
            for i, pl in enumerate(PLs):
                if i == index: continue # why ?
                if (PLs[i] - PLs[index]) < minScore:
                    minScore = PLs[i] - PLs[index] 
            
            if minScore >= cutoffForPLscoreForNonHomozygousAlt:
                ## this is a good call
                if GTs[0] > 0:
                    # select alt
                    alleles[s] = f"{alts[GTs[0]-1]}{alts[GTs[0]-1]}"
                    hasAlt[GTs[0]-1] = 1
                else:
                    alleles[s] = f"{ref}{ref}"
            else:
                alleles[s] = "NN"

        numGoodAlt = 0
        theGoodAlt = -1
        # TODO: fixme
        for i, k in enumerate(hasAlt):
            if k == 1:
                numGoodAlt +=1
                theGoodAlt = i
        if numGoodAlt == 1:
            if theGoodAlt == -1:
                sys.exit("Interal error: the good Alt not found properly")
            if alts[theGoodAlt] not in nucleotide.keys():
                numINDEL +=1

            ## this is a good SNP across the strains. save it for output 
            allAlleles = "\t".join(alleles)
            alt = alts[theGoodAlt]
            # write output here
            ## header 
            # ## LOCAL_IDENTIFIER SS_ID CHROMOSOME ACCESSION_NUM POSITION STRAND ALLELES + strains ...
            D = f"SNP_{chrom}_{pos}"
            outline = f"{D}\t{D}\t{chrom}\t{D}\t{pos}\t.\t{ref}/{alt}\t{allAlleles}\n"
            outputdict[chrom].write(outline)
            numGoodSNP += 1
        else:
            if numGoodAlt == 0:
                numLowQual += 1
            else:
            ## multiple good alternatives found;
                numMultiAlt +=1          

    # close file
    for chrom, ouput in outputdict.items(): output.close()

    print(f"total: {totalVariant}, interested: {totalInterested}, notPass: {numNonPass}, " +\
        f"indels: {numINDEL}, lowQual: {numLowQual}, multiAlt: {numMultiAlt}, good: {numGoodSNP}")


vcf2niehs(snakemake.input['vcf'], snakemake.params['outdir'], 
          snakemake.input['strain'], snakemake.params['chrom'],
          snakemake.params['qual_samtools'], 
          snakemake.params['heterzygote_cutoff'])


### strains included in the VCF file in the exact order
# strains = ("B10", "BTBR", "BUB", "CEJ", "DBA1J", "FVB", "KK", "NON", "NUJ", "NZB", "NZW", "RFJ", "RHJ", "RIIIS", "SJL", 
#   	       "129P2", "129S1", "129S5", "A_J", "AKR", "B_C", "C3H", "C57BL10J", "C57BL6NJ", "C57BRcd", "C57LJ", "C58", 
#   	       "CBA", "DBA", "ILNJ", "LGJ", "LPJ", "MAMy", "MRL", "NOD", "NZO", "PJ", "PLJ", "SEA", "SMJ", "ST", "SWR")

# # inputs and outputs
# chromosome = "X"
# invcf = "/home/fangzq/data/20180101/samtools1/var.chrX.raw.vcf"
# sstrain = "/home/fangzq/data/20180101/samtools1/var.chrX.raw.strain.txt"
# outdir = "/home/fangzq/data/test_convert"

# vcf2niehs(invcf, outdir, sstrain, chromosome, 50,  20)