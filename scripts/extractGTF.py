import gzip

def extract_gtf(gtf, gene_anno, tx2gene, threads):
    # python code
    lines_seen = set()
    tx2gene_seen = set()
    
    if gtf.endswith("gz"):
        gtf2 = gzip.open(gtf, 'rt') # read text mode
    else:
        gtf2 = open(gtf)
    out1 = open(gene_anno, 'w') 
    out2 = open(tx2gene, 'w')
    out1.write("gene_id\tgene_name\tgene_type\tchrom\n")
    out2.write("tx_id\tgene_id\tgene_name\tchrom\n")
    for line in gtf2:
        if line.startswith('#'): continue
        row = line.strip().split('\t')
        #if row[4] != "gene": continue
        chrom = row[0]
        attrs = row[8].split(';')
        attr = dict(item.strip().split(' ')  for item in attrs if item)
        # ignore txid version
        line_1st = attr['gene_id'].strip('\"') +'\t'+ attr['gene_name'].strip('\"') + '\t'
        if 'gene_type' in attr:
            line_1st  = line_1st + attr['gene_type'].strip('\"')
        line_out = line_1st +'\t' +chrom + '\n'

        if line_out not in lines_seen:
            out1.write(line_out)
            lines_seen.add(line_out)
        if 'transcript_id' in attr:
            tx2gene_out = attr['transcript_id'].strip('\"')+'\t'+ attr['gene_id'].strip('\"')+'\t'+ attr['gene_name'].strip('\"')
            tx2gene_out = tx2gene_out + '\t' + chrom +'\n'
            if tx2gene_out not in tx2gene_seen:
                out2.write(tx2gene_out)
                tx2gene_seen.add(tx2gene_out)
    out1.close()
    out2.close()

extract_gtf(snakemake.input[0], snakemake.output['gene_anno'], snakemake.output['tx2gene'], snakemake.threads)
