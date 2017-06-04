def extract_gtf(gtf, gene_anno, tx2gene, threads):
    # python code
    lines_seen = set()
    tx2gene_seen = set()

    out1 = open(gene_anno, 'w') 
    out2 = open(tx2gene, 'w')
    out1.write("gene_id\tgene_name\tgene_type\n")
    out2.write("tx_id\tgene_id\n")
    for line in open(gtf):
        if line.startswith('#'): continue
        attr = dict(item.strip().split(' ')  for item in line.split('\t')[8].strip('\n').split(';') if item)
        line_1st = attr['gene_id'].strip('\"')+'\t'+ attr['gene_name'].strip('\"')+'\t'
        line_2nd = attr['gene_type'].strip('\"')+'\n'
        line_out = line_1st + line_2nd

        if line_out not in lines_seen:
            out1.write(line_out)
            lines_seen.add(line_out)
        if 'transcript_id' in attr:
            tx2gene_out = attr['transcript_id'].strip('\"')+'\t'+ attr['gene_id'].strip('\"')+'\n'
            if tx2gene_out not in tx2gene_seen:
                out2.write(tx2gene_out)
                tx2gene_seen.add(tx2gene_out)
    out1.close()
    out2.close()

extract_gtf(snakemake.input[0], snakemake.output['gene_anno'], snakemake.output['tx2gene'], snakemake.threads)
