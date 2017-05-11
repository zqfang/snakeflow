def merge_stringtie(annotation, filelist, full, tpm, fpkm, threads):
    # python code
    from pandas import read_table, read_excel, concat, ExcelWriter

    anno = read_table(annotation)
    frames = []
    for f in filelist:
        name = f.split("/")[-1].strip(".tab")
        frame = read_table(f, index_col="Gene ID")
        frame = frame.sort_index()
        frame.rename(columns={'Coverage': 'Coverage.'+name,
                              'FPKM': 'FPKM.'+name, 'TPM':'TPM.'+name},
                     inplace=True)
        frames.append(frame)

    #when concat dataframe, differrent dataframes will ordered by their index by default
    result = concat(frames, axis=1)
    df = result.loc[:,~result.columns.duplicated()]
    df_merge = anno.merge(df, left_on='gene_id',right_index=True, how='right')
    df_merge.drop('Gene Name', axis=1, inplace=True)
    df_merge.to_csv(full, index=False)
    col_tpm = ['gene_id','gene_name'] + [item for item in df_merge.columns if item.startswith("TPM.")]
    col_fpkm = ['gene_id','gene_name'] + [item for item in df_merge.columns if item.startswith("FPKM.")]
    df_merge[col_tpm].to_csv(tpm, index=False)
    df_merge[col_fpkm].to_csv(fpkm, index=False)

merge_stringtie(snakemake.input['annotation'], snakemake.input['filelist'], 
                snakemake.output['full'], snakemake.output['tpm'],  snakemake.output['fpkm'],  
                snakemake.threads)
