def barplot(enrichr_res, png, pdf, threads):
    # python code
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    from matplotlib.figure import Figure
    from numpy import log10
    from pandas import read_table

    d = read_table(enrichr_res)
    d['logAP'] = -log10(d['Adjusted P-value']) 
    d = d.sort_values('logAP', ascending=False)
    dd = d.head(10).sort_values('logAP')
    fig = Figure(figsize=(12,6))
    canvas = FigureCanvas(fig)
    ax = fig.add_subplot(111)
    bar = dd.plot.barh(x='Term', y='logAP', color="salmon", alpha=0.75, edgecolor='none',fontsize=32, ax=ax)
    bar.set_xlabel("-log$_{10}$ Adjust P-value", fontsize=32)
    bar.set_ylabel("")
    #bar.set_title("",fontsize=32)
    bar.legend(loc=4)
    fig.savefig(png, bbox_inches='tight')
    fig.savefig(pdf, bbox_inches='tight')

barplot(snakemake.input[0], snakemake.output['png'],snakemake.output['pdf'], snakemake.threads)
