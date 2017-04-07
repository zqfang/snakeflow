def do_something(data_path, out_path, threads, myparam):
        # python code

do_something(snakemake.input[0], snakemake.output[0], snakemake.threads, snakemake.config["myparam"])
