

rule add_position_to_diag_vectors:
    input:
        haplos = rules.subset_on_population.output[0],
        vector = rules.vector.output[0]
    output:
        "{prefix}/1kg/vector/{population}/{window_size}/{chromosome}_with_pos.pq",
    run:
        positions = pd.read_parquet(input.haplos, usecols=0).reset_index().squeeze()
        v = pd.read_parquet(input.vector)
        v.index = positions
        v.to_parquet(output[0])


rule add_position_to_other_vectors:
    input:
        haplos = rules.subset_on_population.output[0],
        vector = rules.normalize_square.output[0]
    output:
        "{prefix}/1kg/vector/{population}/{window_size}/{chromosome}_with_pos.pq",
    run:
        positions = pd.read_parquet(input.haplos, usecols=0).reset_index().squeeze()
        v = pd.read_parquet(input.vector)
        v.index = positions
        v.to_parquet(output[0])


rule split_chromosomes:
    input:
        "data/fourier_ls_{pop}.bed"
    output:
        expand("data/fourier_ls_{{pop}}_{chromosome}.bed", chromosome=chromosomes)
    run:
        outpath_template = "data/fourier_ls_{pop}_{chromosome}.bed"

        gr = pr.read_bed(input[0])

        for c, cdf in gr:
            outpath = outpath_template.format(pop=wildcards.pop, chromosome=c)
            cdf.to_csv(outpath, sep="\t", index=False)

# rule graph_chromosome
