

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

rule find_scores:
    input:
        "{prefix}/1kg/{vector_type}/{population}/{window_size}/{chromosome}.pq"
    output:
        "{prefix}/1kg/vector_summary/{vector_type}/{population}/{window_size}/{chromosome}.txt"
    run:
        s = pd.read_parquet(input[0]).squeeze()

        w = wildcards

        median = s.median()
        mean = s.mean()
        _sum = s.sum()

        measure = {"normalized_square": "window", "vector": "diagonal"}[w.vector_type]

        rowdict = {"WindowSize": w.window_size, "Population": w.population, "Chromosome": w.chromosome, "Measure": measure, "Median": median, "Mean": mean, "Sum": _sum}

        rowdicts = [rowdict]

        df = pd.DataFrame.from_dict(rowdicts)

        df.to_csv(output[0], sep="\t", index=False)
        # with open(output[0], "w+") as o:
        #     o.write(str(median) + "\n")

score_f = "{prefix}/1kg/vector_summary/{vector_type}/{population}/{window_size}/{chromosome}.txt"
rule collect_scores:
    input:
        expand(score_f, prefix=prefix, chromosome=chromosomes, population=pop, window_size=window_sizes, vector_type="vector normalized_square".split())
    output:
        "{prefix}/1kg/vector_summary/summary.txt"
    run:
        df = pd.concat([pd.read_table(f) for f in input])
        df.to_csv(output[0], sep="\t", index=False)





