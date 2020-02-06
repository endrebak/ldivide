

rule variant_samples:
    input:
        rules.fetch_variants.output[0],
        rules.index_variants.output[0]
    output:
        "{prefix}/1kg/{chromosome}_samples.tsv"
    resources:
        instances = 1
    shell:
        "bcftools view -h {input[0]} | tail -1 | tr '\t' '\n' | tail -n +10 > {output[0]}"


rule individuals_in_reference_panel:
    input:
        samples = sample_info,
        in_vcf = rules.variant_samples.output[0]
    output:
        n_ind = "{prefix}/partitions/{population}/{chromosome}/n_individuals.txt",
        samples = "{prefix}/partitions/{population}/{chromosome}/samples.txt"
    run:
        df = pd.read_table(input["samples"])
        in_vcf = pd.read_table(input["in_vcf"], header=None, squeeze=True).to_list()
        print(in_vcf[:5])

        # print(pop_as_list())

        pop = df[(df.Population.isin(populations)) & (df.Sample.isin(in_vcf))]

        pop_size = len(pop.drop_duplicates("Sample"))

        with open(output["n_ind"], "w+") as o:
            o.write(str(pop_size) + "\n")

        pop.Sample.to_frame().to_csv(output["samples"], index=False, sep="\t", header=False)


rule subset_on_population:
    input:
        variants = rules.fetch_variants.output[0],
        samples = rules.individuals_in_reference_panel.output.samples,
        index = rules.index_variants.output[0]
    output:
        "{prefix}/1kg/{population}/{chromosome}.tsv.gz"
    shell:
        'bcftools view --threads 48 --force-samples -S {input.samples} {input.variants} -H | cut -f 10- | tr "|" "\\t" | gzip > {output[0]}'
