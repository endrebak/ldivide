variant_url = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.{chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

rule fetch_variants:
    output:
        "{prefix}/1kg/{chromosome}.vcf.gz"
    resources:
        instances = 1
    run:
        url = variant_url.format(chromosome=wildcards.chromosome)
        shell("axel {url} -q -o {output[0]}")



rule index_variants:
    input:
        rules.fetch_variants.output[0]
    output:
        "{prefix}/1kg/{chromosome}.vcf.gz.tbi"
    shell:
        "tabix {input[0]}"


rule vcf_to_bed:
    input:
        rules.fetch_variants.output
    output:
        "{prefix}/1kg/{chromosome}.bed"
    run:
        import gzip

        rows_to_skip = 0

        with gzip.open(input[0]) as fh:
            for i, l in enumerate(fh, 0):
                rows_to_skip = i
                if not l.decode().startswith("#"):
                    break

        print("rows_to_skip", rows_to_skip)

        df = pd.read_table(input[0], header=None, skiprows=rows_to_skip, usecols=[0, 1], nrows=None)
        print(df.head())
        df.columns = "Chromosome Start".split()
        df.insert(df.shape[1], "End", df.Start + 1)
        df.loc[:, "Chromosome"] = "chr" + df.Chromosome.astype(str)

        df.to_csv(output[0], sep="\t", index=False, header=False)
