variant_url_hg38 = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.{chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
variant_url_hg19 = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.{chromosome}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" 

variants = {"hg38": variant_url_hg38,
            "hg19": variant_url_hg19}

rule fetch_variants:
    output:
        f"{variant_prefix}/1kg/{genome}/{{chromosome}}.vcf.gz"
    resources:
        instances = 1
    run:
        variant_url = variants[genome]
        url = variant_url.format(chromosome=wildcards.chromosome)
        shell("axel {url} -q -o {output[0]}")



rule index_variants:
    input:
        rules.fetch_variants.output[0]
    output:
        f"{variant_prefix}/1kg/{genome}/{{chromosome}}.vcf.gz.tbi"
    shell:
        "tabix {input[0]}"


rule vcf_to_bed:
    input:
        rules.fetch_variants.output[0]
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
