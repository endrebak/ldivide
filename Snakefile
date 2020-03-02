"""
Need to read dataframe in chunks for processing. But how large should they be?
"""

import pandas as pd

from glob import glob

pop = ["CEU", "YRI"]
# pop = "CEU"
genome = "hg19"

def pop_as_list():
    # pop = config["population"]
    if isinstance(pop, str):
        return [pop]
    elif isinstance(pop, list):
        return pop
    else:
        values = []
        for v in list(pop.values())[0]:
            values.append(v)

        return values

populations = pop_as_list()

window_sizes = [2500, 5000, 10000, 15000, 20000]

variant_prefix = "/mnt/work/endrebak/1kg"
sample_info = "data/sample_info.tsv"
chromosomes = ["chr" + str(i) for i in range(1, 23)]

prefix = "/mnt/work/endrebak/ldivide/hg19"

for rule in glob("rules/*.smk"):
    include: rule


to_regex = lambda vs: "|".join([str(v) for v in vs])

wildcard_constraints:
    chromosome = to_regex(chromosomes),
    population = to_regex(populations),


f = "{prefix}/1kg/local_minima/{population}/{chromosome}.pq"
col = "{prefix}/1kg/colvectors/{population}/{chromosome}.pq"
row = "{prefix}/1kg/rowvectors/{population}/{chromosome}.pq"
f = "{prefix}/1kg/normalized_square/{population}/{chromosome}.pq"
f = "{prefix}/1kg/local_minima/{population}/{chromosome}.tsv"
f = "{prefix}/1kg/local_minima/{population}/genome.tsv"
# f = "{prefix}/1kg/local_minima/{population}/genome.tsv"

# f = "{prefix}/1kg/minima/{population}/{window_size}/{chromosome}.tsv"


f = f"{prefix}/1kg/vector_summary/summary.txt"

# f = "{prefix}/1kg/vector_summary/{vector_type}/{population}/{window_size}/{chromosome}.txt"

rule all:
    input:
        f
        # expand("{prefix}/1kg/{chromosome}.vcf.gz", prefix=prefix, chromosome=chromosomes),
        # expand("{prefix}/1kg/{population}/{chromosome}.tsv.gz", prefix=prefix, chromosome=chromosomes, population=pop),
        # expand("{prefix}/1kg/autocorr/{population}/{chromosome}.tsv.gz", prefix=prefix, chromosome=chromosomes, population=pop)
        # expand(f, prefix=prefix, chromosome=chromosomes, population=pop, window_size=window_sizes, vector_type="vector normalized_square".split()),
        # expand(row, prefix=prefix, chromosome=chromosomes, population=pop),
        # expand(col, prefix=prefix, chromosome=chromosomes, population=pop)
        # expand("{prefix}/1kg/vector/{population}/{chromosome}.pq", prefix=prefix, chromosome=chromosomes, population=pop)

