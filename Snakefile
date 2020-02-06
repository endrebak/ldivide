
"""
Need to read dataframe in chunks for processing. But how large should they be?
"""

import pandas as pd

from glob import glob

pop = "CEU"

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


sample_info = "data/sample_info.tsv"

for rule in glob("rules/*.smk"):
    include: rule

prefix = "/mnt/work/endrebak/ldivide"

chromosomes = ["chr" + str(i) for i in range(1, 23)]

to_regex = lambda vs: "|".join([str(v) for v in vs])

wildcard_constraints:
    chromosome = to_regex(chromosomes)


rule all:
    input:
        expand("{prefix}/1kg/{chromosome}.vcf.gz", prefix=prefix, chromosome=chromosomes),
        expand("{prefix}/1kg/{population}/{chromosome}.tsv.gz", prefix=prefix, chromosome=chromosomes, population=pop)

