# ldivide

ldivide efficiently finds independent LD-blocks from genetic data, using only
VCF/BCF files. It is ethnically inclusive in that it relies on no previous data
and has no parameters, hidden or otherwise, which would favor already
well-studied populations.

It is a reimagined implementation of ldetect.

# Rationale

We want LD-blocks to be computed a priori, so results incorporating LD
information between different analyses are comparable. The naive way to do this
is to segment the genome into non-overlapping blocks, which will split regions
in strong LD over multiple blocks.

# Method

ldivide works in two steps:

- it slides a window over the SNPs to compute the summed R-squared statistic of
each SNP to those within the window. This leads to a vector of
association-scores, one per SNP.
- the vector of association scores is smoothed with low-pass filters of
  increasing width until the desired number of breakpoints are found

# Public results

The results from running our pipeline are provided here: [link to come]

# Features

- very fast and memory efficient
- easy to use
- nice visualization of the results
- inclusive; does not use previous data or assumptions which would favor already
  well-studied populations
- few parameters and none hidden;
  - the number of SNPs in each window used to compute association signals
  - the ones used in the signal processing; which are these?

# Comparison to previous work

ldetect (Berisa and Pickrell, 2015) solved this problem by finding blocks where
the SNPs had strong association signals to other SNPs within the block and weak
to those outside. This method incorporated knowledge about genetic distance when
computing association signals. In the supplementaries, we discuss the downsides
of this, but the main points are that the genetic recombination data have
several problems:

* they are based on heuristics and unproven assumptions
* they are computed on old data
* they are not computed for all populations
* they are not reproducible; the data are provided as-is, without any way to
  understand exactly how they were made

Furthermore, ldetect itself had several problems also discussed in the
supplementaries. The most important of these was that

* the provided ldetect code does not work on real world data, even after arduous setup
* the published datasets contain errors
* the code contains logical errors
* the code is slow
* the method relies on genetic recombination data
