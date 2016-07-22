# ALPACA - the ALgebraic PArallel variant CAller

**ALPACA is currently under development and not yet intended for productive use.**

ALPACA is a caller for genomic variants (single nucleotide and small indels) from next-generation sequencing data.
It has two major distinguishing features compared to other variant callers:

* ALPACA incorporates arbitrary filtering of samples against each other into the calling. This is done via an expressive, algebraic query language. It calculates the posterior probability for each locus to not behave like described in the filter query. If that probability is small enough, the locus is called.
* Since the filtering is part of the null hypothesis, controlling the FDR becomes easy and intuitive.

The algebraic query language allows to model calling scenarios in a flexible way, e.g.,

* calling all de-novo mutations of a child: 'child - (mother + father)'
* calling all variants recurrent in at least 3 samples of a group of samples s1,s2,...s5: 's1 x s2 x s3 x s4 x s5 with k = 3'

A complete description of algebraic variant calling can be found in my thesis

> Köster, J. Parallelization, Scalability, and Reproducibility in Next-Generation Sequencing Analysis. PhD-Thesis, TU Dortmund, Germany 2014. ISBN: 978-3737537773.

If you use ALPACA, please cite the thesis for now.

## Installation

The easiest way to install ALPACA is to use the conda package manager (http://conda.pydata.org/docs/install/quick.html). With conda installed, you can issue

    $ conda install -c johanneskoester alpaca

to install or update ALPACA.


## Usage

ALPACA relies on having preprocessed genotype likelihoods for each locus of interest.
For a set of samples, these can be obtained by performing a joint variant calling with e.g.
Samtools, Freebayes or GATK or any other small variant caller that either provides a GL or a
PL field per sample. ALPACA works best if the obtained calls are available in BCF format.

For example, all variants present in sample A, but not in sample B or C while controlling the FDR at 5% can be obtained with

    $ alpaca call --fdr 0.05 'A - (B + C)' < likelihoods.bcf > calls.bcf

with `likelihoods.bcf` being the preprocessed genotype likelihoods for a set of samples consisting of at least sample A, B, and C.

## Author

[Johannes Köster](https://github.com/johanneskoester)

## License

Licensed under the [MIT license](http://opensource.org/licenses/MIT). This project may not be copied, modified, or distributed except according to those terms.
