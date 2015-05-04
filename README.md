# ALPACA - the ALgebraic PArallel variant CAller
ALPACA is a caller for genetic variants from next-generation sequencing data.
It has two major distinguishing features compared to other variant callers:

* ALPACA incorporates arbitrary filtering of samples against each other into the calling. This is done via an expressive, algebraic query language. During calling, the posterior probability for each locus to not behave like described in the filter query. If that probability is small enough, the locus is called.
* Since the filtering becomes part of the null hypothesis, controlling the FDR becomes easy and intuitive.

Alpaca separates the process into four steps.

* preprocessing of each sample into a BCF file,
* merging preprocessed samples into one BCF file containing only relevant loci,
* filtering irrelevant sites
* calling on the merged BCF file.

The separation allows to add samples later without having to redo all the computations. Since most of the work is done during preprocessing the final calling becomes lightweight and can be repeated with different parameters within seconds.
The algebraic query allows to model calling scenarios in a flexible way, e.g.,

* calling all de-novo mutations of a child: 'child - (mother + father)'
* calling all variants recurrent in at least 3 samples of a group of samples s1,s2,...s5: 's1 x s2 x s3 x s4 x s5 with k = 3'

A complete description of algebraic variant calling can be found in my thesis

> Köster, J. Parallelization, Scalability, and Reproducibility in Next-Generation Sequencing Analysis. PhD-Thesis, TU Dortmund, Germany 2014. ISBN: 978-3737537773.

If you use ALPACA, please cite the thesis for now.

## Example usage:

All in one command:

    $ alpaca preprocess A.bam B.bam C.bam | alpaca filter | alpaca call --fdr 0.05 'A - (B + C)' > calls.bcf

Separate preprocessing and merging (this allows to add samples or change queries without redundant computations; alpaca call usually needs a few seconds):

    $ alpaca preprocess A.bam > A.bcf
    $ alpaca preprocess B.bam > B.bcf
    $ alpaca preprocess C.bam > C.bcf
    $ alpaca merge A.bcf B.bcf C.bcf | alpaca filter > all.bcf
    $ alpaca call --fdr 0.05 'A - (B + C)' < all.bcf > calls.bcf