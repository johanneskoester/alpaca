extern crate argparse;
extern crate alpaca;

/*Alpaca is an ALgebraic PArallel variant CAller.
It has two major distinguishing features compared to other variant callers like GATK or MuTect:

* ALPACA incorporates the filtering of samples against each other into the calling. This is done via an expressive, algebraic query language. During calling, the posterior probability for each locus to not behave like described in the filter query. If that probability is small enough, the locus is called.
* Since the filtering becomes part of the null hypothesis, controlling the FDR becomes easy and intuitive.

Alpaca separates the process into three steps.

* preprocessing of each sample into a BCF file,
* merging preprocessed samples into one BCF file containing only relevant loci,
* calling on the merged BCF file.

The separation allows to add samples later without having to redo all the computations. Since most of the work is done during preprocessing the final calling becomes lightweight and can be repeated with different parameters within seconds.
The algebraic query allows to model calling scenarios in a flexible way, e.g.,

* calling all de-novo mutations of a child: 'child - (mother + father)'
* calling all variants recurrent in at least 3 samples of a group of samples s1,s2,...s5: 's1 x s2 x s3 x s4 x s5 with k = 3'
*/



use argparse::{ArgumentParser, Store};

use alpaca::cli;


fn main() {

    let mut query = "".to_string();
    let mut fdr = 0.05;
    let mut heterozygosity = 0.001;
    let mut threads = 1;
    
    {
        let mut ap = ArgumentParser::new();
        ap.set_description(
"Alpaca is an ALgebraic PArallel variant CAller. 
This subcommand calls variants. It takes the merged preprocessed
samples in BCF format from STDIN, and prints the found variants in
BCF format to STDOUT.

Example:

$ alpaca-call --fdr 0.05 'child - (father + mother + brother)' < merged.bcf > calls.bcf");
        ap.refer(&mut query)
          .add_argument("query", Store, "Algebraic query describing the filtering scenario.");
        ap.refer(&mut fdr)
          .add_option(&["--fdr"], Store, "Desired false discovery rate.");
        ap.refer(&mut heterozygosity)
          .add_option(&["--het"], Store, "Expected heterozygosity.");
        ap.refer(&mut threads)
          .add_option(&["--threads", "-t"], Store, "Number of threads to use.");
        ap.parse_args_or_exit();
    }

    cli::call(&query, fdr, threads, heterozygosity);
}
