extern crate argparse;
extern crate alpaca;


use argparse::{ArgumentParser, Store, List};

use alpaca::cli;


fn main() {

    let mut fasta = "".to_string();
    let mut bams: Vec<String> = vec![];
    let mut threads = 1;
    
    {
        let mut ap = ArgumentParser::new();
        ap.set_description(
"Alpaca is an ALgebraic PArallel variant CAller. 
Preprocess loci of one or more samples for calling with ALPACA.
This command takes a FASTA file with the reference genome, one or more
indexed BAM files, calculates genotype likelihoods for each locus
using samtools mpileup, and prints the result to STDOUT in BCF format.

Example:

$ alpaca-preprocess sample.bam > sample.bcf");


        ap.refer(&mut fasta)
          .add_argument("fasta", Store, "FASTA file with reference genome.");        
        ap.refer(&mut bams)
          .add_argument("bam", List, "BAM files to preprocess.");
        ap.refer(&mut threads)
          .add_option(&["--threads", "-t"], Store, "Number of threads to use.");
        ap.parse_args_or_exit();
    }
    cli::preprocess(&fasta, &bams, threads);
}
