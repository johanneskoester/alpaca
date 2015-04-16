extern crate argparse;
extern crate alpaca;


use argparse::{ArgumentParser, List};

use alpaca::cli;


fn main() {

    let mut bcfs: Vec<String> = vec![];
    
    {
        let mut ap = ArgumentParser::new();
        ap.set_description(
"Alpaca is an ALgebraic PArallel variant CAller. 
Merge preprocessed loci of one or more samples for calling
with ALPACA. This command takes two or more BCF files and
prints a merged BCF suitable for calling with alpaca-call
to STDOUT.

Example:

$ alpaca-merge A.bcf B.bcf C.bcf > merged.bcf");
        ap.refer(&mut bcfs)
          .add_argument("bcf", List, "ALPACA-preprocessed BCF files to merge.");
        ap.parse_args_or_exit();
    }
    cli::merge(&bcfs);
}
