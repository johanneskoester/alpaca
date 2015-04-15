
extern crate alpaca;

use alpaca::cli;

#[test]
fn test_preprocess() {
    
    cli::preprocess(&"mpileup.ref.fa", &[&"mpileup.1.bam"], 3);
}
