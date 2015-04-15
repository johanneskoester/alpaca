use std::collections::HashMap;
use std::str;


use htslib::bcf;

use Prob;
use call;

pub mod ast;


peg_file! parser("grammar.rustpeg");


pub fn parse(query: &str, sample_idx: &HashMap<String, usize>, heterozygosity: Prob) -> Box<call::Caller> {
    let expr = parser::expr(query).ok().expect("Error parsing query.");
    expr.caller(sample_idx, heterozygosity)
}


pub fn sample_index(bcf: &bcf::Reader) -> HashMap<String, usize> {
    let mut sample_idx = HashMap::new();
    for (idx, sample) in bcf.header.samples().iter().enumerate() {
        sample_idx.insert(str::from_utf8(sample).ok().expect("Invalid sample name in BCF.").to_string(), idx);
    }

    sample_idx
}
