use std::collections::HashMap;

use Prob;
use call;

pub mod ast;


peg_file! parser("grammar.rustpeg");


pub fn parse(query: &str, sample_idx: &HashMap<String, usize>, heterozygosity: Prob) -> Box<call::Caller> {
    let expr = parser::expr(query).ok().expect("Error parsing query.");
    expr.caller(sample_idx, heterozygosity)
}
