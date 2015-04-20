use std::collections::HashMap;
use std::str;


use htslib::bcf;

use Prob;
use call;

pub mod ast;


peg_file! parser("grammar.rustpeg");


pub fn parse(query: &str, sample_idx: &HashMap<String, usize>, heterozygosity: Prob) -> Box<call::Caller> {
    let expr = parser::expr(query).unwrap();
    expr.caller(sample_idx, heterozygosity)
}


pub fn sample_index(bcf: &bcf::Reader) -> HashMap<String, usize> {
    let mut sample_idx = HashMap::new();
    for (idx, sample) in bcf.header.samples().iter().enumerate() {
        sample_idx.insert(str::from_utf8(sample).ok().expect("Invalid sample name in BCF.").to_string(), idx);
    }

    sample_idx
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_parse() {
        let mut sample_idx = HashMap::new();
        sample_idx.insert("A.0".to_string(), 0);
        sample_idx.insert("B_0".to_string(), 0);
        sample_idx.insert("C0".to_string(), 0);

        parse("A.0 - B_0", &sample_idx, 0.001);
        parse("A.0", &sample_idx, 0.001);
        parse("A.0 - (B_0 + C0)", &sample_idx, 0.001);
        parse("A.0 x B_0 x C0 with k = 2", &sample_idx, 0.001);
    }
}
