use std::collections::HashMap;
use std::str;


use Prob;
use call;

pub mod ast;

mod grammar;

pub struct Parser {
    sample_idx: HashMap<String, usize>,
    heterozygosity: Prob,
    dependency: Prob
}


impl Parser {
    pub fn new(samples: &[&[u8]], heterozygosity: Prob, dependency: Prob) -> Self {
        let mut sample_idx = HashMap::new();
        for (idx, sample) in samples.iter().enumerate() {
            sample_idx.insert(str::from_utf8(sample).ok().expect("Invalid sample name in BCF.").to_string(), idx);
        }


        Parser {
            sample_idx: sample_idx,
            heterozygosity: heterozygosity,
            dependency: dependency
        }
    }

    pub fn parse(&self, query: &str) -> (Box<call::Caller>, Vec<String>) {
        let expr = grammar::expr(query).unwrap();
        let mut samples = Vec::new();
        self.collect_samples(&expr, &mut samples);
        let population = self.sample_idx(&samples);

        (self.caller(&expr, &population, call::Dependency::None), samples)
    }

    fn caller(&self, expr: &Box<ast::Expr>, population: &[usize], dependency: call::Dependency) -> Box<call::Caller> {
        match **expr {
            ast::Expr::SampleUnion(ref samples) => {
                Box::new(call::SampleUnion::new(
                    self.sample_idx(samples),
                    2, // TODO extend to arbitrary ploidy once samtools mpileup has catched up
                    self.heterozygosity,
                    dependency
                ))
            },
            ast::Expr::Diff(ref a, ref b)           => Box::new(call::Diff {
                    left: self.caller(a, population, dependency),
                    right: self.caller(b, population, call::Dependency::GivenVariant(self.dependency))
            }),
            ast::Expr::Union(ref a, ref b)          => Box::new(call::Union {
                left: self.caller(a, population, dependency),
                right: self.caller(b, population, call::Dependency::GivenReference(self.dependency))
            }),
            ast::Expr::RelaxedIntersection(ref exprs, k) => Box::new(call::RelaxedIntersection {
                children: exprs.iter().map(|expr| self.caller(expr, population, dependency)).collect(),
                k: k,
            })
        }
    }

    fn collect_samples(&self, expr: &Box<ast::Expr>, samples: &mut Vec<String>) {
        match **expr {
            ast::Expr::SampleUnion(ref _samples) => samples.extend_from_slice(&_samples),
            ast::Expr::Diff(ref a, ref b) | ast::Expr::Union(ref a, ref b) => {
                self.collect_samples(a, samples);
                self.collect_samples(b, samples);
            },
            ast::Expr::RelaxedIntersection(ref children, _) => {
                for child in children.iter() {
                    self.collect_samples(child, samples);
                }
            }
        }
    }

    fn sample_idx(&self, samples: &[String]) -> Vec<usize> {
        samples.iter().map(|sample| *self.sample_idx.get(sample).expect("Unknown sample name.")).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse() {
        let samples = vec![&b"A.0"[..], &b"B_01"[..], &b"C0"[..]];
        let parser = Parser::new(&samples, 0.001, 0.0);

        parser.parse("A.0 - B_01");
        parser.parse("A.0");
        parser.parse("A.0 - (B_01 + C0)");
        parser.parse("A.0 x B_01 x C0 with k = 2");
    }
}
