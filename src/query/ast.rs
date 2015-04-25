use std::collections::HashMap;

use Prob;
use call;

pub enum Expr {
    SampleUnion(Vec<String>),
    Diff(Box<Expr>, Box<Expr>),
    Union(Box<Expr>, Box<Expr>),
    RelaxedIntersection(Vec<Box<Expr>>, usize)
}


impl Expr {
    pub fn samples(&self) -> Vec<String> {
        let mut samples = vec![];
        self.collect_samples(&mut samples);

        samples
    }

    fn collect_samples(&self, samples: &mut Vec<String>) {
        match self {
            &Expr::SampleUnion(ref _samples) => samples.push_all(_samples),
            &Expr::Diff(ref a, ref b) | &Expr::Union(ref a, ref b) => {
                a.collect_samples(samples);
                b.collect_samples(samples);
            },
            &Expr::RelaxedIntersection(ref children, _) => {
                for child in children.iter() {
                    child.collect_samples(samples);
                }
            }
        }
    }

    pub fn caller(&self, sample_idx: &HashMap<String, usize>, heterozygosity: Prob) -> Box<call::Caller> {
        match self {
            &Expr::SampleUnion(ref samples) => Box::new(call::SampleUnion::new(
                samples.iter().map(|sample| *sample_idx.get(sample).expect("Unknown sample name.")).collect(),
                2, // TODO extend to arbitrary ploidy once samtools mpileup has catched up
                heterozygosity,
            )),
            &Expr::Diff(ref a, ref b)           => Box::new(call::Diff { left: a.caller(sample_idx, heterozygosity), right: b.caller(sample_idx, heterozygosity) }),
            &Expr::Union(ref a, ref b)          => Box::new(call::Union { left: a.caller(sample_idx, heterozygosity), right: b.caller(sample_idx, heterozygosity) }),
            &Expr::RelaxedIntersection(ref exprs, k) => Box::new(call::RelaxedIntersection {
                children: exprs.iter().map(|expr| expr.caller(sample_idx, heterozygosity)).collect(),
                k: k,
            })
        }
    }
}
