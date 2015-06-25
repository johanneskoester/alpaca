use std::collections::HashMap;

use Prob;
use call;


pub trait Expr {
    fn samples(&self) -> Vec<String> {
        let mut samples = vec![];
        self.collect_samples(&mut samples);

        samples
    }

    fn collect_samples(&self, samples: &mut Vec<String>);

    fn caller(&self, sample_idx: &HashMap<String, usize>, heterozygosity: Prob) -> Box<call::Caller>;
}


pub struct SampleUnion {
    pub samples: Vec<String>
}


impl Expr for SampleUnion {
    fn collect_samples(&self, samples: &mut Vec<String>) {
        samples.push_all(&self.samples);
    }

    fn caller(&self, sample_idx: &HashMap<String, usize>, heterozygosity: Prob) -> Box<call::Caller> {
        Box::new(call::SampleUnion::new(
            self.samples.iter().map(|sample| *sample_idx.get(sample).expect("Unknown sample name.")).collect(),
            2, // TODO extend to arbitrary ploidy once samtools mpileup has catched up
            heterozygosity,
        ))
    }
}


pub struct Diff {
    pub left: Box<Expr>,
    pub right: Box<Expr>
}


impl Expr for Diff {
    fn collect_samples(&self, samples: &mut Vec<String>) {
        self.left.collect_samples(samples);
        self.right.collect_samples(samples);
    }

    fn caller(&self, sample_idx: &HashMap<String, usize>, heterozygosity: Prob) -> Box<call::Caller> {
        Box::new(call::Diff { left: self.left.caller(sample_idx, heterozygosity), right: self.right.caller(sample_idx, heterozygosity) })
    }
}


pub struct Union {
    pub left: Box<Expr>,
    pub right: Box<Expr>
}


impl Expr for Union {
    fn collect_samples(&self, samples: &mut Vec<String>) {
        self.left.collect_samples(samples);
        self.right.collect_samples(samples);
    }

    fn caller(&self, sample_idx: &HashMap<String, usize>, heterozygosity: Prob) -> Box<call::Caller> {
        Box::new(call::Union { left: self.left.caller(sample_idx, heterozygosity), right: self.right.caller(sample_idx, heterozygosity) })
    }
}


pub struct RelaxedIntersection {
    pub children: Vec<Box<Expr>>,
    pub k: usize
}


impl Expr for RelaxedIntersection {
    fn collect_samples(&self, samples: &mut Vec<String>) {
        for child in self.children.iter() {
            child.collect_samples(samples);
        }
    }

    fn caller(&self, sample_idx: &HashMap<String, usize>, heterozygosity: Prob) -> Box<call::Caller> {
        Box::new(call::RelaxedIntersection {
            children: self.children.iter().map(|expr| expr.caller(sample_idx, heterozygosity)).collect(),
            k: self.k,
        })
    }
}
