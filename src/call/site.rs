
use htslib::bcf;
use itertools::Itertools;

use utils;
use Prob;


pub struct Site {
    record: bcf::record::Record,
}


impl Site {
    pub fn genotype_likelihoods(&mut self) -> Result<Vec<GenotypeLikelihoods>, bcf::record::TagError> {
        let allele_count = self.record.allele_count() as usize;
        let mut fmt = self.record.format(&b"PL"[..]);
        let pl = try!(fmt.integer());
        Ok(pl.iter().map(|sample_pl| {
            GenotypeLikelihoods {
                likelihoods: sample_pl.iter().map(|&s| s as f64 * utils::PHRED_TO_LOG_FACTOR).collect(),
                allele_count: allele_count,
            }
        }).collect())
    }
}


pub struct GenotypeLikelihoods {
    likelihoods: Vec<Prob>,
    allele_count: usize,
}


impl GenotypeLikelihoods {
    pub fn with_allelefreq(&self, m: usize) -> Vec<Prob> {
        let idx = |j, k| (k * (k + 1) / 2) + j;
        match m {
            0 => vec![0.0],
            1 => (1..self.allele_count).map(|k| self.likelihoods[idx(0, k)]).collect(),
            2 => (1..self.allele_count).cartesian_product(1..self.allele_count)
                                       .map(|(j, k)| self.likelihoods[idx(j, k)]).collect(),
            _ => panic!("Expecting diploid samples.")
        }
    }
}
