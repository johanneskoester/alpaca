
use htslib::bcf;
use itertools::Itertools;

use utils;
use Prob;


pub struct Site {
    record: bcf::record::Record
}


impl Site {
    pub fn new(record: bcf::record::Record) -> Self {
        Site { record: record }
    }

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

    pub fn into_record(mut self, prob: Prob) -> bcf::record::Record {
        {
            let qual = prob * utils::LOG_TO_PHRED_FACTOR;
            self.record.set_qual(qual as f32);

            let likelihoods = self.genotype_likelihoods()
                .ok()
                .expect("Bug: Error reading genotype likelihoods, they should have been read before.");

            let mut fmt = self.record.format(b"GT");
            let mut gtfmt = fmt.integer().ok().expect("Bug: error reading genotype field");

            for (sample, (a, b)) in likelihoods.iter().map(|gl| gl.maximum_likelihood_genotype()).enumerate() {
                // as specified in the BCFv2 docs
                gtfmt[sample][0] = ((a + 1) << 1) as i32;
                gtfmt[sample][1] = ((b + 1) << 1) as i32;
            }
        }        
        self.record
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
            _ => panic!("Bug: Expecting diploid samples.")
        }
    }

    pub fn maximum_likelihood_genotype(&self) -> (u32, u32) {
        // TODO generalize for any ploidy
        let (mut j, mut k) = (0, 0);
        for &l in self.likelihoods.iter() {
            if l == 0.0 {
                return (j, k);
            }
            j += 1;
            if j > k {
                k += 1;
                j = 0;
            }
        }
        panic!("Bug: no likelihood of zero found.");
    }
}
