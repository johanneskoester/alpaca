use htslib::bcf;
use itertools::Itertools;
use bio::stats::logprobs;

use utils;
use LogProb;


pub struct Site {
    pub record: bcf::record::Record
}


impl Site {
    pub fn new(mut record: bcf::record::Record) -> Self {
        // ensure that format data is cached into record (for thread safety)
        record.format(&b"PL"[..]);
        Site { record: record }
    }

    pub fn genotype_likelihoods(&mut self) -> Result<Vec<GenotypeLikelihoods>, bcf::record::TagError> {
        let allele_count = self.record.allele_count() as usize;
        let mut fmt = self.record.format(&b"PL"[..]);
        let pl = try!(fmt.integer());
        Ok(pl.iter().map(|sample_pl| {
            let likelihoods = sample_pl.iter().map(|&s| {
                if s < 0 {
                    None
                }
                else {
                    Some(logprobs::phred_to_log(s as LogProb))
                }
            }).collect();

            GenotypeLikelihoods::new(likelihoods, allele_count)
        }).collect())
    }

    pub fn set_qual(&mut self, prob: LogProb) {
        let qual = logprobs::log_to_phred(prob);
        self.record.set_qual(qual as f32);
    }

    /// Update the record with calling information, this has to happen after translate and subset.
    pub fn calc_genotype(&mut self) {
        let likelihoods = self.genotype_likelihoods()
            .ok()
            .expect("Bug: Error reading genotype likelihoods, they should have been read before.");
        let sample_count = likelihoods.len();

        let mut genotypes = vec![0; sample_count * 2];  // TODO generalize ploidy

        let mut i = 0;
        for (a,b) in likelihoods.iter().map(|lh| lh.maximum_likelihood_genotype()) {
            // as specified in the BCFv2 docs
            genotypes[i] = (a + 1) << 1;
            genotypes[i + 1] = (b + 1) << 1;
            i += 2;
        }
        self.record.push_format_integer(b"GT", &genotypes).ok().expect("Error writing genotype.");
    }
}


pub struct GenotypeLikelihoods {
    likelihoods: Vec<Option<LogProb>>,
    allele_count: usize,
    unknown: bool,
}


impl GenotypeLikelihoods {
    pub fn new(likelihoods: Vec<Option<LogProb>>, allele_count: usize) -> Self {
        let unknown = likelihoods.iter().all(|l| l.is_none());
        GenotypeLikelihoods { likelihoods: likelihoods, allele_count: allele_count, unknown: unknown }
    }

    pub fn with_allelefreq(&self, m: usize) -> Vec<LogProb> {
        let idx = |j, k| (k * (k + 1) / 2) + j;
        match m {
            // in case of no coverage, all log-likelihoods are 0
            _ if self.unknown => vec![0.0],
            // if ref is unknown, we always assume a log-likelihood of 0
            0 => vec![ self.likelihoods[0].expect("Bug: reference likelihood can only be unknown if all likelihoods are unknown.") ],
            // in the other cases, filter out unknown genotypes.
            1 => (1..self.allele_count).filter_map(|k| self.likelihoods[idx(0, k)]).collect(),
            2 => (1..self.allele_count).cartesian_product(1..self.allele_count)
                                       .filter_map(|(j, k)| self.likelihoods[idx(j, k)]).collect(),
            _ => panic!("Bug: Expecting diploid samples.")
        }
    }

    pub fn maximum_likelihood_genotype(&self) -> (i32, i32) {
        // TODO generalize for any ploidy
        if self.unknown {
            return (0, 0);
        }

        let (mut j, mut k) = (0, 0);
        for &l in self.likelihoods.iter() {
            if !l.is_none() && l.unwrap() == 0.0 {
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
