use std::f64;

use bio;

use call::Caller;
use call::site::GenotypeLikelihoods;
use utils;
use {Prob, LogProb};


pub struct SampleUnion {
    samples: Vec<usize>,
    ploidy: usize,
    heterozygosity: LogProb,
    ref_prior: LogProb,
    path_prior: Vec<LogProb>,
}


impl SampleUnion {
    pub fn new(samples: Vec<usize>, ploidy: usize, heterozygosity: Prob) -> Self {
        let ref_prior = (-(1..samples.len() * ploidy + 1)
            .map(|i| 1.0 / i as LogProb)
            .sum::<LogProb>() * heterozygosity)
            .ln_1p();

        let max_allelefreq = (samples.len() * ploidy) as u64;
        let path_prior = (0..max_allelefreq + 1)
            .map(|m| (bio::stats::comb(max_allelefreq, m) as f64).ln())
            .collect();

        SampleUnion {
            samples: samples,
            ploidy: ploidy,
            heterozygosity: heterozygosity.ln(),
            ref_prior: ref_prior,
            path_prior: path_prior,
        }
    }
}


impl Caller for SampleUnion {
    fn call(&self, likelihoods: &[GenotypeLikelihoods]) -> LogProb {
        let (ref_likelihood, marginal) = self.marginal(likelihoods);
        self.prior(0) + ref_likelihood - marginal
    }
}


impl SampleUnion {
    fn prior(&self, m: usize) -> LogProb {
        if m > 0 {
            self.heterozygosity - (m as LogProb).ln()
        }
        else {
            self.ref_prior
        }
    }

    fn allelefreq_likelihood(sample: usize, m: usize, likelihoods: &[GenotypeLikelihoods]) -> LogProb {
        let lh = likelihoods[sample].with_allelefreq(m);
        let prior = (1.0 / lh.len() as f64).ln();
        let aflh = utils::log_prob_sum(&lh) + prior;
        aflh
    }

    fn marginal(&self, likelihoods: &[GenotypeLikelihoods]) -> (LogProb, LogProb) {
        let mut z: Vec<Vec<LogProb>> = utils::matrix(f64::NEG_INFINITY, self.samples.len() + 1, self.ploidy + 1);
        z[0][0] = 0.0;

        let calc_col = |z: &mut Vec<Vec<LogProb>>, k| {
            // the actual index of k in our representation of z
            let k_idx = k % (self.ploidy + 1);

            for j in 1..self.samples.len() + 1 {
                let mut p = vec![];
                for m in 0..self.ploidy + 1 {
                    // the actual index of k - m in our representation of z
                    let km_idx = (if k >= m { k as i32 } else { 0i32 } - m as i32).abs() as usize % (self.ploidy + 1);
                    let lh = Self::allelefreq_likelihood(self.samples[j - 1], m, likelihoods);
                    p.push(z[j-1][km_idx] + lh);
                }
                z[j][k_idx] = utils::log_prob_sum(&p);
            }
            z[self.samples.len()][k_idx]
        };


        // calc column k = 0
        let ref_likelihood = calc_col(&mut z, 0);
        let mut marginal = ref_likelihood - self.path_prior[0] + self.prior(0);
        assert!(marginal <= 0.0, format!("marginal {} > 0", marginal));

        for k in 1..self.samples.len() * self.ploidy + 1 {
            // set z00 to 1 (i.e. 0 in log space) while it is possible to achive the allele freq with the first sample only.
            // afterwards, set it to -inf like the other z0x values
            if k > self.ploidy {
                z[0][0] = f64::NEG_INFINITY;
            }

            let likelihood = calc_col(&mut z, k);

            marginal = utils::log_prob_sum(&[marginal, likelihood - self.path_prior[k] + self.prior(k)]);
            assert!(marginal <= 0.0, format!("marginal {} > 0, {}, {}", marginal, likelihood, self.path_prior[k]));
        }

        (ref_likelihood, marginal)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use call::Caller;
    use utils;
    use LogProb;
    use call::site::GenotypeLikelihoods;

    fn setup(n: usize) -> SampleUnion {
        assert!(n > 0);
        SampleUnion::new((0..n).collect(), 2, 0.001)
    }

    fn ref_likelihoods() -> Vec<GenotypeLikelihoods> {
        vec![ GenotypeLikelihoods::new(vec![Some((1.0f64).ln()), Some((0.1f64).ln()), Some((0.0002f64).ln()), Some((0.0002f64).ln()), Some((0.0002f64).ln()), Some((0.0002f64).ln())], 3) ]
    }

    fn alt_likelihoods() -> Vec<GenotypeLikelihoods> {
        vec![ GenotypeLikelihoods::new(vec![Some(167.0 * utils::PHRED_TO_LOG_FACTOR), Some(0.0), Some(96.0 * utils::PHRED_TO_LOG_FACTOR)], 2) ]
    }

    fn eq(a: LogProb, b: LogProb) -> bool {
        (a - b).abs() < 0.001
    }

    #[test]
    fn test_prior() {
        assert!(eq(setup(1).prior(0).exp(), 0.998));
        // TODO this causes overflow when calculating combinations for path_prior.
        //assert!(eq(setup(200).prior(0).exp(), 0.993));

        let union = setup(1);
        assert!(eq(union.prior(1).exp(), union.heterozygosity.exp()));
        assert!(eq(union.prior(2).exp(), union.heterozygosity.exp() / 2.0));
    }

    #[test]
    fn test_allelefreq_likelihood() {
        let likelihoods = ref_likelihoods();
        assert!(eq(SampleUnion::allelefreq_likelihood(0, 0, &likelihoods).exp(), 1.0));
        assert!(eq(SampleUnion::allelefreq_likelihood(0, 1, &likelihoods).exp(), 0.0501));
        assert!(eq(SampleUnion::allelefreq_likelihood(0, 2, &likelihoods).exp(), 0.0002));
        println!("alt lh {}", SampleUnion::allelefreq_likelihood(0, 1, &likelihoods).exp());
    }

    #[test]
    fn test_marginal() {
        let union = setup(1);
        let likelihoods = ref_likelihoods();

        let (ref_lh, marginal) = union.marginal(&likelihoods);
        println!("marginal {} {}", ref_lh, marginal);
    }

    #[test]
    fn test_call() {
        let union = setup(1);
        let likelihoods = alt_likelihoods();
        let (ref_lh, marginal) = union.marginal(&likelihoods);
        println!("alt marginal {} {}", ref_lh, marginal);

        let posterior = union.call(&likelihoods);

        println!("posterior {} {}", posterior, union.prior(0));
    }
}
