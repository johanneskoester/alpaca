use std::f64;
use std::cmp;

use itertools::Itertools;
use itertools::linspace;
use bio::stats::{logprobs};

use call::{Caller, Dependency};
use call::site::GenotypeLikelihoods;
use utils;
use {Prob, LogProb};


pub struct SampleUnion {
    samples: Vec<usize>,
    ploidy: usize,
    prior: Vec<LogProb>
}


impl SampleUnion {
    pub fn new(samples: Vec<usize>, ploidy: usize, mut heterozygosity: Prob, dependency: Dependency) -> Self {
        // modify heterozygosity according to given dependency
        let het = match dependency {
            Dependency::GivenVariant(dep) if dep > 0.0 => {
                // calc s = sum 1 / i
                let mut s = 0.0;
                for m in 1..samples.len() * ploidy + 1 {
                    s += 1.0 / m as f64;
                }
                // het becomes mixture (by dep) of (het' such that Pr(M=0)=0) and heterozygosity
                dep / s + (1.0 - dep) * heterozygosity
            },
            Dependency::GivenReference(dep) if dep > 0.0 => {
                // set het such that Pr(M=0)=1 if dep = 1
                // dep * 0 + (1 - dep) * heterozygosity
                (1.0 - dep) * heterozygosity
            },
            _ => heterozygosity
        };

        let prior = Self::priors(samples.len(), ploidy, het.ln());

        SampleUnion {
            samples: samples,
            ploidy: ploidy,
            prior: prior
        }
    }

    fn priors(n: usize, ploidy: usize, heterozygosity: LogProb) -> Vec<LogProb> {
        let mut priors = vec![0.0; n * ploidy + 1];
        for m in 1..n * ploidy + 1 {
            priors[m] = heterozygosity - (m as f64).ln();
        }
        // prior for AF=0 is 1 - the rest
        priors[0] = logprobs::ln_1m_exp(logprobs::log_prob_sum(&priors[1..]));
        priors
    }

    fn call_with_prior(&self, allelefreq: usize, likelihoods: &[GenotypeLikelihoods], prior: &[LogProb]) -> LogProb {
        let (allelefreq_likelihoods, marginal) = self.marginal(likelihoods, prior);
        prior[allelefreq] + allelefreq_likelihoods[allelefreq] - marginal
    }
}


impl Caller for SampleUnion {
    fn call(&self, likelihoods: &[GenotypeLikelihoods]) -> LogProb {
        let prob = self.call_with_prior(0, likelihoods, &self.prior);
        debug!("union: {}", prob);
        prob
    }
}


impl SampleUnion {
    fn allelefreq_likelihood(sample: usize, m: usize, likelihoods: &[GenotypeLikelihoods]) -> LogProb {
        let lh = likelihoods[sample].with_allelefreq(m);
        let prior = (1.0 / lh.len() as f64).ln();
        let aflh = logprobs::log_prob_sum(&lh) + prior;
        aflh
    }

    fn marginal(&self, likelihoods: &[GenotypeLikelihoods], prior: &[LogProb]) -> (Vec<LogProb>, LogProb) {
        let mut z: Vec<Vec<LogProb>> = utils::matrix(f64::NEG_INFINITY, self.samples.len() + 1, self.ploidy + 1);
        z[0][0] = 0.0;

        let calc_col = |z: &mut Vec<Vec<LogProb>>, k| {
            // the actual index of k in our representation of z
            let k_idx = k % (self.ploidy + 1);
            // lower limit for sum over z
            let l = cmp::max(0, k as isize - self.ploidy as isize) as usize;

            for j in 1..self.samples.len() + 1 {
                // upper limit for sum over z
                let u = cmp::min(k, (j-1) * self.ploidy);
                if u >= l {
                    // normalize sum over likelihoods by number of summands
                    let path_prior = -((u + 1 - l) as f64).ln();

                    let mut p = vec![];
                    for k_ in l..(u + 1) {
                        let m = k - k_;
                        let lh = Self::allelefreq_likelihood(self.samples[j - 1], m, likelihoods);

                        let km_idx = k_ % (self.ploidy + 1); // TODO fix index
                        p.push(z[j-1][km_idx] + lh + path_prior);
                    }
                    z[j][k_idx] = logprobs::log_prob_sum(&p);
                }
                else {
                    z[j][k_idx] = f64::NEG_INFINITY;
                }
            }
            z[self.samples.len()][k_idx]
        };


        let mut allelefreq_likelihoods = vec![0.0; self.samples.len() * self.ploidy + 1];

        for k in 0..self.samples.len() * self.ploidy + 1 {
            // set z00 to 1 (i.e. 0 in log space) while it is possible to achieve the allele freq with the first sample only.
            // afterwards, set it to -inf like the other z0x values
            if k > self.ploidy {
                z[0][0] = f64::NEG_INFINITY;
            }

            allelefreq_likelihoods[k] = calc_col(&mut z, k);
        }
        let marginal = logprobs::log_prob_sum(&allelefreq_likelihoods.iter().zip(prior).map(|(likelihood, prior)| likelihood + prior).collect_vec());
        assert!(marginal <= 0.00000000001, format!("marginal {} > 0, AFL={:?}, priors={:?}", marginal, allelefreq_likelihoods, prior));

        (allelefreq_likelihoods, marginal.min(0.0))
    }
}


fn het_sum(n: usize, ploidy: usize) -> f64 {
    (1..n * ploidy + 1).map(|i| 1.0 / i as Prob).sum::<Prob>()
}


#[cfg(test)]
mod tests {
    use bio::stats::logprobs;

    use super::*;
    use call::Caller;
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
        vec![ GenotypeLikelihoods::new(vec![Some(logprobs::phred_to_log(167.0)), Some(0.0), Some(logprobs::phred_to_log(96.0))], 2) ]
    }

    fn zero_likelihoods() -> Vec<GenotypeLikelihoods> {
        vec![
            GenotypeLikelihoods::new(vec![Some(0.0), Some(0.0), Some(0.0)], 2),
            GenotypeLikelihoods::new(vec![Some(0.0), Some(0.0), Some(0.0)], 2)
        ]
    }

    fn eq(a: LogProb, b: LogProb) -> bool {
        (a - b).abs() < 0.001
    }

    #[test]
    fn test_prior() {
        assert!(eq(setup(1).prior[0].exp(), 0.998));
        assert!(eq(setup(200).prior[0].exp(), 0.993));

        let union = setup(1);
        assert!(eq(union.prior[1].exp(), 0.001));
        assert!(eq(union.prior[2].exp(), 0.001 / 2.0));
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
        let (allelefreq_likelihoods, marginal) = union.marginal(&likelihoods, &union.prior);
        println!("marginal {} {}", allelefreq_likelihoods[0], marginal);
    }

    #[test]
    fn test_marginal_zero() {
        let union = setup(1);
        let zero_likelihoods = zero_likelihoods();
        let (allelefreq_likelihoods, _) = union.marginal(&zero_likelihoods, &union.prior);
        for p in allelefreq_likelihoods {
            assert!(eq(p, 0.0));
        }
    }

    #[test]
    fn test_marginal_zero_twosample() {
        let union = setup(2);
        let zero_likelihoods = zero_likelihoods();
        let (allelefreq_likelihoods, _) = union.marginal(&zero_likelihoods, &union.prior);
        for p in allelefreq_likelihoods {
            assert!(eq(p, 0.0));
        }
    }

    #[test]
    fn test_call() {
        let union = setup(1);
        let likelihoods = alt_likelihoods();
        let (ref_lh, marginal) = union.marginal(&likelihoods, &union.prior);
        println!("alt marginal {} {}", ref_lh[0], marginal);

        let posterior = union.call(&likelihoods);

        println!("posterior {} {}", posterior, union.prior[0]);
    }
}
