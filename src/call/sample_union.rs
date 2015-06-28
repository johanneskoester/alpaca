use std::f64;
use std::mem;

use itertools::Itertools;
use bio;

use call::Caller;
use call::site::GenotypeLikelihoods;
use utils;
use {Prob, LogProb, EPSILON};


pub struct SampleUnion {
    samples: Vec<usize>,
    ploidy: usize,
    prior: Vec<LogProb>,
    path_prior: Vec<LogProb>,
}


impl SampleUnion {
    pub fn new(samples: Vec<usize>, ploidy: usize, heterozygosity: Prob) -> Self {
        let max_allelefreq = (samples.len() * ploidy) as u64;
        let path_prior = (0..max_allelefreq + 1)
            .map(|m| bio::stats::comb(max_allelefreq, m).ln())
            .collect();
        let prior = Self::priors(samples.len(), ploidy, heterozygosity);

        SampleUnion {
            samples: samples,
            ploidy: ploidy,
            prior: prior,
            path_prior: path_prior,
        }
    }

    fn priors(n: usize, ploidy: usize, heterozygosity: Prob) -> Vec<LogProb> {
        let mut priors = vec![0.0; n * ploidy + 1];
        for m in 1..n * ploidy + 1 {
            priors[m] = (heterozygosity / m as f64).ln();
        }
        // prior for AF=0 is 1 - the rest
        priors[0] = (-utils::log_prob_sum(&priors[1..]).exp()).ln_1p();
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
        let aflh = utils::log_prob_sum(&lh) + prior;
        aflh
    }

    fn marginal(&self, likelihoods: &[GenotypeLikelihoods], prior: &[LogProb]) -> (Vec<LogProb>, LogProb) {
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


        let mut allelefreq_likelihoods = vec![0.0; self.samples.len() * self.ploidy + 1];

        // calc column k = 0
        allelefreq_likelihoods[0] = calc_col(&mut z, 0);
        let mut marginal = allelefreq_likelihoods[0] - self.path_prior[0] + prior[0];
        assert!(marginal <= 0.0, format!("AF={}: marginal {} > 0, AFL={}, PP={}, P={}", 0, marginal, allelefreq_likelihoods[0], self.path_prior[0], prior[0]));

        for k in 1..self.samples.len() * self.ploidy + 1 {
            // set z00 to 1 (i.e. 0 in log space) while it is possible to achive the allele freq with the first sample only.
            // afterwards, set it to -inf like the other z0x values
            if k > self.ploidy {
                z[0][0] = f64::NEG_INFINITY;
            }

            allelefreq_likelihoods[k] = calc_col(&mut z, k);

            let _marginal = utils::log_prob_add(marginal, allelefreq_likelihoods[k] - self.path_prior[k] + prior[k]);
            assert!(_marginal <= 0.0, format!("AF={}: marginal {} > 0, prev_marginal={}, AFL={}, PP={}, P={}", k, _marginal, marginal, allelefreq_likelihoods[k], self.path_prior[k], prior[k]));
            marginal = _marginal;
        }

        (allelefreq_likelihoods, marginal)
    }
}


pub struct DependentSampleUnion {
    population: SampleUnion,
    union: SampleUnion,
    het_sum: f64,
    min_prior: f64,
    max_prior: f64
}


impl DependentSampleUnion {
    pub fn new(population: Vec<usize>, samples: Vec<usize>, ploidy: usize, heterozygosity: Prob) -> Self {
        let filtered_population = population.iter().filter(|s| !samples.contains(s)).cloned().collect_vec();
        let het_sum = (1..filtered_population.len() * ploidy + 1).map(|i| 1.0 / i as Prob).sum::<Prob>().ln();
        let union = SampleUnion::new(samples, ploidy, heterozygosity);

        let mut min_prior = *union.prior.last().unwrap();
        let mut max_prior = union.prior[0];
        if min_prior > max_prior {
            mem::swap(&mut min_prior, &mut max_prior);
        }
        DependentSampleUnion {
            population: SampleUnion::new(filtered_population, ploidy, heterozygosity),
            union: union,
            het_sum: het_sum,
            min_prior: min_prior,
            max_prior: max_prior
        }
    }
}


impl Caller for DependentSampleUnion {
    fn call(&self, likelihoods: &[GenotypeLikelihoods]) -> LogProb {
        if self.population.samples.is_empty() {
            self.union.call(likelihoods)
        }
        else {
            let population_ref = self.population.call(likelihoods).min(self.max_prior).max(self.min_prior);
            let het = (-population_ref.exp()).ln_1p() - self.het_sum;
            let prior = SampleUnion::priors(self.union.samples.len(), self.union.ploidy, het.exp());
            let prob = self.union.call_with_prior(0, likelihoods, &prior);
            debug!("union: {}", prob);
            prob
        }
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

        let (ref_lh, marginal) = union.marginal(&likelihoods, &union.prior);
        println!("marginal {} {}", ref_lh[0], marginal);
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
