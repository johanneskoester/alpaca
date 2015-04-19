use std::f64;

use call::Caller;
use call::site::GenotypeLikelihoods;
use utils;
use Prob;


pub struct SampleUnion {
    samples: Vec<usize>,
    ploidy: usize,
    heterozygosity: Prob,
    ref_prior: Prob,
}


impl SampleUnion {
    pub fn new(samples: Vec<usize>, ploidy: usize, heterozygosity: Prob) -> Self {
        SampleUnion {
            samples: vec![],
            ploidy: ploidy,
            heterozygosity: heterozygosity.ln(),
            ref_prior: (-(1..samples.len() * ploidy + 1).map(|i| 1.0 / i as Prob).sum::<Prob>() * heterozygosity).ln_1p()
        }
    }
}


impl Caller for SampleUnion {
    fn call(&self, likelihoods: &[GenotypeLikelihoods]) -> Prob {
        let (ref_likelihood, marginal) = self.marginal(likelihoods);
        self.prior(0) + ref_likelihood - marginal
    }
}


impl SampleUnion {
    fn prior(&self, m: usize) -> Prob {
        if m > 0 {
            self.heterozygosity - (m as Prob).ln()
        }
        else {
            self.ref_prior
        }
    }

    fn allelefreq_likelihood(sample: usize, m: usize, likelihoods: &[GenotypeLikelihoods]) -> Prob {
        let lh = likelihoods[sample].with_allelefreq(m);
        let prior = (1.0 / lh.len() as f64).ln();
        utils::log_prob_sum(&lh) + prior
    }

    fn marginal(&self, likelihoods: &[GenotypeLikelihoods]) -> (Prob, Prob) {
        let mut z: Vec<Vec<Prob>> = utils::matrix(f64::NEG_INFINITY, self.samples.len() + 1, self.ploidy + 1);
        z[0][0] = 0.0;

        let calc_col = |z: &mut Vec<Vec<Prob>>, k| {
            // the actual index of k in our representation of z
            let k_idx = k % (self.ploidy + 1);

            for j in 1..self.samples.len() + 1 {
                let mut p = vec![];
                for m in 0..self.ploidy + 1 {
                    // the actual index of k - m in our representation of z
                    let km_idx = (if k >= m { k as i32 } else { 0i32 } - m as i32).abs() as usize % (self.ploidy + 1);
                    p.push(z[j-1][km_idx] + Self::allelefreq_likelihood(j - 1, m, likelihoods));
                }
                z[j][k_idx] = utils::log_prob_sum(&p);
            }
            z[self.samples.len()][k_idx]
        };

        
        // calc column k = 0
        let ref_likelihood = calc_col(&mut z, 0);
        let mut marginal = ref_likelihood + self.prior(0);

        for k in 1..self.samples.len() * self.ploidy + 1 {
            // set z00 to 1 (i.e. 0 in log space) while it is possible to achive the allele freq with the first sample only.
            // afterwards, set it to -inf like the other z0x values
            if k > self.ploidy {
                z[0][0] = f64::NEG_INFINITY;
            }

            let likelihood = calc_col(&mut z, k);

            marginal = utils::log_prob_sum(&[marginal, likelihood + self.prior(k)]);
        }

        assert!(marginal <= 0.0);

        (ref_likelihood, marginal)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use Prob;
    use call::site::GenotypeLikelihoods;

    fn setup(n: usize) -> SampleUnion {
        assert!(n > 0);
        SampleUnion::new((0..n).collect(), 2, 0.001)
    }

    fn likelihoods() -> Vec<GenotypeLikelihoods> {
        vec![ GenotypeLikelihoods::new(vec![(1.0f64).ln(), (0.1f64).ln(), (0.0002f64).ln(), (0.0002f64).ln(), (0.0002f64).ln(), (0.0002f64).ln()], 3) ]
    }

    fn eq(a: Prob, b: Prob) -> bool {
        (a - b).abs() < 0.001
    }

    #[test]
    fn test_prior() {
        let epsilon = 0.001;

        assert!(eq(setup(1).prior(0).exp(), 0.998));
        assert!(eq(setup(200).prior(0).exp(), 0.993));

        let union = setup(1);
        assert!(eq(union.prior(1).exp(), union.heterozygosity.exp()));
        assert!(eq(union.prior(2).exp(), union.heterozygosity.exp() / 2.0));
    }

    #[test]
    fn test_allelefreq_likelihood() {
        let likelihoods = likelihoods();
        assert!(eq(SampleUnion::allelefreq_likelihood(0, 0, &likelihoods).exp(), 1.0));
        assert!(eq(SampleUnion::allelefreq_likelihood(0, 1, &likelihoods).exp(), 0.0501));
        assert!(eq(SampleUnion::allelefreq_likelihood(0, 2, &likelihoods).exp(), 0.0002));
    }

    #[test]
    fn test_marginal() {
        let union = setup(1);
        let likelihoods = likelihoods();

        let (ref_lh, marginal) = union.marginal(&likelihoods);
        println!("{} {}", ref_lh, marginal);
    }
}
