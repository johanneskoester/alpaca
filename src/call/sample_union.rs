use std::f64;

use call::Caller;
use call::site::GenotypeLikelihoods;
use utils;
use Prob;


pub struct SampleUnion {
    samples: Vec<usize>,
    ploidy: usize,
    heterozygosity: Prob,
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
            (1.0 - (1..self.samples.len() * self.ploidy).map(|i| 1.0 / i as Prob).sum::<Prob>()).ln()
        }
    }

    fn allelefreq_likelihood(&self, m: usize, sample: usize, likelihoods: &[GenotypeLikelihoods]) -> Prob {
        // TODO uncomment:
        // let prior = (1.0 / utils::binomial(3 + m - 1, m)).ln();
        let prior = 0.1;
        utils::log_prob_sum(&likelihoods[sample].with_allelefreq(m)) + prior
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
                    p.push(z[j-1][km_idx] + self.allelefreq_likelihood(m, j, likelihoods));
                }
                z[j][k_idx] = utils::log_prob_sum(&p);
            }
        };

        // calc column k = 0
        calc_col(&mut z, 0);

        let ref_likelihood = z[self.samples.len()][0];
        let mut marginal = ref_likelihood + self.prior(0);

        for k in 1..self.samples.len() * self.ploidy + 1 {
            // set z00 to 1 (i.e. 0 in log space) while it is possible to achive the allele freq with the first sample only.
            // afterwards, set it to -inf like the other z0x values
            if k > self.ploidy {
                z[0][0] = f64::NEG_INFINITY;
            }
            calc_col(&mut z, k);
            marginal = utils::log_prob_sum(&[marginal, z[self.samples.len()][k]]);
        }

        (ref_likelihood, marginal)
    }
}
