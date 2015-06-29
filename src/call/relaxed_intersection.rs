use std::f64;

use bio::stats::logprobs;

use LogProb;
use call::Caller;
use call::site::GenotypeLikelihoods;
use super::super::utils;


pub struct RelaxedIntersection {
    pub children: Vec<Box<Caller>>,
    pub k: usize,
}


impl Caller for RelaxedIntersection {
    fn call(&self, likelihoods: &[GenotypeLikelihoods]) -> LogProb {
        let prob_ref: Vec<LogProb> = self.children.iter().map(|child| child.call(likelihoods)).collect();

        let mut y = utils::matrix(f64::NEG_INFINITY, self.children.len() + 1, 2);
        y[0][0] = 0.0;

        let mut prob = Vec::with_capacity(self.k);
        prob.push(f64::NEG_INFINITY);
        for k in 0..self.k {
            for j in 1..self.children.len() + 1 {
                let p = prob_ref[j - 1];
                let p_ref = y[j - 1][k % 2] + p;
                let p_var = y[j - 1][(k - 1) % 2] + logprobs::ln_1m_exp(p);
                y[j][k % 2] = logprobs::log_prob_sum(&[p_ref, p_var]);
            }
            prob.push(y[self.children.len()][k % 2]);
        }
        logprobs::log_prob_sum(&prob)
    }
}
