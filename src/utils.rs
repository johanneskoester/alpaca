use std::f64;

use Prob;

pub fn log_prob_sum(probs: &[Prob]) -> Prob {
    let p0 = probs.iter().cloned().fold(f64::NEG_INFINITY, |a, b| a.max(b));
    p0 + (probs.iter().map(|p| (p - p0).exp()).sum::<Prob>() - 1.0).ln_1p()
}


pub fn log_prob_add(p0: Prob, p1: Prob) -> Prob {
    let pmax = p0.max(p1);
    pmax + ((p0 - pmax).exp() + (p1 - pmax).exp() - 1.0).ln_1p()
}


pub fn log_prob_cumsum(probs: &[Prob]) -> Vec<Prob> {
    probs.iter().scan(0.0, |s, p| {
        *s = log_prob_add(*s, *p);
        Some(*s)
    }).collect()
}


pub fn matrix<T: Copy>(v: T, n: usize, m: usize) -> Vec<Vec<T>> {
    let mut matrix = vec![];
    for _ in 0..n {
        matrix.push(vec![v; m]);
    }

    matrix
}
