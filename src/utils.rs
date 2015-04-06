use std::iter::AdditiveIterator;

pub fn log_prob_sum(probs: &[f64]) -> f64 {
    let p0 = probs.iter().max().unwrap();
    p0 + (probs.iter().map(|p| (p - p0).exp()).sum() - 1.0).ln_1p()
}


pub fn log_prob_add(p0: f64, p1: f64) -> f64 {
    let pmax = p0.max(p1);
    pmax + ((p0 - pmax).exp() + (p1 - pmax).exp() - 1.0).ln_1p()
}


pub fn log_prob_cumsum(probs: &[f64]) {
    probs.iter().scan(0, |s, p| {
        *s = log_prob_add(*s, p);
        Some(*s)
    })
}


pub fn matrix<T: Copy>(v: T, n: usize, m: usize) -> Vec<Vec<T>> {
    let mut matrix = vec![];
    for _ in 0..n {
        matrix.push(vec![v; m]);
    }

    matrix
}
