use std::f64;
use std::mem;

use {LogProb, Prob};


pub const LOG_TO_PHRED_FACTOR: f64 = -4.3429448190325175; // -10 * 1 / ln(10)
pub const PHRED_TO_LOG_FACTOR: f64 = -0.23025850929940456; // 1 / (-10 * log10(e))


/// Calcualte the sum of the given probabilities in a numerically stable way.
pub fn log_prob_sum(probs: &[LogProb]) -> LogProb {
    let mut pmax = probs[0];
    let mut imax = 0;
    for (i, &p) in probs[1..].iter().enumerate() {
        if p > pmax {
            pmax = p;
            imax = i;
        }
    }
    if pmax == f64::NEG_INFINITY {
        f64::NEG_INFINITY
    }
    else {
        pmax + (probs.iter().enumerate().filter_map(|(i, p)| if i != imax { Some((p - pmax).exp()) } else { None }).sum::<Prob>()).ln_1p()
    }
}


pub fn log_prob_add(mut p0: LogProb, mut p1: LogProb) -> LogProb {
    if p1 > p0 {
        mem::swap(&mut p0, &mut p1);
    }
    if p0 == f64::NEG_INFINITY {
        f64::NEG_INFINITY
    }
    else {
        p0 + (p1 - p0).exp().ln_1p()
    }
}


pub fn log_prob_cumsum(probs: &[LogProb]) -> Vec<LogProb> {
    probs.iter().scan(f64::NEG_INFINITY, |s, p| {
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
