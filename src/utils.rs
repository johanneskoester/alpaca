use std::f64;

use LogProb;


// TODO check
pub const LOG_TO_PHRED_FACTOR: f64 = -4.3429448190325175; // -10 * 1 / ln(10)
pub const PHRED_TO_LOG_FACTOR: f64 = -0.23025850929940456; // 1 / (-10 * log10(e))


pub fn log_prob_sum(probs: &[LogProb]) -> LogProb {
    let p0 = probs.iter().cloned().fold(f64::NEG_INFINITY, |a, b| a.max(b));
    if p0 == f64::NEG_INFINITY {
        f64::NEG_INFINITY
    }
    else {
        p0 + (probs.iter().map(|p| (p - p0).exp()).sum::<LogProb>()).ln()
    }
}


pub fn log_prob_add(p0: LogProb, p1: LogProb) -> LogProb {
    let pmax = p0.max(p1);
    if pmax == f64::NEG_INFINITY {
        f64::NEG_INFINITY
    }
    else {
        pmax + ((p0 - pmax).exp() + (p1 - pmax).exp()).ln()
    }
}


pub fn log_prob_cumsum(probs: &[LogProb]) -> Vec<LogProb> {
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
