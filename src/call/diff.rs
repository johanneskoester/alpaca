use bio::stats::logprobs;

use call::Caller;
use call::site::GenotypeLikelihoods;
use LogProb;


pub struct Diff {
    pub left: Box<Caller>,
    pub right: Box<Caller>,
}


impl Diff {
    fn diff(p1: LogProb, p2: LogProb) -> LogProb {
        if p2 == 0.0 {
            // increase numerical stability by avoiding unnecessary computation
            p1
        }
        else {
            let mut prob = logprobs::ln_1m_exp(p1);
            debug!("1 - p1 = {}", prob);
            prob += p2;
            debug!("(1 - p1) * p2 = {}", prob);
            prob = logprobs::ln_1m_exp(prob);
            debug!("1 - ((1 - p1) * p2) = {}", prob);
            prob
        }
    }
}


impl Caller for Diff {
    fn call(&self, likelihoods: &[GenotypeLikelihoods]) -> LogProb {
        // 1 - ((1-p_l) * p_r)
        let (left, right) = (self.left.call(likelihoods), self.right.call(likelihoods));
        let prob = Diff::diff(left, right);
        debug!("{} - {} = {}", left, right, prob);
        prob
    }
}


#[cfg(test)]
mod tests {
    use std::f64;

    use super::*;
    use LogProb;
    use EPSILON;


    fn eq(a: LogProb, b: LogProb) -> bool {
        (a - b).abs() < EPSILON
    }

    #[test]
    fn test_call() {
        println!("{}", Diff::diff(f64::NEG_INFINITY, 0.0));
        assert_eq!(Diff::diff(f64::NEG_INFINITY, 0.0), f64::NEG_INFINITY);
        assert!(eq(Diff::diff(f64::NEG_INFINITY, f64::NEG_INFINITY), 0.0));
    }
}
