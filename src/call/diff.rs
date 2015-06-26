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
            ( -( (-p1.exp()).ln_1p() + p2 ).exp() ).ln_1p()
        }
    }
}


impl Caller for Diff {
    fn call(&self, likelihoods: &[GenotypeLikelihoods]) -> LogProb {
        // 1 - ((1-p_l) * p_r)
        info!("{} - {}", self.left.call(likelihoods), self.right.call(likelihoods));
        Diff::diff(self.left.call(likelihoods), self.right.call(likelihoods))
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
/*
    fn test_chr17_pos335() {
        let mut sample_idx = HashMap::new();
        sample_idx.insert("HG00100".to_string(), 0);
        sample_idx.insert("HG00101".to_string(), 0);

        let caller = parse("HG00100 - HG00101", &sample_idx, 0.001);

        let likelihoods = [
            GenotypeLikelihoods::new(
                [14.0,0.0,200.0,38.0,203.0,231.0].iter().map(|p| p * utils::PHRED_TO_LOG_FACTOR).collect(),
                3,
            ),
            GenotypeLikelihoods::new(
                [0.0,0.0,0.0,27,0.0,222].iter().map(|p| p * utils::PHRED_TO_LOG_FACTOR).collect(),
                3,
            ),
        ];

        caller.call(
    }
*/
}
