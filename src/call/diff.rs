use call::Caller;
use call::site::GenotypeLikelihoods;
use Prob;


pub struct Diff {
    pub left: Box<Caller>,
    pub right: Box<Caller>,
}


impl Caller for Diff {
    fn call(&self, likelihoods: &[GenotypeLikelihoods]) -> Prob {
        // 1 - ((1-p_l) * p_r)
        let p = (-self.left.call(likelihoods).exp()).ln_1p() + self.right.call(likelihoods);
        (-p.exp()).ln_1p()
    }
}
