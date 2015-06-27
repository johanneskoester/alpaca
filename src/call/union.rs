use call::Caller;
use call::site::GenotypeLikelihoods;
use LogProb;


pub struct Union {
    pub left: Box<Caller>,
    pub right: Box<Caller>,
}


impl Caller for Union {
    fn call(&self, likelihoods: &[GenotypeLikelihoods]) -> LogProb {
        let (left, right) = (self.left.call(likelihoods), self.right.call(likelihoods));
        debug!("{} - {}", left, right);
        left + right
    }
}
