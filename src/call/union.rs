use call::Caller;
use call::site::GenotypeLikelihoods;


pub struct Union {
    pub left: Box<Caller>,
    pub right: Box<Caller>,
}


impl Caller for Union {
    fn call(&self, likelihoods: &[GenotypeLikelihoods]) -> f64 {
        self.left.call(likelihoods) + self.right.call(likelihoods)
    }
}
