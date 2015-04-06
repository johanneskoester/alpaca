use super::{Caller, GenotypeLikelihoods};
use super::super::utils;


pub struct Union<L: Caller, R: Caller> {
    left: L,
    right: R,
}


impl<L: Caller, R: Caller> Caller for Union<L, R> {
    fn call(&self, likelihoods: &[GenotypeLikelihoods]) -> f64 {
        self.left.call(likelihoods) + self.right.call(likelihoods)
    }
}
