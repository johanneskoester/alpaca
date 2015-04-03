pub trait Caller {
    fn call(&self, likelihoods: &[GenotypeLikelihoods]) -> f64;
}


pub trait GenotypeLikelihoods {
    fn with_allelefreq(&self, m: u16);
}
