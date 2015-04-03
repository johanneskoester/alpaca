use call::Caller;


pub struct RelaxedIntersection {
    children: Vec<Caller>,
}


impl Caller for RelaxedIntersection {
    fn call(&self, likelihoods: &[GenotypeLikelihoods]) -> f64 {
        let mut y = utils::matrix(f64::MIN, self.children.len() + 1, 2);
        y[0][0] = 0;

        let mut p = f64::MIN;
    }
}
