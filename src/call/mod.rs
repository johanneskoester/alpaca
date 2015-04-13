pub mod site;
pub mod diff;
pub mod relaxed_intersection;
pub mod sample_union;
pub mod union;

use Prob;
use call::site::{Site, GenotypeLikelihoods};


pub fn query_probabilities<'a, C: Caller, S: Iterator<Item=&'a mut Site>>(query: C, sites: S) -> Vec<Prob> {
    sites.map(
        |site| query.call(&site.genotype_likelihoods().ok().expect("Error reading genotype likelihoods"))
    ).collect()
}


pub trait Caller {
    fn call(&self, likelihoods: &[GenotypeLikelihoods]) -> Prob;
}
