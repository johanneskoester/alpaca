use simple_parallel;

pub mod site;
pub mod diff;
pub mod relaxed_intersection;
pub mod sample_union;
pub mod union;

pub use call::diff::Diff;
pub use call::relaxed_intersection::RelaxedIntersection;
pub use call::union::Union;
pub use call::sample_union::SampleUnion;

use Prob;
use call::site::{Site, GenotypeLikelihoods};


pub fn query_probabilities<'a, C: Caller, S: Iterator<Item=&'a mut Site> + Send>(query: C, sites: S, pool: &mut simple_parallel::Pool) -> Vec<Prob> {
    let call = |likelihoods: Vec<GenotypeLikelihoods>| query.call(&likelihoods);
    let probs = {    
         pool.map(
            sites.map(|site| site.genotype_likelihoods().ok().expect("Error reading genotype likelihoods")),
            &call
        ).collect()
    };

    probs
}


pub trait Caller: Sync {
    fn call(&self, likelihoods: &[GenotypeLikelihoods]) -> Prob;
}
