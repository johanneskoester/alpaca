use std::path::Path;
use std::convert::AsRef;

use simple_parallel;
use htslib::bcf;

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


pub fn call<P: AsRef<Path>, C: Caller>(path: &P, query: C, fdr: Prob, threads: usize) {
    let reader = bcf::Reader::new(path);
    let mut records = reader.records();

    let mut pool = simple_parallel::Pool::new(threads);
    let call = |likelihoods: &Vec<GenotypeLikelihoods>| query.call(likelihoods);

    let mut buffer = Vec::with_capacity(1000);
    let mut candidates = Vec::with_capacity(1000);

    loop {
        buffer.extend(records.by_ref().take(1000).map(|record| Site::new(record.ok().expect("Error reading BCF."))));
        let likelihoods: Vec<Vec<GenotypeLikelihoods>> = buffer.iter_mut().map(|site| site.genotype_likelihoods()
                                                       .ok()
                                                       .expect("Error reading genotype likelihoods.")
        ).collect();

        candidates.extend(buffer.drain().zip(pool.map(likelihoods.iter(), &call)).filter(|&(_, prob)| prob < fdr));
    }
    // TODO go on with candidates (control FDR, write)
}


pub trait Caller: Sync {
    fn call(&self, likelihoods: &[GenotypeLikelihoods]) -> Prob;
}
