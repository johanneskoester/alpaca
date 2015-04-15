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
use utils;


pub fn call(bcf: &mut bcf::Reader, query: Box<Caller>, fdr: Prob, threads: usize) -> Vec<(Site, Prob)> {
    let mut records = bcf.records();

    let mut pool = simple_parallel::Pool::new(threads);
    let call = |likelihoods: Vec<GenotypeLikelihoods>| query.call(&likelihoods);

    let mut site_buffer = Vec::with_capacity(1000);
    let mut prob_buffer = Vec::with_capacity(1000);
    let mut candidates = Vec::with_capacity(1000);

    loop {
        site_buffer.extend(records.by_ref().take(1000).map(|record| Site::new(record.ok().expect("Error reading BCF."))));

        if site_buffer.is_empty() {
            control_fdr(&mut candidates, fdr);
            return candidates;
        }

        prob_buffer.extend(pool.map(
            site_buffer.iter_mut().map(
                |site| site.genotype_likelihoods().ok().expect("Error reading genotype likelihoods.")
            ), &call)
        );

        candidates.extend(site_buffer.drain().zip(prob_buffer.drain()).filter(|&(_, prob)| prob < fdr));
    }
}


pub fn control_fdr(candidates: &mut Vec<(Site, Prob)>, fdr: Prob) {
    let cmp = |a: &Prob, b: &Prob| a.partial_cmp(b).expect("Bug: NaN probability found.");
    let mut probs: Vec<Prob> = candidates.iter().map(|&(_, prob)| prob).collect();
    probs.sort_by(&cmp);
    let exp_fdr = utils::log_prob_cumsum(&probs);
    let max_prob = match exp_fdr.binary_search_by(|probe| cmp(&probe, &fdr)) {
        Ok(i)           => probs[i],
        Err(i) if i > 0 => probs[i-1],
        _               => 0.0
    };

    candidates.retain(|&(_, prob)| prob <= max_prob);
}


pub trait Caller: Sync {
    fn call(&self, likelihoods: &[GenotypeLikelihoods]) -> Prob;
}
