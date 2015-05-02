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

use LogProb;
use call::site::{Site, GenotypeLikelihoods};
use utils;


pub fn call(bcf: &mut bcf::Reader, query: Box<Caller>, fdr: Option<LogProb>, max_prob: Option<LogProb>, threads: usize) -> Vec<(Site, LogProb)> {
    let mut records = bcf.records();

    let mut pool = simple_parallel::Pool::new(threads);
    let call = |likelihoods: Vec<GenotypeLikelihoods>| query.call(&likelihoods);

    let mut site_buffer = Vec::with_capacity(1000);
    let mut prob_buffer = Vec::with_capacity(1000);
    let mut candidates = Vec::with_capacity(1000);

    let max_prob = match (fdr, max_prob) {
        (Some(fdr), Some(p)) => Some(fdr.min(p)),
        (Some(fdr), None)    => Some(fdr),
        (None, Some(p))      => Some(p),
        _                    => None,
    };

    loop {
        site_buffer.extend(records.by_ref().take(1000).map(|record| Site::new(record.ok().expect("Error reading BCF."))));

        if site_buffer.is_empty() {
            if fdr.is_some() {
                control_fdr(&mut candidates, fdr.unwrap());
            }
            return candidates;
        }

        prob_buffer.extend(pool.map(
            site_buffer.iter_mut().map(
                |site| site.genotype_likelihoods().ok().expect("Error reading genotype likelihoods.")
            ), &call)
        );

        let buffer = site_buffer.drain().zip(prob_buffer.drain());
        match max_prob {
            Some(p) => candidates.extend(buffer.filter(|&(_, prob)| prob < p)),
            None    => candidates.extend(buffer),
        }
    }
}


fn control_fdr(candidates: &mut Vec<(Site, LogProb)>, fdr: LogProb) {
    let cmp = |a: &LogProb, b: &LogProb| a.partial_cmp(b).expect("Bug: NaN probability found.");
    let mut probs: Vec<LogProb> = candidates.iter().map(|&(_, prob)| prob).collect();
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
    fn call(&self, likelihoods: &[GenotypeLikelihoods]) -> LogProb;
}
