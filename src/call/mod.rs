use htslib::bcf;
use itertools::Itertools;
use bio::stats::logprobs;
use cue;

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
use Prob;
use call::site::{Site, GenotypeLikelihoods};


#[derive(Copy)]
#[derive(Clone)]
pub enum Dependency {
    GivenVariant(Prob),
    GivenReference(Prob),
    None
}


pub fn call(bcf: &mut bcf::Reader, query: Box<Caller>, fdr: Option<LogProb>, max_prob: Option<LogProb>, threads: usize) -> Vec<(Site, LogProb)> {
    let buffer_size = 1000000;

    let max_prob = match (fdr, max_prob) {
        (Some(fdr), Some(p)) => fdr.min(p),
        (Some(fdr), None)    => fdr,
        (None, Some(p))      => p,
        _                    => 0.0,
    };

    // read from BCF
    let records = bcf.records().map(|record| Site::new(record.ok().expect("Error reading BCF")));
    let mut candidates = Vec::with_capacity(buffer_size);

    // calculate posterior probabilities in parallel
    cue::pipeline(
        "calling",
        threads,
        records,
        |mut site| {
            let likelihoods = site.genotype_likelihoods().ok().expect("Error reading genotype likelihoods.");
            let posterior = query.call(&likelihoods);
            (site, posterior)
        },
        |(site, prob)| {
            if prob <= max_prob {
                candidates.push((site, prob));
            }
        }
    );

    // control FDR
    if let Some(fdr) = fdr {
        control_fdr(&mut candidates, fdr);
    }
    return candidates;
}


fn control_fdr(candidates: &mut Vec<(Site, LogProb)>, fdr: LogProb) {
    let cmp = |a: &LogProb, b: &LogProb| a.partial_cmp(b).expect("Bug: NaN probability found.");
    let mut probs = candidates.iter().map(|&(_, prob)| prob).collect_vec();
    probs.sort_by(&cmp);
    let exp_fdr = logprobs::cumsum(probs.iter().cloned()).collect_vec();
    let max_prob = match exp_fdr.binary_search_by(|probe| cmp(&probe, &fdr)) {
        Ok(i)           => probs[i],
        Err(i) if i > 0 => probs[i-1],
        _               => 0.0
    };

    candidates.retain(|&(_, prob)| prob <= max_prob);
}


pub trait Caller: Sync + Send {
    fn call(&self, likelihoods: &[GenotypeLikelihoods]) -> LogProb;
}
