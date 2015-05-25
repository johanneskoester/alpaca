use std::io;
use std::io::Write;

use simple_parallel;
use htslib::bcf;
use itertools::Itertools;

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
    let call = |likelihoods: &[Vec<GenotypeLikelihoods>]| {
        likelihoods.iter().map(|lh| query.call(lh)).collect_vec()
    };

    let mut site_buffer = Vec::with_capacity(1000);
    let mut likelihood_buffer = Vec::with_capacity(1000);
    let mut prob_buffer = Vec::with_capacity(1000);
    let mut candidates = Vec::with_capacity(1000);

    let max_prob = match (fdr, max_prob) {
        (Some(fdr), Some(p)) => Some(fdr.min(p)),
        (Some(fdr), None)    => Some(fdr),
        (None, Some(p))      => Some(p),
        _                    => None,
    };

    let mut processed = 0;
    loop {
        site_buffer.extend(records.by_ref().take(1000000).map(|record| Site::new(record.ok().expect("Error reading BCF."))));
        likelihood_buffer.extend(site_buffer.iter_mut().map(|site|
            site.genotype_likelihoods().ok().expect("Error reading genotype likelihoods.")
        ));

        processed += site_buffer.len();

        if site_buffer.is_empty() {
            if fdr.is_some() {
                control_fdr(&mut candidates, fdr.unwrap());
            }
            return candidates;
        }

        for probs in unsafe { pool.map(likelihood_buffer.chunks(1000000 / threads), &call) } {
            prob_buffer.extend(probs);
        }

        likelihood_buffer.clear();
        let buffer = site_buffer.drain(..).zip(prob_buffer.drain(..));
        match max_prob {
            Some(p) => candidates.extend(buffer.filter(|&(_, prob)| prob < p)),
            None    => candidates.extend(buffer),
        }
        writeln!(io::stderr(), "Processed {} records.", processed).ok().expect("Error writing to STDERR.");
    }
}


fn control_fdr(candidates: &mut Vec<(Site, LogProb)>, fdr: LogProb) {
    let cmp = |a: &LogProb, b: &LogProb| a.partial_cmp(b).expect("Bug: NaN probability found.");
    let mut probs = candidates.iter().map(|&(_, prob)| prob).collect_vec();
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
