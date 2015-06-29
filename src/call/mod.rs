use std::io;
use std::io::Write;
use std::thread;
use std::sync::mpsc;

use simple_parallel;
use htslib::bcf;
use itertools::Itertools;
use bio::stats::logprobs;

pub mod site;
pub mod diff;
pub mod relaxed_intersection;
pub mod sample_union;
pub mod union;

pub use call::diff::Diff;
pub use call::relaxed_intersection::RelaxedIntersection;
pub use call::union::Union;
pub use call::sample_union::{SampleUnion, DependentSampleUnion};

use LogProb;
use call::site::{Site, GenotypeLikelihoods};
use utils;


pub fn call(bcf: &mut bcf::Reader, query: Box<Caller>, fdr: Option<LogProb>, max_prob: Option<LogProb>, threads: usize) -> Vec<(Site, LogProb)> {
    let buffer_size = 1000000;

    let max_prob = match (fdr, max_prob) {
        (Some(fdr), Some(p)) => Some(fdr.min(p)),
        (Some(fdr), None)    => Some(fdr),
        (None, Some(p))      => Some(p),
        _                    => None,
    };

    // read chunks from BCF
    let (site_channel_in, site_channel_out) = mpsc::channel();
    let mut records = bcf.records();
    loop {
        let sites = records.by_ref().take(buffer_size).map(|record| Site::new(record.ok().expect("Error reading BCF."))).collect_vec();
        if sites.is_empty() {
            site_channel_in.send(None).ok().expect("Error using channel.");
            break;
        }

        site_channel_in.send(Some(sites)).ok().expect("Error using channel.");
    }

    let processor = thread::spawn(move || {
        let mut pool = simple_parallel::Pool::new(threads);

        let mut candidates = Vec::with_capacity(buffer_size);
        let mut likelihood_buffer = Vec::with_capacity(buffer_size);
        let mut prob_buffer = Vec::with_capacity(buffer_size);

        let call = |likelihoods: &[Vec<GenotypeLikelihoods>]| {
            likelihoods.iter().map(|lh| query.call(lh)).collect_vec()
        };

        let mut processed = 0;
        loop {
            let site_buffer = site_channel_out.recv().ok().expect("Error using channel.");
            match site_buffer {
                None => {
                    return candidates;
                },
                Some(mut site_buffer) => {
                    processed += site_buffer.len();

                    likelihood_buffer.extend(site_buffer.iter_mut().map(|site|
                        site.genotype_likelihoods().ok().expect("Error reading genotype likelihoods.")
                    ));

                    for probs in unsafe { pool.map(likelihood_buffer.chunks(buffer_size / threads), &call) } {
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
        }
    });
    let mut candidates = processor.join().ok().expect("Error joining thread.");
    if fdr.is_some() {
        control_fdr(&mut candidates, fdr.unwrap());
    }
    return candidates;
}


fn control_fdr(candidates: &mut Vec<(Site, LogProb)>, fdr: LogProb) {
    let cmp = |a: &LogProb, b: &LogProb| a.partial_cmp(b).expect("Bug: NaN probability found.");
    let mut probs = candidates.iter().map(|&(_, prob)| prob).collect_vec();
    probs.sort_by(&cmp);
    let exp_fdr = logprobs::log_prob_cumsum(&probs);
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
