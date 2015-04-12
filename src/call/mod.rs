use std::thread;
use std::sync::mpsc::channel;
use std::sync::{Arc, Mutex};
use std::iter;

use Prob;

pub mod diff;
pub mod relaxed_intersection;
pub mod sample_union;
pub mod union;


pub fn query_probabilities<G: GenotypeLikelihoods, S: Site<G>, I: Iterator<Item=S>, C: Caller>(query: C, sites: I, threads: usize) -> Vec<Prob> {
    let (_, exp_site_count) = sites.size_hint();

    let probs = Arc::new(Mutex::new(Vec::with_capacity(exp_site_count.unwrap_or(1000))));
    let (sites_in, sites_out) = channel();

    for _ in 0..threads {
        let sites = sites_out.clone();
        let probs = probs.clone();
        thread::scoped(|| {
            let (i, site) = sites.recv().unwrap();
            let p = query.call(site.genotype_likelihoods());
            let mut probs = probs.lock().unwrap();
            // extend the vector such that i fits in
            probs.extend(iter::repeat(0).take(i + 1 - probs.len()));
            probs[i] = p;
        });
    }

    for (i, site) in sites.enumerate() {
        sites_in.send((i, site));
    }

    *probs.lock().unwrap()
}


pub trait Caller {
    fn call(&self, likelihoods: &[GenotypeLikelihoods]) -> Prob;
}


pub trait Site<G: GenotypeLikelihoods> {
    fn genotype_likelihoods(&self) -> Vec<G>;
}


pub trait GenotypeLikelihoods {
    fn with_allelefreq(&self, m: usize) -> Vec<Prob>;
}
