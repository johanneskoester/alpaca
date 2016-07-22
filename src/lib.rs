extern crate itertools;
extern crate cue;
extern crate rust_htslib as htslib;
extern crate bio;
#[macro_use]
extern crate log;
extern crate fern;

pub type Prob = f64;
pub type LogProb = f64;
pub const EPSILON: Prob = 0.000001;

pub mod call;
pub mod utils;
pub mod query;
pub mod cli;
