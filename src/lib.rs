#![feature(plugin)]
#![plugin(peg_syntax_ext)]
#![feature(step_by)]
#![feature(vec_push_all)]
#![feature(drain)]
#![feature(iter_arith)]


extern crate itertools;
extern crate simple_parallel;
extern crate rust_htslib as htslib;
extern crate tempdir;
extern crate bio;
#[macro_use]
extern crate log;

pub type Prob = f64;
pub type LogProb = f64;
pub const EPSILON: Prob = 0.000001;

pub mod call;
pub mod utils;
pub mod query;
pub mod cli;
