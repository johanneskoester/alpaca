#![feature(core)]


extern crate itertools;
extern crate threadpool;
extern crate rust_htslib as htslib;

pub type Prob = f64;

pub mod call;
pub mod utils;


fn main() {
    println!("Hello, world!");
}
