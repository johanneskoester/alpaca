#![feature(core)]
#![feature(plugin)]
#![plugin(peg_syntax_ext)]
#![feature(collections)]
#![feature(str_char)]

extern crate itertools;
extern crate simple_parallel;
extern crate rust_htslib as htslib;

pub type Prob = f64;

pub mod call;
pub mod utils;
pub mod query;
pub mod io;
pub mod cli;


fn main() {
    println!("Hello, world!");
}
