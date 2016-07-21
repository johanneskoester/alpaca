use std::convert::AsRef;
use std::path::{Path, PathBuf};
use std::process;
use std::fs;
use std::io::Write;
use std::str;

use tempdir;
use itertools::Itertools;
use htslib::bcf;
use simple_parallel;
use bio;
use bio::stats::logprobs;

use LogProb;
use Prob;
use call;
use query;


pub fn filter() {
    filter_cmd()
       .arg("-")
       .status().ok().expect("Failed to execute bcftools view.");
}


pub fn call(
    query: &str,
    fdr: Option<LogProb>,
    max_prob: Option<LogProb>,
    heterozygosity: Prob,
    dependency: Prob,
    threads: usize
) {
    let mut inbcf = bcf::Reader::new(&"-");
    let parser = query::Parser::new(&inbcf.header.samples(), heterozygosity, dependency);
    let (query_caller, samples) = parser.parse(query);

    // create writer
    let mut header = if samples.len() == inbcf.header.sample_count() as usize {
        bcf::Header::with_template(&inbcf.header)
    }
    else {
        bcf::Header::subset_template(
            &inbcf.header, &samples.iter().map(|s| s.as_bytes()).collect_vec()
        ).ok().expect("Unknown sample name.")
    };

    // update header
    header.push_record(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    header.push_record(format!("##alpaca_query={}", query).as_bytes());
    if fdr.is_some() {
        header.push_record(format!("##alpaca_fdr={}", fdr.unwrap()).as_bytes());
    }
    if max_prob.is_some() {
        header.push_record(
            format!(
                "##alpaca_min_qual={}",
                logprobs::log_to_phred(max_prob.unwrap())
            ).as_bytes()
        );
    }
    header.push_record(format!("##alpaca_heterozygosity={}", heterozygosity).as_bytes());

    let mut outbcf = bcf::Writer::new(&"-", &header, false, false);

    // perform the calling
    let mut calls = call::call(&mut inbcf, query_caller, fdr, max_prob, threads);

    for (mut site, prob) in calls.drain(..) {
        outbcf.translate(&mut site.record);
        outbcf.subset(&mut site.record);

        site.set_qual(prob);
        site.calc_genotype();

        // TODO trim alleles causes a segfault
        //site.record.trim_alleles().ok().expect("Error trimming alleles.");
        outbcf.write(&site.record).ok().expect("Error writing calls.");
    }
}
