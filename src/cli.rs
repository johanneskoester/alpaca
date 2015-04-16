use std::convert::AsRef;
use std::path::Path;
use std::process;
use std;
use std::fs;

use tempdir;
use itertools::Itertools;
use htslib::bcf;
use simple_parallel;
use bio;

use Prob;
use call;
use io;
use query;


pub fn preprocess<P: AsRef<Path> + Sync>(fasta: &P, bams: &[P], threads: usize) {
    let fasta_index = bio::io::fasta::Index::with_fasta_file(fasta)
        .ok().expect("Error reading FAI index: FASTA file must be indexed with samtools index.");
    let refs = fasta_index.sequences().into_iter().map(|seq| String::from_utf8(seq.name)
                                      .ok().expect("Invalid reference name.")).collect_vec();

    let bampaths = bams.iter().map(|bam| bam.as_ref().to_str().unwrap()).collect_vec();

    let tmp = tempdir::TempDir::new("alpaca").ok().expect("Failed to create temporary directory.");
    let fifos = (0..refs.len()).map(|i| tmp.path().join(format!("{}.bcf", i)).to_str().unwrap().to_string()).collect_vec();

    for fifo in fifos.iter() {
        if !process::Command::new("mkfifo").arg(fifo).status().ok().expect("Failed to execute mkfifo.").success() {
            panic!("Failed to create fifo.");
        }
    }
    
    let mut pool = simple_parallel::Pool::new(threads);
    let mpileup = |(refname, fifo)| {
        process::Command::new("sh").arg("-c")
            .arg(format!(
                "samtools mpileup -r {} -g -f {} {} > {fifo}",
                refname,
                fasta.as_ref().to_str().expect("Invalid FASTA path."),
                &bampaths.connect(" "),
                fifo=fifo,
            ))
            .status().ok().expect("Failed to execute samtools.")
    };

    let mut concat = process::Command::new("sh").arg("-c").arg(format!(
        "bcftools concat {}",
        fifos.connect(" "),
    )).spawn().ok().expect("Failed to execute bcftools");

    for status in pool.map(refs.into_iter().zip(fifos), &mpileup) {
        if !status.success() {
            panic!("Error during execution of samtools mpileup (see above).");
        }
    }
    concat.wait().ok().expect("Error during execution of bcftools (see above).").success();

    tmp.close().ok().expect("Failed to remove temporary directory");
}


pub fn merge<P: AsRef<Path>>(bcfs: &[P]) {
    let _bcfs: Vec<&str> = bcfs.iter().map(|bcf| bcf.as_ref().to_str().unwrap()).collect();
    let status = process::Command::new("sh").arg("-c").arg(
        format!(
            "bcftools merge -O b {} | bcftools view -O b -e 'PL[0] == 0'",
            &_bcfs.connect(" ")
        )
    ).status().ok().expect("Failed to execute bcftools.");

    if !status.success() {
        panic!("Error during execution of bcftools (see above).");
    }
}


pub fn call(query: &str, fdr: Prob, threads: usize, heterozygosity: Prob) {
    let mut inbcf = bcf::Reader::new(&"-");
    let sample_idx = query::sample_index(&inbcf);
    let query_caller = query::parse(query, &sample_idx, heterozygosity);

    let mut header = bcf::Header::with_template(&inbcf.header);
    header.push_record(b"##FORMAT=<ID=GT,Number=2,Type=Integer,Description=\"Genotype\">"); // TODO generalize ploidy
    let mut outbcf = bcf::Writer::new(&"-", &header);

    let calls = call::call(&mut inbcf, query_caller, fdr, threads);
    io::write(&mut outbcf, calls);
}
