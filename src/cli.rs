use std::convert::AsRef;
use std::path::Path;
use std::ffi::os_str::OsStr;
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

pub fn preprocess<P: AsRef<Path> + AsRef<OsStr> + Sync>(fasta: &P, bams: &[P], threads: usize) {
    let fasta_index = bio::io::fasta::Index::with_fasta_file(fasta)
        .ok().expect("Error reading FAI index: FASTA file must be indexed with samtools index.");
    let seqs = fasta_index.sequences().into_iter().map(|seq| String::from_utf8(seq.name)
                                      .ok().expect("Invalid reference name.")).collect_vec();

    let mut pool = simple_parallel::Pool::new(threads);
    let mpileup = |seq| {
        process::Command::new("samtools").arg("mpileup")
                                         .arg("-r").arg(seq).arg("-f").arg(fasta)
                                         .args(bams)
                                         .stdout(process::Stdio::piped())
                                         .spawn().ok().expect("Failed to execute samtools.")
    };

    let mut writer = None;    
    for process in pool.map(seqs, &mpileup) {
        let mut reader = bcf::Reader::new(&"-");
        
        if writer.is_none() {
            let header = bcf::Header::with_template(reader.header);
            writer = Some(bcf::Writer::new(&"-", header));
        }
        let mut record = bcf::Record::new();
        reader.read(&mut record).ok().expect("Error reading BCF.");
        writer.write(&record).ok().expect("Error writing BCF.");
    }
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
