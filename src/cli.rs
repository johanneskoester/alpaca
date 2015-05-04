use std::convert::AsRef;
use std::path::Path;
use std::process;
use std::fs;
use std::io::Write;
use std::cmp;
use std::str;

use tempdir;
use itertools::Itertools;
use htslib::bcf;
use simple_parallel;
use bio;

use LogProb;
use Prob;
use call;
use query;
use utils;


fn mkfifo<'a>(path: &Path) {
    process::Command::new("mkfifo").arg(path).status().ok().expect("Failed to create FIFO.");
}


fn seqs<P: AsRef<Path>>(fasta: &P) -> Vec<bio::io::fasta::Sequence> {
    let fasta_index = bio::io::fasta::Index::with_fasta_file(fasta)
        .ok().expect("Error reading FAI index: FASTA file must be indexed with samtools index.");
    fasta_index.sequences()
}


pub fn preprocess<P: AsRef<Path> + Sync>(fasta: &P, bams: &[P], threads: usize) {
    let seqs = seqs(fasta);

    let mpileup = |tmp: &Path, region: &str| {
        fs::create_dir(tmp.join(region)).ok().expect("Error creating temporary directory.");
        let fifo_mpileup = tmp.join(region).join("mpileup").with_extension("bcf");
        mkfifo(fifo_mpileup.as_path());
        let fifo_ann = tmp.join(region).join("annotate").with_extension("bcf");
        mkfifo(fifo_ann.as_path());

        let mpileup = process::Command::new("samtools")
            .arg("mpileup")
            .arg("-r").arg(region)
            .arg("-f").arg(fasta.as_ref())
            .arg("-g").arg("-u")
            .arg("-t").arg("DP")
            .arg("-o").arg(&fifo_mpileup)
            .args(&bams.iter().map(|bam| bam.as_ref()).collect_vec())
            .spawn().ok().expect("Failed to execute samtools mpileup.");
        let ann = process::Command::new("bcftools")
            .arg("annotate")
            .arg("-O").arg("b")
            .arg("-o").arg(&fifo_ann)
            .arg("--remove").arg("INFO/INDEL,INFO/IDV,INFO/IMF,INFO/I16,INFO/QS,INFO/DP")
            .arg(fifo_mpileup)
            .spawn().ok().expect("Failed to execute bcftools annotate.");

        (bcf::Reader::new(&fifo_ann), vec![mpileup, ann])
    };

    mapreduce(&seqs, cmp::max(1, threads - threads / 4), mpileup);
}


fn mapreduce<F: Sync>(seqs: &[bio::io::fasta::Sequence], threads: usize, kernel: F) where F: Fn(&Path, &str) -> (bcf::Reader, Vec<process::Child>) {
    let tmp = tempdir::TempDir::new("alpaca").ok().expect("Cannot create temp dir");
    {
        let mut pool = simple_parallel::Pool::new(threads);

        let apply = |region: &String| kernel(tmp.path(), region);

        // create list of regions
        let mut regions: Vec<String> = vec![];
        for seq in seqs {
            for offset in (1..seq.len + 1).step_by(100000) {
                regions.push(format!("{}:{}-{}", str::from_utf8(&seq.name).ok().expect("Invalid sequence name."), offset, offset + 100000));
            }
        }

        let mut readers = pool.map(regions.iter(), &apply);

        let (mut reader, mut procs) = readers.next().unwrap();
        let mut writer = {
            let header = bcf::Header::with_template(&reader.header);
            bcf::Writer::new(&"-", &header, false, false)
        };

        let mut record = bcf::Record::new();
        loop {
            match reader.read(&mut record) {
                Ok(()) => {
                    writer.write(&record).ok().expect("Error writing BCF.");
                },
                Err(bcf::ReadError::Invalid) => panic!("Error reading BCF."),
                Err(bcf::ReadError::NoMoreRecord) => {
                    for mut p in procs {
                        if !p.wait().ok().expect("Failed to retrieve exit status.").success() {
                            panic!("Error executing child process (see above).");
                        }
                    }
                    match readers.next() {
                        Some((r, ps)) => { reader = r; procs = ps },
                        None    => break
                    }                    
                }
            }
        }
    }
    tmp.close().ok().expect("Error removing FIFOs.");
}


pub fn merge<P: AsRef<Path> + Sync>(fasta: &P, bcfs: &[P], threads: usize) {
    let seqs = seqs(fasta);

    let merge = |tmp: &Path, region: &str| {
        fs::create_dir(tmp.join(region)).ok().expect("Error creating temporary directory.");
        let fifo_merge = tmp.join(region).join("mpileup").with_extension("bcf");
        mkfifo(fifo_merge.as_path());
        let fifo_filter = tmp.join(region).join("annotate").with_extension("bcf");
        mkfifo(fifo_filter.as_path());

        let merge = process::Command::new("bcftools")
            .arg("merge")
            .arg("-r").arg(region)
            .arg("-O").arg("u")
            .arg("-o").arg(&fifo_merge)
            .arg("-m").arg("all")
            .args(&bcfs.iter().map(|bcf| bcf.as_ref()).collect_vec())
            .spawn().ok().expect("Failed to execute bcftools merge.");

        let filter = filter_cmd()
            .arg("-o").arg(&fifo_filter)
            .arg(&fifo_merge)
            .spawn().ok().expect("Failed to execute bcftools view.");

        (bcf::Reader::new(&fifo_filter), vec![merge, filter])
    };
    
    mapreduce(&seqs, threads, merge);
}


fn filter_cmd() -> process::Command {
    let mut cmd = process::Command::new("bcftools");
    cmd.arg("view")
       .arg("-O").arg("b")
       .arg("-e").arg("MAX(PL[0]) == 0");

    cmd
}


pub fn filter() {
    filter_cmd()
       .arg("-")
       .status().ok().expect("Failed to execute bcftools view.");
}


pub fn call(query: &str, fdr: Option<LogProb>, max_prob: Option<LogProb>, heterozygosity: Prob, threads: usize) {
    let mut inbcf = bcf::Reader::new(&"-");
    let sample_idx = query::sample_index(&inbcf);
    let (query_caller, samples) = query::parse(query, &sample_idx, heterozygosity);

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
        header.push_record(format!("##alpaca_min_qual={}", max_prob.unwrap() * utils::LOG_TO_PHRED_FACTOR).as_bytes());
    }    
    header.push_record(format!("##alpaca_heterozygosity={}", heterozygosity).as_bytes());

    let mut outbcf = bcf::Writer::new(&"-", &header, false, false);

    // perform the calling
    let mut calls = call::call(&mut inbcf, query_caller, fdr, max_prob, threads);

    for (mut site, prob) in calls.drain(..) {
        outbcf.translate(&mut site.record);
        outbcf.subset(&mut site.record);

        site.update_record(prob);

        // TODO trim alleles causes a segfault
        //site.record.trim_alleles().ok().expect("Error trimming alleles.");
        outbcf.write(&site.record).ok().expect("Error writing calls.");
    }
}
