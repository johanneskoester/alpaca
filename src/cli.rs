use std::convert::AsRef;
use std::path::Path;
use std::process;
use std::fs;
use std::io::Write;
use std::cmp;

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


fn seqs<P: AsRef<Path>>(fasta: &P) -> Vec<String> {
    let fasta_index = bio::io::fasta::Index::with_fasta_file(fasta)
        .ok().expect("Error reading FAI index: FASTA file must be indexed with samtools index.");
    fasta_index.sequences().into_iter().map(|seq| String::from_utf8(seq.name)
                                       .ok().expect("Invalid reference name.")).collect()
}


pub fn preprocess<P: AsRef<Path> + Sync>(fasta: &P, bams: &[P], threads: usize) {
    let seqs = seqs(fasta);

    let mpileup = |tmp: &Path, seq: &str| {
        fs::create_dir(tmp.join(seq)).ok().expect("Error creating temporary directory.");
        let fifo_mpileup = tmp.join(seq).join("mpileup").with_extension("bcf");
        mkfifo(fifo_mpileup.as_path());
        let fifo_ann = tmp.join(seq).join("annotate").with_extension("bcf");
        mkfifo(fifo_ann.as_path());

        let mut mpileup = process::Command::new("samtools")
            .arg("mpileup")
            .arg("-r").arg(seq)
            .arg("-f").arg(fasta.as_ref())
            .arg("-g").arg("-u")
            .arg("-o").arg(&fifo_mpileup)
            .args(&bams.iter().map(|bam| bam.as_ref()).collect_vec())
            .spawn().ok().expect("Failed to execute samtools mpileup.");
        let mut ann = process::Command::new("bcftools")
            .arg("annotate")
            .arg("-O").arg("b")
            .arg("-o").arg(&fifo_ann)
            .arg("--remove").arg("INFO/INDEL,INFO/IDV,INFO/IMF,INFO/I16,INFO/QS")
            .arg(fifo_mpileup)
            .spawn().ok().expect("Failed to execute bcftools annotate.");

        let reader = bcf::Reader::new(&fifo_ann);
        let records = reader.records().map(|rec| rec.ok().expect("Invalid BCF record")).collect_vec();

        if !mpileup.wait().ok().expect("Error retrieving exit status.").success() {
            panic!("Error during execution of samtools mpileup (see above).");
        }
        if !ann.wait().ok().expect("Error retrieving exit status.").success() {
            panic!("Error during execution of bcftools annotate (see above).");
        }

        (reader, records)
    };

    mapreduce(&seqs, threads, mpileup);
}


fn mapreduce<F: Sync>(seqs: &[String], threads: usize, kernel: F) where F: Fn(&Path, &str) -> (bcf::Reader, Vec<bcf::Record>) {
    let tmp = tempdir::TempDir::new("alpaca").ok().expect("Cannot create temp dir");
    {
        let mut pool = simple_parallel::Pool::new(cmp::max(1, threads - threads / 4));

        let apply = |seq: &String| kernel(tmp.path(), seq);

        let mut readers = pool.map(seqs.iter(), &apply);

        let (reader, mut records) = readers.next().unwrap();
        let mut writer = {
            let header = bcf::Header::with_template(&reader.header);
            bcf::Writer::new(&"-", &header, false, false)
        };

        loop {
            for record in records.drain(..) {
                writer.write(&record).ok().expect("Error writing BCF.")
            }
            match readers.next() {
                Some(v) => { records = v.1 },
                None    => break
            }
        }
    }
    tmp.close().ok().expect("Error removing FIFOs.");
}


pub fn merge<P: AsRef<Path> + Sync>(fasta: &P, bcfs: &[P], threads: usize) {
    let seqs = seqs(fasta);

    let merge = |tmp: &Path, seq: &str| {
        fs::create_dir(tmp.join(seq)).ok().expect("Error creating temporary directory.");
        let fifo_merge = tmp.join(seq).join("mpileup").with_extension("bcf");
        mkfifo(fifo_merge.as_path());
        let fifo_filter = tmp.join(seq).join("annotate").with_extension("bcf");
        mkfifo(fifo_filter.as_path());

        let mut merge = process::Command::new("bcftools")
            .arg("merge")
            .arg("-r").arg(seq)
            .arg("-O").arg("u")
            .arg("-o").arg(&fifo_merge)
            .arg("-m").arg("all")
            .args(&bcfs.iter().map(|bcf| bcf.as_ref()).collect_vec())
            .spawn().ok().expect("Failed to execute bcftools merge.");

        let mut filter = filter_cmd()
            .arg("-o").arg(&fifo_filter)
            .arg(&fifo_merge)
            .spawn().ok().expect("Failed to execute bcftools view.");

        let reader = bcf::Reader::new(&fifo_filter);
        let records = reader.records().map(|rec| rec.ok().expect("Invalid BCF record")).collect_vec();

        if !merge.wait().ok().expect("Error retrieving exit status.").success() {
            panic!("Error during execution of bcftools merge (see above).");
        }
        if !filter.wait().ok().expect("Error retrieving exit status.").success() {
            panic!("Error during execution of bcftools view (see above).");
        }

        (reader, records)
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
