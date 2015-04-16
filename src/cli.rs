use std::convert::AsRef;
use std::path::Path;
use std::process;

use tempdir;
use itertools::Itertools;
use htslib::bcf;
use simple_parallel;
use bio;

use Prob;
use call;
use io;
use query;


fn mkfifo<'a>(path: &Path) {
    process::Command::new("mkfifo").arg(path).status().ok().expect("Failed to create FIFO.");
}


pub fn preprocess<P: AsRef<Path> + Sync>(fasta: &P, bams: &[P], threads: usize) {
    let fasta_index = bio::io::fasta::Index::with_fasta_file(fasta)
        .ok().expect("Error reading FAI index: FASTA file must be indexed with samtools index.");
    let seqs = fasta_index.sequences().into_iter().map(|seq| String::from_utf8(seq.name)
                                      .ok().expect("Invalid reference name.")).collect_vec();

    let tmp = tempdir::TempDir::new("alpaca").ok().expect("Cannot create temp dir");
    {
        let mut pool = simple_parallel::Pool::new(threads);
        let mpileup = |seq: &String| {
            let fifo = tmp.path().join(seq).with_extension("bcf");
            mkfifo(fifo.as_path());
            (process::Command::new("samtools").arg("mpileup")
                                             .arg("-r").arg(seq).arg("-f").arg(fasta.as_ref())
                                             .arg("-g").arg("-u").arg("-o").arg(fifo.as_path())
                                             .args(&bams.iter().map(|bam| bam.as_ref()).collect_vec())
                                             .spawn().ok().expect("Failed to execute samtools mpileup."), fifo)
        };

        let mut writer = vec![]; // ugly hack, try Option once it works here.
        for (process, fifo) in pool.map(seqs.iter(), &mpileup) {
            let reader = bcf::Reader::new(&fifo);
            
            if writer.is_empty() {
                let header = bcf::Header::with_template(&reader.header);
                writer.push(bcf::Writer::new(&"-", &header, false, false));
            }

            let mut record = bcf::Record::new();
            loop {
                match reader.read(&mut record) {
                    Ok(()) => writer[0].write(&record).ok().expect("Error writing BCF."),
                    Err(bcf::ReadError::Invalid) => panic!("Invalid BCF record."),
                    Err(bcf::ReadError::NoMoreRecord) => break
                }
            }
            //if !process.wait().ok().expect("Error retrieving exit status.").success() {
            //    panic!("Error during execution of samtools mpileup (see above).");
            //}
        }
    }
    tmp.close().ok().expect("Error removing FIFOs.");
}


pub fn merge<P: AsRef<Path>>(bcfs: &[P]) {
    let tmp = tempdir::TempDir::new("alpaca").ok().expect("Cannot create temp dir");
    let fifo = tmp.path().join("merge").with_extension("bcf");
    mkfifo(fifo.as_path());

    let mut merge = process::Command::new("bcftools").arg("merge")
                                   .arg("-O").arg("b")
                                   .arg("-o").arg(&fifo)
                                   .args(&bcfs.iter().map(|bcf| bcf.as_ref()).collect_vec())
                                   .spawn().ok().expect("Failed to execute bcftools merge.");
    let view = process::Command::new("bcftools").arg("view")
                                   .arg("-O").arg("b")
                                   .arg("-e").arg("PL[0] == 0")
                                   .arg(fifo)
                                   .status().ok().expect("Failed to execute bcftools view.");

    if !merge.wait().ok().expect("Error retrieving exit status.").success() {
        panic!("Error during execution of bcftools merge (see above).");
    }
    if !view.success() {
        panic!("Error during execution of bcftools view (see above).");
    }
    tmp.close().ok().expect("Error removing FIFO.");
}


pub fn call(query: &str, fdr: Prob, threads: usize, heterozygosity: Prob) {
    let mut inbcf = bcf::Reader::new(&"-");
    let sample_idx = query::sample_index(&inbcf);
    let query_caller = query::parse(query, &sample_idx, heterozygosity);

    let mut header = bcf::Header::with_template(&inbcf.header);
    header.push_record(b"##FORMAT=<ID=GT,Number=2,Type=Integer,Description=\"Genotype\">"); // TODO generalize ploidy
    let mut outbcf = bcf::Writer::new(&"-", &header, false, false);

    let calls = call::call(&mut inbcf, query_caller, fdr, threads);
    io::write(&mut outbcf, calls);
}
