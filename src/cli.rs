use std::convert::AsRef;
use std::path::Path;
use std::process;
use std::fs;

use tempdir;
use itertools::Itertools;
use htslib::bcf;
use simple_parallel;
use bio;

use Prob;
use call;
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
            fs::create_dir(tmp.path().join(seq)).ok().expect("Error creating temporary directory.");
            let fifo_mpileup = tmp.path().join(seq).join("mpileup").with_extension("bcf");
            mkfifo(fifo_mpileup.as_path());
            let fifo_ann = tmp.path().join(seq).join("annotate").with_extension("bcf");
            mkfifo(fifo_ann.as_path());

            let mpileup = process::Command::new("samtools").arg("mpileup")
                                             .arg("-r").arg(seq)
                                             .arg("-f").arg(fasta.as_ref())
                                             .arg("-g").arg("-u")
                                             .arg("-o").arg(&fifo_mpileup)
                                             .args(&bams.iter().map(|bam| bam.as_ref()).collect_vec())
                                             .spawn().ok().expect("Failed to execute samtools mpileup.");
            let ann = process::Command::new("bcftools").arg("annotate")
                                   .arg("-O").arg("u")
                                   .arg("-o").arg(&fifo_ann)
                                   .arg("--remove").arg("INFO/INDEL,INFO/IDV,INFO/IMF,INFO/I16,INFO/QS")
                                   .arg(fifo_mpileup)
                                   .spawn().ok().expect("Failed to execute bcftools annotate.");

            (mpileup, ann, fifo_ann)
        };

        let mut writer = vec![]; // ugly hack, try Option once it works here.
        for (mut mpileup, mut ann, fifo) in pool.map(seqs.iter(), &mpileup) {
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
            if !mpileup.wait().ok().expect("Error retrieving exit status.").success() {
                panic!("Error during execution of samtools mpileup (see above).");
            }
            if !ann.wait().ok().expect("Error retrieving exit status.").success() {
                panic!("Error during execution of bcftools annotate (see above).");
            }
        }
    }
    tmp.close().ok().expect("Error removing FIFOs.");
}


pub fn merge<P: AsRef<Path>>(bcfs: &[P]) {
    let tmp = tempdir::TempDir::new("alpaca").ok().expect("Cannot create temp dir");
    let fifo_merge = tmp.path().join("merge").with_extension("bcf");
    mkfifo(fifo_merge.as_path());

    let mut merge = process::Command::new("bcftools").arg("merge")
                                   .arg("-O").arg("u")
                                   .arg("-o").arg(&fifo_merge)
                                   .arg("-m").arg("all")
                                   //.arg("--info-rules").arg("DP:join,VDB:join,RPB:join,MQB:join,BQB:join,MQSB:join,SGB:join,MQ0F:join")
                                   .args(&bcfs.iter().map(|bcf| bcf.as_ref()).collect_vec())
                                   .spawn().ok().expect("Failed to execute bcftools merge.");
    let view = process::Command::new("bcftools").arg("view")
                                   .arg("-O").arg("b")
                                   .arg("-e").arg("MAX(PL[0]) == 0")
                                   .arg(fifo_merge)
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

    header.push_record(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    let mut outbcf = bcf::Writer::new(&"-", &header, false, false);

    let mut calls = call::call(&mut inbcf, query_caller, fdr, threads);

    for (mut site, prob) in calls.drain() {
        outbcf.translate(&mut site.record);
        outbcf.subset(&mut site.record);

        site.update_record(prob);

        //site.record.trim_alleles().ok().expect("Error trimming alleles.");
        outbcf.write(&site.record).ok().expect("Error writing calls.");
    }
}
