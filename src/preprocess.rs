use std::convert::AsRef;
use std::path::Path;
use std::process;
use std::io;


pub fn preprocess<P: AsRef<Path>>(bams: &[P], fasta: P) -> io::Result<process::ExitStatus> {
    process::Command::new("samtools").arg("mpileup")
                                     .arg("-g")
                                     .arg("--fasta-ref").arg(fasta)
                                     .args(bams).status()
}


pub fn merge<P: AsRef<Path>>(bcfs: &[P]) -> io::Result<process::ExitStatus> {
    process::Command::new("sh").arg("-c").arg(
        format!(
            "bcftools merge -O b {} | bcftools view -O b -e 'PL[0] == 0'",
            bcfs.connect(" ")
        )
    ).status()
}
