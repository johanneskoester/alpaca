extern crate argparse;
extern crate alpaca;

use std::str::FromStr;
use std::io::{stdout, stderr};

use argparse::{ArgumentParser, Store, List, StoreTrue};

#[macro_use]
extern crate log;
extern crate fern;

use alpaca::cli;
use alpaca::utils;


#[derive(Debug)]
enum Command {
    Preprocess,
    Merge,
    Filter,
    Call,
    None,
}


impl FromStr for Command {
    type Err = ();

    fn from_str(src: &str) -> Result<Command, ()> {
        return match src {
            "preprocess" => Ok(Command::Preprocess),
            "merge"      => Ok(Command::Merge),
            "filter"     => Ok(Command::Filter),
            "call"       => Ok(Command::Call),
            _            => Err(()),
        };
    }
}


enum OptionalArg<T> {
    Some(T),
    None
}


impl<T> OptionalArg<T> {
    fn into_option(self) -> Option<T> {
        match self {
            OptionalArg::Some(v) => Some(v),
            OptionalArg::None    => None
        }
    }

    fn is_none(&self) -> bool {
        match self {
            &OptionalArg::Some(_) => false,
            &OptionalArg::None    => true
        }
    }
}


impl<T: FromStr> FromStr for OptionalArg<T> {
    type Err = ();

    fn from_str(src: &str) -> Result<OptionalArg<T>, ()> {
        match T::from_str(src) {
            Ok(v) => Ok(OptionalArg::Some(v)),
            _     => Err(())
        }
    }
}


fn main() {
    let logger_config = fern::DispatchConfig {
        format: Box::new(|msg: &str, _: &log::LogLevel, _: &log::LogLocation| {
            msg.to_owned()
        }),
        output: vec![fern::OutputConfig::stderr()],
        level: log::LogLevelFilter::Debug,
    };
    if let Err(e) = fern::init_global_logger(logger_config, log::LogLevelFilter::Debug) {
        panic!("Failed to initialize global logger: {}", e);
    }

    let mut subcommand = Command::None;
    let mut args = vec![];
    {
        let mut ap = ArgumentParser::new();
        ap.set_description(
"ALPACA is a caller for genetic variants from next-generation sequencing data."
        );
        ap.refer(&mut subcommand).required()
            .add_argument("command", Store, "Command to run (preprocess, merge, filter or call)");
        ap.refer(&mut args)
            .add_argument("arguments", List, "Arguments for command");
        ap.stop_on_first_argument(true);
        ap.parse_args_or_exit();
    }

    args.insert(0, format!("subcommand {:?}", subcommand));
    match subcommand {
        Command::Preprocess => preprocess(args),
        Command::Merge      => merge(args),
        Command::Filter     => filter(args),
        Command::Call       => call(args),
        Command::None       => {
            error!("Unknown subcommand.");
            std::process::exit(1);
        }
    }
}


fn preprocess(args: Vec<String>) {
    let mut fasta = "".to_string();
    let mut bams: Vec<String> = vec![];
    let mut threads = 1;
    let mut nobaq = false;

    {
        let mut ap = ArgumentParser::new();
        ap.set_description(
r#"Preprocess loci of one or more samples for calling with ALPACA.
This command takes a FASTA file with the reference genome, one or more
indexed BAM files, calculates genotype likelihoods for each locus
using samtools mpileup, and prints the result to STDOUT in BCF format.
Example: "alpaca preprocess A.bam > A.bcf""#
        );


        ap.refer(&mut fasta)
          .add_argument("fasta", Store, "FASTA file with reference genome.");
        ap.refer(&mut bams)
          .add_argument("bam", List, "BAM files to preprocess.");
        ap.refer(&mut threads)
          .add_option(&["--threads", "-t"], Store, "Number of threads to use.");
        ap.refer(&mut nobaq)
          .add_option(&["--no-BAQ", "-B"], StoreTrue, "Disable BAQ (per-Base Alignment Quality adjustment)");
        parse_args_or_exit(&ap, args);
    }
    cli::preprocess(&fasta, &bams, threads, nobaq);
}


fn merge(args: Vec<String>) {
    let mut fasta = "".to_string();
    let mut bcfs: Vec<String> = vec![];
    let mut threads = 1;

    {
        let mut ap = ArgumentParser::new();
        ap.set_description(
r#"This command merges preprocessed loci of one or more samples for calling
with ALPACA. This command takes two or more BCF files and
prints a merged BCF to STDOUT. For calling, this should be combined
with alpaca-filter, which removes irrelevant sites.
Example: "alpaca merge A.bcf B.bcf C.bcf | alpaca filter > merged.bcf"#
        );

        ap.refer(&mut fasta)
          .add_argument("fasta", Store, "FASTA file with reference genome.");
        ap.refer(&mut bcfs)
          .add_argument("bcf", List, "ALPACA-preprocessed BCF files to merge.");
        ap.refer(&mut threads)
          .add_option(&["--threads", "-t"], Store, "Number of threads to use.");
        parse_args_or_exit(&ap, args);
    }
    cli::merge(&fasta, &bcfs, threads);
}


fn filter(args: Vec<String>) {
    {
        let mut ap = ArgumentParser::new();
        ap.set_description(
r#"Filter preprocessed loci of one or more samples for calling
with ALPACA, removing irrelevant sites.
Example: "alpaca-merge A.bcf B.bcf C.bcf | alpaca-filter > filtered.bcf""#
        );
        parse_args_or_exit(&ap, args);
    }
    cli::filter();
}


fn call(args: Vec<String>) {
    let mut query = "".to_string();
    let mut fdr: OptionalArg<f64> = OptionalArg::None;
    let mut min_qual: OptionalArg<f64> = OptionalArg::None;
    let mut heterozygosity = 0.001;
    let mut dependency = false;
    let mut threads = 1;

    {
        let mut ap = ArgumentParser::new();
        ap.set_description(
r#"Call variants. The command takes the merged preprocessed
samples in BCF format from STDIN, and prints the found variants in
BCF format to STDOUT.
Example: "alpaca call --fdr 0.05 'A - (B + C)' < filtered.bcf > calls.bcf""#
        );
        ap.refer(&mut query)
          .add_argument("query", Store, "Algebraic query describing the filtering scenario.");
        ap.refer(&mut fdr)
          .add_option(&["--fdr"], Store, "Desired false discovery rate.");
        ap.refer(&mut min_qual)
          .add_option(&["--min-qual"], Store, "Minimum variant quality (i.e. PHRED scaled posterior probability for observing the reference genotype given the query).");
        ap.refer(&mut heterozygosity)
          .add_option(&["--het"], Store, "Expected heterozygosity (default 0.001).");
        ap.refer(&mut dependency)
          .add_option(&["--dep", "--dependent-samples"], StoreTrue, "Consider samples to be dependent.");
        ap.refer(&mut threads)
          .add_option(&["--threads", "-t"], Store, "Number of threads to use (default 1).");
        parse_args_or_exit(&ap, args);
    }

    if min_qual.is_none() && fdr.is_none() {
        fdr = OptionalArg::Some(0.05);
    }


    let max_prob = min_qual.into_option().map(|q| q * utils::PHRED_TO_LOG_FACTOR);

    cli::call(&query, fdr.into_option().map(|fdr| fdr.ln()), max_prob, heterozygosity, dependency, threads);
}


fn parse_args_or_exit(ap: &ArgumentParser, args: Vec<String>) {
    match ap.parse(args, &mut stdout(), &mut stderr()) {
        Ok(()) =>  {}
        Err(x) => {
            std::process::exit(x);
        }
    }
}
