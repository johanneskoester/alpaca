extern crate alpaca;

#[macro_use]
extern crate log;
extern crate fern;
#[macro_use]
extern crate clap;
extern crate bio;

use clap::{App,AppSettings};
use bio::stats::logprobs;

use alpaca::cli;



fn main() {
    let yaml = load_yaml!("cli.yaml");
    let matches = App::from_yaml(yaml)
                      .version(env!("CARGO_PKG_VERSION"))
                      .global_settings(&[AppSettings::SubcommandRequired,
                                         AppSettings::ColoredHelp])
                      .get_matches();

    let logger_config = fern::DispatchConfig {
        format: Box::new(|msg: &str, level: &log::LogLevel, _: &log::LogLocation| {
            match level {
                &log::LogLevel::Debug => format!("DEBUG: {}", msg),
                _ => msg.to_owned()
            }
        }),
        output: vec![fern::OutputConfig::stderr()],
        level: log::LogLevelFilter::Debug,
    };

    if let Err(e) = fern::init_global_logger(
        logger_config,
        if matches.is_present("verbose") { log::LogLevelFilter::Debug } else { log::LogLevelFilter::Info }
    ) {
        panic!("Failed to initialize logger: {}", e);
    }

    if let Some(matches) = matches.subcommand_matches("call") {
        let query = matches.value_of("query").unwrap();

        let min_qual = value_t!(matches, "min-qual", f64).ok();

        let mut fdr = value_t!(matches, "fdr", f64).ok();
        if min_qual.is_none() && fdr.is_none() {
            fdr = Some(0.05);
        }
        let fdr = fdr.map(|fdr| fdr.ln());

        let max_prob = min_qual.map(logprobs::phred_to_log);

        let dependency = value_t!(matches, "dependency", f64).unwrap_or(0.0);
        if dependency > 1.0 || dependency < 0.0 {
            panic!("Dependency has to be between 0.0 and 1.0.");
        }

        let heterozygosity = value_t!(matches, "heterozygosity", f64).unwrap_or(0.001);
        if heterozygosity < 0.0 || heterozygosity > 1.0 {
            panic!("Heterozygosity has to be between 0.0 and 1.0.");
        }

        let threads = value_t!(matches, "threads", usize).unwrap_or(1);

        cli::call(&query, fdr, max_prob, heterozygosity, dependency, threads);
    }
}
