//! The pdbiox command line.
//!
//! Argument parsing, output formatting and exit codes over the library. No
//! analysis happens here: anything worth doing is worth doing in the library,
//! where the Rust callers get it too.
//!
//! Results go to standard output and findings go to standard error, always, so
//! a pipeline is never contaminated by warnings.

#![forbid(unsafe_code)]

mod commands;
mod exit;
mod report;

use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;
use std::process::ExitCode;

/// A batteries-included structural bioinformatics engine.
#[derive(Parser, Debug)]
#[command(name = "pdbiox", version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Command,

    /// Output format.
    #[arg(long, global = true, value_enum, default_value_t = OutputFormat::Text)]
    format: OutputFormat,

    /// How much irregularity to tolerate while reading.
    #[arg(long, global = true, value_enum, default_value_t = Tolerance::Permissive)]
    tolerance: Tolerance,

    /// Suppress findings on standard error.
    #[arg(long, global = true)]
    quiet: bool,

    /// Disable colour. Also honoured through the `NO_COLOR` environment variable.
    #[arg(long, global = true)]
    no_color: bool,
}

/// How results are printed.
#[derive(Clone, Copy, PartialEq, Eq, Debug, ValueEnum)]
enum OutputFormat {
    /// Human-readable.
    Text,
    /// Structured, for another program to consume.
    Json,
}

/// How much irregularity a read tolerates.
#[derive(Clone, Copy, PartialEq, Eq, Debug, ValueEnum)]
enum Tolerance {
    /// Anything violating the specification fails the read.
    Strict,
    /// Produce a structure and report what was wrong with it.
    Permissive,
    /// Continue past local errors.
    Recover,
}

impl From<Tolerance> for pdbiox::ParseMode {
    fn from(tolerance: Tolerance) -> Self {
        match tolerance {
            Tolerance::Strict => Self::Strict,
            Tolerance::Permissive => Self::Permissive,
            Tolerance::Recover => Self::Recover,
        }
    }
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Summarise a structure.
    Info {
        /// The file to read.
        input: PathBuf,
        /// Break the summary down per chain.
        #[arg(long)]
        detail: bool,
    },
    /// Convert a structure from one format to another.
    Convert {
        /// The file to read.
        input: PathBuf,
        /// Where to write.
        output: PathBuf,
        /// Rename a chain on the way out, as `FROM=TO`. May be repeated.
        #[arg(long = "chain-map", value_name = "FROM=TO")]
        chain_map: Vec<String>,
        /// Write counts past the field limits in the extended scheme.
        #[arg(long)]
        hybrid36: bool,
    },
    /// Check a structure against the invariants it must satisfy.
    Validate {
        /// The file to read.
        input: PathBuf,
    },
    /// Measure a structure: extent, radius of gyration, centre.
    Measure {
        /// The file to read.
        input: PathBuf,
    },
    /// Compare two structures atom by atom, after fitting one onto the other.
    Rmsd {
        /// The structure to move.
        mobile: PathBuf,
        /// The structure to move it onto.
        reference: PathBuf,
        /// Report the deviation where the structures sit, without fitting.
        #[arg(long)]
        no_fit: bool,
    },
    /// Print the policy an analysis runs under by default.
    Policy,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let color = !cli.no_color && std::env::var_os("NO_COLOR").is_none();
    let context = report::Context {
        format: matches!(cli.format, OutputFormat::Json),
        quiet: cli.quiet,
        color,
        mode: cli.tolerance.into(),
    };

    let exit = match cli.command {
        Command::Info { input, detail } => commands::info(&input, detail, context),
        Command::Convert {
            input,
            output,
            chain_map,
            hybrid36,
        } => commands::convert(&input, &output, &chain_map, hybrid36, context),
        Command::Validate { input } => commands::validate(&input, context),
        Command::Measure { input } => commands::measure(&input, context),
        Command::Rmsd {
            mobile,
            reference,
            no_fit,
        } => commands::rmsd(&mobile, &reference, no_fit, context),
        Command::Policy => commands::policy(context),
    };
    // Every documented code fits a byte; anything else would be a defect here.
    match u8::try_from(exit.code()) {
        Ok(code) => ExitCode::from(code),
        Err(_) => ExitCode::FAILURE,
    }
}
