//! Global output, read and analysis-policy configuration.

use crate::report::{Context, OutputKind};
use crate::{AltlocArgument, AssemblyArgument, ModelArgument, NamespaceArgument};
use clap::{ArgAction, Args, ValueEnum};
use std::io::IsTerminal as _;
use std::path::{Path, PathBuf};

/// Options shared by every command.
#[derive(Args, Debug)]
pub(crate) struct GlobalOptions {
    /// Output format.
    #[arg(long, global = true, value_enum)]
    pub format: Option<OutputFormat>,
    /// Write command results to this path instead of standard output.
    #[arg(long = "result-output", id = "global_output", global = true)]
    pub output: Option<PathBuf>,
    /// Write the effective CLI provenance record to this JSON path.
    #[arg(long, global = true)]
    pub provenance: Option<PathBuf>,
    /// Maximum Rust workers shared by every operation; zero uses pool capacity.
    #[arg(long, global = true, default_value_t = 0)]
    pub workers: usize,
    /// Complete retained-memory budget for the native execution graph.
    #[arg(long, global = true, default_value_t = 100_000_000)]
    pub memory_budget: usize,
    /// Directory for deterministic, bounded native spill files.
    #[arg(long, global = true, requires = "spill_budget")]
    pub spill_directory: Option<PathBuf>,
    /// Maximum spill bytes; requires `--spill-directory`.
    #[arg(long, global = true, default_value_t = 0, requires = "spill_directory")]
    pub spill_budget: u64,
    #[command(flatten)]
    read_behavior: ReadBehavior,
    /// How much irregularity to tolerate while reading.
    #[arg(long, global = true, value_enum, default_value_t = Tolerance::Permissive)]
    pub tolerance: Tolerance,
    /// Load policy overrides from strict TOML or JSON.
    #[arg(long, global = true)]
    pub policy: Option<PathBuf>,
    /// Override assembly: asymmetric-unit, biological:ID or crystal:RADIUS.
    #[arg(long, global = true)]
    pub assembly: Option<AssemblyArgument>,
    /// Override model: first, index:N, all or ensemble.
    #[arg(long = "model", id = "policy_model", global = true)]
    pub model: Option<ModelArgument>,
    /// Override alternate-conformation policy.
    #[arg(long, global = true)]
    pub altloc: Option<AltlocArgument>,
    /// Override identifier namespace: label, auth or explicit.
    #[arg(long, global = true)]
    #[arg(value_enum)]
    pub identifiers: Option<NamespaceArgument>,
    #[command(flatten)]
    display_behavior: DisplayBehavior,
    /// Print the effective policy; repeat for future diagnostic detail.
    #[arg(short, long, global = true, action = ArgAction::Count)]
    pub verbose: u8,
}

#[derive(Args, Debug)]
struct ReadBehavior {
    /// Infer missing element symbols from atom names and record the choice.
    #[arg(long, global = true)]
    infer_elements: bool,
    /// Infer otherwise ambiguous residue boundaries from file order and record it.
    #[arg(long, global = true)]
    infer_residue_boundaries: bool,
}

#[derive(Args, Debug)]
struct DisplayBehavior {
    /// Suppress findings on standard error.
    #[arg(long, global = true)]
    quiet: bool,
    /// Disable colour. Also honoured through the `NO_COLOR` environment variable.
    #[arg(long, global = true)]
    no_color: bool,
}

impl GlobalOptions {
    pub(crate) fn context(self) -> Result<Context, molframe::PolicyConfigError> {
        let mut policy = molframe::AnalysisPolicy::default();
        let mut configured_format = None;
        let mut chemistry = molframe::ChemistryConfiguration::default();
        if let Some(path) = xdg_policy_path().filter(|path| path.is_file()) {
            apply_configuration(&path, &mut policy, &mut configured_format, &mut chemistry)?;
        }
        let local = Path::new("molframe.toml");
        if local.is_file() {
            apply_configuration(local, &mut policy, &mut configured_format, &mut chemistry)?;
        }
        if let Some(path) = self.policy {
            apply_configuration(&path, &mut policy, &mut configured_format, &mut chemistry)?;
        }
        if let Some(value) = self.assembly {
            policy.assembly = value.into();
        }
        if let Some(value) = self.model {
            policy.model = value.into();
        }
        if let Some(value) = self.altloc {
            policy.altloc = value.into();
        }
        if let Some(value) = self.identifiers {
            policy.identifiers = value.into();
        }
        if self.verbose > 0 {
            eprintln!("effective policy:\n{policy}");
        }
        let policy = Box::leak(Box::new(policy));
        let output = self
            .output
            .map(|path| Box::leak(Box::new(path)) as &'static Path);
        let provenance = self
            .provenance
            .map(|path| Box::leak(Box::new(path)) as &'static Path);
        let ccd = chemistry
            .ccd_cache
            .map(|path| Box::leak(Box::new(path)) as &'static Path);
        let ccd_version = chemistry
            .ccd_version
            .map(|version| Box::leak(version.into_boxed_str()) as &'static str);
        let format = match self.format {
            Some(format) => format,
            None => match configured_format {
                Some(format) => format,
                None => default_output_format(),
            },
        };
        let memory = molframe::core::MemoryBudget::new(self.memory_budget).map_err(|error| {
            molframe::PolicyConfigError::InvalidValue {
                field: "memory-budget",
                value: error.to_string(),
            }
        })?;
        let mut execution = molframe::core::ExecutionContext::builder()
            .memory_budget(memory)
            .scratch_policy(molframe::core::ScratchPolicy::new(
                self.memory_budget.min(8_000_000),
            ));
        if self.workers > 0 {
            execution = execution.worker_budget(self.workers);
        }
        if let Some(directory) = self.spill_directory {
            execution = execution.temp_storage_policy(
                molframe::core::TempStoragePolicy::directory(directory, self.spill_budget),
            );
        }
        let execution =
            execution
                .build()
                .map_err(|error| molframe::PolicyConfigError::InvalidValue {
                    field: "execution-context",
                    value: error.to_string(),
                })?;
        let execution = Box::leak(Box::new(execution));
        Ok(Context {
            format: format.into(),
            quiet: self.display_behavior.quiet,
            color: !self.display_behavior.no_color && std::env::var_os("NO_COLOR").is_none(),
            mode: self.tolerance.into(),
            policy,
            output,
            provenance,
            ccd,
            ccd_version,
            execution,
            missing_element_policy: if self.read_behavior.infer_elements {
                molframe::MissingElementPolicy::InferFromAtomName
            } else {
                molframe::MissingElementPolicy::PreserveUnknown
            },
            residue_boundary_policy: if self.read_behavior.infer_residue_boundaries {
                molframe::AmbiguousResidueBoundaryPolicy::InferFromFileOrder
            } else {
                molframe::AmbiguousResidueBoundaryPolicy::Reject
            },
        })
    }
}

fn apply_configuration(
    path: &Path,
    policy: &mut molframe::AnalysisPolicy,
    format: &mut Option<OutputFormat>,
    chemistry: &mut molframe::ChemistryConfiguration,
) -> Result<(), molframe::PolicyConfigError> {
    let configuration = molframe::read_configuration(path)?;
    *policy = configuration.policy.apply_to(policy.clone())?;
    if let Some(value) = configuration.output.format {
        *format = Some(parse_output_format(&value)?);
    }
    if let Some(value) = configuration.chem.ccd_cache {
        chemistry.ccd_cache = Some(value);
    }
    if let Some(value) = configuration.chem.ccd_version {
        chemistry.ccd_version = Some(value);
    }
    Ok(())
}

fn parse_output_format(value: &str) -> Result<OutputFormat, molframe::PolicyConfigError> {
    match value {
        "text" => Ok(OutputFormat::Text),
        "json" => Ok(OutputFormat::Json),
        "jsonl" => Ok(OutputFormat::Jsonl),
        "csv" => Ok(OutputFormat::Csv),
        "tsv" => Ok(OutputFormat::Tsv),
        "arrow" => Ok(OutputFormat::Arrow),
        "parquet" => Ok(OutputFormat::Parquet),
        _ => Err(molframe::PolicyConfigError::InvalidValue {
            field: "output.format",
            value: value.to_owned(),
        }),
    }
}

fn xdg_policy_path() -> Option<PathBuf> {
    std::env::var_os("XDG_CONFIG_HOME")
        .map(PathBuf::from)
        .map(|root| root.join("molframe").join("config.toml"))
}

/// How results are printed.
#[derive(Clone, Copy, PartialEq, Eq, Debug, ValueEnum)]
pub(crate) enum OutputFormat {
    Text,
    Json,
    Jsonl,
    Csv,
    Tsv,
    Arrow,
    Parquet,
}

impl From<OutputFormat> for OutputKind {
    fn from(format: OutputFormat) -> Self {
        match format {
            OutputFormat::Text => Self::Text,
            OutputFormat::Json => Self::Json,
            OutputFormat::Jsonl => Self::JsonLines,
            OutputFormat::Csv => Self::Csv,
            OutputFormat::Tsv => Self::Tsv,
            OutputFormat::Arrow => Self::Arrow,
            OutputFormat::Parquet => Self::Parquet,
        }
    }
}

/// How much irregularity a read tolerates.
#[derive(Clone, Copy, PartialEq, Eq, Debug, ValueEnum)]
pub(crate) enum Tolerance {
    Strict,
    Permissive,
    Recover,
}

impl From<Tolerance> for molframe::ParseMode {
    fn from(tolerance: Tolerance) -> Self {
        match tolerance {
            Tolerance::Strict => Self::Strict,
            Tolerance::Permissive => Self::Permissive,
            Tolerance::Recover => Self::Recover,
        }
    }
}

fn default_output_format() -> OutputFormat {
    if std::io::stdout().is_terminal() {
        OutputFormat::Text
    } else {
        OutputFormat::Json
    }
}
