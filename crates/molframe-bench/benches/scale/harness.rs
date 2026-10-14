//! Case declaration, command line, timing and the sweep loop.

use crate::alloc;
use molframe_bench::Tier;
use molframe_core::ExecutionContext;
use std::time::{Duration, Instant};

/// What a case's size axis counts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Axis {
    /// Atoms; the tier's atom count.
    Atoms,
    /// Frames of a fixed-size system: a hundredth of the tier's atom count.
    Frames,
    /// Whole structures: a hundredth of the tier's atom count.
    Structures,
}

impl Axis {
    pub fn size(self, tier: Tier) -> usize {
        match self {
            Axis::Atoms => tier.atoms(),
            Axis::Frames | Axis::Structures => tier.atoms() / 100,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Axis::Atoms => "atoms",
            Axis::Frames => "frames",
            Axis::Structures => "structures",
        }
    }
}

/// One measurement of one kernel.
#[derive(Clone, Copy, Debug)]
pub struct Sample {
    pub best: Duration,
    pub reps: u32,
    pub peak_bytes: usize,
    pub allocations: usize,
    pub checksum: u64,
}

/// What a case reports back: the sample and the work units it covered.
#[derive(Clone, Copy, Debug)]
pub struct Outcome {
    pub sample: Sample,
    pub items: u64,
}

/// Per-run measurement handle given to every case.
pub struct Run<'a> {
    pub context: ExecutionContext,
    pub settings: &'a Settings,
}

impl Run<'_> {
    /// Times `kernel`, repeating while the repetition budget allows.
    ///
    /// `fingerprint` reduces the kernel's result to a checksum outside the
    /// timed region. Peak heap is measured above the live size at entry, so
    /// inputs built before the call are not charged to it, while the live
    /// result is.
    pub fn measure<T>(
        &self,
        mut kernel: impl FnMut() -> Result<T, String>,
        fingerprint: impl Fn(&T) -> u64,
    ) -> Result<Sample, String> {
        let baseline = alloc::reset_peak();
        let allocations_before = alloc::count();
        let started = Instant::now();
        let first = kernel()?;
        let mut best = started.elapsed();
        let checksum = fingerprint(&first);
        drop(first);
        let mut reps = 1;
        let mut spent = best;
        let mut stable = true;
        while reps < self.settings.max_reps && spent < self.settings.rep_budget {
            let started = Instant::now();
            let again = kernel()?;
            let took = started.elapsed();
            stable &= fingerprint(&again) == checksum;
            best = best.min(took);
            spent += took;
            reps += 1;
        }
        Ok(Sample {
            best,
            reps,
            peak_bytes: alloc::peak().saturating_sub(baseline),
            allocations: (alloc::count() - allocations_before)
                / match usize::try_from(reps) {
                    Ok(reps) => reps,
                    Err(_) => 1,
                },
            checksum: if stable { checksum } else { u64::MAX },
        })
    }
}

/// A declared workload.
pub struct Case {
    pub id: &'static str,
    pub axis: Axis,
    pub tiers: &'static [Tier],
    pub run: fn(&Run<'_>, usize) -> Result<Outcome, String>,
}

/// One row of the report.
#[derive(Clone, Debug)]
pub struct Row {
    pub case: &'static str,
    pub axis: Axis,
    pub tier: Tier,
    pub size: usize,
    pub workers: usize,
    pub outcome: Result<Outcome, String>,
}

/// Command-line settings.
#[derive(Clone, Debug)]
pub struct Settings {
    pub filter: Option<String>,
    pub max_tier: Tier,
    pub workers: Vec<usize>,
    pub budget: Duration,
    pub rep_budget: Duration,
    pub max_reps: u32,
    pub out: Option<String>,
}

impl Settings {
    pub fn from_args(args: impl Iterator<Item = String>) -> Result<Self, String> {
        let mut settings = Self {
            filter: None,
            max_tier: Tier::Medium,
            workers: vec![0],
            budget: Duration::from_secs(60),
            rep_budget: Duration::from_millis(500),
            max_reps: 5,
            out: None,
        };
        let mut args = args.peekable();
        while let Some(flag) = args.next() {
            if flag == "--bench" {
                continue;
            }
            let mut value = || args.next().ok_or_else(|| format!("{flag} needs a value"));
            match flag.as_str() {
                "--filter" => settings.filter = Some(value()?),
                "--max-tier" => {
                    settings.max_tier = Tier::parse(&value()?).map_err(|e| e.to_string())?;
                }
                "--workers" => {
                    settings.workers = value()?
                        .split(',')
                        .map(|w| w.parse::<usize>().map_err(|e| format!("--workers: {e}")))
                        .collect::<Result<_, _>>()?;
                }
                "--budget-secs" => {
                    let secs = value()?
                        .parse::<u64>()
                        .map_err(|e| format!("--budget-secs: {e}"))?;
                    settings.budget = Duration::from_secs(secs);
                }
                "--out" => settings.out = Some(value()?),
                other => return Err(format!("unknown argument `{other}`")),
            }
        }
        Ok(settings)
    }
}

/// Runs every selected case over its tiers and worker counts.
///
/// A tier is skipped when the previous tier's best time, scaled by the size
/// ratio, already exceeds the per-measurement budget (linear growth is the
/// optimistic assumption).
pub fn run(cases: &[&Case], settings: &Settings) -> Vec<Row> {
    let mut rows = Vec::new();
    for &workers in &settings.workers {
        let context = match workers {
            0 => ExecutionContext::default(),
            count => match ExecutionContext::builder().worker_budget(count).build() {
                Ok(context) => context,
                Err(error) => {
                    eprintln!("scale: workers={count}: {error}");
                    continue;
                }
            },
        };
        let effective = context.worker_budget();
        let handle = Run { context, settings };
        for case in cases {
            let mut previous: Option<(usize, Duration)> = None;
            for &tier in case.tiers.iter().filter(|t| **t <= settings.max_tier) {
                let size = case.axis.size(tier);
                if let Some((before, took)) = previous {
                    let ratio = real(size as u64) / real(before as u64);
                    if took.as_secs_f64() * ratio > settings.budget.as_secs_f64() {
                        eprintln!("scale: {} {}: skipped (budget)", case.id, tier.label());
                        continue;
                    }
                }
                eprintln!("scale: {} {} w={effective}", case.id, tier.label());
                let outcome = (case.run)(&handle, size);
                if let Ok(done) = &outcome {
                    previous = Some((size, done.sample.best));
                }
                rows.push(Row {
                    case: case.id,
                    axis: case.axis,
                    tier,
                    size,
                    workers: effective,
                    outcome,
                });
            }
        }
    }
    rows
}

/// Fingerprint for kernels whose result is one real number.
pub fn real_fingerprint(value: f64) -> u64 {
    let mut sum = molframe_bench::Checksum::default();
    sum.real(value);
    sum.value()
}

/// Converts a count to a real for rates and ratios.
pub fn real(count: u64) -> f64 {
    molframe_core::numeric::u64_to_f64(count)
}
