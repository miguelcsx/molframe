//! Tab-separated report with scaling exponents and determinism verdicts.

use crate::harness::{Row, Settings};
use molframe_bench::{Tier, fit_exponent};
use std::collections::BTreeMap;
use std::fmt::Write as _;

pub fn emit(rows: &[Row], settings: &Settings, _tiers: [Tier; 5]) {
    let mut text = String::new();
    let _ = writeln!(
        text,
        "# molframe scale sweep; cpus={} max_tier={}",
        std::thread::available_parallelism().map_or(0, std::num::NonZero::get),
        settings.max_tier.label()
    );
    let _ = writeln!(
        text,
        "case\ttier\taxis\tsize\tworkers\tsecs\titems_per_sec\tpeak_mib\tallocs\treps\tchecksum"
    );
    for row in rows {
        match &row.outcome {
            Ok(done) => {
                let secs = done.sample.best.as_secs_f64();
                let _ = writeln!(
                    text,
                    "{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.4e}\t{:.2}\t{}\t{}\t{:016x}",
                    row.case,
                    row.tier.label(),
                    row.axis.label(),
                    row.size,
                    row.workers,
                    secs,
                    crate::harness::real(done.items) / secs.max(1e-12),
                    crate::harness::real(done.sample.peak_bytes as u64) / 1_048_576.0,
                    done.sample.allocations,
                    done.sample.reps,
                    done.sample.checksum,
                );
            }
            Err(message) => {
                let _ = writeln!(
                    text,
                    "{}\t{}\t{}\t{}\t{}\tERROR {message}",
                    row.case,
                    row.tier.label(),
                    row.axis.label(),
                    row.size,
                    row.workers
                );
            }
        }
    }
    exponents(&mut text, rows);
    determinism(&mut text, rows);
    print!("{text}");
    if let Some(path) = &settings.out
        && let Err(error) = std::fs::write(path, &text)
    {
        eprintln!("scale: cannot write {path}: {error}");
    }
}

fn exponents(text: &mut String, rows: &[Row]) {
    let mut series: BTreeMap<(&str, usize), Vec<(f64, f64)>> = BTreeMap::new();
    for row in rows {
        if let Ok(done) = &row.outcome {
            series.entry((row.case, row.workers)).or_default().push((
                crate::harness::real(row.size as u64),
                done.sample.best.as_secs_f64(),
            ));
        }
    }
    for ((case, workers), samples) in series {
        if let Some(fit) = fit_exponent(&samples) {
            let _ = writeln!(
                text,
                "# exponent {case} workers={workers} k={:.2} r2={:.3}",
                fit.exponent, fit.r_squared
            );
        }
    }
}

fn determinism(text: &mut String, rows: &[Row]) {
    let mut seen: BTreeMap<(&str, usize), (u64, usize)> = BTreeMap::new();
    for row in rows {
        let Ok(done) = &row.outcome else { continue };
        let entry = seen
            .entry((row.case, row.size))
            .or_insert((done.sample.checksum, row.workers));
        if entry.0 != done.sample.checksum {
            let _ = writeln!(
                text,
                "# NONDETERMINISTIC {} size={} workers {} vs {}",
                row.case, row.size, entry.1, row.workers
            );
        }
    }
}
