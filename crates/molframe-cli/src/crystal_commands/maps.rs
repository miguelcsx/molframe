//! Density-map statistics and map-to-map correlation.

use super::super::args::HistogramArguments;
use crate::exit::Exit;
use crate::report::{Context, emit_rows};
use molframe::crystal::DensityMap;
use std::path::Path;

fn read(path: &Path) -> Result<DensityMap, Exit> {
    let bytes = std::fs::read(path).map_err(|error| {
        eprintln!("{} could not be read: {error}", path.display());
        Exit::Failure
    })?;
    DensityMap::from_mrc_bytes(&bytes).map_err(|error| {
        eprintln!("{} is not a usable density map: {error}", path.display());
        Exit::Consistency
    })
}

pub(super) fn statistics(path: &Path, histogram: &HistogramArguments, context: Context) -> Exit {
    let map = match read(path) {
        Ok(map) => map,
        Err(exit) => return exit,
    };
    let stats = match map.statistics() {
        Ok(stats) => stats,
        Err(error) => {
            eprintln!("map statistics failed: {error}");
            return Exit::Consistency;
        }
    };
    let mut rows: Vec<(String, String)> = vec![
        ("count".to_owned(), stats.count.to_string()),
        ("minimum".to_owned(), stats.minimum.to_string()),
        ("maximum".to_owned(), stats.maximum.to_string()),
        ("mean".to_owned(), stats.mean.to_string()),
        ("sigma".to_owned(), stats.sigma.to_string()),
    ];
    if let Err(exit) = add_histogram(&map, histogram, &mut rows) {
        return exit;
    }
    emit_rows(
        context,
        &["property", "value"],
        rows,
        |(property, value)| vec![property, value],
    )
}

/// Appends one `histogram_bin_N` row per bin when a histogram was requested.
fn add_histogram(
    map: &DensityMap,
    histogram: &HistogramArguments,
    rows: &mut Vec<(String, String)>,
) -> Result<(), Exit> {
    let (Some(bins), Some(range)) = (
        histogram.histogram_bins,
        histogram.histogram_range.as_deref(),
    ) else {
        return Ok(());
    };
    let &[minimum, maximum] = range else {
        eprintln!("--histogram-range takes exactly MIN MAX");
        return Err(Exit::Usage);
    };
    let counted = map.histogram(bins, minimum, maximum).map_err(|error| {
        eprintln!("map histogram failed: {error}");
        Exit::Consistency
    })?;
    for (bin, count) in counted.counts.iter().enumerate() {
        rows.push((format!("histogram_bin_{bin}"), count.to_string()));
    }
    Ok(())
}

pub(super) fn correlation(observed: &Path, calculated: &Path, context: Context) -> Exit {
    let (first, second) = match (read(observed), read(calculated)) {
        (Ok(first), Ok(second)) => (first, second),
        (Err(exit), _) | (_, Err(exit)) => return exit,
    };
    match molframe::validation::real_space_map_correlation(&first, &second) {
        Ok(result) => emit_rows(
            context,
            &["property", "value"],
            [
                ("coefficient", result.coefficient.to_string()),
                ("samples", result.sample_count.to_string()),
                ("observed_mean", result.observed_mean.to_string()),
                ("calculated_mean", result.calculated_mean.to_string()),
            ],
            |(property, value)| vec![property.to_owned(), value],
        ),
        Err(error) => {
            eprintln!("map correlation failed: {error}");
            Exit::Consistency
        }
    }
}
