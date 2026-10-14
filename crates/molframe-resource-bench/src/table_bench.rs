//! Allocation probes for the table-driven SASA and clashscore paths.
//!
//! The reference tables are resolved once, before the measured window, so the
//! counts show what a calculation costs once its tables are loaded: per-atom
//! radius lookups are indexed reads and allocate no more than the output.
//! Cost: one table resolution outside the window; inside it, O(atoms).

use super::{ResourceRecord, measure_case};
use std::hint::black_box;

pub(super) fn run_sasa_tables_medium() -> Result<ResourceRecord, String> {
    let structure = molframe_bench::structure(molframe_bench::Sample::Medium);
    let options = molframe::SasaOptions::resolve(
        &molframe::chemistry::TableRequest::profile(),
        molframe_bench::table_store(),
    )
    .map_err(|error| format!("radius table: {error}"))?;
    measure_case("sasa_tables_medium", || {
        let radii = options.atom_radii(&structure);
        let areas = molframe::surface::measure_atoms(
            structure.positions(),
            &radii,
            &options,
            &molframe::ExecutionContext::default(),
        )
        .map_err(|error| format!("SASA failed: {error}"))?;
        Ok(black_box(areas.len() as u64))
    })
}

pub(super) fn run_clashscore_medium() -> Result<ResourceRecord, String> {
    let structure = molframe_bench::structure(molframe_bench::Sample::Medium);
    let options = molframe::validation::ClashScoreOptions::default();
    measure_case("clashscore_medium", || {
        let score = molframe::validation::clash_score(
            &structure,
            &options,
            molframe_bench::table_store(),
            &molframe::ExecutionContext::default(),
        )
        .map_err(|error| format!("clashscore failed: {error}"))?;
        Ok(black_box(score.clashes as u64))
    })
}
