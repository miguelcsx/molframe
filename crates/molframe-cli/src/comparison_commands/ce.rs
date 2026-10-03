//! Combinatorial extension: sequence-order alignment of guide atoms.

use super::command::ComparisonOptions;
use crate::exit::Exit;
use crate::report::Context;
use molframe::compare::{CeOptions, CeSignificanceProfile};

/// Guide-atom coordinates of the selection `query` in one structure.
fn guide(
    structure: &molframe::Structure,
    query: &str,
    context: Context,
    origin: &str,
) -> Result<Vec<[f32; 3]>, Exit> {
    let evaluation =
        crate::commands::select_text(structure, query, context.policy, context.execution).map_err(
            |findings| {
                context.findings(&findings, origin);
                Exit::of(&findings)
            },
        )?;
    let positions = structure.coordinates();
    Ok((&evaluation.selection)
        .into_iter()
        .filter_map(|atom| positions.get(atom as usize).copied())
        .collect())
}

pub(super) fn measure(
    model: &molframe::Structure,
    reference: &molframe::Structure,
    options: &ComparisonOptions<'_>,
    context: Context,
) -> Result<Vec<(&'static str, f64)>, Exit> {
    let extra = options.extra;
    let (
        Some(query),
        Some(window_size),
        Some(max_gap),
        Some(max_paths),
        Some(fragment_similarity_threshold),
        Some(path_similarity_threshold),
        Some(memory_limit_bytes),
    ) = (
        extra.guide.as_deref(),
        extra.ce_window,
        extra.ce_max_gap,
        extra.ce_max_paths,
        extra.ce_fragment_threshold,
        extra.ce_path_threshold,
        extra.ce_memory_limit,
    )
    else {
        eprintln!(
            "ce requires --guide, --ce-window, --ce-max-gap, --ce-max-paths, \
             --ce-fragment-threshold, --ce-path-threshold and --ce-memory-limit"
        );
        return Err(Exit::Usage);
    };
    let ce_options = CeOptions {
        window_size,
        max_gap,
        max_paths,
        fragment_similarity_threshold,
        path_similarity_threshold,
        significance: extra
            .ce_significance
            .then_some(CeSignificanceProfile::OriginalWindowEight),
        memory_limit_bytes,
    };
    let fixed = guide(reference, query, context, "reference")?;
    let mobile = guide(model, query, context, "model")?;
    let alignment = molframe::compare::ce_align(&fixed, &mobile, ce_options).map_err(|error| {
        eprintln!("ce failed: {error}");
        Exit::Consistency
    })?;
    let mut rows = vec![
        ("ce_rmsd", alignment.rmsd),
        ("ce_aligned_atoms", count(alignment.reference_indices.len())),
        ("ce_fragments", count(alignment.fragment_count)),
        ("ce_similarity", alignment.similarity),
    ];
    if let Some(z_score) = alignment.z_score {
        rows.push(("ce_z_score", z_score));
    }
    Ok(rows)
}

#[allow(clippy::cast_precision_loss)]
const fn count(value: usize) -> f64 {
    value as f64
}
