//! Elastic-network normal modes over selected sites.

use crate::args::{ModeTable, NetworkModel, NormalModeArguments};
use crate::commands::{open, select_text};
use crate::exit::Exit;
use crate::report::{Context, emit_rows};
use molframe::analysis::{AnmOptions, GnmOptions};
use molframe::spatial::SpatialBackend;
use molframe_core::selection::AtomSelection;

/// What either network yields, in the shape the tables need.
struct Solved {
    eigenvalues: Vec<f64>,
    sites: Vec<u32>,
    fluctuations: Vec<f64>,
    zero_modes: usize,
}

fn failed(model: &str, error: &dyn std::fmt::Display) -> Exit {
    eprintln!("{model} failed: {error}");
    Exit::Consistency
}

fn gnm(
    args: &NormalModeArguments,
    positions: &[[f32; 3]],
    sites: &AtomSelection,
    context: Context,
) -> Result<Solved, Exit> {
    let Some(reduction) = args.reduction else {
        eprintln!("gnm requires --reduction");
        return Err(Exit::Usage);
    };
    let model = molframe::analysis::gaussian_network_model(
        positions,
        sites,
        GnmOptions {
            contact_distance: args.contact_distance,
            mode_count: args.mode_count,
            zero_mode_tolerance: args.zero_mode_tolerance,
            memory_limit_bytes: args.memory_limit,
            backend: SpatialBackend::Auto,
            reduction: reduction.into(),
        },
        None,
        context.execution,
    )
    .map_err(|error| failed("gnm", &error))?;
    Ok(Solved {
        fluctuations: model.fluctuations(),
        eigenvalues: model.eigenvalues,
        sites: model.sites,
        zero_modes: model.zero_modes,
    })
}

fn anm(
    args: &NormalModeArguments,
    positions: &[[f32; 3]],
    sites: &AtomSelection,
    context: Context,
) -> Result<Solved, Exit> {
    if args.reduction.is_some() {
        eprintln!("--reduction applies only to gnm");
        return Err(Exit::Usage);
    }
    let model = molframe::analysis::anisotropic_network_model(
        positions,
        sites,
        AnmOptions {
            contact_distance: args.contact_distance,
            mode_count: args.mode_count,
            zero_mode_tolerance: args.zero_mode_tolerance,
            memory_limit_bytes: args.memory_limit,
            backend: SpatialBackend::Auto,
        },
        None,
        context.execution,
    )
    .map_err(|error| failed("anm", &error))?;
    Ok(Solved {
        fluctuations: model.fluctuations(),
        eigenvalues: model.eigenvalues,
        sites: model.sites,
        zero_modes: model.zero_modes,
    })
}

pub(crate) fn normal_modes(args: &NormalModeArguments, context: Context) -> Exit {
    let structure = match open(&args.input, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let origin = args.input.display().to_string();
    let sites = match select_text(&structure, &args.sites, context.policy, context.execution) {
        Ok(evaluation) => {
            context.findings(&evaluation.warnings, &origin);
            evaluation.selection
        }
        Err(findings) => {
            context.findings(&findings, &origin);
            return Exit::of(&findings);
        }
    };
    let solved = match args.network {
        NetworkModel::Gnm => gnm(args, structure.coordinates(), &sites, context),
        NetworkModel::Anm => anm(args, structure.coordinates(), &sites, context),
    };
    let solved = match solved {
        Ok(solved) => solved,
        Err(exit) => return exit,
    };
    eprintln!("zero modes: {}", solved.zero_modes);
    match args.table {
        ModeTable::Modes => emit_rows(
            context,
            &["mode", "eigenvalue"],
            solved.eigenvalues.into_iter().enumerate(),
            |(mode, value)| vec![(mode + 1).to_string(), value.to_string()],
        ),
        ModeTable::Fluctuations => emit_rows(
            context,
            &["site", "fluctuation"],
            solved.sites.into_iter().zip(solved.fluctuations),
            |(site, value)| vec![site.to_string(), value.to_string()],
        ),
    }
}
