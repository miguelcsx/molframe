//! Thin rendering over the native interaction kernels.
//!
//! Each command opens the structure, applies the explicit dictionary so donor,
//! acceptor, charge and aromatic roles exist, calls one kernel and prints its
//! table. Nothing is computed here.

use crate::commands::open;
use crate::exit::Exit;
use crate::report::{Context, emit_rows};
use molframe::analysis::{CationPiOptions, HydrogenBondOptions, PiStackingOptions};
use molframe::spatial::SpatialBackend;
use std::path::Path;

/// The structure to analyse and the dictionary that annotates it.
#[derive(Clone, Copy)]
pub(crate) struct Source<'a> {
    pub(crate) input: &'a Path,
    pub(crate) ccd: &'a Path,
    pub(crate) ccd_version: &'a str,
}

fn annotated(source: Source<'_>, context: Context) -> Result<molframe::Structure, Exit> {
    let structure = open(source.input, context)?;
    crate::chemistry::annotate(&structure, source.ccd, source.ccd_version, context)
}

fn failed(what: &str, error: &dyn std::fmt::Display) -> Exit {
    eprintln!("{what} failed: {error}");
    Exit::Usage
}

pub(crate) fn hbonds(source: Source<'_>, options: HydrogenBondOptions, context: Context) -> Exit {
    let structure = match annotated(source, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let table =
        match molframe::analysis::hydrogen_bonds(structure.engine(), options, context.execution) {
            Ok(table) => table,
            Err(error) => return failed("hydrogen-bond search", &error),
        };
    emit_rows(
        context,
        &[
            "donor",
            "hydrogen",
            "acceptor",
            "donor_acceptor_distance_angstrom",
            "hydrogen_acceptor_distance_angstrom",
            "angle_degrees",
        ],
        table.iter(),
        |bond| {
            vec![
                bond.donor.get().to_string(),
                bond.hydrogen.get().to_string(),
                bond.acceptor.get().to_string(),
                bond.donor_acceptor_distance.to_string(),
                bond.hydrogen_acceptor_distance.to_string(),
                bond.angle_degrees.to_string(),
            ]
        },
    )
}

pub(crate) fn salt_bridges(source: Source<'_>, max_distance: f32, context: Context) -> Exit {
    let structure = match annotated(source, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let table = match molframe::analysis::salt_bridges(
        structure.engine(),
        max_distance,
        SpatialBackend::Auto,
        context.execution,
    ) {
        Ok(table) => table,
        Err(error) => return failed("salt-bridge search", &error),
    };
    emit_rows(
        context,
        &["anion", "cation", "distance_angstrom"],
        table.iter(),
        |bridge| {
            vec![
                bridge.anion.get().to_string(),
                bridge.cation.get().to_string(),
                bridge.distance.to_string(),
            ]
        },
    )
}

pub(crate) fn pi_stacking(
    source: Source<'_>,
    options: PiStackingOptions,
    context: Context,
) -> Exit {
    let structure = match annotated(source, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let table = match molframe::analysis::pi_stacking(structure.engine(), options) {
        Ok(table) => table,
        Err(error) => return failed("pi-stacking search", &error),
    };
    emit_rows(
        context,
        &[
            "first_residue",
            "second_residue",
            "centre_distance_angstrom",
            "angle_degrees",
            "kind",
        ],
        table.iter(),
        |stack| {
            vec![
                stack.first.get().to_string(),
                stack.second.get().to_string(),
                stack.centre_distance.to_string(),
                stack.angle.to_string(),
                stacking_name(stack.kind).to_owned(),
            ]
        },
    )
}

const fn stacking_name(kind: molframe::analysis::StackingKind) -> &'static str {
    match kind {
        molframe::analysis::StackingKind::Parallel => "parallel",
        molframe::analysis::StackingKind::TShaped => "t-shaped",
    }
}

pub(crate) fn cation_pi(source: Source<'_>, options: CationPiOptions, context: Context) -> Exit {
    let structure = match annotated(source, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let table = match molframe::analysis::cation_pi(structure.engine(), options) {
        Ok(table) => table,
        Err(error) => return failed("cation-pi search", &error),
    };
    emit_rows(
        context,
        &["cation_residue", "ring_residue", "distance_angstrom"],
        table.iter(),
        |pair| {
            vec![
                pair.cation_residue.get().to_string(),
                pair.ring_residue.get().to_string(),
                pair.distance.to_string(),
            ]
        },
    )
}

pub(crate) fn water_bridges(
    source: Source<'_>,
    options: HydrogenBondOptions,
    context: Context,
) -> Exit {
    let structure = match annotated(source, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let table = match molframe::analysis::water_bridges(
        structure.engine(),
        molframe::analysis::WaterBridgeOptions {
            hydrogen_bonds: options,
        },
        context.execution,
    ) {
        Ok(table) => table,
        Err(error) => return failed("water-bridge search", &error),
    };
    emit_rows(
        context,
        &["water", "first", "second"],
        table.iter(),
        |bridge| {
            vec![
                bridge.water.get().to_string(),
                bridge.first.get().to_string(),
                bridge.second.get().to_string(),
            ]
        },
    )
}
