//! The reference tables shipped inside the crate.
//!
//! Embedded at compile time and decoded once per process. Provenance and
//! licences are recorded here and in `data/README.md`; a result built from
//! these tables reports them through [`TableInfo`].

use std::sync::OnceLock;

use super::reference_data::{RamaGrids, RotamerTables};
use super::reference_tables::{ChiDefinitions, GeometryTargets};
use super::{BuiltinError, Part, ReferenceSet, TableInfo};

const RAMA: &[u8] = include_bytes!("../../data/rama.bin");
const ROTAMER: &[u8] = include_bytes!("../../data/rotamer.bin");
const GEOMETRY: &str = include_str!("../../data/geometry.tsv");
const CHI: &str = include_str!("../../data/chi_atoms.tsv");

const TOP8000: &str = "https://github.com/rlabduke/reference_data";
const TOP8000_COMMIT: &str = "5ee6875fc29eccc3c9dc7cbe705ff9fd1b505d7d";

fn rama_info() -> TableInfo {
    TableInfo::new(
        "top8000-ramachandran",
        TOP8000_COMMIT,
        format!("Richardson Lab Top8000, {TOP8000}"),
        "CC-BY-4.0",
    )
}

fn rotamer_info() -> TableInfo {
    TableInfo::new(
        "top8000-rotamer",
        TOP8000_COMMIT,
        format!("Richardson Lab Top8000, {TOP8000}"),
        "CC-BY-4.0",
    )
}

fn geometry_info() -> TableInfo {
    TableInfo::new(
        "engh-huber-1991",
        "openstructure-stereo_chemical_props",
        "Engh and Huber 1991 via OpenStructure; link terms AlphaFold and CCP4 monomer library",
        "LGPL-3.0 (E&H, CCP4 O-C-N) and Apache-2.0 (AlphaFold link terms)",
    )
}

fn chi_info() -> TableInfo {
    TableInfo::new(
        "chi-atoms",
        "molframe-1",
        "standard chi atom names, molframe-validate/data/chi_atoms.tsv",
        "CC-BY-4.0 (pairs with the Top8000 tables)",
    )
}

fn build() -> Result<ReferenceSet, BuiltinError> {
    Ok(ReferenceSet {
        rama: Part::new(rama_info(), RamaGrids::from_bytes(RAMA)?),
        rotamers: Part::new(rotamer_info(), RotamerTables::from_bytes(ROTAMER)?),
        geometry: Part::new(geometry_info(), GeometryTargets::parse(GEOMETRY)),
        chi: Part::new(chi_info(), ChiDefinitions::parse(CHI)),
    })
}

/// The bundled default reference set, decoded on first use.
///
/// # Errors
///
/// Returns [`BuiltinError`] if an embedded table is malformed.
pub fn bundled_reference() -> Result<&'static ReferenceSet, BuiltinError> {
    static CELL: OnceLock<Result<ReferenceSet, BuiltinError>> = OnceLock::new();
    CELL.get_or_init(build).as_ref().map_err(Clone::clone)
}

#[cfg(test)]
#[path = "bundled_tests.rs"]
mod tests;
