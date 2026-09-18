//! Shared fixture loading and format conversion helpers.

use std::io::Cursor;
use std::sync::OnceLock;

use molframe_core::io::{InputBuffer, Limits, ReadOptions};
use molframe_core::structure::Structure;

use crate::samples::Sample;

/// Wraps bytes already in hand as an [`InputBuffer`] (no decompression).
#[must_use]
pub fn input(bytes: &[u8]) -> InputBuffer {
    InputBuffer::from_bytes(bytes.to_vec())
}

/// Wraps gzip-compressed bytes as an [`InputBuffer`], decompressing them.
///
/// # Panics
///
/// Panics if the bytes are not valid gzip or exceed the default limits.
#[must_use]
pub fn input_gzip(bytes: &[u8]) -> InputBuffer {
    match InputBuffer::from_reader(Cursor::new(bytes.to_vec()), Limits::default()) {
        Ok(buffer) => buffer,
        Err(finding) => panic!("bench fixture gzip decode failed: {finding:?}"),
    }
}

/// Parses a sample into a [`Structure`].
///
/// # Panics
///
/// Panics if the embedded fixture fails to parse.
#[must_use]
pub fn structure(sample: Sample) -> Structure {
    if sample == Sample::Ensemble {
        return match structure_from_pdb(sample) {
            Some(structure) => structure,
            None => panic!("ensemble fixture must ship a parseable PDB file"),
        };
    }

    match sample {
        Sample::Tiny => cached_bcif(&TINY_BCIF, sample),
        Sample::Small => cached_bcif(&SMALL_BCIF, sample),
        Sample::Medium => cached_bcif(&MEDIUM_BCIF, sample),
        Sample::Large => cached_bcif(&LARGE_BCIF, sample),
        Sample::Ensemble => unreachable!("ensemble is loaded from PDB"),
    }
}

/// Parses a sample's legacy PDB into a [`Structure`], where one ships.
///
/// # Panics
///
/// Panics if the embedded fixture fails to parse.
#[must_use]
pub fn structure_from_pdb(sample: Sample) -> Option<Structure> {
    let cache = match sample {
        Sample::Tiny => &TINY_PDB,
        Sample::Small => &SMALL_PDB,
        Sample::Medium => &MEDIUM_PDB,
        Sample::Ensemble => &ENSEMBLE_PDB,
        Sample::Large => return None,
    };
    Some(cache.get_or_init(|| parse_pdb(sample)).clone())
}

fn cached_bcif(cache: &OnceLock<Structure>, sample: Sample) -> Structure {
    cache.get_or_init(|| parse_bcif(sample)).clone()
}

fn parse_bcif(sample: Sample) -> Structure {
    let buffer = input(sample.bcif());
    match molframe_bcif::read(&buffer, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!(
            "bench fixture {} (bcif) failed: {findings:?}",
            sample.label()
        ),
    }
}

fn parse_pdb(sample: Sample) -> Structure {
    let Some(bytes) = sample.pdb() else {
        panic!("sample {} has no PDB fixture", sample.label())
    };
    let buffer = input(bytes);
    match molframe_pdb::read(&buffer, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!(
            "bench fixture {} (pdb) failed: {findings:?}",
            sample.label()
        ),
    }
}

static TINY_BCIF: OnceLock<Structure> = OnceLock::new();
static SMALL_BCIF: OnceLock<Structure> = OnceLock::new();
static MEDIUM_BCIF: OnceLock<Structure> = OnceLock::new();
static LARGE_BCIF: OnceLock<Structure> = OnceLock::new();
static TINY_PDB: OnceLock<Structure> = OnceLock::new();
static SMALL_PDB: OnceLock<Structure> = OnceLock::new();
static MEDIUM_PDB: OnceLock<Structure> = OnceLock::new();
static ENSEMBLE_PDB: OnceLock<Structure> = OnceLock::new();

/// Parses a sample's mmCIF into a [`Structure`], where one ships uncompressed.
///
/// # Panics
///
/// Panics if the embedded fixture fails to parse.
#[must_use]
pub fn structure_from_cif(sample: Sample) -> Option<Structure> {
    let bytes = sample.cif()?;
    let buffer = input(bytes);
    match molframe_cif::read(&buffer, &ReadOptions::new()) {
        Ok((structure, _)) => Some(structure),
        Err(findings) => panic!(
            "bench fixture {} (cif) failed: {findings:?}",
            sample.label()
        ),
    }
}
