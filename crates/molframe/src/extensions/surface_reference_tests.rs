//! SASA against `FreeSASA` 2.2.1 on the committed 1CRN and 4HHB files.
//!
//! The reference totals were computed by `FreeSASA` on the same ATOM records
//! with the same Bondi element radii (haem iron at its Alvarez radius where
//! hetero atoms are included), a 1.4 Å probe and its 960-point
//! Shrake-Rupley and 20-slice Lee-Richards settings, so the comparison isolates
//! the algorithms and the heteroatom rule from the radius table. 4HHB carries
//! haem and 221 waters; both programs leave them out by default.

use super::*;
use crate::prelude::*;
use molframe_chem::{HeteroPolicy, MetalFallback};

const CRAMBIN: &[u8] = include_bytes!("../../../molframe-bench/data/1crn.pdb");
const HEMOGLOBIN: &[u8] = include_bytes!("../../../molframe-bench/data/4hhb.pdb");

fn read(bytes: &[u8]) -> crate::Structure {
    match crate::read_bytes(bytes.to_vec(), Some("e.pdb"), &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture read failed: {findings:?}"),
    }
}

fn total(structure: &crate::Structure, options: &SasaOptions) -> f64 {
    let areas = structure
        .sasa_with(options, &ExecutionContext::default())
        .expect("sasa");
    areas.iter().filter(|area| area.is_finite()).sum()
}

fn within(value: f64, reference: f64, relative: f64) -> bool {
    (value - reference).abs() <= reference * relative
}

fn lee_richards() -> SasaOptions {
    SasaOptions {
        method: SasaMethod::LeeRichards { slices: 20 },
        ..SasaOptions::default()
    }
}

#[test]
fn crambin_matches_freesasa_for_both_methods() {
    let structure = read(CRAMBIN);
    let sr = total(&structure, &SasaOptions::default());
    let lr = total(&structure, &lee_richards());
    assert!(within(sr, 2970.152, 0.005), "SR {sr}");
    assert!(within(lr, 2969.918, 0.005), "LR {lr}");
}

#[test]
fn haemoglobin_default_measures_the_polymer_like_freesasa() {
    let structure = read(HEMOGLOBIN);
    let sr = total(&structure, &SasaOptions::default());
    let lr = total(&structure, &lee_richards());
    assert!(within(sr, 26097.315, 0.005), "SR {sr}");
    assert!(within(lr, 26119.71, 0.005), "LR {lr}");
}

#[test]
fn excluded_heteroatoms_have_nan_area_and_do_not_occlude() {
    let structure = read(HEMOGLOBIN);
    let areas = structure.sasa().expect("sasa");
    let finite = areas.iter().filter(|area| area.is_finite()).count();
    assert_eq!(finite, 4384);
    assert_eq!(areas.len() - finite, structure.atom_count() as usize - 4384);
}

#[test]
fn including_every_heteroatom_matches_freesasa_with_haem_iron_at_alvarez_radius() {
    let structure = read(HEMOGLOBIN);
    let everything = |method| SasaOptions {
        method,
        hetero: HeteroPolicy::ALL,
        ..SasaOptions::default()
    };
    let sr = total(
        &structure,
        &everything(SasaMethod::ShrakeRupley { points: 960 }),
    );
    let lr = total(
        &structure,
        &everything(SasaMethod::LeeRichards { slices: 20 }),
    );
    assert!(within(sr, 24073.036, 0.005), "SR {sr}");
    assert!(within(lr, 24094.362, 0.005), "LR {lr}");
}

#[test]
fn a_metal_without_a_bondi_radius_is_measured_only_when_asked() {
    let structure = read(HEMOGLOBIN);
    let iron = structure
        .atoms()
        .into_iter()
        .find(|atom| atom.element() == Some(molframe_core::Element::IRON))
        .map(|atom| atom.index().as_usize())
        .expect("haem iron");
    let options = |metals| SasaOptions {
        hetero: HeteroPolicy::SOLUTE,
        metals,
        ..SasaOptions::default()
    };
    let ask = |metals| {
        structure
            .sasa_with(&options(metals), &ExecutionContext::default())
            .expect("sasa")[iron]
    };
    assert!(ask(MetalFallback::None).is_nan());
    assert!(ask(MetalFallback::Vdw).is_finite());
}
