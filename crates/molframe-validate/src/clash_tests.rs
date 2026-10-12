use super::{ClashOptions, clashes, clashes_with};
use molframe_chem::RadiusSet;
use molframe_core::bond::{BondProvenance, BondRecord, BondTableBuilder};
use molframe_core::index::AtomIndex;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;
use molframe_core::{BondOrder, ExecutionContext, OptionalSymbol};

use molframe_spatial::SpatialBackend;
use std::fmt::Write;

const HEADER: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

fn structure(source: &str) -> Structure {
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

#[test]
fn two_carbons_deep_inside_each_other_clash() {
    // Bondi carbon radius is 1.70 Å; at 2.0 Å apart they overlap by 1.4 Å.
    let source = format!(
        "{HEADER}\
ATOM 1 C C1 LIG A 1 0 0 0\n\
ATOM 2 C C2 LIG A 2 2 0 0\n"
    );
    let Ok(found) = clashes(
        &structure(&source),
        0.4,
        RadiusSet::Bondi,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    assert_eq!(found.len(), 1);
    let clash = found.iter().next().expect("one clash");
    assert_eq!((clash.first.get(), clash.second.get()), (0, 1));
    assert!(clash.overlap > 0.4);
}

#[test]
fn atoms_that_barely_touch_do_not_clash() {
    // At 3.2 Å the 1.70 Å radii overlap by only 0.2 Å, under the tolerance.
    let source = format!(
        "{HEADER}\
ATOM 1 C C1 LIG A 1 0 0 0\n\
ATOM 2 C C2 LIG A 2 3.2 0 0\n"
    );
    let Ok(found) = clashes(
        &structure(&source),
        0.4,
        RadiusSet::Bondi,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    assert!(found.is_empty());
}

#[test]
fn every_backend_agrees_on_the_clashes() {
    let source = format!(
        "{HEADER}\
ATOM 1 C C1 LIG A 1 0 0 0\n\
ATOM 2 C C2 LIG A 2 2 0 0\n\
ATOM 3 C C3 LIG A 3 2 2 0\n"
    );
    let structure = structure(&source);
    let Ok(baseline) = clashes(
        &structure,
        0.4,
        RadiusSet::Bondi,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    for backend in [SpatialBackend::CellList, SpatialBackend::KdTree] {
        let Ok(other) = clashes(
            &structure,
            0.4,
            RadiusSet::Bondi,
            backend,
            &ExecutionContext::default(),
        ) else {
            panic!("valid");
        };
        assert_eq!(baseline, other, "backend {backend:?} disagreed");
    }
}

#[test]
fn the_clash_list_is_identical_at_every_worker_count() {
    // A dense lattice guarantees many overlapping pairs across several cell
    // blocks, so the parallel path is genuinely exercised.
    // Built through the same mmCIF path the other tests use, at a spacing that
    // guarantees many overlapping pairs across several cell blocks.
    let mut source = String::from(HEADER);
    let mut serial_id = 0;
    for x in 0..6i16 {
        for y in 0..6i16 {
            for z in 0..6i16 {
                serial_id += 1;
                let _ = writeln!(
                    source,
                    "ATOM {serial_id} C C{serial_id} LIG A {serial_id} {:.3} {:.3} {:.3}",
                    f32::from(x) * 1.2,
                    f32::from(y) * 1.2,
                    f32::from(z) * 1.2
                );
            }
        }
    }
    let structure = structure(&source);

    let serial = clashes(
        &structure,
        0.1,
        RadiusSet::Bondi,
        SpatialBackend::CellList,
        &ExecutionContext::default(),
    )
    .expect("serial clash search");
    assert!(!serial.is_empty(), "the fixture must actually clash");

    for workers in [2, 4, 8, 16] {
        let context = ExecutionContext::builder()
            .worker_budget(workers)
            .build()
            .expect("worker context is valid");
        let parallel = clashes(
            &structure,
            0.1,
            RadiusSet::Bondi,
            SpatialBackend::CellList,
            &context,
        )
        .expect("parallel clash search");
        assert_eq!(
            parallel, serial,
            "worker count {workers} changed the result"
        );
    }
}

fn bonded(source: &str, pairs: &[(u32, u32)]) -> Structure {
    let mut data = structure(source).data().clone();
    let mut bonds = BondTableBuilder::new();
    for &(a, b) in pairs {
        bonds.push(BondRecord {
            atom_a: AtomIndex::new(a),
            atom_b: AtomIndex::new(b),
            order: BondOrder::Single,
            provenance: BondProvenance::User,
            kind: OptionalSymbol::NONE,
        });
    }
    data.bonds = bonds.finish();
    Structure::new(data)
}

fn pairs(structure: &Structure, options: ClashOptions) -> Vec<(u32, u32)> {
    let Ok(found) = clashes_with(structure, options, &ExecutionContext::default()) else {
        panic!("valid");
    };
    found
        .iter()
        .map(|clash| (clash.first.get(), clash.second.get()))
        .collect()
}

fn defaults() -> ClashOptions {
    ClashOptions::new(0.4, RadiusSet::Bondi, SpatialBackend::BruteForce)
}

#[test]
fn an_ideal_tetrahedral_carbon_chain_has_no_one_three_clash() {
    // C-C 1.54 Å, C-C-C 109.5°: the 1-3 distance is 2.51 Å, deep inside the
    // summed radii, yet it is fixed by the bonds and must not be reported.
    let source = format!(
        "{HEADER}\
ATOM 1 C C1 LIG A 1 0 0 0\n\
ATOM 2 C C2 LIG A 1 1.54 0 0\n\
ATOM 3 C C3 LIG A 1 2.0539 1.4516 0\n"
    );
    assert!(pairs(&bonded(&source, &[(0, 1), (1, 2)]), defaults()).is_empty());
}

#[test]
fn one_four_pairs_are_excluded_by_default_and_reported_on_request() {
    let source = format!(
        "{HEADER}\
ATOM 1 C C1 LIG A 1 0 0 0\n\
ATOM 2 C C2 LIG A 1 1.54 0 0\n\
ATOM 3 C C3 LIG A 1 2.054 1.452 0\n\
ATOM 4 C C4 LIG A 1 0.857 2.421 0\n"
    );
    // A planar cis chain: bond lengths 1.54 Å and angles 109.5°, so the 1-4
    // distance is 2.57 Å, too long to be perceived as a bond yet overlapping.
    let structure = bonded(&source, &[(0, 1), (1, 2), (2, 3)]);
    assert!(pairs(&structure, defaults()).is_empty());
    let reported = ClashOptions {
        exclude_one_four: false,
        ..defaults()
    };
    assert_eq!(pairs(&structure, reported), vec![(0, 3)]);
}

/// Heavy atoms of residues 1-4 of 1UBQ (Met-Gln-Ile-Phe) exactly as deposited,
/// with no bond table: what a coordinate file gives the reader. The deposited
/// coordinates are the independent reference, since a crystal structure of a
/// well-refined protein has no real overlap between covalent neighbours.
const UBIQUITIN_N_TERMINUS: &str = "\
ATOM 1 N N MET A 1 27.340 24.430 2.614\n\
ATOM 2 C CA MET A 1 26.266 25.413 2.842\n\
ATOM 3 C C MET A 1 26.913 26.639 3.531\n\
ATOM 4 O O MET A 1 27.886 26.463 4.263\n\
ATOM 5 C CB MET A 1 25.112 24.880 3.649\n\
ATOM 6 C CG MET A 1 25.353 24.860 5.134\n\
ATOM 7 S SD MET A 1 23.930 23.959 5.904\n\
ATOM 8 C CE MET A 1 24.447 23.984 7.620\n\
ATOM 9 N N GLN A 2 26.335 27.770 3.258\n\
ATOM 10 C CA GLN A 2 26.850 29.021 3.898\n\
ATOM 11 C C GLN A 2 26.100 29.253 5.202\n\
ATOM 12 O O GLN A 2 24.865 29.024 5.330\n\
ATOM 13 C CB GLN A 2 26.733 30.148 2.905\n\
ATOM 14 C CG GLN A 2 26.882 31.546 3.409\n\
ATOM 15 C CD GLN A 2 26.786 32.562 2.270\n\
ATOM 16 O OE1 GLN A 2 27.783 33.160 1.870\n\
ATOM 17 N NE2 GLN A 2 25.562 32.733 1.806\n\
ATOM 18 N N ILE A 3 26.849 29.656 6.217\n\
ATOM 19 C CA ILE A 3 26.235 30.058 7.497\n\
ATOM 20 C C ILE A 3 26.882 31.428 7.862\n\
ATOM 21 O O ILE A 3 27.906 31.711 7.264\n\
ATOM 22 C CB ILE A 3 26.344 29.050 8.645\n\
ATOM 23 C CG1 ILE A 3 27.810 28.748 8.999\n\
ATOM 24 C CG2 ILE A 3 25.491 27.771 8.287\n\
ATOM 25 C CD1 ILE A 3 27.967 28.087 10.417\n\
ATOM 26 N N PHE A 4 26.214 32.097 8.771\n\
ATOM 27 C CA PHE A 4 26.772 33.436 9.197\n\
ATOM 28 C C PHE A 4 27.151 33.362 10.650\n\
ATOM 29 O O PHE A 4 26.350 32.778 11.395\n\
ATOM 30 C CB PHE A 4 25.695 34.498 8.946\n\
ATOM 31 C CG PHE A 4 25.288 34.609 7.499\n\
ATOM 32 C CD1 PHE A 4 24.147 33.966 7.038\n\
ATOM 33 C CD2 PHE A 4 26.136 35.346 6.640\n\
ATOM 34 C CE1 PHE A 4 23.812 34.031 5.677\n\
ATOM 35 C CE2 PHE A 4 25.810 35.392 5.267\n\
ATOM 36 C CZ PHE A 4 24.620 34.778 4.853\n";

#[test]
fn a_crystal_peptide_without_a_bond_table_has_no_clashes() {
    // Without bond information, 1-2, 1-3 and 1-4 neighbours inside residues, the
    // peptide links between them and the ring carbons of the phenylalanine all
    // sit well inside their summed van der Waals radii.
    let structure = structure(&format!("{HEADER}{UBIQUITIN_N_TERMINUS}"));
    assert!(pairs(&structure, defaults()).is_empty());
}

#[test]
fn without_perceived_bonds_the_same_peptide_is_flooded_with_covalent_pairs() {
    let structure = structure(&format!("{HEADER}{UBIQUITIN_N_TERMINUS}"));
    let raw = ClashOptions {
        perceive_bonds: false,
        ..defaults()
    };
    let found = pairs(&structure, raw);
    // 1-2 pairs such as N-CA (1.46 Å) and the peptide C-N link are in the list.
    assert!(found.contains(&(0, 1)), "{found:?}");
    assert!(found.contains(&(2, 8)), "{found:?}");
    assert!(found.len() > 30, "{}", found.len());
}

#[test]
fn a_phenylalanine_ring_needs_the_one_four_exclusion() {
    // CD1 and CE2 (atoms 32 and 35, zero-based 31 and 34) are 1-4 across the
    // ring and 2.4 Å apart, far inside the 3.4 Å carbon contact distance.
    let structure = structure(&format!("{HEADER}{UBIQUITIN_N_TERMINUS}"));
    let keeps_one_four = ClashOptions {
        exclude_one_four: false,
        ..defaults()
    };
    let found = pairs(&structure, keeps_one_four);
    assert!(found.contains(&(31, 34)), "{found:?}");
    assert!(pairs(&structure, defaults()).is_empty());
}

#[test]
fn a_real_overlap_between_non_neighbours_is_still_reported() {
    // A water oxygen placed 2.0 Å from the Met1 CB, which is 1.4 Å inside the
    // 3.22 Å C...O contact distance, with nothing bonded to it.
    let intruder = "ATOM 37 O O HOH B 5 25.112 24.880 5.649\n";
    let structure = structure(&format!("{HEADER}{UBIQUITIN_N_TERMINUS}{intruder}"));
    let found = pairs(&structure, defaults());
    assert!(found.contains(&(4, 36)), "{found:?}");
    assert!(
        found
            .iter()
            .all(|&(first, second)| first == 36 || second == 36)
    );
}

#[test]
fn a_close_contact_between_unbonded_residues_is_not_mistaken_for_a_bond() {
    // Two residues whose carbons are 1.6 Å apart sit inside the distance window
    // a bond would, but nothing joins them: perception must not erase the very
    // contact this check exists to find.
    let source = format!(
        "{HEADER}\
ATOM 1 C CA ALA A 1 0 0 0\n\
ATOM 2 C CA ALA A 2 1.6 0 0\n"
    );
    assert_eq!(pairs(&structure(&source), defaults()), vec![(0, 1)]);
}

#[test]
fn atoms_in_different_alternate_locations_do_not_clash() {
    let header = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.label_alt_id\n\
_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";
    let run = |a: &str, b: &str| {
        let source =
            format!("{header}ATOM 1 C C1 LIG A 1 {a} 0 0 0\nATOM 2 C C2 LIG A 2 {b} 2 0 0\n");
        pairs(&structure(&source), defaults())
    };
    assert!(run("A", "B").is_empty());
    assert_eq!(run("A", "A"), vec![(0, 1)]);
    assert_eq!(run(".", "B"), vec![(0, 1)]);
}
