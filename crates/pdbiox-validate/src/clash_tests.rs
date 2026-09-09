use super::clashes;
use pdbiox_chem::RadiusSet;
use pdbiox_core::ExecutionContext;
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_core::structure::Structure;
use pdbiox_spatial::SpatialBackend;
use std::fmt::Write;

const HEADER: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

fn structure(source: &str) -> Structure {
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    match pdbiox_cif::read(&input, &ReadOptions::new()) {
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
    assert_eq!((found[0].first.get(), found[0].second.get()), (0, 1));
    assert!(found[0].overlap > 0.4);
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
