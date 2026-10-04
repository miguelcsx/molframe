use super::{chain_atom_ranges, crystal_contact_view, operation_motion};
use crate::{CellTransform, SYMMETRY_EXTENSION, SymmetryExt, SymmetryOperation, lower_symmetry};
use molframe_bench::{Sample, structure_from_cif};
use molframe_cif::parse;
use molframe_core::{Code, ExecutionContext, InputBuffer, ModelIndex, Structure};
use molframe_geom::Rigid;

/// 1CRN as deposited, with the space-group operators its own file carries.
fn crambin() -> Structure {
    let Some(bytes) = Sample::Tiny.cif() else {
        panic!("1CRN ships an mmCIF");
    };
    let input = InputBuffer::from_bytes(bytes.to_vec());
    let document = match parse(&input) {
        Ok((document, _)) => document,
        Err(findings) => panic!("parse failed: {findings:?}"),
    };
    let symmetry = match lower_symmetry(&document) {
        Ok(symmetry) => symmetry,
        Err(findings) => panic!("symmetry failed: {findings:?}"),
    };
    let Some(structure) = structure_from_cif(Sample::Tiny) else {
        panic!("1CRN reads");
    };
    structure.with_extension(SYMMETRY_EXTENSION, symmetry)
}

fn transform_of(structure: &Structure) -> CellTransform {
    let Some(cell) = structure.data().cell else {
        panic!("1CRN has a cell");
    };
    match CellTransform::new(&cell) {
        Ok(transform) => transform,
        Err(finding) => panic!("cell: {finding}"),
    }
}

fn same(left: &Rigid, right: &Rigid) -> bool {
    let (a, b) = (left.apply([1.0, 2.0, 3.0]), right.apply([1.0, 2.0, 3.0]));
    let (c, d) = (left.apply([-4.0, 0.5, 9.0]), right.apply([-4.0, 0.5, 9.0]));
    (0..3).all(|axis| (a[axis] - b[axis]).abs() < 1e-3 && (c[axis] - d[axis]).abs() < 1e-3)
}

#[test]
fn a_symmetry_operation_moves_atoms_as_its_fractional_expression_says() {
    let structure = crambin();
    let transform = transform_of(&structure);
    let Some(symmetry) = structure.symmetry_set() else {
        panic!("1CRN carries operators");
    };
    let point = [5.5_f32, -3.25, 12.0];
    for operation in symmetry.operations() {
        for lattice in [[0, 0, 0], [1, -1, 0], [-2, 0, 1]] {
            let Ok(motion) = operation_motion(&transform, operation, lattice) else {
                panic!("proper operation");
            };
            let fractional = operation.apply_fractional(transform.to_fractional([
                f64::from(point[0]),
                f64::from(point[1]),
                f64::from(point[2]),
            ]));
            let expected = transform.to_cartesian([
                fractional[0] + f64::from(lattice[0]),
                fractional[1] + f64::from(lattice[1]),
                fractional[2] + f64::from(lattice[2]),
            ]);
            let found = motion.apply(point);
            for axis in 0..3 {
                assert!((f64::from(found[axis]) - expected[axis]).abs() < 1e-3);
            }
        }
    }
}

#[test]
fn an_improper_operation_is_refused_rather_than_mirroring_the_molecule() {
    let structure = crambin();
    let transform = transform_of(&structure);
    let Ok(inversion) = SymmetryOperation::parse("inv", "-x,-y,-z") else {
        panic!("a valid expression");
    };
    let refused = operation_motion(&transform, &inversion, [0, 0, 0]);
    assert!(matches!(refused, Err(finding) if finding.code() == Code::E6015));
}

#[test]
fn the_contacts_of_crambin_are_exactly_the_images_with_an_atom_inside_the_radius() {
    let structure = crambin();
    let transform = transform_of(&structure);
    let Some(symmetry) = structure.symmetry_set() else {
        panic!("1CRN carries operators");
    };
    let radius = 4.0_f32;
    let Ok(view) = crystal_contact_view(&structure, radius, &ExecutionContext::default()) else {
        panic!("crambin has crystal contacts");
    };
    let chains = chain_atom_ranges(&structure).len();
    let motions: Vec<Rigid> = view.chains().map(|instance| instance.transform).collect();
    assert_eq!(motions.len() % chains, 0);
    assert!(
        same(&motions[0], &Rigid::IDENTITY),
        "the deposited unit comes first"
    );
    assert!(motions.len() > chains, "crambin touches its symmetry mates");

    // Brute force: every operation and lattice translation, kept when some image atom
    // lies within the radius of some deposited atom.
    let Some(positions) = structure.model_positions(ModelIndex::new(0)) else {
        panic!("model 0");
    };
    let mut expected: Vec<Rigid> = Vec::new();
    for operation in symmetry.operations() {
        for a in -2..=2 {
            for b in -2..=2 {
                for c in -2..=2 {
                    if operation.is_identity() && [a, b, c] == [0, 0, 0] {
                        continue;
                    }
                    let Ok(motion) = operation_motion(&transform, operation, [a, b, c]) else {
                        panic!("proper operation");
                    };
                    let touches = positions.iter().any(|&image| {
                        let moved = motion.apply(image);
                        positions.iter().any(|&atom| {
                            let d: f32 =
                                (0..3).map(|axis| (moved[axis] - atom[axis]).powi(2)).sum();
                            d.sqrt() <= radius
                        })
                    });
                    if touches {
                        expected.push(motion);
                    }
                }
            }
        }
    }
    let found: Vec<&Rigid> = motions.iter().skip(chains).collect();
    assert_eq!(
        found.len(),
        expected.len(),
        "neither missing nor extra images"
    );
    for motion in &expected {
        assert!(
            found.iter().any(|candidate| same(candidate, motion)),
            "a contact the brute force finds is missing"
        );
    }
}

#[test]
fn a_larger_radius_keeps_every_contact_of_a_smaller_one() {
    let structure = crambin();
    let context = ExecutionContext::default();
    let (Ok(near), Ok(far)) = (
        crystal_contact_view(&structure, 3.5, &context),
        crystal_contact_view(&structure, 6.0, &context),
    ) else {
        panic!("crambin has crystal contacts");
    };
    assert!(far.instance_count() >= near.instance_count());
    let far_motions: Vec<Rigid> = far.chains().map(|instance| instance.transform).collect();
    for instance in near.chains() {
        assert!(
            far_motions
                .iter()
                .any(|motion| same(motion, &instance.transform))
        );
    }
}

#[test]
fn the_view_materialises_into_a_structure_whose_copies_sit_where_the_motions_put_them() {
    let structure = crambin();
    let Ok(view) = crystal_contact_view(&structure, 4.0, &ExecutionContext::default()) else {
        panic!("crambin has crystal contacts");
    };
    let Ok(materialised) = view.materialize() else {
        panic!("a valid expansion");
    };
    assert_eq!(
        materialised.atom_count() as usize,
        view.atoms().count(),
        "one atom per generated atom, in generation order"
    );
    let Some(source) = structure.model_positions(ModelIndex::new(0)) else {
        panic!("model 0");
    };
    let Some(placed) = materialised.model_positions(ModelIndex::new(0)) else {
        panic!("model 0");
    };
    for (position, instance) in placed.iter().zip(view.atoms()) {
        let expected = instance.apply(source[instance.source_atom.as_usize()]);
        for axis in 0..3 {
            assert!((position[axis] - expected[axis]).abs() < 1e-3);
        }
    }
}

#[test]
fn a_structure_without_symmetry_operators_has_no_crystal_contacts_to_report() {
    let Some(structure) = structure_from_cif(Sample::Tiny) else {
        panic!("1CRN reads");
    };
    let refused = crystal_contact_view(&structure, 4.0, &ExecutionContext::default());
    assert!(matches!(refused, Err(finding) if finding.code() == Code::E6016));
}
