use super::{PyDihedron, PyHedron, PyInternalAtom, PyInternalCoordinates, place_atom};

#[test]
fn placement_and_geometry_records_are_native_projections() {
    let i = [0.0, 0.0, 0.0];
    let j = [1.0, 0.0, 0.0];
    let k = [1.0, 1.0, 0.0];
    let l = place_atom(i, j, k, 1.5, 1.9, 0.7).expect("non-degenerate placement");
    let hedron = PyHedron::from_points(i, j, k).expect("non-degenerate hedron");
    let dihedron = PyDihedron::from_points(i, j, k, l).expect("non-degenerate dihedron");
    assert!((hedron.first_length - 1.0).abs() < 1.0e-6);
    assert!((dihedron.length - 1.5).abs() < 1.0e-5);
}

#[test]
fn internal_records_convert_typed_indices_at_the_boundary() {
    let value = pdbiox::InternalAtom {
        atom: pdbiox::AtomIndex::new(4),
        references: [
            pdbiox::AtomIndex::new(0),
            pdbiox::AtomIndex::new(1),
            pdbiox::AtomIndex::new(2),
        ],
        coordinate: pdbiox::Dihedron {
            first: pdbiox::Hedron {
                first_length: 1.0,
                angle: 1.1,
                second_length: 1.2,
            },
            angle: 1.3,
            length: 1.4,
            torsion: 1.5,
        },
    };
    let projected: PyInternalAtom = value.into();
    assert_eq!(projected.atom, 4);
    assert_eq!(projected.references, [0, 1, 2]);
}

#[test]
fn internal_coordinate_wrapper_is_owned() {
    let structure = pdbiox::Structure::new(pdbiox::StructureData::empty());
    let result = pdbiox::internal_coordinates(&structure, pdbiox::ModelIndex::new(0));
    let inner = result.expect("empty structure has a deterministic empty frame");
    let type_check: fn(pdbiox::InternalCoordinates) -> PyInternalCoordinates =
        |inner| PyInternalCoordinates { inner };
    let projected = type_check(inner);
    assert!(projected.inner.seeds().is_empty());
}
