use super::{MeasurementOptions, MeasurementValue, measure_constraints};
use crate::{
    AtomSite, ComponentRole, ComponentSpec, Constraint, Motif, NamedConstraint, align_intrinsic,
    map_motif,
};
use pdbiox_core::contract::AnalysisPolicy;
use pdbiox_core::io::{InputBuffer, ReadOptions};

fn fixture() -> pdbiox_core::structure::Structure {
    let source = "data_s\nloop_\n_atom_site.group_PDB\n_atom_site.id\n\
_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 O OG SER A 1 0 0 0\nATOM 2 O OD1 ASP A 2 4 0 0\n\
ATOM 3 O OD2 ASP A 2 2.8 0 0\n";
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    pdbiox_cif::read(&input, &ReadOptions::new()).map_or_else(
        |findings| panic!("fixture failed: {findings:?}"),
        |(structure, _)| structure,
    )
}

fn motif() -> Motif {
    Motif::new(
        [
            (
                "ser".into(),
                ComponentSpec::new(ComponentRole::Residue).component("SER"),
            ),
            (
                "asp".into(),
                ComponentSpec::new(ComponentRole::Residue)
                    .component("ASP")
                    .equivalent(["OD1", "OD2"]),
            ),
        ],
        [NamedConstraint {
            name: "dyad_distance".into(),
            constraint: Constraint::Distance {
                first: AtomSite::new("ser", "OG"),
                second: AtomSite::new("asp", "OD1"),
                target: 2.8,
                tolerance: 0.1,
            },
        }],
    )
    .unwrap_or_else(|error| panic!("motif failed: {error}"))
}

fn options(maximum_alternatives: usize) -> MeasurementOptions {
    MeasurementOptions {
        maximum_alternatives,
        plane_fit: pdbiox_geom::EigenOptions::standard(),
    }
}

#[test]
fn equivalent_atom_choice_is_decomposed_and_deterministic() {
    let structure = fixture();
    let motif = motif();
    let mappings = map_motif(&structure, &motif, None, &AnalysisPolicy::default(), 8)
        .unwrap_or_else(|error| panic!("mapping failed: {error}"));
    let Some(mapping) = mappings.mappings.into_iter().next() else {
        panic!("one complete mapping expected");
    };
    let measured = measure_constraints(&structure, &motif, &align_intrinsic(mapping), options(8))
        .unwrap_or_else(|error| panic!("measurement failed: {error}"));
    let result = &measured.constraints[0];
    assert!(matches!(
        result.value,
        Some(MeasurementValue::Scalar(value)) if (value - 2.8).abs() < 1.0e-6
    ));
    assert_eq!(result.satisfied, Some(true));
    assert_eq!(result.atoms[1], pdbiox_core::index::AtomIndex::new(2));
    assert!(
        measured
            .metrics()
            .get("dyad_distance")
            .is_some_and(|value| (*value - 2.8).abs() < 1.0e-6)
    );
}

#[test]
fn measurement_limit_is_enforced_before_enumeration() {
    let structure = fixture();
    let motif = motif();
    let mappings = map_motif(&structure, &motif, None, &AnalysisPolicy::default(), 8)
        .unwrap_or_else(|error| panic!("mapping failed: {error}"));
    let Some(mapping) = mappings.mappings.into_iter().next() else {
        panic!("one complete mapping expected");
    };
    assert!(
        measure_constraints(&structure, &motif, &align_intrinsic(mapping), options(1)).is_err()
    );
}
