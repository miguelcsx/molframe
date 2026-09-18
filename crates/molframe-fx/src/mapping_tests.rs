use super::map_motif;
use crate::{AtomSite, ComponentRole, ComponentSpec, Constraint, Motif, NamedConstraint};
use molframe_core::contract::AnalysisPolicy;
use molframe_core::io::{InputBuffer, ReadOptions};

#[test]
fn maps_components_and_explicitly_equivalent_atoms_without_name_guessing() {
    let source = "data_s\nloop_\n_atom_site.group_PDB\n_atom_site.id\n\
_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 O OG SER A 1 0 0 0\nATOM 2 O OD2 ASP A 2 2.8 0 0\n";
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let (structure, _) = molframe_cif::read(&input, &ReadOptions::new())
        .unwrap_or_else(|findings| panic!("fixture failed: {findings:?}"));
    let motif = Motif::new(
        [
            (
                "ser".into(),
                ComponentSpec::new(ComponentRole::Residue)
                    .component("SER")
                    .require("OG"),
            ),
            (
                "asp".into(),
                ComponentSpec::new(ComponentRole::Residue)
                    .component("ASP")
                    .equivalent(["OD1", "OD2"]),
            ),
        ],
        [NamedConstraint {
            name: "dyad".into(),
            constraint: Constraint::Distance {
                first: AtomSite::new("ser", "OG"),
                second: AtomSite::new("asp", "OD1"),
                target: 2.8,
                tolerance: 0.4,
            },
        }],
    )
    .unwrap_or_else(|error| panic!("motif failed: {error}"));
    let mappings = map_motif(&structure, &motif, None, &AnalysisPolicy::default(), 16)
        .unwrap_or_else(|error| panic!("mapping failed: {error}"));
    assert_eq!(mappings.mappings.len(), 1);
    assert_eq!(
        mappings.mappings[0]
            .atoms
            .get(&AtomSite::new("asp", "OD1"))
            .map(Vec::as_slice),
        Some([molframe_core::index::AtomIndex::new(1)].as_slice())
    );
}

#[test]
fn roles_are_not_inferred_from_atom_record_types() {
    let source = "data_s\nloop_\n_atom_site.group_PDB\n_atom_site.id\n\
_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 O OG SER A 1 0 0 0\n";
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let (structure, _) = molframe_cif::read(&input, &ReadOptions::new())
        .unwrap_or_else(|findings| panic!("fixture failed: {findings:?}"));
    let motif = Motif::new(
        [(
            "site".into(),
            ComponentSpec::new(ComponentRole::Residue).require("OG"),
        )],
        [],
    )
    .unwrap_or_else(|error| panic!("motif failed: {error}"));

    let mappings = map_motif(&structure, &motif, None, &AnalysisPolicy::default(), 16)
        .unwrap_or_else(|error| panic!("mapping failed: {error}"));

    assert!(mappings.mappings.is_empty());
}
