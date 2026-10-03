use super::*;
use crate::bindings::PyStructure;
use crate::query_messages::QueryError;

const SOURCE: &str = "data_hydrogen
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
HETATM 1 H H UNK A 1 0 0 0
";

#[test]
fn python_hydrogen_selectors_report_unavailable_bonds_as_query_errors() {
    let input = molframe::InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    let (structure, _) = molframe::read_buffer(
        &input,
        Some("hydrogen.cif"),
        &molframe::ReadOptions::new().only_atomic_coords(true),
    )
    .expect("coordinate-only read parses without bond enrichment");
    assert!(!structure.bonds().is_available());
    assert_eq!(structure.atom_count(), 1);
    Python::initialize();
    Python::attach(|py| {
        let structure = Bound::new(py, PyStructure::new(structure)).expect("structure binds");
        let module = PyModule::new(py, "sel").expect("module exists");
        register(&module).expect("selectors register");
        for name in ["polar_hydrogen", "nonpolar_hydrogen"] {
            let query = module
                .getattr(name)
                .expect("selector exists")
                .call0()
                .expect("query constructs");
            let error = query
                .call_method1("select", (&structure,))
                .expect_err("unavailable bonds fail through Python");
            assert!(error.is_instance_of::<QueryError>(py));
            assert!(error.to_string().contains("MOLFRAME-E4003"));
            assert!(error.to_string().contains("bond topology"));
        }
    });
}
