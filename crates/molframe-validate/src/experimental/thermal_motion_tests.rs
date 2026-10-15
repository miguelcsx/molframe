use super::*;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::{AtomSelection, Structure};

fn structure() -> Structure {
    let cif = b"data_b
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
_atom_site.B_iso_or_equiv
ATOM 1 C A1 LIG A 1 0 0 0 10
ATOM 2 C A2 LIG A 1 1 0 0 30
";
    let input = InputBuffer::from_bytes(cif.to_vec());
    molframe_cif::read(&input, &ReadOptions::new()).map_or_else(
        |findings| panic!("fixture failed: {findings:?}"),
        |result| result.0,
    )
}

#[test]
fn distribution_and_declared_tls_are_assessed_without_defaults() {
    let structure = structure();
    let selection = AtomSelection::from_sorted(vec![0, 1]);
    let distribution = b_factor_distribution(&structure, &selection, 0.5)
        .unwrap_or_else(|error| panic!("distribution failed: {error}"));
    assert!((distribution.mean - 20.0).abs() < 1e-12);
    assert_eq!(distribution.outliers.len(), 2);
    let translation = 10.0 / (8.0 * core::f64::consts::PI.powi(2));
    let report = tls_b_factor_consistency(
        &structure,
        &[TlsGroup {
            id: "declared".into(),
            atoms: AtomSelection::from_sorted(vec![0]),
            model: TlsModel {
                origin: [0.0; 3],
                translation: [
                    [translation, 0.0, 0.0],
                    [0.0, translation, 0.0],
                    [0.0, 0.0, translation],
                ],
                libration: [[0.0; 3]; 3],
                screw: [[0.0; 3]; 3],
            },
        }],
        1e-6,
        1e-12,
    )
    .unwrap_or_else(|error| panic!("TLS failed: {error}"));
    assert_eq!(report.assessed, 1);
    assert!(report.flags.is_empty());
}
