use super::internal_coordinates;
use molframe_core::{InputBuffer, ModelIndex, ReadOptions};

const SOURCE: &str = "data_ic\n\
loop_\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n\
_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n\
_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
1 C A1 LIG A 1 0 0 0\n2 C A2 LIG A 1 1 0 0\n3 C A3 LIG A 1 1 1 0\n\
4 C A4 LIG A 1 1 1 1\n5 C A5 LIG A 1 2 1 1\n\
loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n\
_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n\
_struct_conn.ptnr1_label_comp_id\n_struct_conn.ptnr1_label_atom_id\n\
_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n\
_struct_conn.ptnr2_label_comp_id\n_struct_conn.ptnr2_label_atom_id\n\
_struct_conn.pdbx_value_order\n\
1 covale A 1 LIG A1 A 1 LIG A2 SING\n\
2 covale A 1 LIG A2 A 1 LIG A3 SING\n\
3 covale A 1 LIG A3 A 1 LIG A4 SING\n\
4 covale A 1 LIG A4 A 1 LIG A5 SING\n";

#[test]
fn a_bonded_model_round_trips_through_internal_coordinates() {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    let structure = match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    let internal = match internal_coordinates(&structure, ModelIndex::new(0)) {
        Ok(internal) => internal,
        Err(finding) => panic!("conversion failed: {finding}"),
    };
    assert_eq!(internal.seeds().len(), 3);
    assert_eq!(internal.atoms().len(), 2);
    let rebuilt = match internal.rebuild() {
        Ok(rebuilt) => rebuilt,
        Err(finding) => panic!("rebuild failed: {finding}"),
    };
    for (actual, expected) in rebuilt.iter().zip(structure.positions()) {
        let Some(actual) = actual else {
            panic!("coordinate absent")
        };
        assert!(
            actual
                .iter()
                .zip(expected)
                .all(|(a, b)| (*a - *b).abs() < 1e-5)
        );
    }
}

#[test]
fn bat_topology_measures_multiple_frames_and_rebuilds_them() {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    let structure = match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    let internal = internal_coordinates(&structure, ModelIndex::new(0))
        .unwrap_or_else(|error| panic!("internal coordinates failed: {error}"));
    let shifted = structure
        .model_positions(ModelIndex::new(0))
        .unwrap_or_else(|| panic!("model absent"))
        .iter()
        .map(|position| Some([position[0] + 7.0, position[1] - 2.0, position[2] + 1.0]))
        .collect::<Vec<_>>();
    let bat = internal
        .measure_bat(&shifted)
        .unwrap_or_else(|error| panic!("BAT measurement failed: {error}"));
    assert_eq!(bat.coordinates().len(), internal.atoms().len());
    let rebuilt = internal
        .rebuild_bat(&bat)
        .unwrap_or_else(|error| panic!("BAT rebuild failed: {error}"));
    for (observed, expected) in rebuilt.iter().zip(shifted) {
        let observed = observed.unwrap_or_else(|| panic!("rebuilt atom absent"));
        let expected = expected.unwrap_or_else(|| panic!("source atom absent"));
        assert!(
            observed
                .iter()
                .zip(expected)
                .all(|(left, right)| (left - right).abs() < 1.0e-5)
        );
    }
}
