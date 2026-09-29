use molframe_core::io::{InputBuffer, ReadOptions};
use num_traits::ToPrimitive;

fn structure_of(source: &str) -> molframe_core::structure::Structure {
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let (structure, _) = crate::read(&input, &ReadOptions::new()).expect("fixture reads");
    structure
}

fn atom_site_header() -> String {
    String::from(
        "loop_\n\
         _atom_site.group_PDB\n\
         _atom_site.id\n\
         _atom_site.type_symbol\n\
         _atom_site.label_atom_id\n\
         _atom_site.label_comp_id\n\
         _atom_site.label_asym_id\n\
         _atom_site.label_seq_id\n\
         _atom_site.Cartn_x\n\
         _atom_site.Cartn_y\n\
         _atom_site.Cartn_z\n\
         _atom_site.occupancy\n\
         _atom_site.B_iso_or_equiv\n",
    )
}

/// Asserts one stored tensor against expected components, within rounding.
fn assert_tensor(actual: [f32; 6], expected: [f32; 6]) {
    for (actual, expected) in actual.iter().zip(expected) {
        assert!(
            (actual - expected).abs() <= 1e-6 * expected.abs().max(1.0),
            "{actual:?} against {expected:?}"
        );
    }
}

#[test]
fn a_u_form_category_joins_rows_by_atom_site_id() {
    let source = format!(
        "data_u_form\n\
         {}\n\
         ATOM 1 N N ALA A 1 0 0 0 1.0 20.0\n\
         ATOM 2 C CA ALA A 1 1.5 0 0 1.0 20.0\n\
         #\n\
         loop_\n\
         _atom_site_anisotrop.id\n\
         _atom_site_anisotrop.type_symbol\n\
         _atom_site_anisotrop.U[1][1]\n\
         _atom_site_anisotrop.U[2][2]\n\
         _atom_site_anisotrop.U[3][3]\n\
         _atom_site_anisotrop.U[1][2]\n\
         _atom_site_anisotrop.U[1][3]\n\
         _atom_site_anisotrop.U[2][3]\n\
         1 N 0.030 0.041 0.024 0.002 -0.003 0.005\n\
         2 C 0.050 0.061 0.044 0.012 -0.013 0.015\n\
         #\n",
        atom_site_header()
    );
    let structure = structure_of(&source);
    let anisotropy = &structure.data().anisotropy;
    assert!(anisotropy.is_available());
    assert_eq!(anisotropy.len(), 2);
    let Some(record) = anisotropy.get(molframe_core::AnisotropyIndex::new(0)) else {
        panic!("first ellipsoid was lost");
    };
    assert_eq!(record.atom, molframe_core::AtomIndex::new(0));
    assert_tensor(record.u, [0.030, 0.041, 0.024, 0.002, -0.003, 0.005]);
    let Some(second) = anisotropy.get(molframe_core::AnisotropyIndex::new(1)) else {
        panic!("second ellipsoid was lost");
    };
    assert_eq!(second.atom, molframe_core::AtomIndex::new(1));
    assert_tensor(second.u, [0.050, 0.061, 0.044, 0.012, -0.013, 0.015]);
}

#[test]
fn a_reversed_anisotrop_row_still_reaches_its_atom_by_id() {
    let source = format!(
        "data_reversed\n\
         {}\n\
         ATOM 1 N N ALA A 1 0 0 0 1.0 20.0\n\
         ATOM 2 C CA ALA A 1 1.5 0 0 1.0 20.0\n\
         #\n\
         loop_\n\
         _atom_site_anisotrop.id\n\
         _atom_site_anisotrop.U[1][1]\n\
         _atom_site_anisotrop.U[2][2]\n\
         _atom_site_anisotrop.U[3][3]\n\
         _atom_site_anisotrop.U[1][2]\n\
         _atom_site_anisotrop.U[1][3]\n\
         _atom_site_anisotrop.U[2][3]\n\
         2 0.050 0.061 0.044 0.012 -0.013 0.015\n\
         1 0.030 0.041 0.024 0.002 -0.003 0.005\n\
         #\n",
        atom_site_header()
    );
    let structure = structure_of(&source);
    let anisotropy = &structure.data().anisotropy;
    // Rows land in atom order regardless of file order; the id, not the row
    // position, decides which atom each tensor belongs to.
    let Some(first) = anisotropy.get(molframe_core::AnisotropyIndex::new(0)) else {
        panic!("first row was lost");
    };
    assert_eq!(first.atom, molframe_core::AtomIndex::new(0));
    assert!((first.u[0] - 0.030).abs() <= 1e-6);
    let Some(second) = anisotropy.get(molframe_core::AnisotropyIndex::new(1)) else {
        panic!("second row was lost");
    };
    assert_eq!(second.atom, molframe_core::AtomIndex::new(1));
    assert!((second.u[0] - 0.050).abs() <= 1e-6);
}

#[test]
fn a_b_form_category_converts_with_the_eight_pi_squared_factor() {
    let source = format!(
        "data_b_form\n\
         {}\n\
         ATOM 1 N N ALA A 1 0 0 0 1.0 20.0\n\
         #\n\
         loop_\n\
         _atom_site_anisotrop.id\n\
         _atom_site_anisotrop.b[1][1]\n\
         _atom_site_anisotrop.b[2][2]\n\
         _atom_site_anisotrop.b[3][3]\n\
         _atom_site_anisotrop.b[1][2]\n\
         _atom_site_anisotrop.b[1][3]\n\
         _atom_site_anisotrop.b[2][3]\n\
         1 1.0 2.0 3.0 4.0 5.0 6.0\n\
         #\n",
        atom_site_header()
    );
    let structure = structure_of(&source);
    let anisotropy = &structure.data().anisotropy;
    let Some(record) = anisotropy.get(molframe_core::AnisotropyIndex::new(0)) else {
        panic!("converted ellipsoid was lost");
    };
    let factor = 0.012_665_147_955_292_222_f64;
    let expected: [f32; 6] = [1.0_f32, 2.0, 3.0, 4.0, 5.0, 6.0]
        .iter()
        .map(|value| f64::from(*value).mul_add(factor, 0.0))
        .filter_map(|converted| converted.to_f32())
        .collect::<Vec<_>>()
        .try_into()
        .expect("six components");
    assert_tensor(record.u, expected);
    // The conversion scaled rather than passing the values through.
    assert!((record.u[0] - 1.0).abs() > f32::EPSILON);
}

#[test]
fn an_absent_category_leaves_the_table_unavailable() {
    let source = format!(
        "data_no_aniso\n\
         {}\n\
         ATOM 1 N N ALA A 1 0 0 0 1.0 20.0\n\
         #\n",
        atom_site_header()
    );
    let structure = structure_of(&source);
    let anisotropy = &structure.data().anisotropy;
    assert_eq!(anisotropy.len(), 0);
    assert!(!anisotropy.is_available());
}

#[test]
fn the_inline_aniso_spelling_is_accepted() {
    let source = "data_inline\n\
                  loop_\n\
                  _atom_site.group_PDB\n\
                  _atom_site.id\n\
                  _atom_site.type_symbol\n\
                  _atom_site.label_atom_id\n\
                  _atom_site.label_comp_id\n\
                  _atom_site.label_asym_id\n\
                  _atom_site.label_seq_id\n\
                  _atom_site.aniso_U[1][1]\n\
                  _atom_site.aniso_U[2][2]\n\
                  _atom_site.aniso_U[3][3]\n\
                  _atom_site.aniso_U[1][2]\n\
                  _atom_site.aniso_U[1][3]\n\
                  _atom_site.aniso_U[2][3]\n\
                  _atom_site.Cartn_x\n\
                  _atom_site.Cartn_y\n\
                  _atom_site.Cartn_z\n\
                  _atom_site.occupancy\n\
                  _atom_site.B_iso_or_equiv\n\
                  ATOM 1 N N ALA A 1 0.030 0.041 0.024 0.002 -0.003 0.005 0 0 0 1.0 20.0\n\
                  #\n";
    let structure = structure_of(source);
    let anisotropy = &structure.data().anisotropy;
    assert!(anisotropy.is_available());
    assert_eq!(anisotropy.len(), 1);
    let Some(record) = anisotropy.get(molframe_core::AnisotropyIndex::new(0)) else {
        panic!("inline ellipsoid was lost");
    };
    assert_eq!(record.atom, molframe_core::AtomIndex::new(0));
    assert_tensor(record.u, [0.030, 0.041, 0.024, 0.002, -0.003, 0.005]);
}

#[test]
fn rows_named_by_unknown_identifiers_are_left_unattached() {
    let source = format!(
        "data_unknown_id\n\
         {}\n\
         ATOM 1 N N ALA A 1 0 0 0 1.0 20.0\n\
         #\n\
         loop_\n\
         _atom_site_anisotrop.id\n\
         _atom_site_anisotrop.U[1][1]\n\
         _atom_site_anisotrop.U[2][2]\n\
         _atom_site_anisotrop.U[3][3]\n\
         _atom_site_anisotrop.U[1][2]\n\
         _atom_site_anisotrop.U[1][3]\n\
         _atom_site_anisotrop.U[2][3]\n\
         9 0.030 0.041 0.024 0.002 -0.003 0.005\n\
         #\n",
        atom_site_header()
    );
    let structure = structure_of(&source);
    let anisotropy = &structure.data().anisotropy;
    assert!(anisotropy.is_available());
    assert_eq!(anisotropy.len(), 0);
}

#[test]
fn legacy_u_11_item_names_are_accepted() {
    let source = format!(
        "data_legacy_names\n\
         {}\n\
         ATOM 1 N N ALA A 1 0 0 0 1.0 20.0\n\
         #\n\
         loop_\n\
         _atom_site_anisotrop.id\n\
         _atom_site_anisotrop.U_11\n\
         _atom_site_anisotrop.U_22\n\
         _atom_site_anisotrop.U_33\n\
         _atom_site_anisotrop.U_12\n\
         _atom_site_anisotrop.U_13\n\
         _atom_site_anisotrop.U_23\n\
         1 0.030 0.041 0.024 0.002 -0.003 0.005\n\
         #\n",
        atom_site_header()
    );
    let structure = structure_of(&source);
    let anisotropy = &structure.data().anisotropy;
    assert!(anisotropy.is_available());
    let Some(record) = anisotropy.get(molframe_core::AnisotropyIndex::new(0)) else {
        panic!("legacy-spelled ellipsoid was lost");
    };
    assert_tensor(record.u, [0.030, 0.041, 0.024, 0.002, -0.003, 0.005]);
}
