use super::*;

#[cfg(feature = "pdb")]
use crate::PdbOptions;

const DIPEPTIDE: &str = "\
ATOM      1  N   GLY A   1      27.340  24.430   2.614  1.00 10.00           N
ATOM      2  CA  GLY A   1      26.266  25.413   2.842  1.00 11.00           C
END
";

#[cfg(all(feature = "mmcif", feature = "xtal"))]
const MMCIF_XTAL: &str = r"data_x
_cell.length_a 10
_cell.length_b 10
_cell.length_c 10
_cell.angle_alpha 90
_cell.angle_beta 90
_cell.angle_gamma 90
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 C CA GLY A 1 1 2 3
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[3]
1 1 0 0 0 0 1 0 0 0 0 1 0
_pdbx_struct_assembly.id 1
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 1 A
loop_
_struct_ncs_oper.id
_struct_ncs_oper.code
_struct_ncs_oper.matrix[1][1]
_struct_ncs_oper.matrix[1][2]
_struct_ncs_oper.matrix[1][3]
_struct_ncs_oper.vector[1]
_struct_ncs_oper.matrix[2][1]
_struct_ncs_oper.matrix[2][2]
_struct_ncs_oper.matrix[2][3]
_struct_ncs_oper.vector[2]
_struct_ncs_oper.matrix[3][1]
_struct_ncs_oper.matrix[3][2]
_struct_ncs_oper.matrix[3][3]
_struct_ncs_oper.vector[3]
ncs-1 generate 1 0 0 2 0 1 0 0 0 0 1 0
_space_group.IT_number 1
_space_group.name_H-M_alt 'P 1'
_space_group.name_Hall 'P 1'
_space_group_symop.id 1
_space_group_symop.operation_xyz 'x,y,z'
";

#[test]
fn bytes_are_dispatched_to_the_reader_their_content_names() {
    let read = read_bytes(DIPEPTIDE.as_bytes().to_vec(), None, &ReadOptions::new());
    assert_eq!(
        read.map(|(structure, _)| structure.atom_count()).ok(),
        Some(2)
    );
}

#[cfg(feature = "mmcif")]
#[test]
fn structured_text_is_dispatched_to_the_reader_its_content_names() {
    let text = b"data_1ABC\n_entry.id 1ABC\n".to_vec();
    // No coordinates, so the read is refused — but by the reader that owns the
    // format, which is the dispatch this checks.
    let refused = read_bytes(text, None, &ReadOptions::new());
    let findings = refused.err().unwrap_or_default();
    assert_eq!(findings.first().map(Diagnostic::code), Some(Code::E2001));
}

#[cfg(all(feature = "mmcif", feature = "xtal"))]
#[test]
fn mmcif_read_attaches_assemblies_from_the_same_document_parse() {
    use crate::{AssemblyExt, NcsExt, SymmetryExt};

    let structure = match read_bytes(
        MMCIF_XTAL.as_bytes().to_vec(),
        Some("entry.cif"),
        &ReadOptions::new(),
    ) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("read failed: {findings:?}"),
    };

    assert_eq!(
        structure.assembly_set().map(crate::AssemblySet::len),
        Some(1)
    );
    assert!(
        structure
            .assembly("1")
            .is_ok_and(|view| view.instance_count() == 1)
    );
    assert_eq!(structure.ncs_set().map(crate::NcsSet::len), Some(1));
    assert_eq!(
        structure.ncs_generated().map(|view| view.copy_count()),
        Some(1)
    );
    assert_eq!(
        structure
            .symmetry_set()
            .and_then(|set| set.international_number),
        Some(1)
    );
    assert!(
        structure
            .crystal_neighbors(1.0)
            .is_ok_and(|neighbors| neighbors.is_empty())
    );

    #[cfg(feature = "bcif")]
    {
        let input = InputBuffer::from_bytes(MMCIF_XTAL.as_bytes().to_vec());
        let document = match pdbiox_cif::parse(&input) {
            Ok((document, _)) => document,
            Err(findings) => panic!("fixture parse failed: {findings:?}"),
        };
        let bytes = match pdbiox_bcif::write_document(&document) {
            Ok(bytes) => bytes,
            Err(finding) => panic!("fixture BCIF write failed: {finding}"),
        };
        let binary = match read_bytes(bytes, Some("entry.bcif"), &ReadOptions::new()) {
            Ok((structure, _)) => structure,
            Err(findings) => panic!("BCIF read failed: {findings:?}"),
        };
        assert_eq!(binary.assembly_set().map(crate::AssemblySet::len), Some(1));
        assert_eq!(binary.ncs_set().map(crate::NcsSet::len), Some(1));
        assert_eq!(
            binary.symmetry_set().map(|set| set.operations().len()),
            Some(1)
        );
    }
}

#[test]
fn an_unreadable_path_reports_a_finding_rather_than_panicking() {
    assert!(read("no/such/file.pdb").is_err());
}

#[cfg(feature = "chem")]
#[test]
fn component_dictionary_reading_is_explicit_and_returns_a_versioned_provider() {
    use crate::ComponentProvider;

    let ccd = "data_HOH\n_chem_comp.id HOH\n_chem_comp.name WATER\n_chem_comp.type water\n";
    let path = std::env::temp_dir().join(format!("pdbiox-ccd-{}.cif", std::process::id()));
    if let Err(error) = std::fs::write(&path, ccd) {
        panic!("fixture write failed: {error}")
    }
    let read = read_component_dictionary(&path, DictionaryVersion::new("test"));
    let _removed = std::fs::remove_file(path);
    let (provider, findings) = match read {
        Ok(result) => result,
        Err(findings) => panic!("CCD read failed: {findings:?}"),
    };
    assert!(findings.is_empty());
    assert_eq!(provider.version().as_str(), "test");
    assert!(
        provider
            .get("HOH")
            .is_ok_and(|component| component.is_some())
    );
}

#[cfg(feature = "pdb")]
#[test]
fn a_structure_read_through_the_facade_writes_back_through_it() {
    let Ok((structure, _)) = read_bytes(DIPEPTIDE.as_bytes().to_vec(), None, &ReadOptions::new())
    else {
        panic!("expected a structure")
    };
    let written = write_pdb(&structure, &PdbOptions::new());
    assert!(written.unwrap_or_default().contains("ATOM"));
}

#[cfg(all(feature = "pdb", feature = "mmcif", feature = "gzip"))]
#[test]
fn a_generic_write_dispatches_by_name_and_reads_back() {
    let structure = match read_bytes(
        DIPEPTIDE.as_bytes().to_vec(),
        Some("test.pdb"),
        &ReadOptions::new(),
    ) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("read failed: {findings:?}"),
    };
    let path = std::env::temp_dir().join(format!("pdbiox-facade-{}.cif.gz", std::process::id()));
    if let Err(findings) = write(&path, &structure) {
        panic!("write failed: {findings:?}")
    }
    let round_tripped = read(&path);
    let _removed = std::fs::remove_file(path);
    let round_tripped = match round_tripped {
        Ok(structure) => structure,
        Err(findings) => panic!("round trip failed: {findings:?}"),
    };
    assert_eq!(round_tripped.atom_count(), structure.atom_count());
}

#[cfg(feature = "geom")]
#[test]
fn a_rigid_transform_reuses_coordinate_transactions_and_preserves_the_source() {
    let structure = match read_bytes(
        DIPEPTIDE.as_bytes().to_vec(),
        Some("test.pdb"),
        &ReadOptions::new(),
    ) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("read failed: {findings:?}"),
    };
    let moved = match transform(
        &structure,
        &AtomSelection::from_sorted(vec![0]),
        &pdbiox_geom::Rigid::translation([1.0, 2.0, 3.0]),
    ) {
        Ok(structure) => structure,
        Err(findings) => panic!("transform failed: {findings:?}"),
    };

    assert_eq!(structure.generation().get(), 0);
    assert_eq!(moved.generation().get(), 1);
    assert!(
        structure.positions()[1]
            .iter()
            .zip(moved.positions()[1])
            .all(|(original, actual)| (*original - actual).abs() < f32::EPSILON)
    );
    for ((actual, original), shift) in moved.positions()[0]
        .iter()
        .zip(structure.positions()[0])
        .zip([1.0, 2.0, 3.0])
    {
        assert!((*actual - original - shift).abs() < 1.0e-5);
    }
}

#[cfg(feature = "geom")]
#[test]
fn a_transform_rejects_atoms_outside_the_topology_before_editing() {
    let structure = match read_bytes(
        DIPEPTIDE.as_bytes().to_vec(),
        Some("test.pdb"),
        &ReadOptions::new(),
    ) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("read failed: {findings:?}"),
    };
    let result = transform(
        &structure,
        &AtomSelection::from_sorted(vec![2]),
        &pdbiox_geom::Rigid::IDENTITY,
    );
    assert_eq!(
        result
            .err()
            .and_then(|findings| findings.first().map(Diagnostic::code)),
        Some(Code::E6009)
    );
}
