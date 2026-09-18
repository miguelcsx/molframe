use super::*;

#[cfg(feature = "pdb")]
use crate::PdbOptions;

mod namespace_tests;
mod plddt_tests;

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
    let Err(findings) = read_bytes(text, None, &ReadOptions::new()) else {
        panic!("a document with no coordinates should have been refused")
    };
    assert_eq!(findings.first().map(Diagnostic::code), Some(Code::E2001));
}

#[cfg(all(feature = "mmcif", feature = "xtal"))]
#[test]
fn mmcif_read_attaches_assemblies_from_the_same_document_parse() {
    use crate::{AssemblyExt, NcsExt, SymmetryExt};

    let (structure, direct_findings) = match read_bytes(
        MMCIF_XTAL.as_bytes().to_vec(),
        Some("entry.cif"),
        &ReadOptions::new(),
    ) {
        Ok(result) => result,
        Err(findings) => panic!("read failed: {findings:?}"),
    };

    let input = InputBuffer::from_bytes(MMCIF_XTAL.as_bytes().to_vec());
    let options = ReadOptions::new();
    let (document, base, base_findings) =
        molframe_cif::read_with_document(&input, &options).expect("lossless read should succeed");
    let (lossless, lossless_findings) = super::extensions::attach_materialized_cif_metadata(
        &document,
        base,
        base_findings,
        &options,
    )
    .expect("lossless metadata lowering should succeed");
    let difference = crate::structure_difference(
        &structure,
        &lossless,
        crate::StructureDifferenceOptions {
            coordinate_tolerance: 0.0,
        },
    )
    .expect("zero is a valid coordinate tolerance");
    assert!(difference.is_empty(), "difference: {difference:?}");
    assert_eq!(direct_findings, lossless_findings);

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
            .collect_crystal_neighbors(1.0, &molframe_core::ExecutionContext::default())
            .is_ok_and(|neighbors| neighbors.is_empty())
    );
    assert_eq!(
        structure.assembly_set().and_then(|set| set.get("1")),
        lossless.assembly_set().and_then(|set| set.get("1")),
    );
    assert_eq!(
        structure.ncs_set().and_then(|set| set.get("ncs-1")),
        lossless.ncs_set().and_then(|set| set.get("ncs-1")),
    );
    assert_eq!(
        structure
            .symmetry_set()
            .map(molframe_xtal::SymmetrySet::operations),
        lossless
            .symmetry_set()
            .map(molframe_xtal::SymmetrySet::operations),
    );

    #[cfg(feature = "bcif")]
    assert_bcif_roundtrip_matches(&lossless, &lossless_findings, &options);
}

#[cfg(feature = "bcif")]
fn assert_bcif_roundtrip_matches(
    lossless: &crate::Structure,
    lossless_findings: &[Diagnostic],
    options: &ReadOptions,
) {
    use crate::{AssemblyExt, NcsExt, SymmetryExt};

    let input = InputBuffer::from_bytes(MMCIF_XTAL.as_bytes().to_vec());
    let document = match molframe_cif::parse(&input) {
        Ok((document, _)) => document,
        Err(findings) => panic!("fixture parse failed: {findings:?}"),
    };
    let bytes = match molframe_bcif::write_document(&document) {
        Ok(bytes) => bytes,
        Err(finding) => panic!("fixture BCIF write failed: {finding}"),
    };
    let (binary, binary_findings) = match read_bytes(bytes, Some("entry.bcif"), options) {
        Ok(result) => result,
        Err(findings) => panic!("BCIF read failed: {findings:?}"),
    };
    let difference = crate::structure_difference(
        &binary,
        lossless,
        crate::StructureDifferenceOptions {
            coordinate_tolerance: 0.0,
        },
    )
    .expect("zero is a valid coordinate tolerance");
    assert!(difference.is_empty(), "BCIF difference: {difference:?}");
    assert_eq!(binary_findings.as_slice(), lossless_findings);
    assert_eq!(binary.assembly_set().map(crate::AssemblySet::len), Some(1));
    assert_eq!(binary.ncs_set().map(crate::NcsSet::len), Some(1));
    assert_eq!(
        binary.symmetry_set().map(|set| set.operations().len()),
        Some(1)
    );
    assert_eq!(
        binary.assembly_set().and_then(|set| set.get("1")),
        lossless.assembly_set().and_then(|set| set.get("1")),
    );
    assert_eq!(
        binary.ncs_set().and_then(|set| set.get("ncs-1")),
        lossless.ncs_set().and_then(|set| set.get("ncs-1")),
    );
    assert_eq!(
        binary
            .symmetry_set()
            .map(molframe_xtal::SymmetrySet::operations),
        lossless
            .symmetry_set()
            .map(molframe_xtal::SymmetrySet::operations),
    );
}

#[test]
fn an_unreadable_path_reports_a_finding_rather_than_panicking() {
    assert!(read("no/such/file.pdb").is_err());
}

#[cfg(feature = "chem")]
#[test]
fn component_dictionary_reading_is_explicit_and_returns_a_versioned_provider() {
    use crate::ComponentProvider;

    let ccd = "data_HOH\n\
_chem_comp.id HOH\n_chem_comp.name WATER\n_chem_comp.type water\n\
loop_\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n\
_chem_comp_atom.pdbx_aromatic_flag\n_chem_comp_atom.pdbx_leaving_atom_flag\nO O 0 N N\n";
    let path = std::env::temp_dir().join(format!("molframe-ccd-{}.cif", std::process::id()));
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
    let structure = with_entry_id(&structure, "test");
    let path = std::env::temp_dir().join(format!("molframe-facade-{}.cif.gz", std::process::id()));
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

#[cfg(feature = "pdb")]
#[test]
fn mmtf_is_dispatched_by_content_and_suffix() {
    let (structure, _) = read_bytes(
        DIPEPTIDE.as_bytes().to_vec(),
        Some("test.pdb"),
        &ReadOptions::new(),
    )
    .unwrap_or_else(|findings| panic!("fixture read failed: {findings:?}"));
    let mut output = Vec::new();
    let rendered = write_stream(
        &mut output,
        &structure,
        Format::Mmtf,
        OutputOptions::default().memory_limit_bytes,
    );
    assert!(rendered.is_err());
}

#[cfg(feature = "mmcif")]
#[test]
fn pdbml_input_remains_supported_without_an_eager_generic_writer() {
    let (structure, _) = read_bytes(
        DIPEPTIDE.as_bytes().to_vec(),
        Some("test.pdb"),
        &ReadOptions::new(),
    )
    .unwrap_or_else(|findings| panic!("fixture read failed: {findings:?}"));
    let structure = with_entry_id(&structure, "test");
    let canonical = molframe_cif::write_canonical(&structure).expect("canonical CIF render");
    let input = InputBuffer::from_bytes(canonical.into_bytes());
    let (document, _) = molframe_cif::parse(&input).expect("canonical CIF parse");
    let rendered = molframe_cif::write_pdbml(&document)
        .expect("explicit PDBML render")
        .into_bytes();
    let (decoded, findings) =
        read_bytes(rendered, None, &ReadOptions::new()).expect("PDBML content dispatch");
    assert!(findings.is_empty());
    assert_eq!(decoded.atom_count(), structure.atom_count());
    let mut output = Vec::new();
    assert!(
        write_stream(
            &mut output,
            &structure,
            Format::Pdbml,
            OutputOptions::default().memory_limit_bytes,
        )
        .is_err()
    );
    assert!(output.is_empty());
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
        &molframe_geom::Rigid::translation([1.0, 2.0, 3.0]),
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
        &molframe_geom::Rigid::IDENTITY,
    );
    assert_eq!(
        result
            .err()
            .and_then(|findings| findings.first().map(Diagnostic::code)),
        Some(Code::E6009)
    );
}

#[cfg(feature = "modelcif")]
#[test]
fn modelcif_metadata_is_attached_and_written_through_the_facade() {
    use crate::ModelCifExt;

    let source = "data_model\n\
loop_\n_ma_qa_metric.id\n_ma_qa_metric.name\n_ma_qa_metric.type\n_ma_qa_metric.mode\n1 score pLDDT local\n#\n\
loop_\n_ma_qa_metric_local.model_id\n_ma_qa_metric_local.label_asym_id\n_ma_qa_metric_local.label_seq_id\n_ma_qa_metric_local.label_comp_id\n_ma_qa_metric_local.metric_id\n_ma_qa_metric_local.metric_value\n1 A 1 GLY 1 95.0\n#\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 1 0 0 0\n";
    let (structure, _) = read_bytes(
        source.as_bytes().to_vec(),
        Some("model.cif"),
        &ReadOptions::new(),
    )
    .unwrap_or_else(|findings| panic!("read failed: {findings:?}"));
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let (document, _) = molframe_cif::parse(&input)
        .unwrap_or_else(|findings| panic!("lossless parse failed: {findings:?}"));
    let (expected_model, expected_findings) =
        molframe_modelcif::lower(&document).expect("ModelCIF lowering succeeds");
    assert!(expected_findings.is_empty());
    assert_eq!(structure.model_cif(), Some(&expected_model));
    assert_eq!(
        structure
            .confidence()
            .map(|confidence| confidence.plddt().count()),
        Some(1)
    );
    let mut rendered = Vec::new();
    write_mmcif_to(&structure, &mut rendered)
        .unwrap_or_else(|error| panic!("write failed: {error}"));
    assert!(String::from_utf8_lossy(&rendered).contains("_ma_qa_metric_local.metric_value"));

    #[cfg(feature = "bcif")]
    {
        let bytes = molframe_bcif::write_document(&document).expect("fixture BCIF writes");
        let (binary, _) = read_bytes(bytes, Some("model.bcif"), &ReadOptions::new())
            .unwrap_or_else(|findings| panic!("BCIF read failed: {findings:?}"));
        assert_eq!(binary.model_cif(), Some(&expected_model));
        assert_eq!(binary.confidence(), structure.confidence());
    }
}

fn with_entry_id(structure: &Structure, id: &str) -> Structure {
    let mut data = structure.data().clone();
    data.entry.id = Some(id.into());
    Structure::from(data)
}
