use super::*;
use molframe_core::io::{InputBuffer, ReadOptions};
use std::io::{self, Write};

const SOURCE: &str = "data_source\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
_atom_site.pdbx_PDB_model_num\nATOM 44 C CA GLY LONG_CHAIN 9 1 2 3 7\n\
ATOM 45 N N GLY LONG_CHAIN 9 2 2 3 7\n";

#[test]
fn absent_block_identity_refuses_without_an_output_value() {
    let structure = molframe_core::Structure::new(molframe_core::StructureData::empty());
    assert_eq!(
        write_canonical(&structure),
        Err(CifWriteError::MissingBlockId)
    );
}

#[test]
fn explicit_block_identity_preserves_atom_model_and_chain_identifiers() {
    let structure = read(SOURCE);
    let options = CifWriteOptions::new().with_block_id("chosen");
    let Ok(output) = write_canonical_with_options(&structure, &options) else {
        panic!("complete structure should write");
    };
    assert!(output.starts_with("data_chosen\n"), "{output}");
    assert!(output.contains("ATOM 44 C CA . GLY LONG_CHAIN . 9 ? 1.000 2.000 3.000"));
    assert!(output.ends_with(" 7\n#\n"), "{output}");
}

#[test]
fn absent_atom_site_identity_is_never_replaced_by_an_ordinal() {
    let structure = read(&SOURCE.replace("ATOM 44", "ATOM ?"));
    let options = CifWriteOptions::new().with_block_id("chosen");
    assert_eq!(
        write_canonical_with_options(&structure, &options),
        Err(CifWriteError::MissingAtomSiteId { atom: 0 })
    );
}

#[test]
fn anisotropy_writes_its_own_category_joined_by_deposited_id_and_round_trips() {
    let mut data = read(SOURCE).data().clone();
    let mut builder = molframe_core::anisotropy::AnisotropyTableBuilder::new();
    builder.push(molframe_core::AnisotropicDisplacement {
        atom: molframe_core::AtomIndex::new(0),
        u: [0.011, 0.022, 0.033, 0.001, -0.002, 0.003],
    });
    data.anisotropy = builder.finish();
    let structure = molframe_core::Structure::new(data);
    let options = CifWriteOptions::new().with_block_id("chosen");
    let Ok(output) = write_canonical_with_options(&structure, &options) else {
        panic!("structure with ellipsoids should write");
    };
    assert!(output.contains("_atom_site_anisotrop.U[1][1]"), "{output}");
    assert!(
        output.contains("44 C 0.01100 0.02200 0.03300 0.00100 -0.00200 0.003"),
        "{output}"
    );
    // The emitted category reads back as the same tensor on the same atoms.
    let read_back = read(&output);
    let anisotropy = &read_back.data().anisotropy;
    assert!(anisotropy.is_available());
    assert_eq!(anisotropy.len(), 1);
    let record = anisotropy
        .get(molframe_core::AnisotropyIndex::new(0))
        .expect("round-tripped");
    assert_eq!(record.atom, molframe_core::AtomIndex::new(0));
    assert!((record.u[0] - 0.011).abs() <= 1e-6);
    assert!((record.u[5] - 0.003).abs() <= 1e-6);
}

#[test]
fn a_structure_without_anisotropy_writes_no_anisotrop_category() {
    let structure = read(SOURCE);
    let options = CifWriteOptions::new().with_block_id("chosen");
    let Ok(output) = write_canonical_with_options(&structure, &options) else {
        panic!("plain structure should write");
    };
    assert!(!output.contains("_atom_site_anisotrop"), "{output}");
}

#[test]
fn bond_identifiers_require_explicit_generation_permission() {
    let source = format!(
        "{SOURCE}loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n\
         _struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n\
         _struct_conn.ptnr1_label_comp_id\n_struct_conn.ptnr1_label_atom_id\n\
         _struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n\
         _struct_conn.ptnr2_label_comp_id\n_struct_conn.ptnr2_label_atom_id\n\
         c1 covale LONG_CHAIN 9 GLY CA LONG_CHAIN 9 GLY N\n"
    );
    let structure = read(&source);
    let options = CifWriteOptions::new().with_block_id("chosen");
    assert_eq!(
        write_canonical_with_options(&structure, &options),
        Err(CifWriteError::ConnectionIdsNotEnabled)
    );
}

#[derive(Default)]
struct CountingWriter {
    bytes: usize,
    writes: usize,
    largest_write: usize,
}

impl Write for CountingWriter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.bytes += bytes.len();
        self.writes += 1;
        self.largest_write = self.largest_write.max(bytes.len());
        Ok(bytes.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[test]
fn canonical_writer_emits_bounded_fragments_without_a_full_file() {
    let structure = read(SOURCE);
    let options = CifWriteOptions::new().with_block_id("chosen");
    let expected = match write_canonical_with_options(&structure, &options) {
        Ok(output) => output,
        Err(error) => panic!("in-memory write failed: {error}"),
    };
    let mut writer = CountingWriter::default();
    if let Err(error) = write_canonical_to(&structure, &options, &mut writer) {
        panic!("stream write failed: {error}");
    }
    assert_eq!(writer.bytes, expected.len());
    assert!(writer.writes > 1);
    assert!(writer.largest_write < expected.len());
}

#[test]
fn canonical_preflight_refuses_before_touching_the_destination() {
    let structure = molframe_core::Structure::new(molframe_core::StructureData::empty());
    let mut writer = CountingWriter::default();
    assert!(
        matches!(
            write_canonical_to(&structure, &CifWriteOptions::new(), &mut writer),
            Err(CifWriteToError::Projection(CifWriteError::MissingBlockId))
        ),
        "missing identity must remain a projection error"
    );
    assert_eq!(writer.bytes, 0);
    assert_eq!(writer.writes, 0);
}

fn read(source: &str) -> molframe_core::Structure {
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    match crate::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}
