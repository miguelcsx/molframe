use super::*;
use crate::extension_name;
use arrow::array::Array;
use arrow::ffi_stream::ArrowArrayStreamReader;
use pdbiox_core::chunk::{AtomRecord, ChunkBuilder};
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_core::optional::{OptionalI32, OptionalSymbol};
use pdbiox_core::structure::{CoordinateStore, StructureData};
use pdbiox_core::topology::ResidueRecord;
use pdbiox_core::{AltId, Element, Presence};

const SOURCE: &str = "data_a\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA ALA A 1 1 2 3\nATOM 2 N N ALA A 1 4 5 6\n";

fn structure() -> Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    match pdbiox_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

#[test]
fn coordinates_and_plain_columns_alias_the_structure_snapshot() {
    let structure = structure();
    let table = AtomTable::new(&structure);
    let batches = match table.record_batches() {
        Ok(batches) => batches,
        Err(error) => panic!("Arrow export failed: {error}"),
    };
    let Some(batch) = batches.first() else {
        panic!("batch absent")
    };
    let Some(coordinates) = batch
        .column_by_name("coordinates")
        .and_then(|array| array.as_any().downcast_ref::<FixedSizeListArray>())
    else {
        panic!("coordinates absent")
    };
    let Some(values) = coordinates.values().as_any().downcast_ref::<Float32Array>() else {
        panic!("coordinate values have wrong type")
    };
    assert!(std::ptr::eq(
        values.values().as_ptr(),
        structure.positions().as_ptr().cast::<f32>()
    ));
    assert_eq!(coordinates.len(), structure.atom_count() as usize);
    assert_eq!(batch_u32_values(batch, "residue_index"), &[0, 0]);
}

#[test]
fn schema_carries_domain_extensions_and_copy_costs() {
    let table = AtomTable::new(&structure());
    let schema = table.schema();
    let Some(field) = schema.field_with_name("coordinates").ok() else {
        panic!("coordinate field absent")
    };
    assert_eq!(extension_name(field), Some("pdbiox.coordinates3f"));
    assert_eq!(
        field
            .metadata()
            .get("pdbiox:export_cost")
            .map(String::as_str),
        Some("zero-copy")
    );
    assert!(table.arrow_stream().is_ok());
}

#[test]
fn atom_c_stream_yields_one_internal_chunk_at_a_time() {
    let structure = chunked_structure();
    let table = AtomTable::new(&structure);
    let stream = match table.arrow_stream() {
        Ok(stream) => stream,
        Err(error) => panic!("stream creation failed: {error}"),
    };
    let reader = match ArrowArrayStreamReader::try_new(stream.into_ffi()) {
        Ok(reader) => reader,
        Err(error) => panic!("stream import failed: {error}"),
    };
    let rows = reader
        .map(|batch| match batch {
            Ok(batch) => batch.num_rows(),
            Err(error) => panic!("stream pull failed: {error}"),
        })
        .collect::<Vec<_>>();
    assert_eq!(rows, [2, 2, 2]);
}

#[test]
fn pulled_batch_outlives_stream_table_and_structure_handles() {
    let structure = structure();
    let table = AtomTable::new(&structure);
    let stream = match table.arrow_stream() {
        Ok(stream) => stream,
        Err(error) => panic!("stream creation failed: {error}"),
    };
    let mut reader = match ArrowArrayStreamReader::try_new(stream.into_ffi()) {
        Ok(reader) => reader,
        Err(error) => panic!("stream import failed: {error}"),
    };
    let batch = match reader.next() {
        Some(Ok(batch)) => batch,
        Some(Err(error)) => panic!("stream pull failed: {error}"),
        None => panic!("stream batch absent"),
    };
    drop(reader);
    drop(table);
    drop(structure);

    let Some(coordinates) = batch
        .column_by_name("coordinates")
        .and_then(|array| array.as_any().downcast_ref::<FixedSizeListArray>())
    else {
        panic!("coordinates absent")
    };
    let Some(values) = coordinates.values().as_any().downcast_ref::<Float32Array>() else {
        panic!("coordinate values have wrong type")
    };
    assert_eq!(values.values(), &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
}

#[test]
fn every_pulled_batch_remains_live_while_later_batches_are_pulled() {
    let structure = chunked_structure();
    let table = AtomTable::new(&structure);
    let stream = match table.arrow_stream() {
        Ok(stream) => stream,
        Err(error) => panic!("stream creation failed: {error}"),
    };
    let mut reader = match ArrowArrayStreamReader::try_new(stream.into_ffi()) {
        Ok(reader) => reader,
        Err(error) => panic!("stream import failed: {error}"),
    };
    let first = match reader.next() {
        Some(Ok(batch)) => batch,
        Some(Err(error)) => panic!("first stream pull failed: {error}"),
        None => panic!("first stream batch absent"),
    };
    let second = match reader.next() {
        Some(Ok(batch)) => batch,
        Some(Err(error)) => panic!("second stream pull failed: {error}"),
        None => panic!("second stream batch absent"),
    };
    let third = match reader.next() {
        Some(Ok(batch)) => batch,
        Some(Err(error)) => panic!("third stream pull failed: {error}"),
        None => panic!("third stream batch absent"),
    };
    drop(reader);
    drop(table);
    drop(structure);

    assert_eq!(batch_coordinates(&first), &[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]);
    assert_eq!(batch_u32_values(&first, "atom_index"), &[0, 1]);
    assert_eq!(batch_u32_values(&first, "residue_index"), &[0, 1]);
    drop(first);
    assert_eq!(batch_coordinates(&second), &[2.0, 0.0, 0.0, 3.0, 0.0, 0.0]);
    assert_eq!(batch_u32_values(&second, "atom_index"), &[2, 3]);
    assert_eq!(batch_u32_values(&second, "residue_index"), &[2, 3]);
    drop(second);
    assert_eq!(batch_coordinates(&third), &[4.0, 0.0, 0.0, 5.0, 0.0, 0.0]);
    assert_eq!(batch_u32_values(&third, "atom_index"), &[4, 5]);
    assert_eq!(batch_u32_values(&third, "residue_index"), &[4, 5]);
}

fn batch_coordinates(batch: &RecordBatch) -> &[f32] {
    let Some(coordinates) = batch
        .column_by_name("coordinates")
        .and_then(|array| array.as_any().downcast_ref::<FixedSizeListArray>())
    else {
        panic!("coordinates absent")
    };
    let Some(values) = coordinates.values().as_any().downcast_ref::<Float32Array>() else {
        panic!("coordinate values have wrong type")
    };
    values.values()
}

fn batch_u32_values<'a>(batch: &'a RecordBatch, name: &str) -> &'a [u32] {
    let Some(indices) = batch
        .column_by_name(name)
        .and_then(|array| array.as_any().downcast_ref::<UInt32Array>())
    else {
        panic!("{name} values absent")
    };
    indices.values()
}

fn chunked_structure() -> Structure {
    let mut data = StructureData::empty();
    let atom_name = match data.dictionary.intern("CA") {
        Ok(name) => name,
        Err(error) => panic!("atom name failed: {error}"),
    };
    let component = match data.dictionary.intern("ALA") {
        Ok(component) => component,
        Err(error) => panic!("component failed: {error}"),
    };
    let mut builder = ChunkBuilder::with_target(2);
    for atom in 0..6_u32 {
        let coordinate = match u16::try_from(atom) {
            Ok(value) => f32::from(value),
            Err(error) => panic!("coordinate failed: {error}"),
        };
        let residue = match data.topology.residues.push(
            ResidueRecord {
                label_comp_id: component,
                auth_comp_id: OptionalSymbol::NONE,
                label_seq_id: OptionalI32::some(atom.cast_signed() + 1),
                auth_seq_id: OptionalI32::some(atom.cast_signed() + 1),
                ins_code: OptionalSymbol::NONE,
                het: false,
            },
            atom..atom + 1,
        ) {
            Ok(residue) => residue,
            Err(error) => panic!("residue failed: {error}"),
        };
        builder.push(AtomRecord {
            position: Some([coordinate, 0.0, 0.0]),
            element: Element::CARBON,
            atom_name,
            auth_atom_name: OptionalSymbol::NONE,
            alternate_component_id: OptionalSymbol::NONE,
            alt_id: AltId::BLANK,
            residue,
            occupancy: (1.0, Presence::Present),
            b_factor: (10.0, Presence::Present),
            formal_charge: (0, Presence::Inapplicable),
            atom_site_id: atom + 1,
        });
    }
    let (chunks, coordinates) = builder.finish();
    data.chunks = chunks.into();
    data.coords = CoordinateStore::Single(coordinates);
    Structure::new(data)
}
