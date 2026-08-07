use super::*;
use crate::extension_name;
use arrow::array::Array;
use pdbiox_core::io::{InputBuffer, ReadOptions};

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
