use super::*;
use crate::columnar::AtomTable;
use arrow::array::{Array, AsArray};
use arrow::datatypes::Int32Type;
use arrow::record_batch::RecordBatch;
use std::fmt::Write as _;

const ATOMS: usize = 10_000;

/// Two chains (label A/B, author X/Y), two residue names, two atom names,
/// two elements, spread over three internal chunks.
fn structure() -> Structure {
    let mut source = String::from(
        "data_t\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
         _atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n\
         _atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n\
         _atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.B_iso_or_equiv\n\
         _atom_site.auth_asym_id\n",
    );
    for atom in 0..ATOMS {
        let first = atom < ATOMS / 2;
        let (chain, author, comp) = if first {
            ("A", "X", "ALA")
        } else {
            ("B", "Y", "GLY")
        };
        let (element, name) = if atom % 2 == 0 {
            ("C", "CA")
        } else {
            ("N", "N")
        };
        let alt = if atom == 1 { "A" } else { "." };
        let _ = writeln!(
            source,
            "ATOM {} {element} {name} {alt} {comp} {chain} {} 1 2 3 {}.5 {author}",
            atom + 1,
            atom / 2 + 1,
            atom % 100
        );
    }
    crate::test_support::read_cif(&source)
}

fn batches(structure: &Structure) -> Vec<RecordBatch> {
    match AtomTable::new(structure).record_batches() {
        Ok(batches) => batches,
        Err(error) => panic!("export failed: {error}"),
    }
}

/// Decodes one dictionary column of one batch to text, for assertions only.
fn strings(batch: &RecordBatch, name: &str) -> Vec<Option<String>> {
    let Some(column) = batch.column_by_name(name) else {
        panic!("{name} absent")
    };
    let dictionary = column.as_dictionary::<Int32Type>();
    let values = dictionary.values().as_string::<i32>();
    dictionary
        .keys()
        .iter()
        .map(|key| key.map(|key| values.value(key.unsigned_abs() as usize).to_owned()))
        .collect()
}

#[test]
fn categorical_text_columns_are_dictionaries_with_the_right_text() {
    let structure = structure();
    let batches = batches(&structure);
    assert!(batches.len() >= 2, "fixture must span several chunks");
    let mut chain = Vec::new();
    let mut label = Vec::new();
    let mut resname = Vec::new();
    let mut name = Vec::new();
    let mut element = Vec::new();
    let mut altloc = Vec::new();
    for batch in &batches {
        for column in [
            "chain_id",
            "label_asym_id",
            "resname",
            "atom_name",
            "element_symbol",
            "altloc",
        ] {
            assert_eq!(
                batch
                    .schema()
                    .field_with_name(column)
                    .map(|f| f.data_type().clone())
                    .ok(),
                Some(dictionary_type())
            );
        }
        chain.extend(strings(batch, "chain_id"));
        label.extend(strings(batch, "label_asym_id"));
        resname.extend(strings(batch, "resname"));
        name.extend(strings(batch, "atom_name"));
        element.extend(strings(batch, "element_symbol"));
        altloc.extend(strings(batch, "altloc"));
    }
    assert_eq!(chain.len(), ATOMS);
    for atom in 0..ATOMS {
        let first = atom < ATOMS / 2;
        assert_eq!(chain[atom].as_deref(), Some(if first { "X" } else { "Y" }));
        assert_eq!(label[atom].as_deref(), Some(if first { "A" } else { "B" }));
        assert_eq!(
            resname[atom].as_deref(),
            Some(if first { "ALA" } else { "GLY" })
        );
        let even = atom % 2 == 0;
        assert_eq!(name[atom].as_deref(), Some(if even { "CA" } else { "N" }));
        assert_eq!(element[atom].as_deref(), Some(if even { "C" } else { "N" }));
        assert_eq!(
            altloc[atom].as_deref(),
            if atom == 1 { Some("A") } else { None }
        );
    }
}

#[test]
fn the_dictionary_values_are_shared_by_every_batch_not_copied() {
    let batches = batches(&structure());
    assert!(batches.len() >= 2);
    for column in ["chain_id", "resname", "atom_name", "element_symbol"] {
        let pointers = batches
            .iter()
            .filter_map(|batch| batch.column_by_name(column))
            .map(|array| {
                array
                    .as_dictionary::<Int32Type>()
                    .values()
                    .to_data()
                    .buffers()[0]
                    .as_ptr()
            })
            .collect::<Vec<_>>();
        assert_eq!(pointers.len(), batches.len());
        assert!(
            pointers.windows(2).all(|pair| pair[0] == pair[1]),
            "{column}"
        );
    }
}

#[test]
fn dictionaries_hold_only_values_that_occur_and_never_per_atom_strings() {
    let dictionaries = AtomDictionaries::build(&structure());
    assert_eq!(dictionaries.chain.len(), 2);
    assert_eq!(dictionaries.label_asym.len(), 2);
    assert_eq!(dictionaries.resname.len(), 2);
    assert_eq!(dictionaries.atom_name.len(), 2);
    assert_eq!(dictionaries.element.len(), 2);
    assert_eq!(dictionaries.altloc.len(), 1);
    assert!(dictionaries.residue_chain.len() <= ATOMS);
}

#[test]
fn b_factor_and_schema_metadata_ride_along() {
    let structure = structure();
    let table = AtomTable::new(&structure);
    let schema = table.schema();
    let metadata = schema.metadata();
    assert_eq!(
        metadata.get("molframe:schema_version").map(String::as_str),
        Some("2")
    );
    assert_eq!(
        metadata.get("molframe:table").map(String::as_str),
        Some("atoms")
    );
    assert_eq!(
        metadata.get("molframe:atom_count").map(String::as_str),
        Some("10000")
    );
    assert!(
        metadata
            .get("molframe:units")
            .is_some_and(|units| units.contains("angstrom"))
    );
    assert!(metadata.contains_key("molframe:producer"));
    let unit = |name: &str| {
        schema
            .field_with_name(name)
            .ok()
            .and_then(|field| field.metadata().get("molframe:unit").cloned())
    };
    assert_eq!(unit("coordinates").as_deref(), Some("angstrom"));
    assert_eq!(unit("b_factor").as_deref(), Some("angstrom^2"));
    let first = &batches(&structure)[0];
    let b = first
        .column_by_name("b_factor")
        .map(|a| a.as_primitive::<arrow::datatypes::Float32Type>().value(1));
    assert!(b.is_some_and(|b| (b - 1.5).abs() < f32::EPSILON));
}

fn with_plddt(structure: &Structure) -> Structure {
    let mut data = structure.data().clone();
    let entries = (0..ATOMS).map(|atom| {
        let present = if atom == 3 {
            molframe_core::Presence::Unknown
        } else {
            molframe_core::Presence::Present
        };
        (
            match u16::try_from(atom % 100) {
                Ok(value) => f64::from(value) + 0.25,
                Err(_) => 0.25,
            },
            present,
        )
    });
    let column = match molframe_core::AnnotationColumn::from_entries(entries) {
        Ok(column) => column,
        Err(error) => panic!("annotation failed: {error:?}"),
    };
    data.annotations.insert(
        molframe_core::PLDDT_ANNOTATION,
        molframe_core::AtomAnnotation::Real(column),
    );
    Structure::new(data)
}

#[test]
fn plddt_is_exported_only_when_the_structure_carries_it() {
    let plain = structure();
    assert!(
        AtomTable::new(&plain)
            .schema()
            .field_with_name("plddt")
            .is_err()
    );
    let scored = with_plddt(&plain);
    let batches = batches(&scored);
    let column = batches[0].column_by_name("plddt").map(|array| {
        array
            .as_primitive::<arrow::datatypes::Float32Type>()
            .clone()
    });
    let Some(column) = column else {
        panic!("plddt absent")
    };
    assert!((column.value(2) - 2.25).abs() < f32::EPSILON);
    assert!(column.is_null(3), "unknown confidence is null");
    assert_eq!(column.len(), batches[0].num_rows());
}

#[test]
fn files_round_trip_dictionaries_and_merge_metadata_under_caller_keys() {
    let structure = structure();
    let directory = match tempfile::tempdir() {
        Ok(directory) => directory,
        Err(error) => panic!("temporary directory failed: {error}"),
    };
    let metadata = std::collections::BTreeMap::from([
        ("source".to_owned(), "fixture".to_owned()),
        ("molframe:table".to_owned(), "caller-wins".to_owned()),
    ]);
    let ipc = directory.path().join("atoms.arrow");
    let parquet = directory.path().join("atoms.parquet");
    assert!(crate::write_atom_ipc_with_metadata(&ipc, &structure, metadata.clone()).is_ok());
    assert!(crate::write_atom_parquet_with_metadata(&parquet, &structure, metadata).is_ok());

    let Ok(file) = std::fs::File::open(&ipc) else {
        panic!("ipc open failed")
    };
    let Ok(reader) = arrow::ipc::reader::FileReader::try_new(file, None) else {
        panic!("ipc read failed")
    };
    let schema = reader.schema();
    let rows: usize = reader
        .map(|batch| batch.map_or(0, |batch| batch.num_rows()))
        .sum();
    assert_eq!(rows, ATOMS);
    assert_eq!(
        schema
            .field_with_name("chain_id")
            .map(|f| f.data_type().clone())
            .ok(),
        Some(dictionary_type())
    );
    assert_eq!(
        schema.metadata().get("source").map(String::as_str),
        Some("fixture")
    );
    assert_eq!(
        schema.metadata().get("molframe:table").map(String::as_str),
        Some("caller-wins")
    );
    assert_eq!(
        schema
            .metadata()
            .get("molframe:schema_version")
            .map(String::as_str),
        Some("2")
    );

    let Ok(file) = std::fs::File::open(&parquet) else {
        panic!("parquet open failed")
    };
    let Ok(builder) = parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder::try_new(file)
    else {
        panic!("parquet read failed")
    };
    let schema = builder.schema().clone();
    assert_eq!(
        schema
            .field_with_name("resname")
            .map(|f| f.data_type().clone())
            .ok(),
        Some(dictionary_type())
    );
    assert_eq!(
        schema
            .metadata()
            .get("molframe:schema_version")
            .map(String::as_str),
        Some("2")
    );
    let Ok(reader) = builder.build() else {
        panic!("parquet batches failed")
    };
    let first: Vec<RecordBatch> = reader.filter_map(Result::ok).collect();
    let decoded = strings(&first[0], "resname");
    assert_eq!(decoded[0].as_deref(), Some("ALA"));
}
