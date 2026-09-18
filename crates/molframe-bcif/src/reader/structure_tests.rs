use super::*;
use crate::write_document;
use molframe_cif::{CifEventSink, CifScalar, CifValue, DataBlock};
use molframe_core::span::ByteSpan;
use molframe_core::structure::{StructureDifferenceOptions, structure_difference};

const MULTI_MODEL: &str = "\
data_demo
_entry.id DEMO
_struct.title 'direct BinaryCIF'
loop_
_entity.id
_entity.type
1 polymer
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
_space_group.IT_number 1
loop_
_ma_qa_metric.id
_ma_qa_metric.name
_ma_qa_metric.type
_ma_qa_metric.mode
1 score pLDDT local
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
_atom_site.future_annotation
ATOM 1 N N GLY A 1 1 0 0 0 1 same
ATOM 2 C CA GLY A 1 1 1 0 0 1 same
ATOM 1 N N GLY A 1 1 0 1 0 2 same
ATOM 2 C CA GLY A 1 1 1 1 0 2 same
";

fn document() -> Document {
    let mut document = Document::new();
    let mut block = DataBlock::new("test");
    let category = block.category_mut("entry", ByteSpan::default());
    category.column_mut("id").push(
        CifValue::Text("1ABC".into()),
        molframe_cif::lexer::Quoting::Bare,
    );
    document.push(block);
    document
}

#[test]
fn opening_a_document_keeps_columns_lazy_until_requested() {
    let bytes = write_document(&document()).expect("document writes");
    let input = InputBuffer::from_bytes(bytes);
    let binary = read_document(&input, Limits::default()).expect("document opens");
    assert_eq!(binary.block_count(), 1);
    let entry = binary
        .category(0, "entry")
        .expect("category decodes")
        .expect("entry exists");
    assert_eq!(entry.text("id", 0), Some("1ABC"));
}

#[test]
fn direct_and_lossless_binary_reads_are_semantically_identical() {
    let input = binary_input(MULTI_MODEL);
    let options = ReadOptions::new();
    let (direct, direct_findings) = read(&input, &options).expect("direct read succeeds");
    let (_, lossless, lossless_findings) =
        read_with_document(&input, &options).expect("lossless read succeeds");
    let difference = structure_difference(
        &direct,
        &lossless,
        StructureDifferenceOptions {
            coordinate_tolerance: 0.0,
        },
    )
    .expect("zero is a valid tolerance");

    assert!(difference.is_empty(), "difference: {difference:?}");
    assert_eq!(direct_findings, lossless_findings);
    assert_eq!(direct.data().entry.id.as_deref(), Some("DEMO"));
    assert_eq!(direct.model_count(), 2);
}

#[test]
fn an_unconsumed_atom_annotation_still_selects_ragged_storage() {
    let changed = MULTI_MODEL.replacen(
        "2 C CA GLY A 1 1 1 1 0 2 same",
        "2 C CA GLY A 1 1 1 1 0 2 changed",
        1,
    );
    let input = binary_input(&changed);
    let options = ReadOptions::new();
    let (direct, direct_findings) = read(&input, &options).expect("direct read succeeds");
    let (_, lossless, lossless_findings) =
        read_with_document(&input, &options).expect("lossless read succeeds");
    let difference = structure_difference(
        &direct,
        &lossless,
        StructureDifferenceOptions {
            coordinate_tolerance: 0.0,
        },
    )
    .expect("zero is a valid tolerance");

    assert!(direct.ragged_models().is_some());
    assert!(difference.is_empty(), "difference: {difference:?}");
    assert_eq!(direct_findings, lossless_findings);
}

#[test]
fn projected_core_xtal_and_modelcif_categories_match_the_lossless_document() {
    let input = binary_input(MULTI_MODEL);
    let options = ReadOptions::new();
    let (projected, _, _) = read_with_metadata(&input, &options, keep_extension_category)
        .expect("projected read succeeds");
    let (lossless, _, _) = read_with_document(&input, &options).expect("lossless read succeeds");
    let projected_block = projected.first_block().expect("projected block exists");
    let lossless_block = lossless.first_block().expect("lossless block exists");

    assert!(projected_block.category("atom_site").is_none());
    for name in [
        "entry",
        "struct",
        "entity",
        "struct_asym",
        "space_group",
        "ma_qa_metric",
    ] {
        let actual = projected_block.category(name).expect("category projected");
        let expected = lossless_block.category(name).expect("category decoded");
        assert_eq!(
            actual.items().collect::<Vec<_>>(),
            expected.items().collect::<Vec<_>>()
        );
        assert_eq!(actual.row_count(), expected.row_count());
        for item in expected.items() {
            for row in 0..expected.row_count() {
                assert_eq!(actual.value(item, row), expected.value(item, row));
            }
        }
    }
}

#[derive(Default)]
struct ModelProbe {
    cells: Vec<(String, String, CifValue)>,
}

impl CifEventSink for ModelProbe {
    type Output = Self;

    fn block(&mut self, _name: &str) {}

    fn accepts_category(&self, category: &str) -> bool {
        category.starts_with("ma_")
    }

    fn value(&mut self, category: &str, item: &str, value: CifScalar<'_>, _span: ByteSpan) {
        let owned = match value {
            CifScalar::Inapplicable => CifValue::Inapplicable,
            CifScalar::Unknown => CifValue::Unknown,
            CifScalar::Text(value) => CifValue::Text(value.into()),
            CifScalar::Integer(value) => CifValue::Integer(value),
            CifScalar::Float(value) => CifValue::Float(value),
        };
        self.cells
            .push((category.to_owned(), item.to_owned(), owned));
    }

    fn finish(self) -> Self::Output {
        self
    }
}

#[test]
fn compact_projection_receives_exact_binary_cells_without_a_modelcif_dom() {
    let input = binary_input(MULTI_MODEL);
    let options = ReadOptions::new();
    let (metadata, _, projected, findings) =
        read_with_projection(&input, &options, |_| false, ModelProbe::default())
            .expect("combined BinaryCIF projection succeeds");
    let binary = read_document(&input, Limits::default()).expect("lossless binary opens");
    let lossless = binary.to_document().expect("lossless binary decodes");
    let category = lossless
        .first_block()
        .and_then(|block| block.category("ma_qa_metric"))
        .expect("fixture has ModelCIF metadata");
    let mut expected = Vec::new();
    for item in category.items() {
        for row in 0..category.row_count() {
            let value = category
                .value(item, row)
                .expect("complete BinaryCIF column")
                .clone();
            expected.push((category.name().to_owned(), item.to_owned(), value));
        }
    }

    assert!(findings.is_empty());
    assert_eq!(projected.cells, expected);
    assert!(
        metadata
            .first_block()
            .and_then(|block| block.category("ma_qa_metric"))
            .is_none()
    );
}

fn binary_input(text: &str) -> InputBuffer {
    let source = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let (document, findings) = molframe_cif::parse(&source).expect("fixture parses");
    assert!(findings.is_empty(), "fixture findings: {findings:?}");
    let bytes = write_document(&document).expect("fixture encodes");
    InputBuffer::from_bytes(bytes)
}

fn keep_extension_category(category: &str) -> bool {
    category.starts_with("ma_")
        || matches!(
            category,
            "pdbx_struct_oper_list"
                | "pdbx_struct_assembly"
                | "pdbx_struct_assembly_gen"
                | "struct_ncs_oper"
                | "space_group"
                | "symmetry"
                | "space_group_symop"
                | "symmetry_equiv"
        )
}
