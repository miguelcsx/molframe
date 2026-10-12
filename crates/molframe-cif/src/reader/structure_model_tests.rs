//! Multi-model reads: dense frames, ragged ensembles, first-model reads.

use super::*;

#[derive(Default)]
struct ProjectionProbe {
    atom_values: usize,
    extension_values: usize,
}

impl crate::CifEventSink for ProjectionProbe {
    type Output = Self;

    fn block(&mut self, _name: &str) {}

    fn value(
        &mut self,
        category: &str,
        _item: &str,
        _value: crate::CifScalar<'_>,
        _span: molframe_core::span::ByteSpan,
    ) {
        if category == "atom_site" {
            self.atom_values += 1;
        } else if category == "ma_probe" {
            self.extension_values += 1;
        }
    }

    fn finish(self) -> Self::Output {
        self
    }
}

#[test]
fn structure_and_external_projection_share_one_parser_event_stream() {
    let text = "data_x\n_ma_probe.id 7\nloop_\n\
_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n\
_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n\
_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
1 C CA GLY A 1 0 0 0\n";
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let (_, structure, probe, findings) = read_with_projection(
        &input,
        &ReadOptions::new(),
        |_| false,
        ProjectionProbe::default(),
    )
    .expect("combined direct read succeeds");
    assert!(findings.is_empty());
    assert_eq!(structure.atom_count(), 1);
    assert_eq!(probe.atom_values, 9);
    assert_eq!(probe.extension_values, 1);
}

#[test]
fn dense_models_keep_numbers_and_identity_changes_become_ragged() {
    let text = "\
data_x
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
_atom_site.pdbx_PDB_model_num
1 N N GLY A 1 0 0 0 5
2 C CA GLY A 1 1 0 0 5
1 N N GLY A 1 0 1 0 9
2 C CA GLY A 1 1 1 0 9
#
";
    let (structure, direct_findings) = parse_structure(text);
    let dense_input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let (_, dense_lossless, lossless_findings) =
        read_with_document(&dense_input, &ReadOptions::new())
            .expect("lossless dense read should succeed");
    let dense_difference = structure_difference(
        &structure,
        &dense_lossless,
        StructureDifferenceOptions {
            coordinate_tolerance: 0.0,
        },
    )
    .expect("zero is a valid coordinate tolerance");
    assert!(
        dense_difference.is_empty(),
        "difference: {dense_difference:?}"
    );
    assert_eq!(direct_findings, lossless_findings);
    let numbers: Vec<_> = structure
        .data()
        .models()
        .filter_map(molframe_core::structure::ModelRef::number)
        .collect();
    assert_eq!(numbers, [5, 9]);
    assert_eq!(structure.model_count(), 2);

    let mismatched = text.replacen("2 C CA GLY A 1 1 1 0 9", "2 O O GLY A 1 1 1 0 9", 1);
    let input = InputBuffer::from_bytes(mismatched.into_bytes());
    let options = ReadOptions::new().mode(ParseMode::Recover);
    let (ragged, direct_findings) = match read(&input, &options) {
        Ok(result) => result,
        Err(findings) => panic!("ragged read failed: {findings:?}"),
    };
    let (_, ragged_lossless, lossless_findings) =
        read_with_document(&input, &options).expect("lossless ragged read should succeed");
    let ragged_difference = structure_difference(
        &ragged,
        &ragged_lossless,
        StructureDifferenceOptions {
            coordinate_tolerance: 0.0,
        },
    )
    .expect("zero is a valid coordinate tolerance");
    assert!(
        ragged_difference.is_empty(),
        "difference: {ragged_difference:?}"
    );
    assert_eq!(direct_findings, lossless_findings);
    let Some(models) = ragged.ragged_models() else {
        panic!("identity-changing models must be ragged")
    };
    assert_eq!(models.len(), 2);
    let names: Vec<Vec<_>> = models
        .iter()
        .map(|model| model.data().atoms().filter_map(AtomRef::name).collect())
        .collect();
    assert_eq!(names, [["N", "CA"], ["N", "O"]]);
    let numbers: Vec<_> = ragged
        .data()
        .models()
        .filter_map(molframe_core::structure::ModelRef::number)
        .collect();
    assert_eq!(numbers, [5, 9]);
}

#[test]
fn an_unknown_per_atom_annotation_change_still_makes_models_ragged() {
    let text = "\
data_x
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
_atom_site.pdbx_PDB_model_num
_atom_site.future_annotation
1 N N GLY A 1 0 0 0 1 same
1 N N GLY A 1 1 1 1 2 changed
#
";
    let (structure, _) = parse_structure(text);
    assert!(structure.ragged_models().is_some());
}
#[test]
fn reading_only_the_first_cif_model_does_not_append_an_empty_frame() {
    let text = "\
data_x
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
_atom_site.pdbx_PDB_model_num
1 N N GLY A 1 0 0 0 4
1 N N GLY A 1 1 1 1 8
#
";
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let options = ReadOptions::new().only_first_model(true);
    let (structure, _) = match read(&input, &options) {
        Ok(result) => result,
        Err(findings) => panic!("read failed: {findings:?}"),
    };
    assert_eq!(structure.model_count(), 1);
    let number = structure
        .data()
        .models()
        .next()
        .and_then(molframe_core::structure::ModelRef::number);
    assert_eq!(number, Some(4));

    let (_, lossless, _) =
        read_with_document(&input, &options).expect("lossless first-model read should succeed");
    let difference = structure_difference(
        &structure,
        &lossless,
        StructureDifferenceOptions {
            coordinate_tolerance: 0.0,
        },
    )
    .expect("zero is a valid coordinate tolerance");
    assert!(difference.is_empty(), "difference: {difference:?}");
}

const PER_MODEL_VALUES: &str = "\
data_x
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
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
1 N N GLY A 1 0 0 0 1.00 10.00 1
1 N N GLY A 1 1 1 1 0.50 20.00 2
1 N N GLY A 1 2 2 2 0.25 30.00 3
#
";

#[test]
fn per_model_occupancy_and_b_factor_stay_dense_and_are_reported_once() {
    let (structure, findings) = parse_structure(PER_MODEL_VALUES);
    assert!(structure.ragged_models().is_none());
    assert_eq!(structure.model_count(), 3);
    let atom = structure.data().atoms().next().expect("one atom");
    assert_eq!(atom.occupancy(), Some(1.0));
    assert_eq!(atom.b_factor(), Some(10.0));
    let noted = findings
        .iter()
        .filter(|finding| finding.code() == molframe_core::diagnostic::Code::W3015)
        .count();
    assert_eq!(noted, 1);
}

#[test]
fn the_lossless_read_reports_per_model_values_as_the_direct_read_does() {
    let count_w3015 = |findings: &[molframe_core::diagnostic::Diagnostic]| {
        findings
            .iter()
            .filter(|finding| finding.code() == molframe_core::diagnostic::Code::W3015)
            .count()
    };
    let (direct, direct_findings) = parse_structure(PER_MODEL_VALUES);
    let input = InputBuffer::from_bytes(PER_MODEL_VALUES.as_bytes().to_vec());
    let (_, lossless, lossless_findings) = read_with_document(&input, &ReadOptions::new())
        .expect("lossless read of per-model values succeeds");
    assert_eq!(count_w3015(&direct_findings), 1);
    assert_eq!(count_w3015(&lossless_findings), 1);
    assert!(lossless.ragged_models().is_none());
    assert_eq!(lossless.model_count(), direct.model_count());
    let atom = lossless.data().atoms().next().expect("one atom");
    assert_eq!(atom.occupancy(), Some(1.0));
    assert_eq!(atom.b_factor(), Some(10.0));
}

#[test]
fn models_with_the_same_values_report_nothing_on_the_lossless_read() {
    let same = PER_MODEL_VALUES
        .replace("0.50 20.00", "1.00 10.00")
        .replace("0.25 30.00", "1.00 10.00");
    let input = InputBuffer::from_bytes(same.into_bytes());
    let (_, _, findings) = read_with_document(&input, &ReadOptions::new()).expect("reads");
    assert!(
        findings
            .iter()
            .all(|finding| finding.code() != molframe_core::diagnostic::Code::W3015)
    );
}
