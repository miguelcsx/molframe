//! A real `AlphaFold` model joined to its published PAE survives a cif and a
//! `BinaryCIF` round trip with every confidence cell intact.

use molframe::formats::modelcif::{ModelCifExt, PaeMatrix, PaePolicy, read_pae_json};
use molframe::{ReadOptions, Structure};
use std::fmt::Write as _;
use std::path::Path;

const DATA: &str = "/tmp/rt-data/af";
const PAIRWISE: &str = "loop_\n_ma_qa_metric_local_pairwise.ordinal_id\n_ma_qa_metric_local_pairwise.model_id\n_ma_qa_metric_local_pairwise.label_asym_id_1\n_ma_qa_metric_local_pairwise.label_seq_id_1\n_ma_qa_metric_local_pairwise.label_comp_id_1\n_ma_qa_metric_local_pairwise.label_asym_id_2\n_ma_qa_metric_local_pairwise.label_seq_id_2\n_ma_qa_metric_local_pairwise.label_comp_id_2\n_ma_qa_metric_local_pairwise.metric_id\n_ma_qa_metric_local_pairwise.metric_value\n";

/// Residue names in sequence order, taken from the file's own pLDDT table.
fn residues(text: &str) -> Vec<String> {
    text.lines()
        .filter(|line| line.starts_with("A ") && line.split_whitespace().count() == 7)
        .filter_map(|line| line.split_whitespace().nth(1).map(str::to_owned))
        .collect()
}

/// Joins the model file and its PAE JSON into one `ModelCIF` document.
fn synthesise(text: &str, matrix: &PaeMatrix) -> String {
    let marker = "2 local  pLDDT 1 pLDDT \n";
    assert!(text.contains(marker), "metric table layout changed");
    let mut joined = text.replacen(marker, &format!("{marker}3 local-pairwise PAE 1 PAE \n"), 1);
    let names = residues(text);
    assert_eq!(names.len(), matrix.size());
    joined.push_str(PAIRWISE);
    for row in 0..matrix.size() {
        for column in 0..matrix.size() {
            let value = matrix.get_f64(row, column).unwrap_or(f64::NAN);
            let _ = writeln!(
                joined,
                "{} 1 A {} {} A {} {} 3 {value}",
                row * matrix.size() + column + 1,
                row + 1,
                names[row],
                column + 1,
                names[column],
            );
        }
    }
    joined.push_str("#\n");
    joined
}

fn read(bytes: Vec<u8>, name: &str) -> Structure {
    match molframe::read_bytes(bytes, Some(name), &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(failure) => panic!("{name}: {failure:?}"),
    }
}

fn assert_same(structure: &Structure, expected: &PaeMatrix, label: &str) {
    let dense = structure
        .model_cif()
        .and_then(|model| model.dense_pae())
        .unwrap_or_else(|| panic!("{label}: no dense PAE"));
    let got = dense.matrix();
    assert_eq!(got.size(), expected.size(), "{label}");
    for row in 0..expected.size() {
        for column in 0..expected.size() {
            assert_eq!(
                got.get_f64(row, column),
                expected.get_f64(row, column),
                "{label} ({row},{column})"
            );
        }
    }
}

fn check(accession: &str) {
    let base = Path::new(DATA);
    let model = base.join(format!("AF-{accession}-F1-model_v6.cif"));
    let json = base.join(format!("AF-{accession}-F1-predicted_aligned_error_v6.json"));
    let (Ok(text), Ok(json)) = (std::fs::read_to_string(model), std::fs::read(json)) else {
        return;
    };
    let matrix = read_pae_json(&json, PaePolicy::new()).expect("pae json");
    let joined = synthesise(&text, &matrix);
    let structure = read(joined.clone().into_bytes(), "joined.cif");
    assert_same(&structure, &matrix, "cif read");

    let cif = molframe::write_mmcif(&structure).expect("write cif");
    assert_same(&read(cif.into_bytes(), "again.cif"), &matrix, "cif write");

    let bcif = molframe::write_bcif(&structure).expect("write bcif");
    let from_bcif = read(bcif.clone(), "again.bcif");
    assert_same(&from_bcif, &matrix, "bcif");
    assert!(
        from_bcif
            .engine()
            .annotations()
            .get(molframe::PLDDT_ANNOTATION)
            .is_some(),
        "pLDDT lost through bcif"
    );
    let second = molframe::write_bcif(&from_bcif).expect("rewrite bcif");
    assert_eq!(bcif, second, "bcif output must be stable");
}

#[test]
fn hemoglobin_alpha_round_trips() {
    check("P69905");
}

#[test]
fn p53_round_trips() {
    check("P04637");
}

#[test]
fn egfr_round_trips() {
    check("P00533");
}
