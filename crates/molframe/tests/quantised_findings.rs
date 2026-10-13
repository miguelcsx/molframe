//! A PAE matrix stored under an approximate policy is reported by both
//! CIF-family writers, and an exact one is not.

use molframe::ReadOptions;
use molframe::formats::modelcif::{ModelCifExt, PaeKind, PaeLayout, PaeMatrix};
use molframe_cif::CifWriteOptions;
use std::fmt::Write as _;

const SIZE: usize = 4;

fn document(fraction: &str) -> String {
    let mut text = String::from(
        "data_q\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.pdbx_PDB_model_num\n",
    );
    for index in 0..SIZE {
        let _ = writeln!(
            text,
            "ATOM {} C CA GLY A {} {index}.0 0.0 0.0 1",
            index + 1,
            index + 1
        );
    }
    text.push_str("#\nloop_\n_ma_qa_metric.id\n_ma_qa_metric.name\n_ma_qa_metric.type\n_ma_qa_metric.mode\n_ma_qa_metric.software_group_id\n2 PAE PAE local-pairwise 1\n#\nloop_\n");
    for item in [
        "ordinal_id",
        "model_id",
        "label_asym_id_1",
        "label_seq_id_1",
        "label_comp_id_1",
        "label_asym_id_2",
        "label_seq_id_2",
        "label_comp_id_2",
        "metric_id",
        "metric_value",
    ] {
        let _ = writeln!(text, "_ma_qa_metric_local_pairwise.{item}");
    }
    for row in 0..SIZE {
        for column in 0..SIZE {
            let _ = writeln!(
                text,
                "{} 1 A {} GLY A {} GLY 2 {}{fraction}",
                row * SIZE + column + 1,
                row + 1,
                column + 1,
                (row + column) % 9
            );
        }
    }
    text.push_str("#\n");
    text
}

fn structure(fraction: &str, approximate: bool) -> molframe::Structure {
    let text = document(fraction);
    let (structure, _) =
        molframe::read_bytes(text.into_bytes(), Some("q.cif"), &ReadOptions::new()).expect("read");
    let mut model = structure.model_cif().expect("model").clone();
    if approximate {
        let layout = PaeLayout {
            kind: PaeKind::U8,
            denominator: 1,
            decimals: None,
        };
        let mut matrix = PaeMatrix::zeroed(SIZE, layout, 1 << 20).expect("matrix");
        for row in 0..SIZE {
            for column in 0..SIZE {
                matrix.set_f64(row, column, 1.5).expect("cell");
            }
        }
        assert!(!matrix.is_exact());
        model.replace_pae_matrix(matrix).expect("replace");
    }
    structure
        .engine()
        .with_extension(molframe::formats::modelcif::MODEL_CIF_EXTENSION, model)
        .into()
}

fn mentions_quantised(structure: &molframe::Structure) -> (bool, bool) {
    let (_, cif) =
        molframe::write_mmcif_with_findings(structure, &CifWriteOptions::new()).expect("cif");
    let (_, bcif) =
        molframe::write_bcif_with_findings(structure, &CifWriteOptions::new()).expect("bcif");
    let has = |findings: &molframe::Findings| {
        findings
            .to_vec()
            .iter()
            .any(|finding| finding.message().contains("quantised"))
    };
    (has(&cif), has(&bcif))
}

#[test]
fn an_approximate_matrix_is_reported_by_both_writers() {
    let approximate = structure(".5", true);
    assert_eq!(mentions_quantised(&approximate), (true, true));
}

#[test]
fn an_exact_matrix_is_not_reported() {
    let exact = structure(".5", false);
    assert_eq!(mentions_quantised(&exact), (false, false));
}
