//! A dense PAE matrix is written to `BinaryCIF` column by column and read
//! back exactly, within a debug-build time budget.

use molframe::ReadOptions;
use molframe::formats::modelcif::ModelCifExt;
use std::fmt::Write as _;
use std::time::{Duration, Instant};

const SIZE: usize = 160;

fn document() -> String {
    let mut text = String::from(
        "data_scale\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.pdbx_PDB_model_num\n",
    );
    for index in 0..SIZE {
        let _ = writeln!(
            text,
            "ATOM {} C CA GLY A {} {}.0 0.0 0.0 1",
            index + 1,
            index + 1,
            index
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
                "{} 1 A {} GLY A {} GLY 2 {}.25",
                row * SIZE + column + 1,
                row + 1,
                column + 1,
                (row * 7 + column * 3) % 31
            );
        }
    }
    text.push_str("#\n");
    text
}

#[test]
fn dense_pae_survives_a_direct_bcif_write_quickly() {
    let (structure, _) = molframe::read_bytes(
        document().into_bytes(),
        Some("scale.cif"),
        &ReadOptions::new(),
    )
    .expect("read");
    let dense = structure
        .model_cif()
        .and_then(|model| model.dense_pae())
        .expect("dense PAE");
    assert_eq!(dense.matrix().size(), SIZE);

    let started = Instant::now();
    let bytes = molframe::write_bcif(&structure).expect("write bcif");
    assert!(
        started.elapsed() < Duration::from_secs(5),
        "bcif write took {:?}",
        started.elapsed()
    );

    let (again, _) =
        molframe::read_bytes(bytes, Some("scale.bcif"), &ReadOptions::new()).expect("reread");
    let back = again
        .model_cif()
        .and_then(|model| model.dense_pae())
        .expect("dense PAE after bcif");
    for row in 0..SIZE {
        for column in 0..SIZE {
            assert_eq!(
                back.matrix().get_f64(row, column),
                dense.matrix().get_f64(row, column),
                "({row},{column})"
            );
        }
    }
}
