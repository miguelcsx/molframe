//! Corrupt, truncated and hostile PDB text must error, never panic.

use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_pdb::{read, read_pdbqt, read_pqr};
use proptest::prelude::*;

const SAMPLE: &[u8] = include_bytes!("../../molframe-bench/data/4hhb.pdb");

fn exercise(bytes: Vec<u8>) {
    let input = InputBuffer::from_bytes(bytes);
    let options = ReadOptions::default();
    let _ = read(&input, &options);
    let _ = read_pqr(&input, &options);
    let _ = read_pdbqt(&input, &options);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    #[test]
    fn truncation_never_panics(cut in 0..SAMPLE.len()) {
        exercise(SAMPLE[..cut].to_vec());
    }

    #[test]
    fn byte_flips_never_panic(
        edits in proptest::collection::vec((0..SAMPLE.len(), any::<u8>()), 1..48),
    ) {
        let mut bytes = SAMPLE.to_vec();
        for (at, value) in edits {
            bytes[at] = value;
        }
        exercise(bytes);
    }

    #[test]
    fn arbitrary_lines_never_panic(
        lines in proptest::collection::vec("[A-Z0-9 .+-]{0,90}", 0..40),
    ) {
        exercise(lines.join("\n").into_bytes());
    }

    #[test]
    fn record_prefixed_noise_never_panics(
        record in prop::sample::select(vec![
            "ATOM  ", "HETATM", "CONECT", "MODEL ", "ENDMDL", "TER   ", "CRYST1",
            "ANISOU", "SCALE1", "MTRIX1", "REMARK", "HELIX ", "SHEET ", "SSBOND",
        ]),
        tail in "[ -~]{0,80}",
        repeat in 1usize..30,
    ) {
        exercise(format!("{record}{tail}\n").repeat(repeat).into_bytes());
    }
}

#[test]
fn pathological_shapes_error_or_parse() {
    let shapes = [
        "MODEL     99999999999\n".repeat(1000),
        format!("CONECT{}\n", "99999".repeat(40)),
        "ENDMDL\n".repeat(100_000),
        format!("ATOM  {}\n", "9".repeat(5000)),
        "\u{feff}ATOM  \u{1f600}\u{1f600}\u{1f600}\n".to_owned(),
        "\r\r\r\n\n\0\0\0".to_owned(),
    ];
    for shape in shapes {
        exercise(shape.into_bytes());
    }
}
