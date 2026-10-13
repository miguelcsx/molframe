//! Corrupt, truncated and hostile mmCIF text must error, never panic.

use molframe_cif::{parse, read};
use molframe_core::io::{InputBuffer, ReadOptions};
use proptest::prelude::*;

const SAMPLE: &[u8] = include_bytes!("../../molframe-bench/data/1crn.cif");

fn exercise(bytes: Vec<u8>) {
    let input = InputBuffer::from_bytes(bytes);
    let _ = parse(&input);
    let _ = read(&input, &ReadOptions::default());
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(150))]

    #[test]
    fn truncation_never_panics(cut in 0..SAMPLE.len()) {
        exercise(SAMPLE[..cut].to_vec());
    }

    #[test]
    fn byte_flips_never_panic(
        edits in proptest::collection::vec((0..SAMPLE.len(), any::<u8>()), 1..32),
    ) {
        let mut bytes = SAMPLE.to_vec();
        for (at, value) in edits {
            bytes[at] = value;
        }
        exercise(bytes);
    }

    #[test]
    fn spliced_slices_never_panic(a in 0..SAMPLE.len(), len in 0usize..4096, b in 0..SAMPLE.len()) {
        let mut bytes = SAMPLE[..a].to_vec();
        bytes.extend_from_slice(&SAMPLE[b..SAMPLE.len().min(b + len)]);
        exercise(bytes);
    }

    #[test]
    fn arbitrary_text_never_panics(text in "[ -~\n\t'\";#_.\\[\\]?-]{0,600}") {
        exercise(text.into_bytes());
    }
}

#[test]
fn pathological_shapes_error_or_parse() {
    let shapes = [
        "data_x\nloop_\n".to_owned(),
        format!("data_x\n_a.b '{}", "q".repeat(1 << 20)),
        format!("data_x\n;{}", "line\n".repeat(100_000)),
        format!("data_x\nloop_\n{}", "_a.b\n".repeat(50_000)),
        format!("{}\n", "data_x\n".repeat(50_000)),
        format!("data_x\nloop_\n_a.b\n{}", "1 ".repeat(200_000)),
        "data_x\n_atom_site.Cartn_x 1e99999999\n".to_owned(),
    ];
    for shape in shapes {
        exercise(shape.into_bytes());
    }
}
