//! Corrupt, truncated and hostile crystallographic inputs must error, never
//! panic or allocate what they merely claim.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss
)]

use molframe_cif::parse;
use molframe_core::io::InputBuffer;
use molframe_xtal::{
    DensityMap, MrcBlockOptions, MrcBlockReader, OperExpression, SymmetryOperation,
    lower_assemblies, lower_ncs, lower_structure_factor_cif, lower_symmetry, read_cube, read_dx,
    read_mtz,
};
use std::io::Cursor;

fn mrc_header(dims: [i32; 3], mode: i32, extra: usize) -> Vec<u8> {
    let mut bytes = vec![0_u8; 1024 + extra];
    for (slot, value) in dims.iter().enumerate() {
        bytes[slot * 4..slot * 4 + 4].copy_from_slice(&value.to_le_bytes());
    }
    bytes[12..16].copy_from_slice(&mode.to_le_bytes());
    for (slot, value) in dims.iter().enumerate() {
        bytes[28 + slot * 4..32 + slot * 4].copy_from_slice(&value.to_le_bytes());
    }
    for slot in 0..3 {
        bytes[40 + slot * 4..44 + slot * 4].copy_from_slice(&1.0_f32.to_le_bytes());
        bytes[64 + slot * 4..68 + slot * 4].copy_from_slice(&(slot as i32 + 1).to_le_bytes());
        bytes[52 + slot * 4..56 + slot * 4].copy_from_slice(&90.0_f32.to_le_bytes());
    }
    bytes[208..212].copy_from_slice(b"MAP ");
    bytes[212..216].copy_from_slice(&[0x44, 0x41, 0, 0]);
    bytes
}

#[test]
fn an_mrc_declaring_giant_dimensions_errors_without_allocating() {
    for dims in [
        [i32::MAX; 3],
        [1 << 20, 1 << 20, 1 << 20],
        [65_536, 65_536, 2],
    ] {
        let bytes = mrc_header(dims, 2, 64);
        assert!(DensityMap::from_mrc_bytes(&bytes).is_err());
        let options = MrcBlockOptions::default().with_memory_limit(1 << 20);
        let _ = MrcBlockReader::new(Cursor::new(bytes), options);
    }
}

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn below(&mut self, bound: u64) -> u64 {
        self.next() % bound.max(1)
    }

    fn text(&mut self, alphabet: &[u8], max: u64) -> String {
        let len = self.below(max);
        (0..len)
            .map(|_| char::from(alphabet[self.below(alphabet.len() as u64) as usize]))
            .collect()
    }
}

const PRINTABLE: &[u8] = b" abcxyz0123456789,()-:+*/'.;_\n#?";

#[test]
fn mrc_headers_never_panic() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    for _ in 0..600 {
        let dims = [0, 1, 2].map(|_| rng.next() as i32);
        let mode = rng.below(22) as i32 - 2;
        let mut bytes = mrc_header(dims, mode, rng.below(128) as usize);
        for _ in 0..rng.below(12) {
            let at = rng.below(1100) as usize;
            if at < bytes.len() {
                bytes[at] = rng.next() as u8;
            }
        }
        let _ = DensityMap::from_mrc_bytes(&bytes);
        let options = MrcBlockOptions::default().with_memory_limit(1 << 20);
        let _ = MrcBlockReader::new(Cursor::new(bytes), options);
    }
}

#[test]
fn mtz_bytes_never_panic() {
    let mut rng = Rng(0x1234_5678_9abc_def1);
    for _ in 0..600 {
        let mut bytes = b"MTZ ".to_vec();
        let word = rng.below(402) as i32 - 1;
        bytes.extend_from_slice(&word.to_le_bytes());
        bytes.extend_from_slice(&[0x44, 0x41, 0, 0]);
        for _ in 0..rng.below(1500) {
            bytes.push(rng.next() as u8);
        }
        let _ = read_mtz(&bytes);
    }
}

#[test]
fn grid_and_expression_text_never_panics() {
    let mut rng = Rng(0xdead_beef_cafe_f00d);
    for _ in 0..1500 {
        let text = rng.text(PRINTABLE, 400);
        let _ = read_cube(text.as_bytes());
        let _ = read_dx(text.as_bytes());
        let _ = OperExpression::parse(&text);
        let _ = OperExpression::parse_with_limit(&text, 64);
        let _ = SymmetryOperation::parse("1", &text);
    }
}

#[test]
fn cif_lowerers_never_panic() {
    let categories = [
        "pdbx_struct_assembly_gen",
        "pdbx_struct_oper_list",
        "struct_ncs_oper",
        "symmetry",
        "space_group_symop",
        "refln",
        "cell",
    ];
    let mut rng = Rng(0x0bad_cafe_1357_9bdf);
    for round in 0..600_usize {
        let category = categories[round % categories.len()];
        let mut rows = String::new();
        for _ in 0..rng.below(12) {
            rows.push_str(&rng.text(PRINTABLE, 60));
            rows.push('\n');
        }
        let text = format!(
            "data_x\nloop_\n_{category}.id\n_{category}.details\n_{category}.oper_expression\n_{category}.asym_id_list\n{rows}\n"
        );
        let input = InputBuffer::from_bytes(text.into_bytes());
        if let Ok((document, _)) = parse(&input) {
            let _ = lower_assemblies(&document);
            let _ = lower_ncs(&document);
            let _ = lower_symmetry(&document);
            let _ = lower_structure_factor_cif(&document);
        }
    }
}

#[test]
fn giant_operator_expressions_stay_within_the_instance_limit() {
    let expression = "(1-60)(1-60)(1-60)(1-60)(1-60)(1-60)";
    assert!(OperExpression::parse(expression).is_err());
    let range = "(1-4294967295)";
    assert!(OperExpression::parse(range).is_err());
}

#[test]
fn cube_and_dx_headers_with_giant_counts_error() {
    let cube = "c\nc\n999999999 0 0 0\n999999999 1 0 0\n999999999 0 1 0\n999999999 0 0 1\n";
    assert!(read_cube(cube.as_bytes()).is_err());
    let dx = "object 1 class gridpositions counts 999999999 999999999 999999999\norigin 0 0 0\ndelta 1 0 0\ndelta 0 1 0\ndelta 0 0 1\nobject 2 class array type double rank 0 items 99999999999 data follows\n1\n";
    assert!(read_dx(dx.as_bytes()).is_err());
}
