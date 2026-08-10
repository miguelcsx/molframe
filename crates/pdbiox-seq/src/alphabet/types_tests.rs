use super::{Alphabet, CustomAlphabet, DNA, PROTEIN, Sequence};

#[test]
fn dna_encodes_to_dense_indices() {
    let Some(codes) = DNA.encode(b"ACGT") else {
        panic!("all bases are in the DNA alphabet");
    };
    assert_eq!(codes, vec![0, 1, 2, 3]);
}

#[test]
fn encoding_then_decoding_is_the_identity() {
    let sequence = b"MKVLA";
    let Some(codes) = PROTEIN.encode(sequence) else {
        panic!("all residues are standard");
    };
    let Some(back) = PROTEIN.decode(&codes) else {
        panic!("codes are in range");
    };
    assert_eq!(back, sequence);
}

#[test]
fn a_symbol_outside_the_alphabet_fails_to_encode() {
    assert!(DNA.encode(b"ACXT").is_none());
}

#[test]
fn a_code_out_of_range_fails_to_decode() {
    assert!(DNA.decode(&[0, 9]).is_none());
}

#[test]
fn the_alphabets_have_their_expected_sizes() {
    assert_eq!(DNA.len(), 4);
    assert_eq!(PROTEIN.len(), 27);
}

#[test]
fn a_typed_sequence_encodes_once_and_decodes_losslessly() {
    let Ok(sequence) = Sequence::new(PROTEIN, b"MZ-U") else {
        panic!("extended protein symbols are valid");
    };
    assert_eq!(sequence.symbols(), b"MZ-U");
    assert_eq!(sequence.codes().len(), 4);
}

#[test]
fn custom_alphabets_refuse_duplicate_symbols() {
    assert!(CustomAlphabet::new(&b"AAB"[..]).is_err());
}
