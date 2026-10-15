//! Typed sequence alphabets.

mod types;

pub use types::{
    Alphabet, AlphabetError, Code, CustomAlphabet, DNA, DnaAlphabet, NUCLEOTIDE,
    NucleotideAlphabet, PROTEIN, ProteinAlphabet, RNA, RnaAlphabet, Sequence, SequenceError,
};
