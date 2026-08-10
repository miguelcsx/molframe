//! Typed alphabets and sequences encoded once into dense integer codes.

use std::fmt::{Display, Formatter};

/// Dense symbol code used by sequence kernels.
pub type Code = u8;

/// An ordered sequence alphabet.
pub trait Alphabet: Clone {
    /// Symbols in code order.
    fn symbols(&self) -> &[u8];

    /// Dense code for one symbol.
    fn encode_symbol(&self, symbol: u8) -> Option<Code> {
        self.symbols()
            .iter()
            .position(|candidate| *candidate == symbol)
            .and_then(|position| Code::try_from(position).ok())
    }

    /// Symbol represented by one dense code.
    fn decode_symbol(&self, code: Code) -> Option<u8> {
        self.symbols().get(usize::from(code)).copied()
    }

    /// Encodes all symbols, refusing the first value outside the alphabet.
    fn encode(&self, sequence: &[u8]) -> Option<Vec<Code>> {
        sequence
            .iter()
            .map(|symbol| self.encode_symbol(*symbol))
            .collect()
    }

    /// Decodes all codes, refusing the first value outside the alphabet.
    fn decode(&self, codes: &[Code]) -> Option<Vec<u8>> {
        codes.iter().map(|code| self.decode_symbol(*code)).collect()
    }

    /// Number of symbols.
    fn len(&self) -> usize {
        self.symbols().len()
    }

    /// Whether the alphabet has no symbols.
    fn is_empty(&self) -> bool {
        self.symbols().is_empty()
    }
}

macro_rules! static_alphabet {
    ($name:ident, $symbols:expr) => {
        #[doc = concat!("The `", stringify!($name), "` built-in alphabet.")]
        #[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
        pub struct $name;

        impl Alphabet for $name {
            fn symbols(&self) -> &[u8] {
                $symbols
            }
        }
    };
}

static_alphabet!(ProteinAlphabet, b"ACDEFGHIKLMNPQRSTVWYBZXUO-*" as &[u8]);
static_alphabet!(DnaAlphabet, b"ACGT" as &[u8]);
static_alphabet!(RnaAlphabet, b"ACGU" as &[u8]);
static_alphabet!(NucleotideAlphabet, b"ACGTURYSWKMBDHVN-" as &[u8]);

/// Protein alphabet value for concise call sites.
pub const PROTEIN: ProteinAlphabet = ProteinAlphabet;
/// DNA alphabet value for concise call sites.
pub const DNA: DnaAlphabet = DnaAlphabet;
/// RNA alphabet value for concise call sites.
pub const RNA: RnaAlphabet = RnaAlphabet;
/// DNA/RNA ambiguity alphabet value for concise call sites.
pub const NUCLEOTIDE: NucleotideAlphabet = NucleotideAlphabet;

/// A validated runtime alphabet for reduced or structural representations.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CustomAlphabet {
    symbols: Box<[u8]>,
}

impl CustomAlphabet {
    /// Builds a non-empty alphabet with at most 256 unique symbols.
    ///
    /// # Errors
    ///
    /// Returns [`AlphabetError`] for an empty, oversized, or duplicate set.
    pub fn new(symbols: impl Into<Box<[u8]>>) -> Result<Self, AlphabetError> {
        let symbols = symbols.into();
        if symbols.is_empty() {
            return Err(AlphabetError::Empty);
        }
        if symbols.len() > usize::from(Code::MAX) + 1 {
            return Err(AlphabetError::TooManySymbols {
                count: symbols.len(),
            });
        }
        let mut present = [false; 256];
        for symbol in &symbols {
            let slot = &mut present[usize::from(*symbol)];
            if *slot {
                return Err(AlphabetError::DuplicateSymbol { symbol: *symbol });
            }
            *slot = true;
        }
        Ok(Self { symbols })
    }
}

impl Alphabet for CustomAlphabet {
    fn symbols(&self) -> &[u8] {
        &self.symbols
    }
}

/// A sequence whose symbols are validated and densely encoded at construction.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Sequence<A: Alphabet> {
    alphabet: A,
    codes: Vec<Code>,
}

impl<A: Alphabet> Sequence<A> {
    /// Validates and encodes raw symbols exactly once.
    ///
    /// # Errors
    ///
    /// Returns [`SequenceError`] at the first symbol outside the alphabet.
    pub fn new(alphabet: A, symbols: &[u8]) -> Result<Self, SequenceError> {
        let mut codes = Vec::with_capacity(symbols.len());
        for (position, symbol) in symbols.iter().copied().enumerate() {
            let Some(code) = alphabet.encode_symbol(symbol) else {
                return Err(SequenceError::UnknownSymbol { position, symbol });
            };
            codes.push(code);
        }
        Ok(Self { alphabet, codes })
    }

    /// Validates already encoded codes without re-encoding them.
    ///
    /// # Errors
    ///
    /// Returns [`SequenceError`] at the first code outside the alphabet.
    pub fn from_codes(alphabet: A, codes: Vec<Code>) -> Result<Self, SequenceError> {
        for (position, code) in codes.iter().copied().enumerate() {
            if alphabet.decode_symbol(code).is_none() {
                return Err(SequenceError::UnknownCode { position, code });
            }
        }
        Ok(Self { alphabet, codes })
    }

    /// Encoded storage borrowed without allocation.
    #[must_use]
    pub fn codes(&self) -> &[Code] {
        &self.codes
    }

    /// Alphabet that defines the code mapping.
    #[must_use]
    pub const fn alphabet(&self) -> &A {
        &self.alphabet
    }

    /// Number of encoded symbols.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.codes.len()
    }

    /// Whether the sequence carries no symbols.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.codes.is_empty()
    }

    /// Decodes to display/file symbols.
    #[must_use]
    pub fn symbols(&self) -> Vec<u8> {
        self.codes
            .iter()
            .map(|code| self.alphabet.symbols()[usize::from(*code)])
            .collect()
    }
}

/// Invalid alphabet definition.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AlphabetError {
    /// No symbols were supplied.
    Empty,
    /// A byte code cannot address the requested symbol count.
    TooManySymbols {
        /// Requested count.
        count: usize,
    },
    /// One symbol appears more than once.
    DuplicateSymbol {
        /// Repeated byte.
        symbol: u8,
    },
}

/// Invalid raw or encoded sequence content.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SequenceError {
    /// A raw symbol is outside the alphabet.
    UnknownSymbol {
        /// Zero-based sequence position.
        position: usize,
        /// Rejected byte.
        symbol: u8,
    },
    /// A dense code is outside the alphabet.
    UnknownCode {
        /// Zero-based sequence position.
        position: usize,
        /// Rejected code.
        code: Code,
    },
}

impl Display for AlphabetError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Empty => formatter.write_str("alphabet must contain at least one symbol"),
            Self::TooManySymbols { count } => write!(formatter, "alphabet has {count} symbols"),
            Self::DuplicateSymbol { symbol } => write!(formatter, "duplicate symbol {symbol}"),
        }
    }
}

impl Display for SequenceError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnknownSymbol { position, symbol } => {
                write!(
                    formatter,
                    "symbol {symbol} at position {position} is outside the alphabet"
                )
            }
            Self::UnknownCode { position, code } => {
                write!(
                    formatter,
                    "code {code} at position {position} is outside the alphabet"
                )
            }
        }
    }
}

impl std::error::Error for AlphabetError {}
impl std::error::Error for SequenceError {}

#[cfg(test)]
#[path = "types_tests.rs"]
mod tests;
