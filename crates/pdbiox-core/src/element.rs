//! Chemical element identity, and the fallback for recovering one from an atom
//! name when a file does not declare it.
//!
//! An element is its atomic number in one byte. That fits the seven bits the
//! per-atom column budgets and, more usefully, makes the element its own ordinal
//! in the 128-bit set each chunk carries — so "does this chunk contain zinc?" is
//! one bit test against a summary rather than a scan of every row in it.
//!
//! Element *properties* — radii, mass, electronegativity — are deliberately not
//! here. They are reference data with a version and a provenance, and they
//! belong with the rest of the chemistry.

use std::fmt;

/// Symbols in atomic-number order. Index zero is the unknown element.
const SYMBOLS: [&str; 119] = [
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
    "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
];

/// A symbol reduced to a fixed two-byte key: first letter upper-case, second
/// letter lower-case, zero where there is no second letter.
type Key = [u8; 2];

/// Keys in atomic-number order, derived from [`SYMBOLS`] at compile time.
///
/// Lookup then compares two bytes per candidate with no case folding and no
/// allocation. The scan is linear over 118 entries and runs once per distinct
/// symbol in a file, because parsed symbols are interned.
const KEYS: [Key; 119] = build_keys();

const fn build_keys() -> [Key; 119] {
    let mut keys = [[0u8; 2]; 119];
    let mut z = 1;
    while z < SYMBOLS.len() {
        let bytes = SYMBOLS[z].as_bytes();
        keys[z] = if bytes.len() == 1 {
            [bytes[0], 0]
        } else {
            [bytes[0], bytes[1]]
        };
        z += 1;
    }
    keys
}

/// Reduces up to two alphabetic bytes to a lookup key, or `None` if they are not
/// a plausible symbol.
const fn key_of(bytes: &[u8]) -> Option<Key> {
    match bytes {
        [a] if a.is_ascii_alphabetic() => Some([a.to_ascii_uppercase(), 0]),
        [a, b] if a.is_ascii_alphabetic() && b.is_ascii_alphabetic() => {
            Some([a.to_ascii_uppercase(), b.to_ascii_lowercase()])
        }
        _ => None,
    }
}

/// A chemical element, represented by its atomic number.
///
/// [`Element::UNKNOWN`] means no element could be determined. That is a
/// different statement from an element being absent from a selection, and the
/// two are never conflated.
///
/// # Examples
///
/// ```
/// use pdbiox_core::Element;
///
/// assert_eq!(Element::from_symbol("FE"), Some(Element::IRON));
/// assert_eq!(Element::IRON.symbol(), "Fe");
/// assert_eq!(Element::from_symbol("Unobtainium"), None);
/// ```
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
#[repr(transparent)]
pub struct Element(u8);

impl Element {
    /// No element could be determined.
    pub const UNKNOWN: Self = Self(0);
    /// Hydrogen. Deuterium and tritium resolve here too.
    pub const HYDROGEN: Self = Self(1);
    /// Carbon.
    pub const CARBON: Self = Self(6);
    /// Nitrogen.
    pub const NITROGEN: Self = Self(7);
    /// Oxygen.
    pub const OXYGEN: Self = Self(8);
    /// Phosphorus.
    pub const PHOSPHORUS: Self = Self(15);
    /// Sulfur.
    pub const SULFUR: Self = Self(16);
    /// Calcium — what `CA` means when the atom name starts in the first column.
    pub const CALCIUM: Self = Self(20);
    /// Iron.
    pub const IRON: Self = Self(26);
    /// Zinc.
    pub const ZINC: Self = Self(30);
    /// Selenium, as in selenomethionine.
    pub const SELENIUM: Self = Self(34);

    /// The largest atomic number this type represents.
    pub const MAX_ATOMIC_NUMBER: u8 = 118;

    /// Creates an element from an atomic number, or [`Element::UNKNOWN`] when the
    /// number is past the end of the periodic table.
    #[must_use]
    pub const fn from_atomic_number(z: u8) -> Self {
        if (z as usize) < SYMBOLS.len() {
            Self(z)
        } else {
            Self::UNKNOWN
        }
    }

    /// Returns the atomic number, which doubles as this element's bit position
    /// in a chunk's element set.
    #[must_use]
    pub const fn atomic_number(self) -> u8 {
        self.0
    }

    /// Returns the capitalised symbol, or the empty string when unknown.
    #[must_use]
    pub fn symbol(self) -> &'static str {
        SYMBOLS[self.0 as usize]
    }

    /// Returns true when no element could be determined.
    #[must_use]
    pub const fn is_unknown(self) -> bool {
        self.0 == 0
    }

    /// Returns true for hydrogen.
    ///
    /// Hydrogen presence is summarised per chunk and is what a policy that
    /// excludes or includes hydrogens tests.
    #[must_use]
    pub const fn is_hydrogen(self) -> bool {
        self.0 == 1
    }

    /// Parses a symbol, ignoring surrounding whitespace and letter case.
    ///
    /// `D` and `T` are hydrogen isotopes, not elements, and resolve to
    /// [`Element::HYDROGEN`]; neutron structures spell deuterium that way.
    #[must_use]
    pub fn from_symbol(text: &str) -> Option<Self> {
        let text = text.trim();
        if text.eq_ignore_ascii_case("D") || text.eq_ignore_ascii_case("T") {
            return Some(Self::HYDROGEN);
        }
        let key = key_of(text.as_bytes())?;
        let z = KEYS.iter().position(|candidate| *candidate == key)?;
        Some(Self(z as u8))
    }

    /// Recovers an element from a fixed-column atom-name field.
    ///
    /// `field` must be the four-column name **untrimmed**, because the leading
    /// space is the only thing that distinguishes the two readings:
    ///
    /// - a name that does not start in the first column has a one-letter
    ///   element, taken from its first letter;
    /// - a name that starts in the first column has a two-letter element when its
    ///   first two characters name one, and a one-letter element otherwise.
    ///
    /// This is why `" CA "` in an amino acid is carbon and `"CA  "` in an ion is
    /// calcium. It is only a fallback — a declared element always wins, and
    /// falling back is recorded as a warning, because the rule is a convention
    /// that files are free to violate.
    #[must_use]
    pub fn infer_from_pdb_atom_name(field: &str) -> Self {
        let bytes = field.as_bytes();
        if bytes.first().is_some_and(u8::is_ascii_alphabetic)
            && let Some(two) = bytes.get(..2)
            && let Some(key) = key_of(two)
            && let Some(z) = KEYS.iter().position(|candidate| *candidate == key)
        {
            return Self(z as u8);
        }
        Self::first_letter_element(bytes)
    }

    /// Recovers an element from a name with no column convention.
    ///
    /// Tries the leading two letters, then the leading one. Without a leading
    /// space to disambiguate, `CA` reads as calcium; formats that lack the
    /// convention make the element column mandatory for exactly this reason.
    #[must_use]
    pub fn infer_from_name(name: &str) -> Self {
        let bytes = name.as_bytes();
        let letters = bytes.iter().take_while(|b| b.is_ascii_alphabetic()).count();
        if letters >= 2
            && let Some(two) = bytes.get(..2)
            && let Some(key) = key_of(two)
            && let Some(z) = KEYS.iter().position(|candidate| *candidate == key)
        {
            return Self(z as u8);
        }
        Self::first_letter_element(bytes)
    }

    fn first_letter_element(bytes: &[u8]) -> Self {
        let Some(letter) = bytes.iter().copied().find(u8::is_ascii_alphabetic) else {
            return Self::UNKNOWN;
        };
        let Some(key) = key_of(&[letter]) else {
            return Self::UNKNOWN;
        };
        match KEYS.iter().position(|candidate| *candidate == key) {
            Some(z) => Self(z as u8),
            None => Self::UNKNOWN,
        }
    }
}

/// An element's atomic number is also its bit position in the per-chunk element
/// set, so the periodic table must fit that set for the summary to be usable.
const _: () = assert!((Element::MAX_ATOMIC_NUMBER as usize) < 128);

impl fmt::Debug for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_unknown() {
            f.write_str("Element(unknown)")
        } else {
            write!(f, "Element({})", self.symbol())
        }
    }
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.symbol())
    }
}

#[cfg(test)]
#[path = "element_tests.rs"]
mod tests;
