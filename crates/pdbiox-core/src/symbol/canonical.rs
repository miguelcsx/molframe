//! The compiled-in half of the identifier dictionary.
//!
//! These are the strings that recur in every structure ever deposited: element
//! symbols, the standard residue codes, and the atom names of the standard
//! components. Giving them fixed identifiers has a consequence beyond saving a
//! little memory — two structures agree on what `CA` means without consulting
//! each other, so comparing them is integer work, and a query compiled against
//! one structure can be replayed against another.
//!
//! The list is **append-only**. Reassigning an identifier would silently change
//! the meaning of a serialised query or a cached plan, so new strings go on the
//! end and nothing is ever removed or reordered.

use hashbrown::HashTable;
use std::hash::{BuildHasher, RandomState};
use std::sync::OnceLock;

/// Strings with fixed identifiers, in identifier order. Never reorder this.
pub(super) static CANONICAL: &[&str] = &[
    // Element symbols as files spell them, upper-case.
    "H", "HE", "LI", "BE", "B", "C", "N", "O", "F", "NE", "NA", "MG", "AL", "SI", "P", "S", "CL",
    "AR", "K", "CA", "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS",
    "SE", "BR", "KR", "RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN",
    "SN", "SB", "TE", "I", "XE", "CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB",
    "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG", "TL",
    "PB", "BI", "PO", "AT", "RN", "FR", "RA", "AC", "TH", "PA", "U", "NP", "PU", "AM", "CM", "BK",
    "CF", "ES", "FM", "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG", "CN", "NH",
    "FL", "MC", "LV", "TS", "OG",
    // Standard amino acids, the common non-standard ones, and water.
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE", "SEC", "PYL", "ASX", "GLX", "UNK",
    "HOH", "WAT", "DOD", // Nucleotides, ribo and deoxy.
    "DA", "DC", "DG", "DT", "DU", "DI", "A", "G", "T", "UNL",
    // Protein backbone and side-chain atom names. `N`, `CA`, `C`, `O`, `CD`,
    // `CE`, `NE`, `OG` and `SG` are already declared above as element symbols:
    // one string is one identifier, and which reading applies is decided by the
    // column it was read from, not by the dictionary.
    "CB", "CG", "CG1", "CG2", "OG1", "CD1", "CD2", "ND1", "ND2", "OD1", "OD2", "SD", "CE1", "CE2",
    "CE3", "NE1", "NE2", "OE1", "OE2", "CZ", "CZ2", "CZ3", "NZ", "OH", "NH1", "NH2", "CH2", "OXT",
    "AD1", "AD2", "AE1", "AE2", // Nucleic acid atom names.
    "OP1", "OP2", "OP3", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N1", "C2",
    "O2", "N3", "C4", "N4", "O4", "C5", "C6", "N6", "O6", "N7", "C8", "N9", "N2", "C7", "O1P",
    "O2P", "O3P",
    // Hydrogen names common enough to be worth a fixed identifier. `HG` is
    // already declared above as mercury's symbol.
    "H1", "H2", "H3", "HA", "HA2", "HA3", "HB", "HB1", "HB2", "HB3", "HG1", "HG2", "HG3", "HD1",
    "HD2", "HD3", "HE1", "HE2", "HE3", "HZ", "HZ1", "HZ2", "HZ3", "HH", "HH2", "HN", "HXT",
    // Alternate-location labels and single-character chain identifiers.
    "D", "E", "L", "M", "Q", "R", "X", "Z",
];

/// Lazily built lookup from string to canonical identifier.
///
/// The table cannot be a sorted array searched by comparison, because that would
/// force the list into alphabetical order and identifiers must follow
/// declaration order to stay stable. It is built once for the process.
static LOOKUP: OnceLock<Lookup> = OnceLock::new();

struct Lookup {
    table: HashTable<u32>,
    hasher: RandomState,
}

fn lookup() -> &'static Lookup {
    LOOKUP.get_or_init(|| {
        let hasher = RandomState::new();
        let mut table = HashTable::with_capacity(CANONICAL.len());
        for (ordinal, text) in CANONICAL.iter().enumerate() {
            let hash = hasher.hash_one(*text);
            // A duplicate would make one of the two identifiers unreachable; the
            // test below proves there are none, so the later entry is dropped.
            table
                .entry(
                    hash,
                    |&other| CANONICAL[other as usize] == *text,
                    |&other| hasher.hash_one(CANONICAL[other as usize]),
                )
                .or_insert(ordinal as u32);
        }
        Lookup { table, hasher }
    })
}

/// Returns the canonical identifier for `text`, if it has one.
pub(super) fn ordinal_of(text: &str) -> Option<u32> {
    let lookup = lookup();
    let hash = lookup.hasher.hash_one(text);
    lookup
        .table
        .find(hash, |&ordinal| CANONICAL[ordinal as usize] == text)
        .copied()
}

/// Returns the string a canonical identifier names.
pub(super) fn text_of(ordinal: u32) -> Option<&'static str> {
    CANONICAL.get(ordinal as usize).copied()
}

#[cfg(test)]
#[path = "canonical_tests.rs"]
mod tests;
