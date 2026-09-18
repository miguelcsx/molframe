//! Interned identifiers.
//!
//! No kernel compares strings. Every textual identifier a file carries — atom
//! name, residue code, chain label — is interned once at parse time and becomes
//! a 32-bit integer, so `name CA and resname ALA` reduces to two integer
//! comparisons per atom instead of two string comparisons.
//!
//! Identifiers come from two dictionaries. The canonical one is compiled in and
//! its identifiers are the same in every structure and every release. The local
//! one belongs to a single structure and covers everything else: ligand codes,
//! chain labels, unusual atom names.
//!
//! Local strings live in one arena rather than in one allocation each. A file
//! with a hundred thousand atoms has a few hundred distinct identifiers, and
//! interning them costs a few hundred bytes copied once, not one allocation per
//! atom.

mod altloc;
mod canonical;
mod interner;

pub use altloc::*;
pub use interner::*;
