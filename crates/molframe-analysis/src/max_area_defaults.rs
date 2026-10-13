//! The built-in maximum residue areas, as data.
//!
//! This is the one place the numbers live; the algorithm that divides an area
//! by a maximum never names them. A caller who wants other values supplies a
//! [`MaxAreaTable`](crate::MaxAreaTable) and leaves this module alone.

use molframe_chem::TableIdentity;

/// The table's name as a report records it.
pub const MAX_RESIDUE_AREA_REFERENCE: &str =
    "Tien et al. 2013, theoretical maximum ASA (Gly-X-Gly extended, all-atom)";

/// Where the table comes from.
pub const MAX_RESIDUE_AREA_CITATION: &str = "Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO (2013) Maximum allowed solvent accessibilities of residues in proteins. PLoS ONE 8(11): e80635";

/// Identity of the built-in table. The values are published measurements of a
/// physical quantity; they are reproduced here with the citation.
pub const TIEN_2013: TableIdentity = TableIdentity::builtin(
    "tien-2013-theoretical",
    "2013",
    MAX_RESIDUE_AREA_CITATION,
    "values from a published article, reproduced with citation",
);

/// Theoretical maximum accessible area of each standard amino acid, Å².
pub const TIEN_2013_AREAS: &[(&str, f64)] = &[
    ("ALA", 129.0),
    ("ARG", 274.0),
    ("ASN", 195.0),
    ("ASP", 193.0),
    ("CYS", 167.0),
    ("GLN", 225.0),
    ("GLU", 223.0),
    ("GLY", 104.0),
    ("HIS", 224.0),
    ("ILE", 197.0),
    ("LEU", 201.0),
    ("LYS", 236.0),
    ("MET", 224.0),
    ("PHE", 240.0),
    ("PRO", 159.0),
    ("SER", 155.0),
    ("THR", 172.0),
    ("TRP", 285.0),
    ("TYR", 263.0),
    ("VAL", 174.0),
];
