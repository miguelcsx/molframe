//! Catalogued operations of the `formats` namespace.

use super::{Capability, Domain};
use molframe::Cost;

const D: Domain = Domain("formats");

pub(super) const ROWS: &[Capability] = &[
    D.row("to_mmcif", "structure", "text", "cpu,memory"),
    D.row("to_pdb", "structure", "text", "cpu,memory"),
    D.row("to_bcif", "structure", "bytes", "cpu,memory"),
    D.row("write", "structure, path", "file", "cpu,memory"),
    D.row("preserve_mmcif", "path, path", "file", "cpu,memory")
        .cost(Cost::Copy),
    Domain("modelcif")
        .row("pae", "structure", "PaeMatrix", "cpu")
        .cost(Cost::Borrow)
        .feature("formats"),
    Domain("modelcif")
        .row("plddt_band_colors", "none", "uint8[4,3]", "cpu")
        .feature("formats"),
    Domain("modelcif")
        .row("plddt_band_names", "none", "str[]", "cpu")
        .feature("formats"),
    Domain("modelcif")
        .row("plddt_bands", "structure", "uint8[]", "cpu")
        .feature("formats"),
    Domain("modelcif")
        .row("read_pae_json", "path", "PaeMatrix", "cpu,memory")
        .cost(Cost::Decode)
        .feature("formats"),
];
