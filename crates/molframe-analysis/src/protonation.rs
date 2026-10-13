//! Maps protonation and disulphide variants of the standard amino acids to
//! their parent residue, so per-residue tables are keyed once.
//!
//! The variants are the AMBER and CHARMM names that molecular-dynamics
//! preparation writes into coordinate files: `HID HIE HIP` and `HSD HSE HSP`
//! (histidine), `CYX CYM` (cysteine), `ASH` (aspartate), `GLH` (glutamate) and
//! `LYN` (lysine). Lookup is a single match, constant time.

/// The parent residue of a protonation variant, or the name itself.
#[must_use]
pub(crate) fn parent_residue(name: &str) -> &str {
    match name {
        "HID" | "HIE" | "HIP" | "HSD" | "HSE" | "HSP" => "HIS",
        "CYX" | "CYM" => "CYS",
        "ASH" => "ASP",
        "GLH" => "GLU",
        "LYN" => "LYS",
        other => other,
    }
}

#[cfg(test)]
#[path = "protonation_tests.rs"]
mod tests;
