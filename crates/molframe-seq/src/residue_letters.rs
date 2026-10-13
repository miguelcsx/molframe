//! Residue and component names as the one-letter codes a sequence is written in.

/// The one-letter code of a component name: the twenty amino acids, the common variants a
/// deposited sequence uses, and the nucleotides. `X` for anything else (an unknown or modified
/// residue). Case and surrounding whitespace do not matter.
#[must_use]
pub fn one_letter(name: &str) -> char {
    match name.trim().to_ascii_uppercase().as_str() {
        "ALA" | "A" | "DA" => 'A',
        "ARG" => 'R',
        "ASN" => 'N',
        "ASP" | "ASH" => 'D',
        "CYS" | "CYX" | "CYM" | "CSS" | "C" | "DC" => 'C',
        "GLN" => 'Q',
        "GLU" | "GLH" => 'E',
        "GLY" | "G" | "DG" => 'G',
        "HIS" | "HID" | "HIE" | "HIP" | "HSD" | "HSE" | "HSP" => 'H',
        "ILE" => 'I',
        "LEU" => 'L',
        "LYS" | "LYN" => 'K',
        "MET" | "MSE" => 'M',
        "PHE" => 'F',
        "PRO" => 'P',
        "SER" => 'S',
        "THR" | "T" | "DT" => 'T',
        "TRP" => 'W',
        "TYR" => 'Y',
        "VAL" => 'V',
        "SEC" | "U" | "DU" => 'U',
        "PYL" => 'O',
        _ => 'X',
    }
}

#[cfg(test)]
#[path = "residue_letters_tests.rs"]
mod tests;
