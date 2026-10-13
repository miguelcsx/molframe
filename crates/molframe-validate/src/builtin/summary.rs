//! Per-structure counts over the per-residue built-in tables.
//!
//! The tables hold one code per residue; the summaries fold them into the
//! counts a report needs, so callers never re-implement the class meanings.

use super::nucleic_scoring::{NucleicTable, PUCKER_C2_ENDO, PUCKER_C3_ENDO, PUCKER_OTHER};
use super::protein_conformation::{NOT_EVALUATED, ProteinTable};

/// Counts of Ramachandran, rotamer and peptide-bond classes.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ProteinSummary {
    /// Residues with a Ramachandran class.
    pub rama_evaluated: usize,
    /// Residues in the favoured region.
    pub rama_favoured: usize,
    /// Residues in the allowed region.
    pub rama_allowed: usize,
    /// Residues in the outlier region.
    pub rama_outliers: usize,
    /// Residues with a rotamer class.
    pub rotamer_evaluated: usize,
    /// Side chains in a favoured rotamer.
    pub rotamer_favoured: usize,
    /// Side chains in an allowed rotamer.
    pub rotamer_allowed: usize,
    /// Side chains in an outlier rotamer.
    pub rotamer_outliers: usize,
    /// Peptide bonds with an omega class.
    pub omega_evaluated: usize,
    /// Cis peptide bonds.
    pub cis_peptides: usize,
    /// Peptide bonds twisted out of plane.
    pub twisted_peptides: usize,
}

/// Counts of sugar pucker and suite-sanity classes.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct NucleicSummary {
    /// Residues with a pucker class.
    pub pucker_evaluated: usize,
    /// Sugars in C3'-endo.
    pub c3_endo: usize,
    /// Sugars in C2'-endo.
    pub c2_endo: usize,
    /// Sugars in neither pucker.
    pub other_pucker: usize,
    /// Residues with a suite sanity verdict.
    pub suite_evaluated: usize,
    /// Residues whose suite was flagged.
    pub suite_flagged: usize,
}

/// Percentage of `part` in `whole`, `NaN` when `whole` is zero.
#[must_use]
pub fn percent(part: usize, whole: usize) -> f64 {
    if whole == 0 {
        return f64::NAN;
    }
    #[allow(clippy::cast_precision_loss)]
    let ratio = part as f64 / whole as f64;
    100.0 * ratio
}

impl ProteinTable {
    /// Fold the per-residue classes into counts.
    #[must_use]
    pub fn summary(&self) -> ProteinSummary {
        let mut out = ProteinSummary::default();
        for &code in self.rama_class.iter().filter(|&&c| c != NOT_EVALUATED) {
            out.rama_evaluated += 1;
            match code {
                0 => out.rama_outliers += 1,
                1 => out.rama_allowed += 1,
                _ => out.rama_favoured += 1,
            }
        }
        for &code in self.rotamer_class.iter().filter(|&&c| c != NOT_EVALUATED) {
            out.rotamer_evaluated += 1;
            match code {
                0 => out.rotamer_outliers += 1,
                1 => out.rotamer_allowed += 1,
                _ => out.rotamer_favoured += 1,
            }
        }
        for &code in self.omega_class.iter().filter(|&&c| c != NOT_EVALUATED) {
            out.omega_evaluated += 1;
            match code {
                1 => out.cis_peptides += 1,
                2 => out.twisted_peptides += 1,
                _ => {}
            }
        }
        out
    }
}

impl NucleicTable {
    /// Fold the per-residue classes into counts.
    #[must_use]
    pub fn summary(&self) -> NucleicSummary {
        let mut out = NucleicSummary::default();
        for &code in self.pucker_class.iter().filter(|&&c| c != NOT_EVALUATED) {
            out.pucker_evaluated += 1;
            match code {
                PUCKER_C3_ENDO => out.c3_endo += 1,
                PUCKER_C2_ENDO => out.c2_endo += 1,
                PUCKER_OTHER => out.other_pucker += 1,
                _ => {}
            }
        }
        for &code in self.suite_class.iter().filter(|&&c| c != NOT_EVALUATED) {
            out.suite_evaluated += 1;
            if code == 1 {
                out.suite_flagged += 1;
            }
        }
        out
    }
}

#[cfg(test)]
#[path = "summary_tests.rs"]
mod tests;
