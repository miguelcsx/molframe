//! Parameters and results of the residue-aware lDDT.

use crate::CompareError;
use crate::EmptyLddtPolicy;
use molframe_core::index::ResidueIndex;

/// How chemically interchangeable atoms are named in the model.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SymmetryPolicy {
    /// Score atoms under the names the model carries.
    Disabled,
    /// Rename interchangeable atoms of each residue to whichever naming
    /// preserves the most reference distances.
    ///
    /// Interchangeability comes from the bond graph of the residue's component
    /// with bond orders ignored, because a dictionary may mark only one of two
    /// carboxylate oxygens as double bonded. A swap is admitted only about a
    /// planar centre, so the two methyls of valine or leucine, which are
    /// distinguishable by chirality, are never exchanged.
    Resolve {
        /// Most graph automorphisms enumerated for one residue. A residue that
        /// exceeds it keeps the model's naming and is reported.
        automorphism_limit: usize,
    },
}

/// Definition of a residue-aware lDDT.
#[derive(Clone, Debug, PartialEq)]
pub struct ResidueLddtOptions {
    /// Largest reference distance that is tested, in ångström.
    pub inclusion_radius: f64,
    /// Absolute-error tolerances, ascending; a pair scores the fraction it meets.
    pub tolerances: Box<[f64]>,
    /// Pairs of residues of one chain at most this many residues apart are not
    /// tested. Zero excludes only the pairs inside a single residue.
    pub sequence_separation: u32,
    /// Naming of interchangeable atoms.
    pub symmetry: SymmetryPolicy,
    /// Behaviour when the reference has no testable pair.
    pub empty_policy: EmptyLddtPolicy,
}

impl ResidueLddtOptions {
    /// The conventional profile: 15 Å, 0.5/1/2/4 Å tolerances, symmetry resolved.
    #[must_use]
    pub fn standard() -> Self {
        Self {
            inclusion_radius: 15.0,
            tolerances: Box::new([0.5, 1.0, 2.0, 4.0]),
            sequence_separation: 0,
            symmetry: SymmetryPolicy::Resolve {
                automorphism_limit: 64,
            },
            empty_policy: EmptyLddtPolicy::Perfect,
        }
    }

    pub(super) fn validate(&self) -> Result<(), CompareError> {
        let tolerances = !self.tolerances.is_empty()
            && self
                .tolerances
                .iter()
                .all(|value| value.is_finite() && *value > 0.0)
            && self.tolerances.windows(2).all(|pair| pair[0] < pair[1]);
        if !self.inclusion_radius.is_finite() || self.inclusion_radius <= 0.0 || !tolerances {
            return Err(CompareError::InvalidScoreInput);
        }
        Ok(())
    }
}

/// The tests that involve one reference residue.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ResidueScore {
    /// The reference residue.
    pub reference: ResidueIndex,
    /// The model residue it was compared with; `None` when the model has none.
    pub model: Option<ResidueIndex>,
    /// Reference distances from this residue's atoms that were tested.
    pub considered: u64,
    /// Tolerance tests passed, summed over every tolerance.
    pub preserved: u64,
    /// Fraction of tolerance tests passed; `None` when nothing was tested.
    pub score: Option<f64>,
}

/// A residue-aware lDDT.
#[derive(Clone, Debug, PartialEq)]
pub struct ResidueLddt {
    /// Fraction of tolerance tests passed over every tested distance.
    pub score: f64,
    /// Reference distances tested, each pair counted once.
    pub considered: u64,
    /// One entry per reference residue, in the order supplied.
    pub residues: Vec<ResidueScore>,
    /// Residues whose interchangeable atoms were renamed.
    pub renamed_residues: usize,
    /// Residues whose symmetry exceeded the automorphism limit and was skipped.
    pub unresolved_symmetry: usize,
}
