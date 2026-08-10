use crate::MappedMotif;
use pdbiox_geom::Rigid;

/// How a mapped motif was placed before measurement.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AlignmentKind {
    /// All constraints are internal and therefore rigid-transform invariant.
    NotRequired,
    /// A caller supplied the rigid transform for a shared reference frame.
    CallerSupplied,
}

/// A mapping paired with an explicit coordinate-frame decision.
#[derive(Clone, Debug, PartialEq)]
pub struct AlignedMotif {
    /// Chemistry-aware component and atom mapping.
    pub mapping: MappedMotif,
    /// Transform applied to positions before measurement.
    pub transform: Rigid,
    /// Why that transform was selected.
    pub kind: AlignmentKind,
}

/// Records that no fit is required for rigid-transform-invariant constraints.
#[must_use]
pub fn align_intrinsic(mapping: MappedMotif) -> AlignedMotif {
    AlignedMotif {
        mapping,
        transform: Rigid::IDENTITY,
        kind: AlignmentKind::NotRequired,
    }
}

/// Applies a caller-computed rigid transform without guessing correspondences.
#[must_use]
pub fn align_with_transform(mapping: MappedMotif, transform: Rigid) -> AlignedMotif {
    AlignedMotif {
        mapping,
        transform,
        kind: AlignmentKind::CallerSupplied,
    }
}

#[cfg(test)]
#[path = "alignment_tests.rs"]
mod tests;
