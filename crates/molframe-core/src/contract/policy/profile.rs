//! The policy itself, and the named profile it defaults to.
//!
//! "We used the defaults" is not a statement anyone can reproduce. "We used
//! molframe-default-1.0" is, which is why the default is a named, versioned
//! profile recorded exactly as an explicit policy would be.

use super::fields::{
    AlignmentPolicy, AltlocPolicy, AssemblyChoice, ContactDefinition, EquivalencePolicy,
    HydrogenPolicy, MissingPolicy, ModelChoice, Namespace, PeriodicPolicy, Precision, ProfileId,
    RadiiSet, SymmetryPolicy, Tolerance,
};
use super::fingerprint::{Fingerprint, Fnv1a};
use std::fmt;
use std::hash::{Hash, Hasher};

/// Every semantic decision an analysis depends on.
///
/// # Examples
///
/// ```
/// use molframe_core::contract::{AnalysisPolicy, Namespace, ProfileId};
///
/// let policy = AnalysisPolicy::default();
/// assert_eq!(policy.profile(), Some(ProfileId::DEFAULT));
/// assert_eq!(policy.identifiers, Namespace::Auth);
///
/// let explicit = AnalysisPolicy::default().with_identifiers(Namespace::Explicit);
/// assert_eq!(explicit.profile(), None, "no longer the named profile");
/// assert_ne!(explicit.fingerprint(), policy.fingerprint());
/// ```
#[derive(Clone, PartialEq, Debug, Default)]
pub struct AnalysisPolicy {
    /// Which arrangement of the contents to analyse.
    pub assembly: AssemblyChoice,
    /// Which model.
    pub model: ModelChoice,
    /// How alternate conformations resolve.
    pub altloc: AltlocPolicy,
    /// Which identifier namespace an unqualified selector means.
    pub identifiers: Namespace,
    /// What to do about atoms that were not modelled.
    pub missing_atoms: MissingPolicy,
    /// How hydrogens enter.
    pub hydrogens: HydrogenPolicy,
    /// How atoms are matched between structures.
    pub atom_equivalence: EquivalencePolicy,
    /// Whether symmetry copies participate.
    pub symmetry: SymmetryPolicy,
    /// How structures are put into correspondence.
    pub alignment: AlignmentPolicy,
    /// The width accumulation runs at.
    pub precision: Precision,
    /// How periodic boundaries are handled.
    pub periodic: PeriodicPolicy,
    /// Which radii set.
    pub vdw_radii: RadiiSet,
    /// What counts as a contact.
    pub contact_def: ContactDefinition,
    /// How close two floating-point values must be to count as equal.
    pub float_tolerance: Tolerance,
}

impl AnalysisPolicy {
    /// The named profile this policy is, if it has not been altered.
    ///
    /// A result recording a profile name is reproducible by someone who has only
    /// the name. A result recording "the defaults" is not.
    #[must_use]
    pub fn profile(&self) -> Option<ProfileId> {
        (*self == Self::default()).then_some(ProfileId::DEFAULT)
    }

    /// A stable fingerprint of every field.
    ///
    /// Two policies that differ anywhere differ here, and the same policy always
    /// fingerprints the same way, within and across runs. It is a fingerprint
    /// for telling policies apart, not a cryptographic commitment.
    #[must_use]
    pub fn fingerprint(&self) -> Fingerprint {
        let mut hasher = Fnv1a::new();
        self.assembly.hash(&mut hasher);
        self.model.hash(&mut hasher);
        self.altloc.hash(&mut hasher);
        self.identifiers.hash(&mut hasher);
        self.missing_atoms.hash(&mut hasher);
        self.hydrogens.hash(&mut hasher);
        self.atom_equivalence.hash(&mut hasher);
        self.symmetry.hash(&mut hasher);
        self.alignment.hash(&mut hasher);
        self.precision.hash(&mut hasher);
        self.periodic.hash(&mut hasher);
        self.vdw_radii.hash(&mut hasher);
        match self.contact_def {
            ContactDefinition::DistanceCutoff { tolerance } => {
                (0u8, tolerance.to_bits()).hash(&mut hasher);
            }
            ContactDefinition::SurfaceBased { probe } => {
                (1u8, probe.to_bits()).hash(&mut hasher);
            }
        }
        (
            self.float_tolerance.relative.to_bits(),
            self.float_tolerance.absolute.to_bits(),
        )
            .hash(&mut hasher);
        Fingerprint::from_hash(hasher.finish())
    }

    /// Returns this policy with a different identifier namespace.
    #[must_use]
    pub fn with_identifiers(mut self, namespace: Namespace) -> Self {
        self.identifiers = namespace;
        self
    }

    /// Returns this policy with a different alternate-conformation rule.
    #[must_use]
    pub fn with_altloc(mut self, altloc: AltlocPolicy) -> Self {
        self.altloc = altloc;
        self
    }

    /// Returns this policy with a different model choice.
    #[must_use]
    pub fn with_model(mut self, model: ModelChoice) -> Self {
        self.model = model;
        self
    }

    /// Returns this policy with a different assembly choice.
    #[must_use]
    pub fn with_assembly(mut self, assembly: AssemblyChoice) -> Self {
        self.assembly = assembly;
        self
    }

    /// Returns this policy with a different rule for unmodelled atoms.
    #[must_use]
    pub fn with_missing_atoms(mut self, missing: MissingPolicy) -> Self {
        self.missing_atoms = missing;
        self
    }
}

impl fmt::Display for AnalysisPolicy {
    /// Writes the policy as one field per line, which is what the command line
    /// prints and what a written record contains.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.profile() {
            Some(profile) => writeln!(f, "{profile}")?,
            None => writeln!(f, "(modified from {})", ProfileId::DEFAULT)?,
        }
        writeln!(f, "  assembly         = {:?}", self.assembly)?;
        writeln!(f, "  model            = {:?}", self.model)?;
        writeln!(f, "  altloc           = {:?}", self.altloc)?;
        writeln!(f, "  identifiers      = {:?}", self.identifiers)?;
        writeln!(f, "  missing_atoms    = {:?}", self.missing_atoms)?;
        writeln!(f, "  hydrogens        = {:?}", self.hydrogens)?;
        writeln!(f, "  atom_equivalence = {:?}", self.atom_equivalence)?;
        writeln!(f, "  symmetry         = {:?}", self.symmetry)?;
        writeln!(f, "  alignment        = {:?}", self.alignment)?;
        writeln!(f, "  precision        = {:?}", self.precision)?;
        writeln!(f, "  periodic         = {:?}", self.periodic)?;
        writeln!(f, "  vdw_radii        = {:?}", self.vdw_radii)?;
        writeln!(f, "  contact_def      = {:?}", self.contact_def)?;
        write!(f, "  float_tolerance  = {:?}", self.float_tolerance)
    }
}

#[cfg(test)]
#[path = "profile_tests.rs"]
mod tests;
