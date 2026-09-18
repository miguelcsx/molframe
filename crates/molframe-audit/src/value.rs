use molframe_core::contract::{
    AlignmentPolicy, AltlocPolicy, AssemblyChoice, ContactDefinition, EquivalencePolicy,
    HydrogenPolicy, MissingPolicy, ModelChoice, Namespace, PeriodicPolicy, Precision, RadiiSet,
    SymmetryPolicy, Tolerance,
};

/// One typed value that can replace a field in an analysis policy.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum PolicyValue {
    /// An assembly choice.
    Assembly(AssemblyChoice),
    /// A model choice.
    Model(ModelChoice),
    /// An alternate-location rule.
    Altloc(AltlocPolicy),
    /// An identifier namespace.
    Identifiers(Namespace),
    /// A missing-atom rule.
    MissingAtoms(MissingPolicy),
    /// A hydrogen rule.
    Hydrogens(HydrogenPolicy),
    /// An atom-equivalence rule.
    AtomEquivalence(EquivalencePolicy),
    /// A symmetry rule.
    Symmetry(SymmetryPolicy),
    /// An alignment rule.
    Alignment(AlignmentPolicy),
    /// An accumulation precision.
    Precision(Precision),
    /// A periodic-boundary rule.
    Periodic(PeriodicPolicy),
    /// A van der Waals radius set.
    VdwRadii(RadiiSet),
    /// A contact definition.
    ContactDef(ContactDefinition),
    /// A floating-point tolerance.
    FloatTolerance(Tolerance),
}
