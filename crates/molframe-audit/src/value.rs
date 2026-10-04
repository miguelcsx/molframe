use molframe_core::contract::{
    AlignmentPolicy, AltlocPolicy, AssemblyChoice, ContactDefinition, EquivalencePolicy,
    HydrogenPolicy, MissingPolicy, ModelChoice, Namespace, PeriodicPolicy, PolicyField,
    PolicyParseError, Precision, RadiiSet, SymmetryPolicy, Tolerance,
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

impl PolicyValue {
    /// Reads the value of `field` written as a word of a written policy.
    ///
    /// The same parser as the policy's own is used, so `"first"`, `"biological:1"` and
    /// `"crystal:8.0"` mean here what they mean there. `float_tolerance` is written
    /// `"relative,absolute"`.
    ///
    /// # Errors
    ///
    /// Returns the policy parse error naming the field, the word and the accepted set.
    pub fn named(field: PolicyField, word: &str) -> Result<Self, PolicyParseError> {
        fn parse<T: std::str::FromStr<Err = PolicyParseError>>(
            word: &str,
            wrap: fn(T) -> PolicyValue,
        ) -> Result<PolicyValue, PolicyParseError> {
            word.parse().map(wrap)
        }
        match field {
            PolicyField::Assembly => parse(word, Self::Assembly),
            PolicyField::Model => parse(word, Self::Model),
            PolicyField::Altloc => parse(word, Self::Altloc),
            PolicyField::Identifiers => parse(word, Self::Identifiers),
            PolicyField::MissingAtoms => parse(word, Self::MissingAtoms),
            PolicyField::Hydrogens => parse(word, Self::Hydrogens),
            PolicyField::AtomEquivalence => parse(word, Self::AtomEquivalence),
            PolicyField::Symmetry => parse(word, Self::Symmetry),
            PolicyField::Alignment => parse(word, Self::Alignment),
            PolicyField::Precision => parse(word, Self::Precision),
            PolicyField::Periodic => parse(word, Self::Periodic),
            PolicyField::VdwRadii => parse(word, Self::VdwRadii),
            PolicyField::ContactDef => parse(word, Self::ContactDef),
            PolicyField::FloatTolerance => tolerance(word).map(Self::FloatTolerance),
            _ => Err(PolicyParseError::new(
                "policy field",
                field.name(),
                "a field this version of the audit knows",
            )),
        }
    }
}

fn tolerance(word: &str) -> Result<Tolerance, PolicyParseError> {
    let refuse = || PolicyParseError::new("float_tolerance", word, "<relative>,<absolute>");
    let (relative, absolute) = word.split_once(',').ok_or_else(refuse)?;
    let parse = |text: &str| {
        text.trim()
            .parse::<f64>()
            .ok()
            .filter(|v| v.is_finite() && *v >= 0.0)
    };
    match (parse(relative), parse(absolute)) {
        (Some(relative), Some(absolute)) => Ok(Tolerance { relative, absolute }),
        _ => Err(refuse()),
    }
}
