use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

/// A component-local atom addressed as `component.atom`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AtomSite {
    /// Motif component name.
    pub component: Box<str>,
    /// Chemical atom name.
    pub atom: Box<str>,
}

impl AtomSite {
    /// Creates an atom reference without string parsing.
    #[must_use]
    pub fn new(component: impl Into<Box<str>>, atom: impl Into<Box<str>>) -> Self {
        Self {
            component: component.into(),
            atom: atom.into(),
        }
    }
}

/// Chemical role a component must fulfil.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum ComponentRole {
    /// A polymer residue.
    #[default]
    Residue,
    /// A discrete ligand.
    Ligand,
    /// A cofactor distinguished from an ordinary ligand by the specification.
    Cofactor,
    /// A monoatomic metal centre.
    Metal,
}

/// Declarative chemical requirements for one motif component.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ComponentSpec {
    /// Accepted component identifiers; empty accepts any identifier of the role.
    pub component_ids: BTreeSet<Box<str>>,
    /// Chemical role.
    pub role: ComponentRole,
    /// Atoms that must be present.
    pub required_atoms: BTreeSet<Box<str>>,
    /// Explicit interchangeable atom-name groups.
    pub equivalent_atoms: Vec<BTreeSet<Box<str>>>,
}

impl ComponentSpec {
    /// Creates a component specification of one role.
    #[must_use]
    pub const fn new(role: ComponentRole) -> Self {
        Self {
            component_ids: BTreeSet::new(),
            role,
            required_atoms: BTreeSet::new(),
            equivalent_atoms: Vec::new(),
        }
    }

    /// Adds an accepted chemical component identifier.
    #[must_use]
    pub fn component(mut self, id: impl Into<Box<str>>) -> Self {
        self.component_ids.insert(id.into());
        self
    }

    /// Adds a required atom.
    #[must_use]
    pub fn require(mut self, atom: impl Into<Box<str>>) -> Self {
        self.required_atoms.insert(atom.into());
        self
    }

    /// Declares atom names interchangeable for this motif.
    #[must_use]
    pub fn equivalent(mut self, atoms: impl IntoIterator<Item = impl Into<Box<str>>>) -> Self {
        self.equivalent_atoms
            .push(atoms.into_iter().map(Into::into).collect());
        self
    }
}

/// One geometric or chemical constraint.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum Constraint {
    /// Distance in ångström.
    Distance {
        /// First atom.
        first: AtomSite,
        /// Second atom.
        second: AtomSite,
        /// Desired distance.
        target: f64,
        /// Allowed absolute deviation.
        tolerance: f64,
    },
    /// Angle in degrees.
    Angle {
        /// Three atoms, with the vertex in the middle.
        atoms: [AtomSite; 3],
        /// Desired angle.
        target: f64,
        /// Allowed absolute deviation.
        tolerance: f64,
    },
    /// Dihedral in degrees.
    Dihedral {
        /// Four ordered atoms.
        atoms: [AtomSite; 4],
        /// Desired signed angle.
        target: f64,
        /// Allowed circular deviation.
        tolerance: f64,
    },
    /// Required tetrahedral handedness.
    Chirality {
        /// Centre followed by three ordered neighbours.
        atoms: [AtomSite; 4],
        /// Positive or negative scalar-triple sign.
        positive: bool,
    },
    /// Maximum deviation from a common plane.
    Planarity {
        /// At least three atoms.
        atoms: Vec<AtomSite>,
        /// Maximum RMS deviation in ångström.
        tolerance: f64,
    },
    /// Number of partners within a distance of a centre.
    Coordination {
        /// Coordination centre.
        centre: AtomSite,
        /// Potential coordinating atoms.
        partners: Vec<AtomSite>,
        /// Required partner count.
        count: usize,
        /// Maximum centre-partner distance.
        max_distance: f64,
    },
    /// Minimum allowed separation.
    StericExclusion {
        /// First atom.
        first: AtomSite,
        /// Second atom.
        second: AtomSite,
        /// Minimum distance.
        min_distance: f64,
    },
}

impl Constraint {
    pub(crate) fn sites(&self) -> Vec<&AtomSite> {
        match self {
            Self::Distance { first, second, .. } | Self::StericExclusion { first, second, .. } => {
                vec![first, second]
            }
            Self::Angle { atoms, .. } => atoms.iter().collect(),
            Self::Dihedral { atoms, .. } | Self::Chirality { atoms, .. } => atoms.iter().collect(),
            Self::Planarity { atoms, .. } => atoms.iter().collect(),
            Self::Coordination {
                centre, partners, ..
            } => std::iter::once(centre).chain(partners).collect(),
        }
    }
}

/// Stable name paired with a constraint.
#[derive(Clone, Debug, PartialEq)]
pub struct NamedConstraint {
    /// Name used in decomposed output and verdict profiles.
    pub name: Box<str>,
    /// Constraint definition.
    pub constraint: Constraint,
}

/// Immutable declarative functional motif.
#[derive(Clone, Debug, PartialEq)]
pub struct Motif {
    components: BTreeMap<Box<str>, ComponentSpec>,
    constraints: Vec<NamedConstraint>,
}

impl Motif {
    /// Validates component names, constraint names and every atom reference.
    ///
    /// # Errors
    ///
    /// Returns [`MotifError`] for duplicates, unknown components or malformed
    /// constraints.
    pub fn new(
        components: impl IntoIterator<Item = (Box<str>, ComponentSpec)>,
        constraints: impl IntoIterator<Item = NamedConstraint>,
    ) -> Result<Self, MotifError> {
        let mut component_map = BTreeMap::new();
        for (name, component) in components {
            if component_map.insert(name.clone(), component).is_some() {
                return Err(MotifError::DuplicateComponent(name));
            }
        }
        let constraints: Vec<_> = constraints.into_iter().collect();
        let mut names = BTreeSet::new();
        for named in &constraints {
            if !names.insert(named.name.clone()) {
                return Err(MotifError::DuplicateConstraint(named.name.clone()));
            }
            validate_constraint(named, &component_map)?;
        }
        Ok(Self {
            components: component_map,
            constraints,
        })
    }

    /// Components in stable lexical order.
    #[must_use]
    pub const fn components(&self) -> &BTreeMap<Box<str>, ComponentSpec> {
        &self.components
    }

    /// Constraints in declaration order.
    #[must_use]
    pub fn constraints(&self) -> &[NamedConstraint] {
        &self.constraints
    }
}

fn validate_constraint(
    named: &NamedConstraint,
    components: &BTreeMap<Box<str>, ComponentSpec>,
) -> Result<(), MotifError> {
    for site in named.constraint.sites() {
        if !components.contains_key(site.component.as_ref()) {
            return Err(MotifError::UnknownComponent(site.component.clone()));
        }
    }
    let valid = match &named.constraint {
        Constraint::Distance {
            target, tolerance, ..
        }
        | Constraint::Angle {
            target, tolerance, ..
        }
        | Constraint::Dihedral {
            target, tolerance, ..
        } => target.is_finite() && tolerance.is_finite() && *tolerance >= 0.0,
        Constraint::Chirality { .. } => true,
        Constraint::Planarity { atoms, tolerance } => {
            atoms.len() >= 3 && tolerance.is_finite() && *tolerance >= 0.0
        }
        Constraint::Coordination {
            partners,
            count,
            max_distance,
            ..
        } => *count <= partners.len() && max_distance.is_finite() && *max_distance >= 0.0,
        Constraint::StericExclusion { min_distance, .. } => {
            min_distance.is_finite() && *min_distance >= 0.0
        }
    };
    if valid {
        Ok(())
    } else {
        Err(MotifError::InvalidConstraint(named.name.clone()))
    }
}

/// Invalid declarative motif.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum MotifError {
    /// A component name occurs twice.
    DuplicateComponent(Box<str>),
    /// A constraint name occurs twice.
    DuplicateConstraint(Box<str>),
    /// An atom refers to a component absent from the motif.
    UnknownComponent(Box<str>),
    /// Numeric or cardinality fields are invalid.
    InvalidConstraint(Box<str>),
}

impl fmt::Display for MotifError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DuplicateComponent(name) => write!(formatter, "duplicate component {name}"),
            Self::DuplicateConstraint(name) => write!(formatter, "duplicate constraint {name}"),
            Self::UnknownComponent(name) => write!(formatter, "unknown component {name}"),
            Self::InvalidConstraint(name) => write!(formatter, "invalid constraint {name}"),
        }
    }
}

impl std::error::Error for MotifError {}

#[cfg(test)]
#[path = "specification_tests.rs"]
mod tests;
