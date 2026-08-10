//! Side-chain torsions assessed against explicit versioned reference profiles.
//!
//! The profile names component-local torsion paths and empirical distributions;
//! the component provider verifies that each path exists in CCD connectivity.
//! No residue table, atom-name priority, angular well or cutoff is built into
//! this kernel. The caller selects every scientific reference and threshold.

use crate::{ReferenceError, ReferenceLibrary};
use pdbiox_chem::{Component, ComponentProvider};
use pdbiox_core::contract::DictionaryVersion;
use pdbiox_core::index::ResidueIndex;
use pdbiox_core::structure::{ResidueRef, Structure};
use pdbiox_core::{AnalysisPolicy, Code, Diagnostic};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

/// One component-local χ definition and its empirical distribution.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RotamerDefinition {
    /// Exact CCD component identifier.
    pub component_id: Box<str>,
    /// One-based torsion ordinal.
    pub chi_index: u8,
    /// Ordered component atom identifiers defining the dihedral.
    pub atoms: [Box<str>; 4],
    /// Scalar distribution name inside the selected reference set.
    pub distribution: Box<str>,
}

/// Inspectable component-to-reference mapping with stable identity.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RotamerProfile {
    id: Box<str>,
    version: Box<str>,
    definitions: Vec<RotamerDefinition>,
}

impl RotamerProfile {
    /// Builds a deterministic profile and rejects incomplete or duplicate keys.
    ///
    /// # Errors
    ///
    /// Returns [`RotamerError::InvalidProfile`] for invalid identity, torsions,
    /// atom identifiers, distribution names or duplicate component/χ pairs.
    pub fn new(
        id: impl Into<Box<str>>,
        version: impl Into<Box<str>>,
        definitions: impl IntoIterator<Item = RotamerDefinition>,
    ) -> Result<Self, RotamerError> {
        let id = id.into();
        let version = version.into();
        let mut definitions = definitions.into_iter().collect::<Vec<_>>();
        definitions.sort_by(|left, right| {
            left.component_id
                .cmp(&right.component_id)
                .then(left.chi_index.cmp(&right.chi_index))
        });
        let invalid = id.trim().is_empty()
            || version.trim().is_empty()
            || definitions.is_empty()
            || definitions.iter().any(invalid_definition)
            || definitions.windows(2).any(|pair| {
                pair[0].component_id == pair[1].component_id
                    && pair[0].chi_index == pair[1].chi_index
            });
        if invalid {
            Err(RotamerError::InvalidProfile)
        } else {
            Ok(Self {
                id,
                version,
                definitions,
            })
        }
    }

    /// Stable profile family identifier.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Exact profile release.
    #[must_use]
    pub fn version(&self) -> &str {
        &self.version
    }
}

/// Explicit decision policy applied to empirical bin probabilities.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RotamerOptions {
    /// Bins below this normalized probability are outliers.
    pub minimum_probability: f64,
}

/// One torsion below the selected empirical support threshold.
#[derive(Clone, Debug, PartialEq)]
pub struct RotamerFlag {
    /// Residue carrying the torsion.
    pub residue: ResidueIndex,
    /// One-based χ ordinal from the selected profile.
    pub chi_index: u8,
    /// Measured torsion in degrees.
    pub chi_degrees: f64,
    /// Normalized mass of its empirical histogram bin.
    pub probability: f64,
    /// Distribution that produced the probability.
    pub distribution: Box<str>,
}

/// Rotamer findings with all reference identities needed for reproduction.
#[derive(Clone, Debug)]
pub struct RotamerReport {
    /// Sorted actionable torsion findings.
    pub flags: Vec<RotamerFlag>,
    /// Unknown-component and altloc findings.
    pub findings: Vec<Diagnostic>,
    /// Exact CCD provider version.
    pub dictionary_version: DictionaryVersion,
    /// Selected rotamer profile identity.
    pub profile_id: Box<str>,
    /// Selected rotamer profile release.
    pub profile_version: Box<str>,
    /// Selected empirical reference-set identity.
    pub reference_id: Box<str>,
    /// Selected empirical reference-set release.
    pub reference_version: Box<str>,
    /// Probability policy used for classification.
    pub options: RotamerOptions,
    /// Profile-defined residue torsions selected for validation.
    pub intended: usize,
    /// Selected torsions whose coordinates and references were assessable.
    pub assessed: usize,
}

/// Failure to construct or apply an explicit rotamer validation policy.
#[derive(Clone, Debug, thiserror::Error)]
#[non_exhaustive]
pub enum RotamerError {
    /// Profile target counts exceed public coverage counters.
    #[error("rotamer target count exceeds u32 coverage limits")]
    CoverageOverflow,
    /// Profile identity or one of its definitions is invalid.
    #[error("rotamer profile must be named, versioned, non-empty and uniquely keyed")]
    InvalidProfile,
    /// Probability threshold is outside the normalized range.
    #[error("rotamer minimum probability must be finite and in [0, 1]")]
    InvalidOptions,
    /// Profile path does not exist in the selected CCD component graph.
    #[error("rotamer profile path for component {component} chi {chi_index} is not CCD-connected")]
    InvalidPath {
        /// Component carrying the invalid path.
        component: Box<str>,
        /// One-based torsion ordinal.
        chi_index: u8,
    },
    /// Component provider failed.
    #[error("component provider failed: {0:?}")]
    Provider(Diagnostic),
    /// Empirical reference lookup or observation failed.
    #[error(transparent)]
    Reference(#[from] ReferenceError),
}

/// Assesses profile-defined χ torsions against named empirical distributions.
///
/// Runs in `O(residues × definitions-per-component)` with provider lookup
/// caching. No compatibility profile or probability threshold is inferred.
///
/// # Errors
///
/// Returns [`RotamerError`] for invalid policy, provider failure, a profile path
/// inconsistent with CCD, or an unavailable/incompatible reference distribution.
pub fn rotamer_outliers(
    structure: &Structure,
    provider: &dyn ComponentProvider,
    policy: &AnalysisPolicy,
    references: &ReferenceLibrary,
    profile: &RotamerProfile,
    options: RotamerOptions,
) -> Result<RotamerReport, RotamerError> {
    validate_options(options)?;
    let selected = structure.resolve_altlocs(policy);
    let mut findings = selected.warnings;
    let mut flags = Vec::new();
    let mut cache: BTreeMap<Box<str>, Arc<Component>> = BTreeMap::new();
    let mut unresolved = BTreeSet::new();
    let mut intended = 0usize;
    let mut assessed = 0usize;
    for residue in structure.data().residues() {
        let Some(component_id) = residue.name() else {
            continue;
        };
        let definitions: Vec<_> = definitions(profile, component_id).collect();
        if definitions.is_empty() {
            continue;
        }
        intended += definitions.len();
        let component = match cache.get(component_id) {
            Some(component) => Some(component.clone()),
            None => provider.get(component_id).map_err(RotamerError::Provider)?,
        };
        let Some(component) = component else {
            unresolved.insert(Box::<str>::from(component_id));
            continue;
        };
        cache
            .entry(Box::<str>::from(component_id))
            .or_insert_with(|| component.clone());
        for definition in definitions {
            validate_path(&component, definition)?;
            let Some(angle) = torsion(residue, &selected.value, definition) else {
                continue;
            };
            let assessment = references.assess_scalar(&definition.distribution, angle)?;
            assessed += 1;
            if assessment.probability < options.minimum_probability {
                flags.push(RotamerFlag {
                    residue: residue.index(),
                    chi_index: definition.chi_index,
                    chi_degrees: angle,
                    probability: assessment.probability,
                    distribution: definition.distribution.clone(),
                });
            }
        }
    }
    findings.extend(
        unresolved
            .into_iter()
            .map(|component| Diagnostic::new(Code::W3201).with_context("component", component)),
    );
    flags.sort_by_key(|flag| (flag.residue.get(), flag.chi_index));
    Ok(RotamerReport {
        flags,
        findings,
        dictionary_version: provider.version().clone(),
        profile_id: profile.id.clone(),
        profile_version: profile.version.clone(),
        reference_id: references.id().into(),
        reference_version: references.version().into(),
        options,
        intended,
        assessed,
    })
}

fn definitions<'a>(
    profile: &'a RotamerProfile,
    component_id: &'a str,
) -> impl Iterator<Item = &'a RotamerDefinition> {
    profile
        .definitions
        .iter()
        .filter(move |definition| definition.component_id.as_ref() == component_id)
}

fn torsion(
    residue: ResidueRef<'_>,
    selected: &pdbiox_core::AtomSelection,
    definition: &RotamerDefinition,
) -> Option<f64> {
    let positions = definition.atoms.each_ref().map(|name| {
        let atom = residue.atom(name)?;
        selected
            .contains(atom.index().get())
            .then(|| atom.position())
            .flatten()
    });
    let [Some(first), Some(second), Some(third), Some(fourth)] = positions else {
        return None;
    };
    pdbiox_geom::dihedral(first, second, third, fourth).map(pdbiox_geom::degrees)
}

fn validate_path(
    component: &Component,
    definition: &RotamerDefinition,
) -> Result<(), RotamerError> {
    let connected = definition.atoms.windows(2).all(|pair| {
        component.bonds.iter().any(|bond| {
            (bond.atom_a == pair[0] && bond.atom_b == pair[1])
                || (bond.atom_a == pair[1] && bond.atom_b == pair[0])
        })
    });
    let atoms_exist = definition
        .atoms
        .iter()
        .all(|name| component.atom(name).is_some());
    if connected && atoms_exist {
        Ok(())
    } else {
        Err(RotamerError::InvalidPath {
            component: definition.component_id.clone(),
            chi_index: definition.chi_index,
        })
    }
}

fn invalid_definition(definition: &RotamerDefinition) -> bool {
    definition.component_id.trim().is_empty()
        || definition.chi_index == 0
        || definition.distribution.trim().is_empty()
        || definition.atoms.iter().any(|name| name.trim().is_empty())
        || definition.atoms.iter().collect::<BTreeSet<_>>().len() != 4
}

fn validate_options(options: RotamerOptions) -> Result<(), RotamerError> {
    if options.minimum_probability.is_finite() && (0.0..=1.0).contains(&options.minimum_probability)
    {
        Ok(())
    } else {
        Err(RotamerError::InvalidOptions)
    }
}

#[cfg(test)]
#[path = "rotamer_tests.rs"]
mod tests;
