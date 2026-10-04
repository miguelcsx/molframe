//! The system a governed analysis runs over.
//!
//! `AnalysisPolicy::assembly` chooses the arrangement of the contents: the
//! deposited asymmetric unit, a biological assembly the file defines, or the unit
//! with the crystal contacts around it. This module turns that choice into a
//! structure and, for every atom of it, where the atom came from: which atom of
//! the input and under which rigid motion. That map is what lets the same
//! per-frame coordinates feed copies, so a trajectory of the unit drives its
//! assembly.

use molframe_core::contract::{
    AnalysisPolicy, AssemblyChoice, Assumption, AssumptionSource, Impact, PolicyField,
};
use molframe_core::{Code, Diagnostic, ExecutionContext, Structure};
use molframe_geom::Rigid;
use molframe_xtal::{AssemblyExt, AssemblyView, crystal_contact_view};

/// Where one atom of the analysed system came from.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Placement {
    /// The atom of the input structure this atom is a copy of.
    pub(crate) source: usize,
    /// The motion applied to the input coordinates; `None` for the atom itself.
    pub(crate) motion: Option<Rigid>,
}

/// What was built, for the record.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum SystemKind {
    /// The input as deposited.
    AsymmetricUnit,
    /// A biological assembly, by identifier, with its chain count.
    Biological { id: Box<str>, chains: usize },
    /// The unit and its crystal contacts, with the radius and the chain count.
    Crystal { radius: f32, chains: usize },
}

/// The structure to analyse and the origin of each of its atoms.
#[derive(Debug)]
pub(crate) struct AnalysisSystem {
    pub(crate) structure: Structure,
    pub(crate) placements: Vec<Placement>,
    pub(crate) kind: SystemKind,
}

impl AnalysisSystem {
    /// Builds the system the policy asks for.
    ///
    /// # Errors
    ///
    /// Returns the findings that explain why it cannot be built: a contradictory
    /// policy, a missing assembly, absent symmetry operators or cell, an improper
    /// operation, or a failed expansion.
    pub(crate) fn build(
        template: &Structure,
        policy: &AnalysisPolicy,
        context: &ExecutionContext,
    ) -> Result<Self, Vec<Diagnostic>> {
        policy
            .check_consistency()
            .map_err(|finding| vec![finding])?;
        match &policy.assembly {
            AssemblyChoice::AsymmetricUnit => Ok(Self::unit(template)),
            AssemblyChoice::Biological(id) => {
                let view = template.assembly(id).map_err(|finding| vec![finding])?;
                Self::from_view(
                    &view,
                    SystemKind::Biological {
                        id: id.clone(),
                        chains: view.instance_count(),
                    },
                )
            }
            AssemblyChoice::Crystal { radius } => {
                let view = crystal_contact_view(template, *radius, context)
                    .map_err(|finding| vec![finding])?;
                Self::from_view(
                    &view,
                    SystemKind::Crystal {
                        radius: *radius,
                        chains: view.instance_count(),
                    },
                )
            }
            _ => Err(vec![
                Diagnostic::new(Code::E6103).with_context("field", "assembly"),
            ]),
        }
    }

    fn unit(template: &Structure) -> Self {
        Self {
            structure: template.clone(),
            placements: (0..template.atom_count() as usize)
                .map(|source| Placement {
                    source,
                    motion: None,
                })
                .collect(),
            kind: SystemKind::AsymmetricUnit,
        }
    }

    fn from_view(view: &AssemblyView, kind: SystemKind) -> Result<Self, Vec<Diagnostic>> {
        let structure = view.materialize()?;
        let placements: Vec<Placement> = view
            .atoms()
            .map(|atom| Placement {
                source: atom.source_atom.as_usize(),
                motion: Some(atom.transform),
            })
            .collect();
        if placements.len() != structure.atom_count() as usize {
            return Err(vec![Diagnostic::new(Code::E9001).with_context(
                "reason",
                "the materialised atoms do not match the generated placements",
            )]);
        }
        Ok(Self {
            structure,
            placements,
            kind,
        })
    }

    /// Whether atoms of the system are copies of one another.
    pub(crate) const fn is_replicated(&self) -> bool {
        !matches!(self.kind, SystemKind::AsymmetricUnit)
    }

    /// The decision this system records, when it is not the deposited unit.
    pub(crate) fn assumption(&self) -> Option<Assumption> {
        let value = match &self.kind {
            SystemKind::AsymmetricUnit => return None,
            SystemKind::Biological { id, chains } => {
                format!("biological assembly {id}: {chains} chain instances")
            }
            SystemKind::Crystal { radius, chains } => {
                format!("crystal contacts within {radius} A: {chains} chain instances")
            }
        };
        Some(Assumption::new(
            PolicyField::Assembly,
            value,
            AssumptionSource::Explicit,
            Impact::Unmeasured,
        ))
    }
}

#[cfg(test)]
#[path = "system_tests.rs"]
mod tests;
