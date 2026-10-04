//! Contacts whose definition and radii come from the policy.
//!
//! A distance cutoff and a shared-surface test are different answers to "are these
//! atoms in contact", and which radii the distance is measured against changes the
//! answer again. `contacts_kernel` takes one cutoff and ignores both decisions; this
//! kernel reads them, so varying them in an audit varies what the analysis does.

use super::common::{backend, descriptor, float};
use super::{FrameKernelResult, StructureKernel, structure_kernel};
use crate::{Contact, ContactTable, SurfaceContactOptions, atom_contacts, surface_contacts};
use molframe_chem::{RadiusSet, atom_radii};
use molframe_core::contract::{AnalysisPolicy, ContactDefinition, Coverage, PolicyField, Quality};
use molframe_core::{ExecutionContext, Structure};
use molframe_spatial::{SpatialBackend, SpatialError};
use molframe_surface::SasaError;

/// The parameters of the surface-based definition that the policy does not carry.
///
/// A policy says a contact is `surface:<probe>`; how finely the surface is sampled
/// and how large an exposed patch must be are properties of the algorithm.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceSampling {
    /// Extra separation accepted beyond the two radii, in ångström.
    pub tolerance: f32,
    /// Surface samples per square ångström.
    pub density: f32,
    /// Minimum exposed patch on both atoms, in square ångström.
    pub minimum_area: f32,
}

/// Why a policy-defined contact analysis could not be evaluated.
#[derive(Debug, thiserror::Error)]
pub enum DefinitionError {
    /// The pair search refused its inputs.
    #[error(transparent)]
    Spatial(#[from] SpatialError),
    /// The surface test refused its inputs.
    #[error(transparent)]
    Surface(#[from] SasaError),
    /// The policy names a radius set this version does not have.
    #[error("the policy's radius set {0:?} is not one this kernel knows")]
    UnknownRadiusSet(Box<str>),
    /// A surface contact cannot be decided for an atom with no radius.
    #[error(
        "{atoms} atoms have an element with no radius in the {set} set, and a surface contact \
         cannot be decided without one"
    )]
    UnknownRadius {
        /// How many atoms lack a radius.
        atoms: usize,
        /// The radius set consulted.
        set: &'static str,
    },
}

/// Governed contacts under the policy's `contact_def` and `vdw_radii`.
///
/// With `distance:<tolerance>`, two atoms are in contact when their distance is at
/// most the sum of their radii plus the tolerance. An atom whose element has no
/// radius in the set cannot be decided and is reported as missing coverage, so the
/// missing-atoms policy governs what happens. With `surface:<probe>`, two atoms are
/// in contact when their expanded surfaces touch and both have an exposed patch.
#[must_use]
pub fn definition_contacts_kernel(
    spatial: SpatialBackend,
    surface: SurfaceSampling,
) -> impl StructureKernel<Output = ContactTable, Error = DefinitionError> {
    structure_kernel(
        descriptor("contacts-by-definition")
            .reading(&[PolicyField::ContactDef, PolicyField::VdwRadii])
            .with_parameter("spatial_backend", backend(spatial))
            .with_parameter("surface_tolerance", float(surface.tolerance))
            .with_parameter("surface_density", float(surface.density))
            .with_parameter("surface_minimum_area", float(surface.minimum_area)),
        move |structure: &Structure, policy: &AnalysisPolicy, context: &ExecutionContext| {
            evaluate(structure, policy, spatial, surface, context)
        },
    )
}

fn evaluate(
    structure: &Structure,
    policy: &AnalysisPolicy,
    spatial: SpatialBackend,
    surface: SurfaceSampling,
    context: &ExecutionContext,
) -> Result<FrameKernelResult<ContactTable>, DefinitionError> {
    let set: RadiusSet = policy
        .vdw_radii
        .name()
        .parse()
        .map_err(|_| DefinitionError::UnknownRadiusSet(policy.vdw_radii.name().into()))?;
    let radii = atom_radii(structure, set);
    let unknown = radii.iter().filter(|radius| !radius.is_finite()).count();
    let atoms = structure.atom_count();
    let table = match policy.contact_def {
        ContactDefinition::DistanceCutoff { tolerance } => {
            by_distance(structure, &radii, tolerance, spatial, context)?
        }
        ContactDefinition::SurfaceBased { probe } => {
            if unknown > 0 {
                return Err(DefinitionError::UnknownRadius {
                    atoms: unknown,
                    set: set.name(),
                });
            }
            let options = SurfaceContactOptions {
                tolerance: surface.tolerance,
                probe,
                surface_density: surface.density,
                minimum_area: surface.minimum_area,
                backend: spatial,
            };
            surface_contacts(structure, &radii, options, context)?
                .into_iter()
                .collect::<ContactTable>()
        }
        // A definition this version does not know cannot be silently treated as a distance.
        _ => return Err(DefinitionError::UnknownRadiusSet("contact_def".into())),
    };
    // More atoms than `u32` is not representable as a structure, so the count always fits.
    let missing = match u32::try_from(unknown) {
        Ok(count) => count.min(atoms),
        Err(_) => atoms,
    };
    let coverage = Coverage {
        intended: atoms,
        used: atoms - missing,
        missing,
        ambiguous: 0,
    };
    Ok(FrameKernelResult::governed(
        table,
        Quality::Complete,
        coverage,
    ))
}

fn by_distance(
    structure: &Structure,
    radii: &[f32],
    tolerance: f32,
    spatial: SpatialBackend,
    context: &ExecutionContext,
) -> Result<ContactTable, DefinitionError> {
    let widest = radii
        .iter()
        .copied()
        .filter(|radius| radius.is_finite())
        .fold(0.0_f32, f32::max);
    let candidates = atom_contacts(structure, 2.0 * widest + tolerance, spatial, context)?;
    Ok(candidates
        .iter()
        .filter(
            |Contact {
                 first,
                 second,
                 distance,
             }| {
                match (radii.get(first.as_usize()), radii.get(second.as_usize())) {
                    (Some(a), Some(b)) if a.is_finite() && b.is_finite() => {
                        *distance <= a + b + tolerance
                    }
                    _ => false,
                }
            },
        )
        .collect())
}

#[cfg(test)]
#[path = "definition_tests.rs"]
mod tests;
