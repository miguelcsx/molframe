//! Deterministic native plan execution.

pub(super) mod inputs;
pub(super) mod value;

#[cfg(feature = "compare")]
use super::comparison;
use super::geometry;
use super::physical;
use super::requests::ContactsRequest;
use super::spatial;
use super::spatial_cache::SpatialContext;
use super::structure;
#[cfg(feature = "surface")]
use super::surface;
#[cfg(feature = "traj")]
use super::trajectory;
use crate::QueryStructure;
use inputs::{CoordinateInput, FrameInput, IndexInput, PlanInput};
use molframe_analysis::Contact;
use molframe_core::contract::{Analysis, AnalysisPolicy, Coverage, Status};
use molframe_core::diagnostic::Findings;
use molframe_core::structure::Structure;
use molframe_query::Groups;
use molframe_spatial::{PeriodicBox, SpatialBackend, SpatialSearchOptions, StructureSpatial};
use std::collections::BTreeMap;
use std::fmt;
use value::{ExecutionPlanError, PlanOperation, PlanResult, PlanResultEntry, PlanValue};

/// A deterministic native operation plan.
#[derive(Clone, Debug, Default)]
pub struct Plan {
    operations: BTreeMap<Box<str>, PlanOperation>,
}

impl Plan {
    /// Starts an empty plan.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            operations: BTreeMap::new(),
        }
    }

    /// Adds one typed operation.
    ///
    /// Every family's request converts into a [`PlanOperation`], so a plan is
    /// written in one voice regardless of family:
    /// `plan.add("contacts", ContactsRequest::new(..)?)`.
    ///
    /// # Errors
    ///
    /// Returns [`ExecutionPlanError::DuplicateId`] when the stable id is already used.
    pub fn add(
        &mut self,
        id: impl Into<Box<str>>,
        operation: impl Into<PlanOperation>,
    ) -> Result<(), ExecutionPlanError> {
        let id = id.into();
        if self.operations.contains_key(&id) {
            return Err(ExecutionPlanError::DuplicateId(id));
        }
        self.operations.insert(id, operation.into());
        Ok(())
    }

    /// Number of operation nodes.
    #[must_use]
    pub fn len(&self) -> usize {
        self.operations.len()
    }

    /// Whether the plan has no operation nodes.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.operations.is_empty()
    }

    /// Executes all nodes in Rust without materialising borrowed arrays.
    ///
    /// A default, non-disordered structure gets one shared spatial resolver,
    /// so compatible contacts reuse bounded indices.
    ///
    /// # Errors
    ///
    /// Returns a typed setup or native-kernel error without publishing partial
    /// results.
    pub fn execute(
        &self,
        input: PlanInput<'_>,
        context: &molframe_core::ExecutionContext,
    ) -> Result<PlanResult, ExecutionPlanError> {
        let has_contacts = self.operations.values().any(|operation| {
            matches!(
                operation,
                PlanOperation::Contacts(_) | PlanOperation::Selection(_)
            )
        });
        let default_policy = AnalysisPolicy::default();
        let reusable_spatial =
            reusable_spatial(input.structure, has_contacts, &default_policy, context)?;
        let mut coordinate_contexts = coordinate_contexts(&self.operations, input.arrays)?;
        let mut entries = Vec::with_capacity(self.operations.len());
        for (id, operation) in &self.operations {
            let value = Self::execute_one(
                id,
                operation,
                input,
                reusable_spatial.as_ref(),
                &mut coordinate_contexts,
                context,
            )?;
            entries.push(PlanResultEntry {
                id: id.clone(),
                value,
            });
        }
        let structure_index_count = reusable_spatial
            .as_ref()
            .map_or(0, StructureSpatial::cached_index_count);
        let cached_index_count = coordinate_contexts
            .iter()
            .map(SpatialContext::cached_index_count)
            .fold(structure_index_count, usize::saturating_add);
        Ok(PlanResult {
            entries,
            cached_index_count,
        })
    }

    fn execute_one(
        id: &str,
        operation: &PlanOperation,
        input: PlanInput<'_>,
        spatial: Option<&StructureSpatial<'_>>,
        coordinate_contexts: &mut [SpatialContext<'_>],
        context: &molframe_core::ExecutionContext,
    ) -> Result<PlanValue, ExecutionPlanError> {
        match operation {
            PlanOperation::Selection(request) => {
                let structure = required_structure(id, input.structure)?;
                let groups = Groups::new();
                let evaluation = if request.policy() == &AnalysisPolicy::default()
                    && let Some(spatial) = spatial
                {
                    request
                        .query()
                        .evaluate(structure, request.policy(), &groups, Some(spatial))
                        .map_err(Findings::from)
                } else {
                    structure.select_with_options(
                        request.query(),
                        request.policy(),
                        &groups,
                        context,
                    )
                };
                evaluation
                    .map(|value| PlanValue::Selection(Box::new(value)))
                    .map_err(ExecutionPlanError::Selection)
            }
            PlanOperation::Contacts(request) => {
                let structure = required_structure(id, input.structure)?;
                let spatial = spatial.filter(|_| request.policy() == &AnalysisPolicy::default());
                Ok(PlanValue::Contacts(Box::new(execute_contacts(
                    request, structure, spatial, context,
                )?)))
            }
            PlanOperation::Rmsd(request) => {
                let mobile = array_slot(id, input.arrays, request.mobile().slot())?;
                let reference = array_slot(id, input.arrays, request.reference().slot())?;
                Ok(PlanValue::Rmsd(
                    molframe_geom::rmsd(mobile.positions, reference.positions)
                        .map_err(ExecutionPlanError::Rmsd)?,
                ))
            }
            PlanOperation::Structure(request) => {
                let structure = required_structure(id, input.structure)?;
                Ok(PlanValue::Structure(Box::new(structure::execute(
                    request, structure, context,
                )?)))
            }
            PlanOperation::Physical(request) => {
                let structure = required_structure(id, input.structure)?;
                Ok(PlanValue::Physical(Box::new(physical::execute(
                    id, request, structure, input, context,
                )?)))
            }
            PlanOperation::Geometry(request) => Ok(PlanValue::Geometry(Box::new(
                geometry::execute(id, request, input)?,
            ))),
            PlanOperation::Spatial(request) => Ok(PlanValue::Spatial(Box::new(spatial::execute(
                id,
                request,
                input,
                coordinate_contexts
                    .iter_mut()
                    .find(|context| context.matches(request)),
                context,
            )?))),
            #[cfg(feature = "surface")]
            PlanOperation::Surface(request) => Ok(PlanValue::Surface(Box::new(surface::execute(
                id, request, input, context,
            )?))),
            PlanOperation::BondInference(options) => {
                let structure = required_structure(id, input.structure)?;
                crate::infer_bonds(structure, *options, context)
                    .map(|report| PlanValue::BondInference(Box::new(report)))
                    .map_err(ExecutionPlanError::Chemistry)
            }
            #[cfg(feature = "compare")]
            PlanOperation::Comparison(request) => {
                let mobile = array_slot(id, input.arrays, request.mobile().slot())?;
                let reference = array_slot(id, input.arrays, request.reference().slot())?;
                Ok(PlanValue::Comparison(comparison::execute(
                    request,
                    mobile.positions,
                    reference.positions,
                    context,
                )?))
            }
            #[cfg(feature = "traj")]
            PlanOperation::Trajectory(request) => {
                let frames = frame_slot(id, input.frames, request.frame_slot())?;
                let atoms = match request.atom_slot() {
                    Some(slot) => index_slot(id, input.indices, slot)?,
                    None => &[],
                };
                Ok(PlanValue::Trajectory(Box::new(trajectory::execute(
                    request, frames, atoms,
                )?)))
            }
        }
    }
}

fn coordinate_contexts<'a>(
    operations: &BTreeMap<Box<str>, PlanOperation>,
    arrays: &'a [CoordinateInput<'a>],
) -> Result<Vec<SpatialContext<'a>>, ExecutionPlanError> {
    let mut contexts: Vec<SpatialContext<'a>> = Vec::new();
    for (id, operation) in operations {
        let PlanOperation::Spatial(request) = operation else {
            continue;
        };
        if contexts.iter().any(|context| context.matches(request)) {
            continue;
        }
        let (slot, options, periodic) = spatial_key(request);
        let input = array_slot(id, arrays, slot)?;
        contexts.push(SpatialContext::new(
            slot,
            input.positions,
            options,
            periodic,
        )?);
    }
    Ok(contexts)
}

fn spatial_key(
    request: &super::spatial::SpatialRequest,
) -> (usize, SpatialSearchOptions, Option<PeriodicBox>) {
    match request {
        super::spatial::SpatialRequest::NeighborPairs {
            positions,
            options,
            periodic,
            ..
        }
        | super::spatial::SpatialRequest::AtomsWithin {
            positions,
            options,
            periodic,
            ..
        } => (*positions, *options, *periodic),
    }
}

fn reusable_spatial<'a>(
    structure: Option<&'a Structure>,
    needed: bool,
    policy: &AnalysisPolicy,
    context: &'a molframe_core::ExecutionContext,
) -> Result<Option<StructureSpatial<'a>>, ExecutionPlanError> {
    let Some(structure) = structure.filter(|_| needed) else {
        return Ok(None);
    };
    let resolution = structure.resolve_altlocs(policy);
    if resolution.status != Status::Complete || resolution.coverage.ambiguous != 0 {
        return Ok(None);
    }
    StructureSpatial::new_with_options(
        structure,
        policy,
        SpatialSearchOptions::with_backend(SpatialBackend::Auto),
        context,
    )
    .map(Some)
    .map_err(ExecutionPlanError::Spatial)
}

fn required_structure<'a>(
    operation: &str,
    structure: Option<&'a Structure>,
) -> Result<&'a Structure, ExecutionPlanError> {
    structure.ok_or_else(|| ExecutionPlanError::MissingInput {
        operation: operation.into(),
        input: "structure",
    })
}

fn array_slot<'a>(
    operation: &str,
    arrays: &'a [CoordinateInput<'a>],
    slot: usize,
) -> Result<&'a CoordinateInput<'a>, ExecutionPlanError> {
    arrays
        .get(slot)
        .ok_or_else(|| ExecutionPlanError::ArraySlot {
            operation: operation.into(),
            slot,
        })
}

fn frame_slot<'a>(
    operation: &str,
    frames: &'a [FrameInput<'a>],
    slot: usize,
) -> Result<&'a FrameInput<'a>, ExecutionPlanError> {
    frames
        .get(slot)
        .ok_or_else(|| ExecutionPlanError::FrameSlot {
            operation: operation.into(),
            slot,
        })
}

fn index_slot<'a>(
    operation: &str,
    indices: &'a [IndexInput<'a>],
    slot: usize,
) -> Result<&'a [usize], ExecutionPlanError> {
    indices
        .get(slot)
        .map(|input| input.indices)
        .ok_or_else(|| ExecutionPlanError::IndexSlot {
            operation: operation.into(),
            slot,
        })
}

fn execute_contacts(
    request: &ContactsRequest,
    structure: &Structure,
    spatial: Option<&StructureSpatial<'_>>,
    context: &molframe_core::ExecutionContext,
) -> Result<Analysis<Vec<Contact>>, ExecutionPlanError> {
    let groups = Groups::new();
    let (left, right, contacts) = if let Some(spatial) = spatial {
        let left = request
            .left_query()
            .evaluate(structure, request.policy(), &groups, Some(spatial))
            .map_err(|findings| ExecutionPlanError::Governed(format!("{findings:?}").into()))?;
        let right = request
            .right_query()
            .evaluate(structure, request.policy(), &groups, Some(spatial))
            .map_err(|findings| ExecutionPlanError::Governed(format!("{findings:?}").into()))?;
        let contacts = molframe_analysis::atom_contacts_between_with_spatial(
            structure,
            &left.selection,
            &right.selection,
            request.cutoff(),
            request.backend(),
            spatial,
        )
        .map_err(ExecutionPlanError::Spatial)?;
        (left, right, contacts)
    } else {
        let left = structure
            .select_with_options(request.left_query(), request.policy(), &groups, context)
            .map_err(|findings| ExecutionPlanError::Governed(format!("{findings:?}").into()))?;
        let right = structure
            .select_with_options(request.right_query(), request.policy(), &groups, context)
            .map_err(|findings| ExecutionPlanError::Governed(format!("{findings:?}").into()))?;
        let contacts = molframe_analysis::atom_contacts_between(
            structure,
            &left.selection,
            &right.selection,
            request.cutoff(),
            request.backend(),
            context,
        )
        .map_err(|error| ExecutionPlanError::Governed(error.to_string().into()))?;
        (left, right, contacts)
    };
    let mut analysis = Analysis::complete(
        contacts,
        Coverage::complete(structure.atom_count()),
        request.policy(),
    );
    analysis.warnings.extend(left.warnings);
    analysis.warnings.extend(right.warnings);
    Ok(analysis)
}

impl fmt::Display for Plan {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("Plan")
            .field("operation_count", &self.operations.len())
            .field(
                "operation_ids",
                &self
                    .operations
                    .keys()
                    .map(AsRef::as_ref)
                    .collect::<Vec<_>>(),
            )
            .finish()
    }
}
