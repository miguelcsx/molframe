//! Reusable structure-bound queries evaluated against changing frame coordinates.

use crate::Timestep;
use pdbiox_core::ExecutionContext;
use pdbiox_core::contract::AnalysisPolicy;
use pdbiox_core::diagnostic::{Diagnostic, Strictness};
use pdbiox_core::index::ModelIndex;
use pdbiox_core::structure::{Structure, validate};
use pdbiox_query::{Evaluation, Groups, PhysicalQuery, Query};
use pdbiox_spatial::{SpatialBackend, StructureSpatial};

/// Failure while projecting one frame onto a fixed topology and query plan.
#[derive(Clone, Debug, thiserror::Error)]
#[non_exhaustive]
pub enum UpdatingSelectionError {
    /// Coordinate copy-on-write could not fit the shared execution account.
    #[error("frame coordinate copy exceeds the execution memory budget: {0}")]
    Memory(#[from] pdbiox_core::MemoryBudgetError),
    /// The fixed topology atom count cannot be represented on this platform.
    #[error("topology atom count exceeds the platform index range")]
    AtomCountOverflow,
    /// Frame and topology atom counts differ.
    #[error("updating selection expected {expected} atoms, found {found}")]
    AtomCountMismatch {
        /// Topology atom count.
        expected: usize,
        /// Frame atom count.
        found: usize,
    },
    /// Coordinate publication violated a structure invariant.
    #[error("frame coordinates are invalid: {0:?}")]
    Coordinates(Vec<Diagnostic>),
    /// The spatial resolver could not honour the policy or cell.
    #[error("frame spatial resolver failed: {0}")]
    Spatial(Diagnostic),
    /// Query evaluation failed for this frame.
    #[error("frame query failed: {0:?}")]
    Query(Vec<Diagnostic>),
}

/// One compiled query rebound only to coordinates that change each frame.
#[derive(Clone, Debug)]
pub struct UpdatingSelection {
    topology: Structure,
    plan: PhysicalQuery,
    policy: AnalysisPolicy,
    groups: Groups,
    backend: SpatialBackend,
}

impl UpdatingSelection {
    /// Binds symbols and static predicates once against a fixed topology.
    #[must_use]
    pub fn new(
        topology: Structure,
        query: &Query,
        policy: AnalysisPolicy,
        groups: Groups,
        backend: SpatialBackend,
    ) -> Self {
        let plan = query.plan(&topology, &policy);
        Self {
            topology,
            plan,
            policy,
            groups,
            backend,
        }
    }

    /// Evaluates static and geometric predicates against one timestep.
    ///
    /// # Errors
    ///
    /// Refuses atom-count mismatches, invalid coordinates, unavailable periodic
    /// metadata and query evaluation failures.
    pub fn evaluate(
        &self,
        timestep: &Timestep,
        context: &ExecutionContext,
    ) -> Result<Evaluation, UpdatingSelectionError> {
        let expected = usize::try_from(self.topology.atom_count())
            .map_err(|_| UpdatingSelectionError::AtomCountOverflow)?;
        if timestep.positions.len() != expected {
            return Err(UpdatingSelectionError::AtomCountMismatch {
                expected,
                found: timestep.positions.len(),
            });
        }
        let mut editor = self.topology.edit_coordinates(context)?;
        let positions = editor
            .positions_mut(ModelIndex::new(0))
            .ok_or_else(|| UpdatingSelectionError::Coordinates(Vec::new()))?;
        positions.copy_from_slice(&timestep.positions);
        let structure = editor
            .commit()
            .map_err(UpdatingSelectionError::Coordinates)?;
        let structure = with_frame_cell(structure, timestep.cell)?;
        let spatial = StructureSpatial::new(&structure, &self.policy, self.backend, context)
            .map_err(UpdatingSelectionError::Spatial)?;
        self.plan
            .evaluate(&structure, &self.policy, &self.groups, Some(&spatial))
            .map_err(UpdatingSelectionError::Query)
    }
}

fn with_frame_cell(
    structure: Structure,
    cell: Option<pdbiox_core::structure::UnitCell>,
) -> Result<Structure, UpdatingSelectionError> {
    let Some(cell) = cell else {
        return Ok(structure);
    };
    let mut data = structure.data().clone();
    data.cell = Some(cell);
    let findings = validate(&data);
    if findings
        .iter()
        .any(|finding| finding.is_error(Strictness::Medium))
    {
        return Err(UpdatingSelectionError::Coordinates(findings));
    }
    Ok(Structure::new(data))
}

#[cfg(test)]
#[path = "selection_tests.rs"]
mod tests;
