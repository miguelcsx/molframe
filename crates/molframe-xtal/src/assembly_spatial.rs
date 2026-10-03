//! Spatial search over lazy biological-assembly transforms.

use crate::AssemblyView;
use molframe_core::selection::AtomSelection;
use molframe_core::{AtomIndex, Diagnostic, ExecutionContext, InstanceId, ModelIndex};
use molframe_spatial::{SpatialBackend, pairs_within};

/// One unique unordered pair between generated assembly atoms.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AssemblyNeighbor {
    /// Stable chain instance containing the first atom.
    pub first_instance: InstanceId,
    /// Deposited atom represented by the first assembly atom.
    pub first_atom: AtomIndex,
    /// Stable chain instance containing the second atom.
    pub second_instance: InstanceId,
    /// Deposited atom represented by the second assembly atom.
    pub second_atom: AtomIndex,
    /// Squared Cartesian distance in ångström².
    pub distance_squared: f32,
}

/// A geometrically perceived covalent link between distinct chain instances.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AssemblyBond {
    /// Stable first chain instance identifier.
    pub first_instance: InstanceId,
    /// Deposited first atom index.
    pub first_atom: AtomIndex,
    /// Stable second chain instance identifier.
    pub second_instance: InstanceId,
    /// Deposited second atom index.
    pub second_atom: AtomIndex,
    /// Perceived order; a distance alone establishes only a single bond.
    pub order: molframe_core::BondOrder,
}

impl AssemblyView {
    /// Finds cross-instance covalent links with the canonical chemistry
    /// predicate over the shared assembly spatial search. No topology or
    /// expanded structure is copied; search coordinates are temporary.
    ///
    /// Metal coordination and same-instance pairs are excluded. Results use
    /// stable chain instance identifiers, not operator-list indices.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic for an unavailable model or spatial search failure.
    pub fn covalent_links(
        &self,
        model: ModelIndex,
        context: &ExecutionContext,
    ) -> Result<Vec<AssemblyBond>, Diagnostic> {
        let pairs = self.neighbors(model, 2.6, SpatialBackend::Auto, context)?;
        Ok(pairs
            .into_iter()
            .filter_map(|pair| {
                if pair.first_instance == pair.second_instance {
                    return None;
                }
                let first = self.source().atom(pair.first_atom)?;
                let second = self.source().atom(pair.second_atom)?;
                molframe_chem::covalent_pair(first, second, pair.distance_squared).then_some(
                    AssemblyBond {
                        first_instance: pair.first_instance,
                        first_atom: pair.first_atom,
                        second_instance: pair.second_instance,
                        second_atom: pair.second_atom,
                        order: molframe_core::BondOrder::Single,
                    },
                )
            })
            .collect())
    }
}

impl AssemblyView {
    /// Searches transformed assembly atoms without materialising a structure.
    ///
    /// Results are ordered by the assembly's stable virtual-atom order.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic for an invalid cutoff or unavailable model.
    pub fn neighbors(
        &self,
        model: ModelIndex,
        cutoff: f32,
        backend: SpatialBackend,
        context: &ExecutionContext,
    ) -> Result<Vec<AssemblyNeighbor>, Diagnostic> {
        let mut atoms = Vec::new();
        let mut positions = Vec::new();
        for (instance, position) in self.positions(model) {
            let Some(position) = position else {
                continue;
            };
            if !position.iter().all(|value| value.is_finite()) {
                continue;
            }
            atoms.push((instance.instance_id, instance.source_atom));
            positions.push(position);
        }
        if self.source().model_positions(model).is_none() {
            return Err(molframe_core::Diagnostic::new(molframe_core::Code::E6003)
                .with_context("model", model.to_string()));
        }
        let count = u32::try_from(positions.len())
            .map_err(|_| molframe_core::Diagnostic::new(molframe_core::Code::E1901))?;
        let selection = AtomSelection::All(count);
        // Materialising is inherent here: the result is one neighbour record per
        // pair, so streaming the intermediate would remove a constant factor and
        // not the growth. Bounding it means changing what this returns.
        let pairs = pairs_within(
            &positions, &selection, &selection, cutoff, backend, None, context,
        )
        .map_err(molframe_spatial::SpatialError::into_diagnostic)?;
        pairs
            .into_iter()
            .map(|pair| {
                let first = atoms
                    .get(pair.first as usize)
                    .copied()
                    .ok_or_else(invariant)?;
                let second = atoms
                    .get(pair.second as usize)
                    .copied()
                    .ok_or_else(invariant)?;
                Ok(AssemblyNeighbor {
                    first_instance: first.0,
                    first_atom: first.1,
                    second_instance: second.0,
                    second_atom: second.1,
                    distance_squared: pair.distance_squared,
                })
            })
            .collect()
    }
}

fn invariant() -> Diagnostic {
    Diagnostic::new(molframe_core::Code::E9001)
}

#[cfg(test)]
#[path = "assembly_spatial_tests.rs"]
mod tests;
