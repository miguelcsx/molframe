//! Typed projection of residue confidence onto existing atom rows.

use crate::{LocalMetric, ModelCif, ModelCifError};
use pdbiox_core::annotation::AnnotationColumn;
use pdbiox_core::column::Presence;
use pdbiox_core::{ModelIndex, Structure};
use std::collections::BTreeMap;

/// Projects residue-scoped pLDDT to one atom-aligned confidence column.
///
/// Values are copied once from compact `ModelCIF` metadata into immutable native
/// property storage. Every atom in a residue receives the same confidence.
/// Polymer residues omitted from the metric table are `Unknown`; rows without
/// a label sequence identity are `Inapplicable`. The values remain confidence
/// scores and are never written to displacement or uncertainty columns.
///
/// The projection costs `O(m log m + r log m + a)`, where `m`, `r`, and `a`
/// are pLDDT rows, residues, and atoms. Its retained storage is one value plus
/// two validity bits per atom in the mixed-presence case.
///
/// # Errors
///
/// Rejects model ambiguity, duplicate residue metrics, or an atom count that
/// exceeds the annotation column domain.
pub(crate) fn lower_plddt_annotation(
    structure: &Structure,
    model: &ModelCif,
) -> Result<Option<AnnotationColumn<f64>>, ModelCifError> {
    let mut metrics = BTreeMap::new();
    let mut model_id: Option<&str> = None;
    for metric in model.confidence().plddt() {
        match model_id {
            Some(first) if first != metric.model_id => {
                return Err(ModelCifError::MultiplePlddtModels {
                    first: first.into(),
                    second: metric.model_id.into(),
                });
            }
            None => model_id = Some(metric.model_id),
            Some(_) => {}
        }
        let key = (metric.chain_id, metric.sequence_id);
        if metrics.insert(key, metric).is_some() {
            return Err(ModelCifError::DuplicatePlddtResidue {
                chain: metric.chain_id.into(),
                sequence: metric.sequence_id,
            });
        }
    }
    if metrics.is_empty() {
        return Ok(None);
    }
    if structure.model_count() != 1 {
        return Err(ModelCifError::MultipleStructureModels {
            count: structure.model_count(),
        });
    }

    let mut entries = Vec::new();
    entries
        .try_reserve_exact(structure.atom_count() as usize)
        .map_err(|_| ModelCifError::Capacity)?;
    let Some(structure_model) = structure.model(ModelIndex::new(0)) else {
        return AnnotationColumn::from_entries(entries)
            .map(Some)
            .map_err(|_| ModelCifError::Capacity);
    };
    for chain in structure_model.chains() {
        let chain_id = chain.label();
        for residue in chain.residues() {
            let sequence = residue.label_seq_id().map(i64::from);
            let metric = chain_id
                .zip(sequence)
                .and_then(|key| metrics.get(&key).copied())
                .filter(|value| component_matches(*value, residue.name()));
            let (value, presence) = match (sequence, metric) {
                (_, Some(metric)) => (metric.value, Presence::Present),
                (Some(_), None) => (0.0, Presence::Unknown),
                (None, None) => (0.0, Presence::Inapplicable),
            };
            entries.extend(residue.atoms().map(|_| (value, presence)));
        }
    }
    AnnotationColumn::from_entries(entries)
        .map(Some)
        .map_err(|_| ModelCifError::Capacity)
}

impl ModelCif {
    /// Projects this model's residue pLDDT onto an existing atom topology.
    ///
    /// # Errors
    ///
    /// Rejects model ambiguity, duplicate residue metrics, or capacity overflow.
    pub fn plddt_annotation(
        &self,
        structure: &Structure,
    ) -> Result<Option<AnnotationColumn<f64>>, ModelCifError> {
        lower_plddt_annotation(structure, self)
    }
}

fn component_matches(metric: LocalMetric<'_>, residue: Option<&str>) -> bool {
    match metric.component_id {
        Some(component) => residue == Some(component),
        None => true,
    }
}

#[cfg(test)]
#[path = "plddt_tests.rs"]
mod tests;
