//! A dense PAE matrix together with the residues and metric it describes.

use super::pae::PaeMatrix;
use super::residues::ResidueIndex;
use crate::PairwiseMetric;

/// Residue-pair aligned error for one model, stored as a matrix.
///
/// Row `i` and column `j` are the residues at positions `i` and `j` of
/// [`Self::residues`]. This replaces the pair-per-row `ModelCIF` table, whose
/// seven columns cost far more per pair than one matrix cell.
#[derive(Clone, Debug, PartialEq)]
pub struct DensePae {
    pub(crate) model_id: Box<str>,
    pub(crate) metric_id: Box<str>,
    pub(crate) residues: ResidueIndex,
    pub(crate) items: Vec<Box<str>>,
    pub(crate) slot: usize,
    pub(crate) matrix: PaeMatrix,
}

impl DensePae {
    /// Identifier of the model the errors belong to.
    #[must_use]
    pub fn model_id(&self) -> &str {
        &self.model_id
    }

    /// Identifier of the QA metric definition.
    #[must_use]
    pub fn metric_id(&self) -> &str {
        &self.metric_id
    }

    /// Residues labelling the matrix rows and columns.
    #[must_use]
    pub const fn residues(&self) -> &ResidueIndex {
        &self.residues
    }

    /// The error matrix.
    #[must_use]
    pub const fn matrix(&self) -> &PaeMatrix {
        &self.matrix
    }

    /// Items of the serialised table in file order.
    #[must_use]
    pub fn items(&self) -> &[Box<str>] {
        &self.items
    }

    /// Every ordered residue pair in row-major order, as the table would list them.
    pub(crate) fn pairs(&self) -> impl Iterator<Item = PairwiseMetric<'_>> {
        let size = self.matrix.size();
        (0..size * size).filter_map(move |index| {
            let first = self.residues.get(index / size)?;
            let second = self.residues.get(index % size)?;
            Some(PairwiseMetric {
                model_id: &self.model_id,
                first_chain_id: first.chain,
                first_sequence_id: first.sequence,
                second_chain_id: second.chain,
                second_sequence_id: second.sequence,
                metric_id: &self.metric_id,
                value: self.matrix.get_f64(index / size, index % size)?,
            })
        })
    }
}
