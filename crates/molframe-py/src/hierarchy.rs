//! Allocation-free Python navigation over the facade hierarchy.

use pyo3::prelude::*;

mod atom;
mod chain;
mod model;
mod residue;

pub(crate) use atom::{PyAtom, PyAtoms};
pub(crate) use chain::{PyChain, PyChains};
pub(crate) use model::{PyModel, PyModels};
pub(crate) use residue::{PyResidue, PyResidueSelection, PyResidues};

pub(super) fn position(index: isize, len: usize) -> PyResult<usize> {
    let len_signed = isize::try_from(len)
        .map_err(|_| pyo3::exceptions::PyOverflowError::new_err("collection is too large"))?;
    let normalized = if index < 0 { len_signed + index } else { index };
    if normalized < 0 || normalized >= len_signed {
        return Err(crate::error::index("hierarchy index is out of range"));
    }
    usize::try_from(normalized).map_err(|_| crate::error::index("invalid hierarchy index"))
}

pub(super) fn contiguous_range<I>(mut indices: I) -> (u32, usize)
where
    I: Iterator<Item = u32>,
{
    match indices.next() {
        Some(first) => (first, 1 + indices.count()),
        None => (0, 0),
    }
}
