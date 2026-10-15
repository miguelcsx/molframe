//! Lossless, order-preserving CIF document model.

mod model;
mod number;

pub(crate) use model::CifValueRef;
pub use model::{Category, CifValue, Column, DataBlock, Document};
