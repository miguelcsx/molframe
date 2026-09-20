//! The per-format convenience verbs: the readers and writers a caller names
//! by format, which the facade's dispatch does not need in line.
//!
//! The file ladder (`read`, `write`, `write_with_options`) lives in the parent
//! and streams through the same kernels these wrap. Each format's verbs live
//! in their own file, gated by that format's feature, so a lean build never
//! carries another format's imports.

#[cfg(feature = "bcif")]
mod bcif;
#[cfg(feature = "chemistry")]
mod chem;
#[cfg(feature = "mmcif")]
mod mmcif;
#[cfg(feature = "pdb")]
mod pdb;
#[cfg(feature = "geometry")]
mod transform;

#[cfg(feature = "bcif")]
pub use bcif::{write_bcif, write_bcif_with_options};
#[cfg(feature = "chemistry")]
pub use chem::read_component_dictionary;
#[cfg(feature = "mmcif")]
pub(crate) use mmcif::{cif_write_findings, write_mmcif_to_with_options};
#[cfg(feature = "mmcif")]
pub use mmcif::{read_document, write_mmcif, write_mmcif_with_options};
#[cfg(feature = "pdb")]
pub use pdb::write_pdb;
#[cfg(feature = "geometry")]
pub use transform::transform;
