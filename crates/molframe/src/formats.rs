//! The structure formats, grouped by the crate that reads and writes them.
//!
//! The read and write verbs at the crate root ([`crate::read`], [`crate::write`],
//! [`crate::write_mmcif`], …) dispatch across whichever of these are linked;
//! this module is where each format's own types live — the document model, the
//! per-family write options, the header records — for a caller that needs more
//! than the dispatched verbs give it.

#[cfg(feature = "mmcif")]
pub use molframe_cif as cif;

#[cfg(feature = "bcif")]
pub use molframe_bcif as bcif;

#[cfg(feature = "pdb")]
pub use molframe_pdb as pdb;

#[cfg(feature = "modelcif")]
pub use molframe_modelcif as modelcif;
