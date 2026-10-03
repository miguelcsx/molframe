//! Readers for the reference inputs that validation checks take as explicit data.
//!
//! The validation kernels refuse to carry a built-in scientific reference, so a
//! reference distribution, a rotamer profile, a set of plane restraints or a set
//! of TLS groups arrives from the caller. These functions read those inputs
//! from JSON or TOML (chosen by suffix) into the kernels' own types and refuse
//! anything the kernels would refuse, naming the field at fault. Each schema is
//! documented on its reader.

mod error;
mod library;
mod rotamer;
mod selections;

pub use error::ValidationInputError;
pub use library::read_reference_library;
pub use rotamer::read_rotamer_profile;
pub use selections::{read_plane_restraints, read_tls_groups};
