//! Where the structure an analysis was built from came from.

use molframe_core::contract::{ContentDigest, Provenance, SourceRef};
use molframe_core::structure::Structure;

/// Where the structure was read from, when the reader said.
pub(super) fn input_source(template: &Structure) -> SourceRef {
    match &template.data().entry.input_name {
        Some(name) => SourceRef::path(name.as_ref()),
        None => SourceRef::Memory,
    }
}

pub(super) fn with_digest(provenance: Provenance, digest: Option<ContentDigest>) -> Provenance {
    match digest {
        Some(digest) => provenance.with_input_digest(digest),
        None => provenance,
    }
}
