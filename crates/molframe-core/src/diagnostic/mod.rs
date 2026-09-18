//! Findings about data, accumulated and returned alongside success.
//!
//! A reader hands back the structure *and* the list of everything wrong with it.
//! Stopping at the first problem makes real archive files unusable, and throwing
//! warnings away on success makes their problems invisible — which is worse,
//! because the caller then has a plausible result and no reason to doubt it.
//!
//! Findings are ordered before they are returned: worst first, then by position
//! in the file, then by code. Two runs over the same input therefore produce the
//! same list in the same order, which is what makes a result hashable and a
//! regression test meaningful.

mod code;
mod finding;
mod findings;
mod registry;
mod render;

pub use code::{Class, Code, Kind, Severity, Strictness};
pub use finding::{ContextItem, Diagnostic, Diagnostics};
pub use findings::Findings;
pub use render::Rendered;
