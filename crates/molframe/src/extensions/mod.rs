//! The ergonomic layer: domain operations reachable as methods on a handle.
//!
//! Each submodule here declares one extension trait over [`Structure`], and
//! each trait method is a forwarder — the body it calls is the free function in
//! the facet's own module, reached by static dispatch through
//! [`Structure::engine`]. There is no vtable, no wrapper value and no second
//! implementation: `structure.atom_contacts(4.5)` and
//! `molframe::analysis::atom_contacts(structure.engine(), 4.5, backend, ctx)`
//! execute the same code, and a caller who wants the long form writes it.
//!
//! # Why the methods carry the kernel's name
//!
//! A method name is the kernel function's name, exactly, and the `_with`
//! suffix marks the form that takes the backend and the [`ExecutionContext`]
//! explicitly. Inventing a shorter alias for a method — `contacts` for
//! `atom_contacts` — would put two names on one concept, which is the drift
//! this façade exists to remove. `use molframe::prelude::*` brings the traits
//! in; the method and the function then coexist without ambiguity, because a
//! method call is not a path.
//!
//! # Which domains have a trait, and which do not
//!
//! A trait exists where the domain's kernels take a structure first:
//! [`AnalysisExt`], [`ValidationExt`] and [`CompareExt`]. Geometry, surface and
//! the spatial planners are coordinate-first — their kernels take
//! `positions: &[[f32; 3]]` and the radii or selections that go with it — so a
//! structure-level method would not be a forwarder. It would be new façade
//! logic deciding what to gather and what to do about an atom whose element
//! has no radius in the chosen set, and that is a design decision with an error
//! vocabulary of its own, not a rewrite of a call. Those kernels stay free
//! functions over coordinates, reached with
//! [`Selection::to_coordinates`](crate::Selection::to_coordinates) and
//! [`Structure::coordinates`](crate::Structure::coordinates) as the two
//! borrow-or-gather doors.
//!
//! [`ExecutionContext`]: crate::ExecutionContext

#[cfg(feature = "analysis")]
mod analysis;
#[cfg(feature = "analysis")]
mod analysis_impl;
#[cfg(feature = "compare")]
mod compare;
#[cfg(feature = "validation")]
mod validation;

#[cfg(feature = "analysis")]
pub use analysis::AnalysisExt;
#[cfg(feature = "compare")]
pub use compare::CompareExt;
#[cfg(feature = "validation")]
pub use validation::ValidationExt;

#[cfg(test)]
#[path = "extension_tests.rs"]
mod tests;
