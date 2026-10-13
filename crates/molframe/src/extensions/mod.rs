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
//! [`AnalysisExt`], [`ValidationExt`] and [`CompareExt`]. [`SurfaceExt`] is the
//! one deliberate exception: surface kernels are coordinate-first, so its
//! methods are not forwarders. They gather the coordinates of one model, take
//! radii from a [`RadiusSet`](molframe_chem::RadiusSet), and give an atom with
//! no radius a `NaN` area that does not occlude its neighbours. That policy is
//! stated once, in `surface`, instead of being rediscovered by every caller.
//! Geometry and the spatial planners stay free functions over coordinates,
//! reached with [`Selection::to_coordinates`](crate::Selection::to_coordinates)
//! and [`Structure::coordinates`](crate::Structure::coordinates).
//!
//! [`ExecutionContext`]: crate::ExecutionContext

#[cfg(feature = "analysis")]
mod analysis;
#[cfg(feature = "analysis")]
mod analysis_impl;
#[cfg(feature = "compare")]
mod compare;
#[cfg(feature = "compare")]
mod compare_impl;
#[cfg(all(feature = "surface", feature = "chemistry"))]
mod surface;
#[cfg(feature = "validation")]
mod validation;

#[cfg(feature = "analysis")]
pub use analysis::AnalysisExt;
#[cfg(feature = "compare")]
pub use compare::CompareExt;
#[cfg(all(feature = "surface", feature = "chemistry"))]
pub use surface::{SasaMethod, SasaOptions, SasaReport, SurfaceExt, SurfaceExtError};
#[cfg(feature = "validation")]
pub use validation::ValidationExt;

#[cfg(test)]
#[path = "extension_tests.rs"]
mod tests;

#[cfg(all(test, feature = "surface", feature = "chemistry"))]
#[path = "surface_tests.rs"]
mod surface_tests;

#[cfg(all(test, feature = "surface", feature = "chemistry"))]
#[path = "surface_reference_tests.rs"]
mod surface_reference_tests;
