//! The extension traits a [`Structure`] handle carries.
//!
//! The crystal and `ModelCIF` metadata traits are implemented by the crates that
//! own the metadata, on the engine's snapshot. Forwarding them onto the
//! facade's own handle keeps a caller on the curated type: `structure.assembly("1")`
//! resolves by static dispatch to the same body it always did, one reference
//! narrower, and nothing has to reach for [`Structure::engine`] to write the
//! ordinary call.
//!
//! Each impl here is a forwarder, not a second implementation — the trait has
//! exactly one body, in the crate that defines it.

#[cfg(any(feature = "crystal", feature = "modelcif"))]
use super::handle::Structure;

#[cfg(feature = "crystal")]
use molframe_core::index::ModelIndex;
#[cfg(feature = "crystal")]
use molframe_core::{ExecutionContext, diagnostic::Diagnostic};

#[cfg(feature = "crystal")]
impl molframe_xtal::AssemblyExt for Structure {
    fn assembly_set(&self) -> Option<&molframe_xtal::AssemblySet> {
        molframe_xtal::AssemblyExt::assembly_set(self.engine())
    }

    fn assembly(&self, id: &str) -> Result<molframe_xtal::AssemblyView, Diagnostic> {
        molframe_xtal::AssemblyExt::assembly(self.engine(), id)
    }
}

#[cfg(feature = "crystal")]
impl molframe_xtal::NcsExt for Structure {
    fn ncs_set(&self) -> Option<&molframe_xtal::NcsSet> {
        molframe_xtal::NcsExt::ncs_set(self.engine())
    }

    fn ncs_generated(&self) -> Option<molframe_xtal::NcsView<'_>> {
        molframe_xtal::NcsExt::ncs_generated(self.engine())
    }
}

#[cfg(feature = "crystal")]
impl molframe_xtal::SymmetryExt for Structure {
    fn symmetry_set(&self) -> Option<&molframe_xtal::SymmetrySet> {
        molframe_xtal::SymmetryExt::symmetry_set(self.engine())
    }

    fn collect_crystal_neighbors(
        &self,
        cutoff: f64,
        context: &ExecutionContext,
    ) -> Result<Vec<molframe_xtal::CrystalNeighbor>, Diagnostic> {
        molframe_xtal::SymmetryExt::collect_crystal_neighbors(self.engine(), cutoff, context)
    }

    fn collect_crystal_neighbors_in_model(
        &self,
        model: ModelIndex,
        cutoff: f64,
        context: &ExecutionContext,
    ) -> Result<Vec<molframe_xtal::CrystalNeighbor>, Diagnostic> {
        molframe_xtal::SymmetryExt::collect_crystal_neighbors_in_model(
            self.engine(),
            model,
            cutoff,
            context,
        )
    }
}

#[cfg(feature = "modelcif")]
impl molframe_modelcif::ModelCifExt for Structure {
    fn model_cif(&self) -> Option<&molframe_modelcif::ModelCif> {
        molframe_modelcif::ModelCifExt::model_cif(self.engine())
    }
}
