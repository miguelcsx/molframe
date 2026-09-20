use super::analysis::AnalysisExt;
use crate::ExecutionContext;
use crate::structure::{Selection, Structure};
use molframe_chem::ComponentProvider;
use molframe_spatial::{SpatialBackend, SpatialError};

impl AnalysisExt for Structure {
    fn atom_contacts(&self, cutoff: f32) -> Result<molframe_analysis::ContactTable, SpatialError> {
        self.atom_contacts_with(cutoff, SpatialBackend::Auto, &ExecutionContext::default())
    }

    fn atom_contacts_with(
        &self,
        cutoff: f32,
        backend: SpatialBackend,
        context: &ExecutionContext,
    ) -> Result<molframe_analysis::ContactTable, SpatialError> {
        molframe_analysis::atom_contacts(self.engine(), cutoff, backend, context)
    }

    fn atom_contacts_between(
        &self,
        left: &Selection,
        right: &Selection,
        cutoff: f32,
    ) -> Result<molframe_analysis::ContactTable, SpatialError> {
        self.atom_contacts_between_with(
            left,
            right,
            cutoff,
            SpatialBackend::Auto,
            &ExecutionContext::default(),
        )
    }

    fn atom_contacts_between_with(
        &self,
        left: &Selection,
        right: &Selection,
        cutoff: f32,
        backend: SpatialBackend,
        context: &ExecutionContext,
    ) -> Result<molframe_analysis::ContactTable, SpatialError> {
        molframe_analysis::atom_contacts_between(
            self.engine(),
            left.atom_selection(),
            right.atom_selection(),
            cutoff,
            backend,
            context,
        )
    }

    fn visit_atom_contacts<F: FnMut(molframe_analysis::Contact)>(
        &self,
        cutoff: f32,
        emit: F,
    ) -> Result<(), SpatialError> {
        self.visit_atom_contacts_with(
            cutoff,
            SpatialBackend::Auto,
            &ExecutionContext::default(),
            emit,
        )
    }

    fn visit_atom_contacts_with<F: FnMut(molframe_analysis::Contact)>(
        &self,
        cutoff: f32,
        backend: SpatialBackend,
        context: &ExecutionContext,
        emit: F,
    ) -> Result<(), SpatialError> {
        molframe_analysis::visit_atom_contacts(self.engine(), cutoff, backend, context, emit)
    }

    fn residue_contact_map(
        &self,
        cutoff: f32,
        min_separation: u32,
    ) -> Result<molframe_analysis::ContactMap, SpatialError> {
        self.residue_contact_map_with(
            cutoff,
            min_separation,
            SpatialBackend::Auto,
            &ExecutionContext::default(),
        )
    }

    fn residue_contact_map_with(
        &self,
        cutoff: f32,
        min_separation: u32,
        backend: SpatialBackend,
        context: &ExecutionContext,
    ) -> Result<molframe_analysis::ContactMap, SpatialError> {
        molframe_analysis::residue_contact_map(
            self.engine(),
            cutoff,
            min_separation,
            backend,
            context,
        )
    }

    fn hydrogen_bonds(
        &self,
        options: molframe_analysis::HydrogenBondOptions,
    ) -> Result<molframe_analysis::HydrogenBondTable, molframe_analysis::HydrogenBondError> {
        self.hydrogen_bonds_with(options, &ExecutionContext::default())
    }

    fn hydrogen_bonds_with(
        &self,
        options: molframe_analysis::HydrogenBondOptions,
        context: &ExecutionContext,
    ) -> Result<molframe_analysis::HydrogenBondTable, molframe_analysis::HydrogenBondError> {
        molframe_analysis::hydrogen_bonds(self.engine(), options, context)
    }

    fn water_bridges(
        &self,
        options: molframe_analysis::WaterBridgeOptions,
    ) -> Result<molframe_analysis::WaterBridgeTable, molframe_analysis::HydrogenBondError> {
        self.water_bridges_with(options, &ExecutionContext::default())
    }

    fn water_bridges_with(
        &self,
        options: molframe_analysis::WaterBridgeOptions,
        context: &ExecutionContext,
    ) -> Result<molframe_analysis::WaterBridgeTable, molframe_analysis::HydrogenBondError> {
        molframe_analysis::water_bridges(self.engine(), options, context)
    }

    fn salt_bridges(
        &self,
        max_distance: f32,
    ) -> Result<molframe_analysis::SaltBridgeTable, SpatialError> {
        self.salt_bridges_with(
            max_distance,
            SpatialBackend::Auto,
            &ExecutionContext::default(),
        )
    }

    fn salt_bridges_with(
        &self,
        max_distance: f32,
        backend: SpatialBackend,
        context: &ExecutionContext,
    ) -> Result<molframe_analysis::SaltBridgeTable, SpatialError> {
        molframe_analysis::salt_bridges(self.engine(), max_distance, backend, context)
    }

    fn pi_stacking(
        &self,
        options: molframe_analysis::PiStackingOptions,
    ) -> Result<molframe_analysis::PiStackingTable, molframe_analysis::PiStackingError> {
        molframe_analysis::pi_stacking(self.engine(), options)
    }

    fn cation_pi(
        &self,
        options: molframe_analysis::CationPiOptions,
    ) -> Result<molframe_analysis::CationPiTable, molframe_analysis::CationPiError> {
        molframe_analysis::cation_pi(self.engine(), options)
    }

    fn base_pairs(
        &self,
        provider: &dyn ComponentProvider,
        options: molframe_analysis::BasePairOptions,
    ) -> Result<Vec<molframe_analysis::BasePair>, molframe_analysis::BasePairError> {
        molframe_analysis::base_pairs(
            self.engine(),
            provider,
            options,
            &ExecutionContext::default(),
        )
    }

    fn base_pairs_with(
        &self,
        provider: &dyn ComponentProvider,
        options: molframe_analysis::BasePairOptions,
        context: &ExecutionContext,
    ) -> Result<Vec<molframe_analysis::BasePair>, molframe_analysis::BasePairError> {
        molframe_analysis::base_pairs(self.engine(), provider, options, context)
    }

    fn secondary_structure(
        &self,
        options: &molframe_analysis::DsspOptions,
    ) -> Result<molframe_analysis::SseTable, molframe_analysis::DsspError> {
        molframe_analysis::secondary_structure(self.engine(), options)
    }

    fn chain_interface(
        &self,
        first: &str,
        second: &str,
        cutoff: f32,
    ) -> Result<Vec<molframe_core::index::ResidueIndex>, SpatialError> {
        self.chain_interface_with(
            first,
            second,
            cutoff,
            SpatialBackend::Auto,
            &ExecutionContext::default(),
        )
    }

    fn chain_interface_with(
        &self,
        first: &str,
        second: &str,
        cutoff: f32,
        backend: SpatialBackend,
        context: &ExecutionContext,
    ) -> Result<Vec<molframe_core::index::ResidueIndex>, SpatialError> {
        molframe_analysis::chain_interface(self.engine(), first, second, cutoff, backend, context)
    }

    fn native_contact_fraction(
        &self,
        target: &Structure,
        cutoff: f32,
        tolerance: f32,
    ) -> Result<molframe_analysis::NativeContacts, molframe_analysis::NativeError> {
        self.native_contact_fraction_with(
            target,
            cutoff,
            tolerance,
            SpatialBackend::Auto,
            &ExecutionContext::default(),
        )
    }

    fn native_contact_fraction_with(
        &self,
        target: &Structure,
        cutoff: f32,
        tolerance: f32,
        backend: SpatialBackend,
        context: &ExecutionContext,
    ) -> Result<molframe_analysis::NativeContacts, molframe_analysis::NativeError> {
        molframe_analysis::native_contact_fraction(
            self.engine(),
            target.engine(),
            cutoff,
            tolerance,
            backend,
            context,
        )
    }

    fn half_sphere_exposure(
        &self,
        radius: f32,
    ) -> Result<Vec<molframe_analysis::HalfSphereExposure>, molframe_analysis::HseError> {
        self.half_sphere_exposure_with(radius, SpatialBackend::Auto, &ExecutionContext::default())
    }

    fn half_sphere_exposure_with(
        &self,
        radius: f32,
        backend: SpatialBackend,
        context: &ExecutionContext,
    ) -> Result<Vec<molframe_analysis::HalfSphereExposure>, molframe_analysis::HseError> {
        molframe_analysis::half_sphere_exposure(self.engine(), radius, backend, context)
    }

    fn nucleic_torsions(
        &self,
    ) -> Result<Vec<molframe_analysis::NucleicTorsions>, molframe_analysis::NucleicTorsionError>
    {
        molframe_analysis::nucleic_torsions(self.engine())
    }

    fn surface_contacts(
        &self,
        radii: &[f32],
        options: molframe_analysis::SurfaceContactOptions,
    ) -> Result<Vec<molframe_analysis::Contact>, molframe_surface::SasaError> {
        molframe_analysis::surface_contacts(
            self.engine(),
            radii,
            options,
            &ExecutionContext::default(),
        )
    }

    fn surface_contacts_with(
        &self,
        radii: &[f32],
        options: molframe_analysis::SurfaceContactOptions,
        context: &ExecutionContext,
    ) -> Result<Vec<molframe_analysis::Contact>, molframe_surface::SasaError> {
        molframe_analysis::surface_contacts(self.engine(), radii, options, context)
    }
}
