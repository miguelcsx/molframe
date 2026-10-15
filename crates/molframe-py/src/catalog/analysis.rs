//! Catalogued operations of the `analysis` namespace.

use super::{Capability, Domain};

const D: Domain = Domain("analysis");

pub(super) const ROWS: &[Capability] = &[
    D.row("hydrogen_bond_candidates", "structure,start,end,options", "Table", "cpu,spatial,memory"),
    D.row("heavy_atom_hydrogen_bonds", "structure,max_distance,min_distance,min_acceptor_angle,backend", "Table", "cpu,spatial,memory"),
    D.row("dssp", "structure", "Table", "cpu,spatial,memory"),
    D.row("contact_potential", "structure,charges,grid", "ScalarGrid", "cpu,spatial,memory"),
    D.row("atom_contacts", "structure|selection", "ContactTable", "cpu,spatial,memory").policy().workflow(),
    D.row("contacts", "structure", "Analysis", "cpu,spatial,memory").policy(),
    D.row("hydrogen_bonds", "structure", "Analysis", "cpu,spatial,memory").policy(),
    D.row("salt_bridges", "structure", "Analysis", "cpu,spatial,memory").policy(),
    D.row("contact_map", "structure,cutoff,min_separation,backend", "Analysis<Table>", "cpu,spatial,memory").policy(),
    D.row("pi_stacking", "structure,max_centre_distance,max_parallel_angle,min_perpendicular_angle", "Analysis<Table>", "cpu,spatial,memory").policy(),
    D.row("cation_pi", "structure,max_distance,max_face_angle", "Analysis<Table>", "cpu,spatial,memory").policy(),
    D.row("water_bridges", "structure,max_distance,min_angle,backend", "Analysis<Table>", "cpu,spatial,memory").policy(),
    D.row("chain_interface", "structure,first_chain,second_chain,cutoff,backend", "Analysis<Table>", "cpu,spatial,memory").policy(),
    D.row("half_sphere_exposure", "structure,radius,backend", "Analysis<Table>", "cpu,spatial,memory").policy(),
    D.row("nucleic_torsions", "structure", "Analysis<Table>", "cpu,spatial,memory").policy(),
    D.row("native_contacts", "reference,target,cutoff,tolerance,match,mapping,backend", "Analysis<dict>", "cpu,spatial,memory").policy(),
    D.row("leaflets", "structure,sites,connection_distance,backend", "Analysis<Table>", "cpu,spatial,memory").policy(),
    D.row("radial_distribution", "structure,first,second,minimum_distance,maximum_distance,bins,volume,backend", "Analysis<Table>", "cpu,spatial,memory").policy(),
    D.row("coordination_numbers", "structure,first,second,minimum_distance,maximum_distance,backend", "Analysis<Table>", "cpu,spatial,memory").policy(),
    D.row("linear_density", "structure,axis,minimum,maximum,bins,weights", "Analysis<Table>", "cpu,spatial,memory").policy(),
    D.row("density_map", "structure,origin,spacing,shape,weights", "Analysis<DensityGrid>", "cpu,spatial,memory").policy(),
    D.row("pore_profile", "structure,axis_origin,axis_direction,start,end,samples,search_radius,grid_spacing,probe_radius,radii,memory_limit_bytes", "Analysis<Table>", "cpu,spatial,memory").policy(),
    D.row("surface_contacts", "structure,tolerance,probe,surface_density,minimum_area,radii,backend", "Analysis<Table>", "cpu,spatial,memory").policy(),
    D.row("sasa", "structure,probe,method,n_points,hetero,metals", "Analysis<ndarray>", "cpu,memory").policy(),
    D.row("residue_sasa", "structure,probe,method,n_points,hetero,metals", "Analysis<ResidueAreas>", "cpu,memory").policy(),
    D.row("contacts_by_definition", "structure,backend,surface_tolerance,surface_density,surface_minimum_area", "Analysis<Table>", "cpu,spatial,memory").policy(),
];
