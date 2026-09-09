use super::{SurfaceContactOptions, surface_contacts};
use pdbiox_core::ExecutionContext;
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_spatial::SpatialBackend;

fn structure(xs: &[f32]) -> pdbiox_core::Structure {
    let mut source = "data_s\nloop_\n_atom_site.group_PDB\n_atom_site.id\n\
_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
        .to_owned();
    for (index, x) in xs.iter().enumerate() {
        use std::fmt::Write as _;
        let _ = writeln!(source, "ATOM {} C C LIG A {} {x} 0 0", index + 1, index + 1);
    }
    let input = InputBuffer::from_bytes(source.into_bytes());
    pdbiox_cif::read(&input, &ReadOptions::new())
        .unwrap_or_else(|findings| panic!("fixture failed: {findings:?}"))
        .0
}

#[test]
fn exposed_facing_patches_support_a_contact() {
    let structure = structure(&[0.0, 3.2]);
    let contacts = surface_contacts(
        &structure,
        &[1.7, 1.7],
        SurfaceContactOptions {
            tolerance: 0.5,
            probe: 1.4,
            surface_density: 20.0,
            minimum_area: 0.2,
            backend: SpatialBackend::BruteForce,
        },
        &ExecutionContext::default(),
    )
    .unwrap_or_else(|error| panic!("surface failed: {error}"));
    assert_eq!(contacts.len(), 1);
}

#[test]
fn an_occluding_atom_removes_the_facing_surface_support() {
    let structure = structure(&[0.0, 3.2, 1.6]);
    let contacts = surface_contacts(
        &structure,
        &[1.7, 1.7, 2.2],
        SurfaceContactOptions {
            tolerance: 0.5,
            probe: 1.4,
            surface_density: 20.0,
            minimum_area: 0.2,
            backend: SpatialBackend::BruteForce,
        },
        &ExecutionContext::default(),
    )
    .unwrap_or_else(|error| panic!("surface failed: {error}"));
    assert!(
        !contacts
            .iter()
            .any(|contact| contact.first.get() == 0 && contact.second.get() == 1)
    );
}
