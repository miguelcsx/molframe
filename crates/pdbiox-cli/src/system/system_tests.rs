use pdbiox::traj::DmsSystem;

#[test]
fn empty_system_has_aligned_default_arrays() {
    let system = DmsSystem::default();
    assert_eq!(
        system.topology.particles.len(),
        system.frame.positions.len()
    );
    assert_eq!(system.frame.positions.len(), system.frame.velocities.len());
}
