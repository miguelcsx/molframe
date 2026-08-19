use super::*;

#[test]
fn particle_binding_preserves_absence() {
    let particle = PyDmsParticle {
        value: DmsParticle {
            atomic_number: Some(6),
            name: Some("C1".into()),
            ..DmsParticle::default()
        },
    };
    assert_eq!(particle.value.atomic_number, Some(6));
    assert_eq!(particle.value.name.as_deref(), Some("C1"));
    assert_eq!(particle.value.charge, None);
}
