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
    assert_eq!(particle.atomic_number(), Some(6));
    assert_eq!(particle.name(), Some("C1"));
    assert_eq!(particle.charge(), None);
}
