use super::*;

#[test]
fn none_policy_variants_have_python_safe_constructors() {
    assert_eq!(PyEquivalencePolicy::none(), PyEquivalencePolicy::None);
    assert_eq!(PySymmetryPolicy::none(), PySymmetryPolicy::None);
    assert_eq!(PyPeriodicPolicy::none(), PyPeriodicPolicy::None);
}

#[test]
fn policy_space_projection_preserves_native_cartesian_cost() {
    let dimension = PyPolicyDimension::hydrogens(vec![
        PyHydrogenPolicy::ExplicitOnly,
        PyHydrogenPolicy::Exclude,
    ]);
    let space = PyPolicySpace::new(None).vary(&dimension);
    assert_eq!(space.cost().expect("valid dimension"), 2);
    let plan = space.plan().expect("bounded plan");
    assert_eq!(plan.cost(), 2);
    assert_eq!(plan.fields(), vec![PyPolicyField::Hydrogens]);
}
