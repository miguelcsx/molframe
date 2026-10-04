use super::*;

#[test]
fn every_name_round_trips_and_underscores_are_read_as_hyphens() {
    for name in NodeFeature::NAMES {
        assert_eq!(
            name.parse::<NodeFeature>().map(NodeFeature::name),
            Ok(*name)
        );
    }
    for name in EdgeFeature::NAMES {
        assert_eq!(
            name.parse::<EdgeFeature>().map(EdgeFeature::name),
            Ok(*name)
        );
    }
    assert_eq!("formal_charge".parse(), Ok(NodeFeature::FormalCharge));
    assert_eq!("residues".parse(), Ok(NodeLevel::Residues));
    assert_eq!("symmetric".parse(), Ok(EdgeDirection::Symmetric));
    let refused = "molecules".parse::<NodeLevel>().expect_err("not a level");
    assert_eq!(refused.field, "nodes");
}

#[test]
fn an_edge_algorithm_takes_exactly_its_own_parameter() {
    assert_eq!(
        EdgeKind::from_parts("bonds", None, None),
        Ok(EdgeKind::Bonds)
    );
    assert_eq!(
        EdgeKind::from_parts("contacts", Some(4.5), None),
        Ok(EdgeKind::Contacts { cutoff: 4.5 })
    );
    assert_eq!(
        EdgeKind::from_parts("k_nearest", None, Some(8)),
        Ok(EdgeKind::KNearest { neighbors: 8 })
    );
    assert!(EdgeKind::from_parts("contacts", None, None).is_err());
    assert!(EdgeKind::from_parts("radius", Some(5.0), Some(3)).is_err());
    assert!(EdgeKind::from_parts("bonds", Some(2.0), None).is_err());
    assert!(EdgeKind::from_parts("delaunay", None, None).is_err());
}
