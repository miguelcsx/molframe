use super::*;

#[test]
fn bond_topology_payload_kind_round_trips_without_an_adapter_default() {
    let python = PyPayloadKind::from(molframe::PayloadKind::BondTopology);

    assert_eq!(python, PyPayloadKind::BondTopology);
    assert_eq!(
        molframe::PayloadKind::from(python),
        molframe::PayloadKind::BondTopology
    );
}
