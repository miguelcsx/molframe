use super::*;
use crate::PointMatch;

fn mapping() -> PointMapping {
    PointMapping::new(
        [
            PointMatch {
                reference: 0,
                model: 0,
            },
            PointMatch {
                reference: 1,
                model: 1,
            },
        ],
        2,
        2,
    )
    .expect("explicit mapping is valid")
}

#[test]
fn typed_regions_share_the_explicit_mapping_measurement() {
    let reference = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let model = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
    let interface = interface_rmsd(
        &reference,
        &model,
        &mapping(),
        ComparisonAlignment::NotRequired,
    );
    let pocket = pocket_rmsd(
        &reference,
        &model,
        &mapping(),
        ComparisonAlignment::NotRequired,
    );
    assert_eq!(interface.0.to_bits(), pocket.0.to_bits());
    assert!((interface.0 - 2.0_f64.sqrt()).abs() < 1.0e-12);
}
