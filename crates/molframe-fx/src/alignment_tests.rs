use super::{AlignmentKind, align_intrinsic, align_with_transform};
use crate::MappedMotif;
use molframe_geom::Rigid;
use std::collections::BTreeMap;

#[test]
fn alignment_origin_is_explicit() {
    let intrinsic = align_intrinsic(MappedMotif {
        components: BTreeMap::default(),
        atoms: BTreeMap::default(),
    });
    assert_eq!(intrinsic.kind, AlignmentKind::NotRequired);

    let moved = align_with_transform(intrinsic.mapping, Rigid::translation([1.0, 2.0, 3.0]));
    assert_eq!(moved.kind, AlignmentKind::CallerSupplied);
    assert!(
        moved
            .transform
            .translation
            .iter()
            .zip([1.0, 2.0, 3.0])
            .all(|(observed, expected)| (observed - expected).abs() < f64::EPSILON)
    );
}
