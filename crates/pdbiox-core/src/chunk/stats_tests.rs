use super::*;

#[test]
fn an_element_set_answers_membership_without_a_search() {
    let mask: ElementMask = [Element::CARBON, Element::NITROGEN, Element::ZINC]
        .into_iter()
        .collect();
    assert!(mask.contains(Element::ZINC));
    assert!(!mask.contains(Element::IRON));
    assert_eq!(mask.len(), 3);
    assert!(!mask.is_empty());
}

#[test]
fn a_chunk_without_zinc_is_rejected_before_any_row_is_read() {
    let stats = AtomChunkStats {
        elements: [Element::CARBON, Element::OXYGEN].into_iter().collect(),
        ..AtomChunkStats::default()
    };
    assert!(stats.excludes_element(Element::ZINC));
    assert!(!stats.excludes_element(Element::CARBON));
}

#[test]
fn a_chunk_too_far_from_a_region_is_rejected_by_comparing_boxes() {
    let stats = AtomChunkStats {
        bounds: Aabb {
            min: [0.0; 3],
            max: [1.0; 3],
        },
        ..AtomChunkStats::default()
    };
    let far = Aabb {
        min: [100.0; 3],
        max: [101.0; 3],
    };
    assert!(stats.excludes_region(&far, 5.0));
    assert!(!stats.excludes_region(&far, 200.0));
}

#[test]
fn a_residue_range_outside_the_chunk_is_rejected_from_the_summary() {
    let stats = AtomChunkStats {
        residue_min: 10,
        residue_max: 20,
        ..AtomChunkStats::default()
    };
    assert!(stats.excludes_residues(0, 9));
    assert!(stats.excludes_residues(21, 100));
    assert!(!stats.excludes_residues(15, 15));
    assert!(!stats.excludes_residues(0, 100));
}

#[test]
fn extremes_ignore_values_that_are_not_numbers() {
    let mut extremes = Extremes::default();
    extremes.observe(f32::NAN);
    assert_eq!(extremes.min(), None);
    extremes.observe(30.0);
    extremes.observe(10.0);
    extremes.observe(f32::INFINITY);
    assert_eq!(extremes.min(), Some(10.0));
    assert_eq!(extremes.max(), Some(30.0));
}

#[test]
fn a_temperature_factor_filter_skips_a_chunk_whose_range_cannot_satisfy_it() {
    let mut extremes = Extremes::default();
    extremes.observe(5.0);
    extremes.observe(20.0);
    assert!(extremes.excludes_above(30.0));
    assert!(!extremes.excludes_above(15.0));
    assert!(extremes.excludes_below(1.0));
    assert!(!extremes.excludes_below(10.0));
}
