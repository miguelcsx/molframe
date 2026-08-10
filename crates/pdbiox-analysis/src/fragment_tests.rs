use super::{FragmentMappingError, FragmentReference, map_fragments};

fn reference(id: &str, y: f32) -> FragmentReference {
    FragmentReference {
        id: id.into(),
        coordinates: vec![[0.0, y, 0.0], [1.0, y, 0.0], [1.0, 1.0 + y, 0.0]],
    }
}

#[test]
fn rigidly_translated_window_maps_to_the_reference() {
    let trace = vec![
        Some([10.0, 3.0, 0.0]),
        Some([11.0, 3.0, 0.0]),
        Some([11.0, 4.0, 0.0]),
    ];
    let Ok(matches) = map_fragments(&trace, &[reference("turn", 0.0)], 1e-5) else {
        panic!("valid fragment mapping");
    };
    assert_eq!(matches.len(), 1);
    assert_eq!(matches[0].fragment_id.as_ref(), "turn");
}

#[test]
fn missing_trace_sites_skip_only_affected_windows() {
    let trace = vec![
        Some([0.0, 0.0, 0.0]),
        Some([1.0, 0.0, 0.0]),
        Some([1.0, 1.0, 0.0]),
        None,
    ];
    let Ok(matches) = map_fragments(&trace, &[reference("turn", 0.0)], 1e-5) else {
        panic!("valid partial trace");
    };
    assert_eq!(matches.len(), 1);
}

#[test]
fn fragment_library_is_explicit_and_unique() {
    let duplicate = [reference("same", 0.0), reference("same", 1.0)];
    assert!(map_fragments(&[], &duplicate, 1.0).is_err());
}

#[test]
fn degenerate_complete_window_is_an_error_not_an_omission() {
    let line = FragmentReference {
        id: "line".into(),
        coordinates: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
    };
    let trace = line
        .coordinates
        .iter()
        .copied()
        .map(Some)
        .collect::<Vec<_>>();
    assert!(matches!(
        map_fragments(&trace, &[line], 1.0),
        Err(FragmentMappingError::Superpose(_))
    ));
}
