use crate::dataset::DatasetFilter;
use std::collections::BTreeMap;
use std::path::PathBuf;

use super::*;

#[test]
fn manifest_filter_and_lazy_loader_do_not_open_other_entries() {
    let dataset = Dataset::new(vec![
        entry("a", "AAAA", "x", "2020-01-01", 1.5),
        entry("b", "CCCC", "y", "2021-01-01", 3.0),
    ])
    .expect("valid dataset");
    let selected = dataset
        .filter(&DatasetFilter {
            resolution_below: Some(2.0),
            method: Some("X-RAY".into()),
            ..DatasetFilter::default()
        })
        .expect("valid filter");
    assert_eq!(selected.len(), 1);
    let loaded = selected
        .load_with(0, |item| Ok::<_, std::io::Error>(item.id.clone()))
        .expect("lazy load");
    assert_eq!(loaded.as_ref(), "a");
}

fn entry(id: &str, sequence: &str, cluster: &str, date: &str, resolution: f32) -> ManifestEntry {
    ManifestEntry {
        id: id.into(),
        path: PathBuf::from(format!("{id}.cif")),
        atom_count: 10,
        resolution: Some(resolution),
        method: Some("X-RAY".into()),
        deposition_date: Some(date.into()),
        sequence: Some(sequence.into()),
        structure_cluster: Some(cluster.into()),
        tags: vec!["test".into()],
        statistics: BTreeMap::new(),
    }
}
