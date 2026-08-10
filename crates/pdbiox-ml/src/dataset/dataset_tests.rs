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

#[test]
fn sequence_identity_keeps_connected_sequences_together() {
    let dataset = Dataset::new(vec![
        entry("a", "AAAA", "x", "2020-01-01", 1.0),
        entry("b", "AAAT", "y", "2020-02-01", 1.0),
        entry("c", "CCCC", "z", "2020-03-01", 1.0),
        entry("d", "GGGG", "w", "2020-04-01", 1.0),
    ])
    .expect("valid dataset");
    let split = dataset
        .split(&SplitOptions {
            strategy: SplitStrategy::SequenceIdentity { threshold: 0.75 },
            ratios: SplitRatios::new(0.5, 0.25, 0.25).expect("valid ratios"),
        })
        .expect("sequence split");
    let partition = |id: &str| {
        [&split.train, &split.validation, &split.test]
            .iter()
            .position(|part| part.entries().any(|entry| entry.id.as_ref() == id))
    };
    assert_eq!(partition("a"), partition("b"));
    assert!(split.warnings.is_empty());
}

#[test]
fn structural_clusters_are_atomic_and_random_split_warns() {
    let dataset = Dataset::new(vec![
        entry("a", "AAAA", "same", "2020-01-01", 1.0),
        entry("b", "CCCC", "same", "2020-02-01", 1.0),
        entry("c", "GGGG", "other", "2020-03-01", 1.0),
    ])
    .expect("valid dataset");
    let ratios = SplitRatios::new(0.5, 0.25, 0.25).expect("valid ratios");
    let clustered = dataset
        .split(&SplitOptions {
            strategy: SplitStrategy::StructuralCluster,
            ratios,
        })
        .expect("cluster split");
    let together = [&clustered.train, &clustered.validation, &clustered.test]
        .iter()
        .any(|part| {
            let ids = part
                .entries()
                .map(|entry| entry.id.as_ref())
                .collect::<Vec<_>>();
            ids.contains(&"a") && ids.contains(&"b")
        });
    assert!(together);
    let random = dataset
        .split(&SplitOptions {
            strategy: SplitStrategy::Random { seed: 17 },
            ratios,
        })
        .expect("random split");
    assert_eq!(random.warnings, vec![DatasetWarning::RandomSplitMayLeak]);
}

#[test]
fn temporal_split_is_chronological() {
    let dataset = Dataset::new(vec![
        entry("new", "AAAA", "x", "2022-01-01", 1.0),
        entry("old", "CCCC", "y", "2020-01-01", 1.0),
        entry("middle", "GGGG", "z", "2021-01-01", 1.0),
        entry("latest", "TTTT", "w", "2023-01-01", 1.0),
    ])
    .expect("valid dataset");
    let split = dataset
        .split(&SplitOptions {
            strategy: SplitStrategy::Temporal,
            ratios: SplitRatios::new(0.5, 0.25, 0.25).expect("valid ratios"),
        })
        .expect("temporal split");
    assert_eq!(
        split
            .train
            .entries()
            .map(|entry| entry.id.as_ref())
            .collect::<Vec<_>>(),
        vec!["old", "middle"]
    );
    assert_eq!(
        split.test.entries().next().map(|entry| entry.id.as_ref()),
        Some("latest")
    );
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
