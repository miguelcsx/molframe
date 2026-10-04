use super::*;

fn close(left: f64, right: f64) -> bool {
    (left - right).abs() < 1e-12
}

#[test]
fn set_overlap_is_one_minus_the_jaccard_index() {
    let metric = SetOverlap::new(|set: &BTreeSet<u32>| set.clone());
    let a: BTreeSet<u32> = [1, 2, 3].into();
    let b: BTreeSet<u32> = [2, 3, 4].into();
    assert!(close(metric.distance(&a, &b), 0.5));
    assert!(close(metric.distance(&a, &a), 0.0));
    assert!(close(
        metric.distance(&BTreeSet::new(), &BTreeSet::new()),
        0.0
    ));
    assert!(close(metric.distance(&a, &[9].into()), 1.0));
}

#[test]
fn scalar_error_is_absolute_or_relative_and_refuses_to_compare_non_finite_values() {
    let absolute = ScalarError::new(|value: &f64| *value, ScalarMode::Absolute);
    let relative = ScalarError::new(|value: &f64| *value, ScalarMode::Relative);
    assert!(close(absolute.distance(&3.0, &5.0), 2.0));
    assert!(close(relative.distance(&3.0, &5.0), 0.4));
    assert!(close(relative.distance(&0.0, &0.0), 0.0));
    assert!(close(relative.distance(&-4.0, &4.0), 2.0));
    assert!(close(absolute.distance(&f64::NAN, &f64::NAN), 0.0));
    assert!(absolute.distance(&f64::NAN, &1.0).is_infinite());
    assert!(relative.distance(&f64::INFINITY, &1.0).is_infinite());
}

#[test]
fn vector_difference_is_rms_or_correlation_distance() {
    let rms = VectorDifference::new(|value: &Vec<f64>| value.clone(), VectorMode::Rms);
    let correlation =
        VectorDifference::new(|value: &Vec<f64>| value.clone(), VectorMode::Correlation);
    let a = vec![1.0, 2.0, 3.0];
    assert!(close(rms.distance(&a, &vec![2.0, 3.0, 4.0]), 1.0));
    // A shifted and scaled copy has the same shape.
    assert!(close(correlation.distance(&a, &vec![3.0, 5.0, 7.0]), 0.0));
    // An inverted profile is as far as standardised vectors go.
    assert!(close(correlation.distance(&a, &vec![3.0, 2.0, 1.0]), 2.0));
    assert!(close(
        correlation.distance(&vec![1.0; 3], &a),
        std::f64::consts::SQRT_2
    ));
    assert!(close(
        correlation.distance(&vec![1.0; 3], &vec![5.0; 3]),
        0.0
    ));
    assert!(rms.distance(&a, &vec![1.0]).is_infinite());
    assert!(rms.distance(&a, &vec![1.0, f64::NAN, 3.0]).is_infinite());
}

#[test]
fn ranking_distance_counts_discordant_pairs_and_half_penalises_ties() {
    let metric = RankingDistance::new(|ranking: &Vec<char>| ranking.clone());
    let abc = vec!['a', 'b', 'c'];
    assert!(close(metric.distance(&abc, &abc), 0.0));
    assert!(close(metric.distance(&abc, &vec!['c', 'b', 'a']), 1.0));
    // One adjacent swap flips one of three pairs.
    assert!(close(
        metric.distance(&abc, &vec!['b', 'a', 'c']),
        1.0 / 3.0
    ));
    // `a` alone: b and c are tied below it, a pair the other ranking orders.
    assert!(close(metric.distance(&abc, &vec!['a']), 0.5 / 3.0));
    // `b` alone puts b above a, against the other ranking's order.
    assert!(close(metric.distance(&vec!['a', 'b'], &vec!['b']), 1.0));
    assert!(close(metric.distance(&Vec::new(), &Vec::new()), 0.0));
}

#[test]
fn graph_difference_separates_nodes_from_edges_and_ignores_edge_direction() {
    let mut first = Graph::new();
    first.add_edge(1, 2);
    first.add_edge(2, 3);
    let mut second = Graph::new();
    second.add_edge(2, 1);
    second.add_edge(3, 4);
    assert!(first.edges().contains(&(1, 2)));
    assert!(close(first.node_difference(&second), 1.0 - 3.0 / 4.0));
    assert!(close(first.edge_difference(&second), 1.0 - 1.0 / 3.0));

    let both = GraphDifference::new(|graph: &Graph<u32>| graph.clone(), GraphPart::Both);
    // Shared identities: nodes 1, 2, 3 and edge (1, 2) of seven in the union.
    assert!(close(both.distance(&first, &second), 1.0 - 4.0 / 7.0));
    let nodes = GraphDifference::new(|graph: &Graph<u32>| graph.clone(), GraphPart::Nodes);
    assert_eq!(nodes.name(), "node-difference");
    assert!(close(nodes.distance(&first, &first), 0.0));
}

#[test]
fn a_conclusion_flip_is_zero_or_one_whatever_the_size_of_the_change() {
    let acceptable = CategoricalFlip::new(|score: &f64| *score >= 0.23);
    assert!(close(acceptable.distance(&0.64, &0.61), 0.0));
    assert!(close(acceptable.distance(&0.24, &0.22), 1.0));
    assert_eq!(acceptable.name(), "conclusion-flip");
}
