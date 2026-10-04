//! The comparisons an audit can make, over what Python hands back.
//!
//! A run's value is reduced to the datum the chosen metric compares; this module
//! owns that reduction and the distance, over the metrics the audit crate defines.

use molframe::audit::{
    CategoricalFlip, Graph, GraphDifference, GraphPart, OutcomeMetric, RankingDistance,
    ScalarError, ScalarMode, SetOverlap, VectorDifference, VectorMode,
};
use pyo3::prelude::*;
use std::collections::BTreeSet;

/// What the answer of one run was reduced to, for comparison.
pub(crate) enum Projected {
    Set(BTreeSet<String>),
    Number(f64),
    Vector(Vec<f64>),
    Ranking(Vec<String>),
    Category(String),
    Graph(Graph<String>),
}

/// How two projected answers are compared.
#[derive(Clone, Copy)]
pub(crate) enum MetricKind {
    Set,
    Absolute,
    Relative,
    Rms,
    Correlation,
    Ranking,
    Graph(GraphPart),
    Flip,
}

const METRICS: &str = "set, absolute, relative, rms, correlation, ranking, graph, graph-nodes, \
                       graph-edges or flip";

impl MetricKind {
    pub(crate) fn parse(name: &str) -> PyResult<Self> {
        Ok(match name.replace('_', "-").as_str() {
            "set" => Self::Set,
            "absolute" => Self::Absolute,
            "relative" => Self::Relative,
            "rms" => Self::Rms,
            "correlation" => Self::Correlation,
            "ranking" => Self::Ranking,
            "graph" => Self::Graph(GraphPart::Both),
            "graph-nodes" => Self::Graph(GraphPart::Nodes),
            "graph-edges" => Self::Graph(GraphPart::Edges),
            "flip" => Self::Flip,
            other => {
                return Err(crate::error::value(format!(
                    "metric {other:?} is not one of {METRICS}"
                )));
            }
        })
    }

    /// Reduces a Python datum to what this metric compares.
    pub(crate) fn project(self, value: &Bound<'_, PyAny>) -> PyResult<Projected> {
        let items = |datum: &Bound<'_, PyAny>| -> PyResult<Vec<String>> {
            datum
                .try_iter()?
                .map(|item| {
                    item.and_then(|item| item.repr())
                        .map(|repr| repr.to_string())
                })
                .collect()
        };
        Ok(match self {
            Self::Set => Projected::Set(items(value)?.into_iter().collect()),
            Self::Absolute | Self::Relative => Projected::Number(value.extract()?),
            Self::Rms | Self::Correlation => Projected::Vector(value.extract()?),
            Self::Ranking => Projected::Ranking(items(value)?),
            Self::Flip => Projected::Category(value.repr()?.to_string()),
            Self::Graph(_) => {
                let (nodes, edges): (Bound<'_, PyAny>, Bound<'_, PyAny>) = value.extract()?;
                let mut graph = Graph::new();
                for node in items(&nodes)? {
                    graph.add_node(node);
                }
                for edge in edges.try_iter()? {
                    let (first, second): (Bound<'_, PyAny>, Bound<'_, PyAny>) = edge?.extract()?;
                    graph.add_edge(first.repr()?.to_string(), second.repr()?.to_string());
                }
                Projected::Graph(graph)
            }
        })
    }
}

impl OutcomeMetric<Projected> for MetricKind {
    fn name(&self) -> &'static str {
        match self {
            Self::Set => "jaccard-distance",
            Self::Absolute => "absolute-error",
            Self::Relative => "relative-error",
            Self::Rms => "rms-difference",
            Self::Correlation => "correlation-distance",
            Self::Ranking => "kendall-tau-distance",
            Self::Graph(GraphPart::Both) => "graph-difference",
            Self::Graph(GraphPart::Nodes) => "node-difference",
            Self::Graph(GraphPart::Edges) => "edge-difference",
            Self::Flip => "conclusion-flip",
        }
    }

    fn distance(&self, first: &Projected, second: &Projected) -> f64 {
        use Projected as P;
        match (self, first, second) {
            (Self::Set, P::Set(a), P::Set(b)) => {
                SetOverlap::new(|set: &BTreeSet<String>| set.clone()).distance(a, b)
            }
            (Self::Absolute, P::Number(a), P::Number(b)) => {
                ScalarError::new(|value: &f64| *value, ScalarMode::Absolute).distance(a, b)
            }
            (Self::Relative, P::Number(a), P::Number(b)) => {
                ScalarError::new(|value: &f64| *value, ScalarMode::Relative).distance(a, b)
            }
            (Self::Rms, P::Vector(a), P::Vector(b)) => {
                VectorDifference::new(|value: &Vec<f64>| value.clone(), VectorMode::Rms)
                    .distance(a, b)
            }
            (Self::Correlation, P::Vector(a), P::Vector(b)) => {
                VectorDifference::new(|value: &Vec<f64>| value.clone(), VectorMode::Correlation)
                    .distance(a, b)
            }
            (Self::Ranking, P::Ranking(a), P::Ranking(b)) => {
                RankingDistance::new(|value: &Vec<String>| value.clone()).distance(a, b)
            }
            (Self::Graph(part), P::Graph(a), P::Graph(b)) => {
                GraphDifference::new(|graph: &Graph<String>| graph.clone(), *part).distance(a, b)
            }
            (Self::Flip, P::Category(a), P::Category(b)) => {
                CategoricalFlip::new(|category: &String| category.clone()).distance(a, b)
            }
            // Every run of one audit is projected by the same metric, so this is unreachable;
            // an incomparable pair is as far apart as anything can be.
            _ => f64::INFINITY,
        }
    }
}
