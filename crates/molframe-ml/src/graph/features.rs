//! Typed feature materialisation with explicit absence handling.

use super::edges::Edge;
use super::nodes::Nodes;
use super::{EdgeFeature, FeatureMatrix, GraphError, MissingFeaturePolicy, NodeFeature, NodeLevel};
use molframe_core::structure::AtomRef;
use molframe_core::{
    AtomAnnotation, BondOrder, FORMAL_CHARGE_ANNOTATION, PARTIAL_CHARGE_ANNOTATION, Structure,
};
use std::ops::Range;

use crate::numeric::{f64_to_f32, i64_to_f64, u32_to_f32};

const AROMATIC_BOND_ORDER: f32 = 1.5;
const SINGLE_BOND_ORDER: f32 = 1.0;
const DOUBLE_BOND_ORDER: f32 = 2.0;
const TRIPLE_BOND_ORDER: f32 = 3.0;
const QUADRUPLE_BOND_ORDER: f32 = 4.0;

pub(super) fn nodes(
    structure: &Structure,
    nodes: &Nodes,
    level: NodeLevel,
    features: &[NodeFeature],
    missing: MissingFeaturePolicy,
) -> Result<FeatureMatrix<NodeFeature>, GraphError> {
    let capacity = nodes
        .positions
        .len()
        .checked_mul(features.len())
        .ok_or(GraphError::ResourceLimit)?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(capacity)
        .map_err(|_| GraphError::ResourceLimit)?;
    for row in 0..nodes.positions.len() {
        for feature in features {
            let value = node_value(structure, nodes, level, row, *feature)?;
            values.push(required(value, node_label(*feature), row, missing)?);
        }
    }
    Ok(FeatureMatrix {
        values: values.into_boxed_slice(),
        rows: nodes.positions.len(),
        columns: features.len(),
        features: features.into(),
    })
}

pub(super) fn edges(
    edges: &[Edge],
    features: &[EdgeFeature],
    missing: MissingFeaturePolicy,
) -> Result<FeatureMatrix<EdgeFeature>, GraphError> {
    let capacity = edges
        .len()
        .checked_mul(features.len())
        .ok_or(GraphError::ResourceLimit)?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(capacity)
        .map_err(|_| GraphError::ResourceLimit)?;
    for (row, edge) in edges.iter().enumerate() {
        for feature in features {
            let value = match feature {
                EdgeFeature::Distance => Some(edge.distance_squared.sqrt()),
                EdgeFeature::BondOrder => edge.bond_order.and_then(bond_order),
            };
            values.push(required(value, edge_label(*feature), row, missing)?);
        }
    }
    Ok(FeatureMatrix {
        values: values.into_boxed_slice(),
        rows: edges.len(),
        columns: features.len(),
        features: features.into(),
    })
}

fn node_value(
    structure: &Structure,
    nodes: &Nodes,
    level: NodeLevel,
    row: usize,
    feature: NodeFeature,
) -> Result<Option<f32>, GraphError> {
    let Some(atoms) = nodes.atom_ranges.get(row).cloned() else {
        return Err(GraphError::IndexOverflow);
    };
    match (level, feature) {
        (NodeLevel::Atoms, NodeFeature::Element) => Ok(structure
            .data()
            .atom(molframe_core::AtomIndex::new(atoms.start))
            .and_then(AtomRef::element)
            .map(|element| f32::from(element.atomic_number()))),
        (NodeLevel::Atoms, NodeFeature::AtomCount)
        | (NodeLevel::Residues, NodeFeature::Element) => Err(GraphError::UnsupportedFeature {
            feature: node_label(feature),
        }),
        (_, NodeFeature::FormalCharge) => Ok(annotation_sum(
            structure,
            atoms,
            FORMAL_CHARGE_ANNOTATION,
            AnnotationKind::Integer,
        )),
        (_, NodeFeature::PartialCharge) => Ok(annotation_sum(
            structure,
            atoms,
            PARTIAL_CHARGE_ANNOTATION,
            AnnotationKind::Real,
        )),
        (_, NodeFeature::BFactor) => atom_mean(structure, atoms, Scalar::BFactor),
        (_, NodeFeature::Occupancy) => atom_mean(structure, atoms, Scalar::Occupancy),
        (NodeLevel::Residues, NodeFeature::AtomCount) => {
            let count = atoms
                .end
                .checked_sub(atoms.start)
                .ok_or(GraphError::IndexOverflow)?;
            Ok(Some(u32_to_f32(count)))
        }
        (_, NodeFeature::PositionX) => position(nodes, row, 0),
        (_, NodeFeature::PositionY) => position(nodes, row, 1),
        (_, NodeFeature::PositionZ) => position(nodes, row, 2),
    }
}

#[derive(Clone, Copy)]
enum AnnotationKind {
    Integer,
    Real,
}

fn annotation_sum(
    structure: &Structure,
    atoms: Range<u32>,
    name: &str,
    kind: AnnotationKind,
) -> Option<f32> {
    let annotation = structure.annotations().get(name)?;
    let mut sum = 0.0f64;
    for atom in atoms {
        let value = match (kind, annotation) {
            (AnnotationKind::Integer, AtomAnnotation::Integer(column)) => column
                .get(atom)
                .filter(|(_, presence)| presence.is_present())
                .map(|(value, _)| i64_to_f64(value)),
            (AnnotationKind::Real, AtomAnnotation::Real(column)) => column
                .get(atom)
                .filter(|(_, presence)| presence.is_present())
                .map(|(value, _)| value),
            _ => return None,
        };
        let value = value?;
        sum += value;
    }
    sum.is_finite().then_some(f64_to_f32(sum))
}

#[derive(Clone, Copy)]
enum Scalar {
    BFactor,
    Occupancy,
}

fn atom_mean(
    structure: &Structure,
    atoms: Range<u32>,
    scalar: Scalar,
) -> Result<Option<f32>, GraphError> {
    let mut sum = 0.0f64;
    let mut count = 0u32;
    for atom in atoms {
        let Some(atom) = structure.data().atom(molframe_core::AtomIndex::new(atom)) else {
            return Err(GraphError::IndexOverflow);
        };
        let value = match scalar {
            Scalar::BFactor => atom.b_factor(),
            Scalar::Occupancy => atom.occupancy(),
        };
        let Some(value) = value else {
            return Ok(None);
        };
        sum += f64::from(value);
        count += 1;
    }
    if count == 0 || !sum.is_finite() {
        Ok(None)
    } else {
        Ok(Some(f64_to_f32(sum / f64::from(count))))
    }
}

fn position(nodes: &Nodes, row: usize, axis: usize) -> Result<Option<f32>, GraphError> {
    if !nodes
        .position_valid
        .get(row)
        .copied()
        .is_some_and(|valid| valid)
    {
        return Ok(None);
    }
    nodes
        .positions
        .get(row)
        .and_then(|position| position.get(axis))
        .copied()
        .map(Some)
        .ok_or(GraphError::IndexOverflow)
}

fn required(
    value: Option<f32>,
    feature: &'static str,
    row: usize,
    policy: MissingFeaturePolicy,
) -> Result<f32, GraphError> {
    match value {
        Some(value) if value.is_finite() => Ok(value),
        Some(_) | None => match policy {
            MissingFeaturePolicy::Error => Err(GraphError::MissingFeature { feature, row }),
            MissingFeaturePolicy::Fill(value) if value.is_finite() => Ok(value),
            MissingFeaturePolicy::Fill(_) => Err(GraphError::InvalidParameter),
        },
    }
}

fn bond_order(order: BondOrder) -> Option<f32> {
    match order {
        BondOrder::Single | BondOrder::Polymeric => Some(SINGLE_BOND_ORDER),
        BondOrder::Double => Some(DOUBLE_BOND_ORDER),
        BondOrder::Triple => Some(TRIPLE_BOND_ORDER),
        BondOrder::Quadruple => Some(QUADRUPLE_BOND_ORDER),
        BondOrder::Aromatic => Some(AROMATIC_BOND_ORDER),
        BondOrder::Unknown => None,
    }
}

pub(super) const fn node_label(feature: NodeFeature) -> &'static str {
    match feature {
        NodeFeature::Element => "element",
        NodeFeature::FormalCharge => "formal_charge",
        NodeFeature::PartialCharge => "partial_charge",
        NodeFeature::BFactor => "b_factor",
        NodeFeature::Occupancy => "occupancy",
        NodeFeature::AtomCount => "atom_count",
        NodeFeature::PositionX => "position_x",
        NodeFeature::PositionY => "position_y",
        NodeFeature::PositionZ => "position_z",
    }
}

pub(super) const fn edge_label(feature: EdgeFeature) -> &'static str {
    match feature {
        EdgeFeature::Distance => "distance",
        EdgeFeature::BondOrder => "bond_order",
    }
}
