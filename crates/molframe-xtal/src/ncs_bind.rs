//! Binding NCS operators to the chains they relate.
//!
//! An operator table alone does not say which molecules it maps. Declared
//! pairings name the domains; when the file declares none, the operator is
//! tested against every pair of chains with the same residue sequence and kept
//! where the mapped backbone lands on the other chain within a tolerance.

use crate::{
    AssemblyDef, AssemblySet, AssemblyView, Generator, NcsCode, NcsDomain, NcsOperator, NcsSet,
    OperExpression, Operator,
};
use molframe_core::{ChainIndex, Code, Diagnostic, Structure};
use std::collections::HashMap;

const NCS_ASSEMBLY: &str = "ncs";

/// Largest chain-pair count a geometric search may test.
pub const DEFAULT_NCS_PAIR_LIMIT: usize = 4_000_000;

/// Limits of a binding search.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NcsBindOptions {
    /// Largest backbone RMSD, in ångström, at which a chain pair is bound.
    pub max_rmsd: f64,
    /// Whether `generate` operators are tested too; their copies are absent
    /// from the coordinates, so they bind nothing unless a file is unusual.
    pub include_generate: bool,
    /// Ceiling on operator-by-chain-pair tests.
    pub pair_limit: usize,
}

impl Default for NcsBindOptions {
    fn default() -> Self {
        Self {
            max_rmsd: 2.0,
            include_generate: false,
            pair_limit: DEFAULT_NCS_PAIR_LIMIT,
        }
    }
}

/// One operator bound to the chain it maps and the chain it lands on.
#[derive(Clone, Debug, PartialEq)]
pub struct NcsChainPair {
    /// Operator identifier.
    pub operator: Box<str>,
    /// Chain the operator is applied to.
    pub source: ChainIndex,
    /// Chain the operator produces.
    pub target: ChainIndex,
    /// Backbone RMSD of the mapped source against the target, when the two
    /// have the same number of representative atoms.
    pub rmsd: Option<f64>,
    /// Number of atom pairs the RMSD covers.
    pub atoms: usize,
    /// Whether the file declared the pairing rather than geometry finding it.
    pub declared: bool,
}

/// A chain with the representative atoms inside a domain's limits.
type DomainChain<'a> = (&'a Backbone, Vec<[f32; 3]>);

/// One chain's representative atoms in residue order.
struct Backbone {
    chain: ChainIndex,
    label: Box<str>,
    names: Vec<Box<str>>,
    seq_ids: Vec<Option<i32>>,
    positions: Vec<[f32; 3]>,
}

/// Binds every operator of `set` to the chain pairs it relates.
///
/// # Errors
///
/// Returns [`Code::E6014`] when the declared pairings name a chain the
/// structure lacks, and [`Code::E6015`] when the search would exceed the
/// configured pair limit.
pub fn bind_ncs(
    structure: &Structure,
    set: &NcsSet,
    options: &NcsBindOptions,
) -> Result<Vec<NcsChainPair>, Diagnostic> {
    let backbones = backbones(structure);
    if !set.domains().pairings.is_empty() {
        return declared_pairs(&backbones, set);
    }
    let operators: Vec<&NcsOperator> = set
        .operators()
        .filter(|operator| options.include_generate || operator.code == NcsCode::Given)
        .collect();
    geometric_pairs(&backbones, &operators, options)
}

fn backbones(structure: &Structure) -> Vec<Backbone> {
    let mut out = Vec::new();
    let Some(model) = structure.data().models().next() else {
        return out;
    };
    for chain in model.chains() {
        let label = match chain.label() {
            Some(label) => label,
            None => "",
        };
        let mut backbone = Backbone {
            chain: chain.index(),
            label: label.into(),
            names: Vec::new(),
            seq_ids: Vec::new(),
            positions: Vec::new(),
        };
        for residue in chain.residues() {
            let atom = residue.atom("CA").or_else(|| residue.atom("C4'"));
            if let (Some(atom), Some(name)) = (atom, residue.name())
                && let Some(position) = atom.position()
            {
                backbone.names.push(name.into());
                backbone.seq_ids.push(residue.label_seq_id());
                backbone.positions.push(position);
            }
        }
        if !backbone.positions.is_empty() {
            out.push(backbone);
        }
    }
    out
}

fn rmsd(operator: &NcsOperator, source: &[[f32; 3]], target: &[[f32; 3]]) -> Option<f64> {
    if source.len() != target.len() || source.is_empty() {
        return None;
    }
    let sum: f64 = source
        .iter()
        .zip(target)
        .map(|(from, to)| {
            let moved = operator.transform.apply(*from);
            (0..3)
                .map(|axis| f64::from(moved[axis] - to[axis]).powi(2))
                .sum::<f64>()
        })
        .sum();
    let count = u32::try_from(source.len()).ok()?;
    Some((sum / f64::from(count)).sqrt())
}

fn geometric_pairs(
    backbones: &[Backbone],
    operators: &[&NcsOperator],
    options: &NcsBindOptions,
) -> Result<Vec<NcsChainPair>, Diagnostic> {
    let mut groups: HashMap<&[Box<str>], Vec<&Backbone>> = HashMap::new();
    for backbone in backbones {
        groups.entry(&backbone.names).or_default().push(backbone);
    }
    let mut tests = 0_usize;
    for members in groups.values() {
        tests = members
            .len()
            .checked_mul(members.len())
            .and_then(|pairs| pairs.checked_mul(operators.len()))
            .and_then(|count| tests.checked_add(count))
            .filter(|count| *count <= options.pair_limit)
            .ok_or_else(|| {
                Diagnostic::new(Code::E6015).with_context("limit", options.pair_limit.to_string())
            })?;
    }
    let mut pairs = Vec::new();
    for operator in operators {
        for backbone in backbones {
            let Some(members) = groups.get(backbone.names.as_slice()) else {
                continue;
            };
            for other in members {
                if other.chain == backbone.chain {
                    continue;
                }
                let Some(value) = rmsd(operator, &backbone.positions, &other.positions) else {
                    continue;
                };
                if value <= options.max_rmsd {
                    pairs.push(NcsChainPair {
                        operator: operator.id.clone(),
                        source: backbone.chain,
                        target: other.chain,
                        rmsd: Some(value),
                        atoms: backbone.positions.len(),
                        declared: false,
                    });
                }
            }
        }
    }
    Ok(pairs)
}

/// Chain-by-chain mapping of two declared domains.
fn declared_pairs(backbones: &[Backbone], set: &NcsSet) -> Result<Vec<NcsChainPair>, Diagnostic> {
    let mut pairs = Vec::new();
    for pairing in &set.domains().pairings {
        let (Some(first), Some(second), Some(operator)) = (
            set.domains().domain(&pairing.domain_1),
            set.domains().domain(&pairing.domain_2),
            set.get(&pairing.operator),
        ) else {
            continue;
        };
        let from = domain_chains(backbones, first)?;
        let to = domain_chains(backbones, second)?;
        for ((source, source_points), (target, target_points)) in from.iter().zip(&to) {
            let value = rmsd(operator, source_points, target_points);
            pairs.push(NcsChainPair {
                operator: operator.id.clone(),
                source: source.chain,
                target: target.chain,
                rmsd: value,
                atoms: if value.is_some() {
                    source_points.len()
                } else {
                    0
                },
                declared: true,
            });
        }
    }
    Ok(pairs)
}

/// The chains a domain's limits touch, each with the representative atoms
/// inside the limits.
fn domain_chains<'a>(
    backbones: &'a [Backbone],
    domain: &NcsDomain,
) -> Result<Vec<DomainChain<'a>>, Diagnostic> {
    let mut chains: Vec<DomainChain<'a>> = Vec::new();
    for limit in &domain.limits {
        let backbone = backbones
            .iter()
            .find(|candidate| *candidate.label == *limit.begin.asym_id)
            .ok_or_else(|| {
                Diagnostic::new(Code::E6014)
                    .with_context("domain", domain.id.to_string())
                    .with_context("label_asym_id", limit.begin.asym_id.to_string())
            })?;
        let (begin, end) = (limit.begin.seq_id, limit.end.seq_id);
        let inside = backbone
            .seq_ids
            .iter()
            .zip(&backbone.positions)
            .filter(|(seq, _)| match (seq, begin, end) {
                (Some(seq), Some(begin), Some(end)) => (begin..=end).contains(seq),
                _ => true,
            })
            .map(|(_, position)| *position);
        match chains
            .iter_mut()
            .find(|(known, _)| known.chain == backbone.chain)
        {
            Some((_, points)) => points.extend(inside),
            None => chains.push((backbone, inside.collect())),
        }
    }
    Ok(chains)
}

/// The structure with every `generate` operator applied to all of its chains.
///
/// The deposited coordinates come first, unchanged, followed by one copy of
/// the whole model per operator.
///
/// # Errors
///
/// Returns a diagnostic when an operator is not a proper rotation, the copy
/// count exceeds the instance limit, or materialisation fails.
pub fn expand_ncs(structure: &Structure, set: &NcsSet) -> Result<Structure, Vec<Diagnostic>> {
    let generators: Vec<&NcsOperator> = set.generators().collect();
    if generators.is_empty() {
        return Ok(structure.clone());
    }
    let mut operators = vec![Operator::from_matrix(
        "0",
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [0.0; 3],
    )];
    for (position, operator) in generators.iter().enumerate() {
        operators.push(Operator::from_matrix(
            (position + 1).to_string(),
            operator.transform.matrix,
            operator.transform.translation,
        ));
    }
    let expression = (0..operators.len())
        .map(|id| id.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let asym_ids: Box<[Box<str>]> = structure
        .data()
        .models()
        .next()
        .into_iter()
        .flat_map(molframe_core::structure::ModelRef::chains)
        .filter_map(|chain| chain.label().map(Box::from))
        .collect();
    let definition = AssemblyDef {
        id: NCS_ASSEMBLY.into(),
        details: None,
        method: None,
        oligomeric: None,
        generators: vec![Generator {
            oper_expression: OperExpression::parse(&format!("({expression})"))
                .map_err(|error| vec![error])?,
            asym_ids,
        }],
    };
    let assemblies =
        AssemblySet::from_definitions([definition], operators).map_err(|error| vec![error])?;
    AssemblyView::new(structure, &assemblies, NCS_ASSEMBLY)
        .map_err(|error| vec![error])?
        .materialize()
}

#[cfg(test)]
#[path = "ncs_bind_tests.rs"]
mod tests;
