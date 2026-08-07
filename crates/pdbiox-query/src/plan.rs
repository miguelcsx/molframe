//! Storage-independent logical plans and structure-bound physical plans.

use crate::ast::{Column, Expr};
use crate::glob::Glob;
use pdbiox_core::contract::AnalysisPolicy;
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::structure::Structure;
use pdbiox_core::symbol::SymbolId;
use std::collections::BTreeSet;

/// An immutable storage-independent selection plan.
#[derive(Clone, Debug, PartialEq)]
pub struct LogicalPlan(pub(crate) Expr);

/// A selection plan bound to one structure dictionary and policy.
///
/// Glob patterns have become stable symbol sets and constant branches have
/// been removed. Reuse this value for repeated frames sharing the topology.
#[derive(Clone, Debug)]
pub struct PhysicalQuery {
    pub(crate) expr: PhysicalExpr,
    pub(crate) warnings: Vec<Diagnostic>,
}

#[derive(Clone, Debug)]
pub(crate) enum PhysicalExpr {
    All,
    None,
    And(Box<Self>, Box<Self>),
    Or(Box<Self>, Box<Self>),
    Not(Box<Self>),
    ResolvedMembership {
        column: Column,
        symbols: BTreeSet<SymbolId>,
    },
    Logical(Expr),
}

pub(crate) fn bind(
    logical: &LogicalPlan,
    structure: &Structure,
    policy: &AnalysisPolicy,
    mut warnings: Vec<Diagnostic>,
) -> PhysicalQuery {
    let expr = lower(&logical.0, structure, policy, &mut warnings);
    warnings.sort_by_key(Diagnostic::code);
    PhysicalQuery { expr, warnings }
}

fn lower(
    expr: &Expr,
    structure: &Structure,
    policy: &AnalysisPolicy,
    warnings: &mut Vec<Diagnostic>,
) -> PhysicalExpr {
    match expr {
        Expr::All => PhysicalExpr::All,
        Expr::None => PhysicalExpr::None,
        Expr::And(left, right) => conjunction(
            lower(left, structure, policy, warnings),
            lower(right, structure, policy, warnings),
        ),
        Expr::Or(left, right) => disjunction(
            lower(left, structure, policy, warnings),
            lower(right, structure, policy, warnings),
        ),
        Expr::Not(target) => match lower(target, structure, policy, warnings) {
            PhysicalExpr::All => PhysicalExpr::None,
            PhysicalExpr::None => PhysicalExpr::All,
            target => PhysicalExpr::Not(Box::new(target)),
        },
        Expr::Membership { column, values }
            if symbol_column(*column, policy)
                && !(matches!(column, Column::AlternateLocation)
                    && values
                        .iter()
                        .any(|value| value.eq_ignore_ascii_case("none"))) =>
        {
            let column = crate::predicate::resolved_column(*column, policy);
            let symbols = resolve_symbols(structure, values);
            if symbols.is_empty() {
                warnings.push(Diagnostic::new(Code::W4003));
                PhysicalExpr::None
            } else {
                PhysicalExpr::ResolvedMembership { column, symbols }
            }
        }
        other => PhysicalExpr::Logical(other.clone()),
    }
}

fn conjunction(left: PhysicalExpr, right: PhysicalExpr) -> PhysicalExpr {
    match (left, right) {
        (PhysicalExpr::None, _) | (_, PhysicalExpr::None) => PhysicalExpr::None,
        (PhysicalExpr::All, right) => right,
        (left, PhysicalExpr::All) => left,
        (left, right) if cost(&right) < cost(&left) => {
            PhysicalExpr::And(Box::new(right), Box::new(left))
        }
        (left, right) => PhysicalExpr::And(Box::new(left), Box::new(right)),
    }
}

fn disjunction(left: PhysicalExpr, right: PhysicalExpr) -> PhysicalExpr {
    match (left, right) {
        (PhysicalExpr::All, _) | (_, PhysicalExpr::All) => PhysicalExpr::All,
        (PhysicalExpr::None, right) => right,
        (left, PhysicalExpr::None) => left,
        (left, right) => PhysicalExpr::Or(Box::new(left), Box::new(right)),
    }
}

fn cost(expr: &PhysicalExpr) -> u8 {
    match expr {
        PhysicalExpr::None | PhysicalExpr::All => 0,
        PhysicalExpr::ResolvedMembership { .. } => 1,
        PhysicalExpr::And(_, _) | PhysicalExpr::Or(_, _) | PhysicalExpr::Not(_) => 2,
        PhysicalExpr::Logical(Expr::Geometric(_)) => 5,
        PhysicalExpr::Logical(_) => 3,
    }
}

fn symbol_column(column: Column, policy: &AnalysisPolicy) -> bool {
    if policy.identifiers == pdbiox_core::contract::Namespace::Explicit
        && matches!(
            column,
            Column::Chain | Column::ResidueName | Column::AtomName
        )
    {
        return false;
    }
    matches!(
        crate::predicate::resolved_column(column, policy),
        Column::LabelChain
            | Column::AuthChain
            | Column::LabelResidueName
            | Column::AuthResidueName
            | Column::LabelAtomName
            | Column::AuthAtomName
            | Column::AlternateLocation
            | Column::Entity
            | Column::Element
            | Column::SegmentId
            | Column::InsertionCode
    )
}

fn resolve_symbols(structure: &Structure, values: &[Box<str>]) -> BTreeSet<SymbolId> {
    let patterns: Vec<Glob> = values.iter().map(|value| Glob::new(value)).collect();
    structure
        .data()
        .dictionary
        .iter()
        .filter_map(|(symbol, text)| {
            patterns
                .iter()
                .any(|pattern| pattern.matches(text))
                .then_some(symbol)
        })
        .collect()
}

#[cfg(test)]
#[path = "plan_tests.rs"]
mod tests;
