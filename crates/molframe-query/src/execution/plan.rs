//! Storage-independent logical plans and structure-bound physical plans.

use crate::ast::{Column, Expr};
use crate::glob::Glob;
use molframe_core::contract::AnalysisPolicy;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::structure::Structure;
use molframe_core::symbol::SymbolId;
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

/// Binds an expression directly without first cloning it into `LogicalPlan`.
///
/// This path allows `Query::plan` to eliminate one complete AST clone while
/// retaining exactly the same physical-plan representation.
pub(crate) fn bind_expr(
    expr: &Expr,
    structure: &Structure,
    policy: &AnalysisPolicy,
    mut warnings: Vec<Diagnostic>,
) -> PhysicalQuery {
    let expr = lower(expr, structure, policy, &mut warnings);

    warnings.sort_by_key(Diagnostic::code);

    PhysicalQuery { expr, warnings }
}

/// Lowers one logical expression into a structure-bound physical expression.
///
/// Constant branches are folded and supported symbol memberships are resolved
/// against the structure dictionary exactly once.
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

/// Folds and cost-orders one binary conjunction.
///
/// The right operand is evaluated inside the left operand's result, so operands
/// are exchanged only when the right-hand side has strictly lower estimated
/// cost *and* both read nothing but the atom they test. Exchanging an operand
/// that reads its surroundings (`within`, `byres`, `bonded`, ...) would change
/// what it sees and therefore the answer: `within 4 of X and not X` selects the
/// shell around `X`, while its mirror would test `X` against a universe that
/// already excludes it.
fn conjunction(left: PhysicalExpr, right: PhysicalExpr) -> PhysicalExpr {
    match (left, right) {
        (PhysicalExpr::None, _) | (_, PhysicalExpr::None) => PhysicalExpr::None,
        (PhysicalExpr::All, right) => right,
        (left, PhysicalExpr::All) => left,
        (left, right) if cost(&right) < cost(&left) && atom_local(&left) && atom_local(&right) => {
            PhysicalExpr::And(Box::new(right), Box::new(left))
        }
        (left, right) => PhysicalExpr::And(Box::new(left), Box::new(right)),
    }
}

/// Constant-folds one binary disjunction.
fn disjunction(left: PhysicalExpr, right: PhysicalExpr) -> PhysicalExpr {
    match (left, right) {
        (PhysicalExpr::All, _) | (_, PhysicalExpr::All) => PhysicalExpr::All,
        (PhysicalExpr::None, right) => right,
        (left, PhysicalExpr::None) => left,
        (left, right) => PhysicalExpr::Or(Box::new(left), Box::new(right)),
    }
}

/// Whether the expression's answer for one atom is independent of the universe
/// it is evaluated in, so its position in a conjunction cannot change results.
fn atom_local(expr: &PhysicalExpr) -> bool {
    match expr {
        PhysicalExpr::All | PhysicalExpr::None | PhysicalExpr::ResolvedMembership { .. } => true,
        PhysicalExpr::And(left, right) | PhysicalExpr::Or(left, right) => {
            atom_local(left) && atom_local(right)
        }
        PhysicalExpr::Not(target) => atom_local(target),
        PhysicalExpr::Logical(logical) => logical_atom_local(logical),
    }
}

fn logical_atom_local(expr: &Expr) -> bool {
    match expr {
        Expr::All
        | Expr::None
        | Expr::Comparison { .. }
        | Expr::Membership { .. }
        | Expr::Group(_)
        | Expr::Atom { .. } => true,
        Expr::And(left, right) | Expr::Or(left, right) => {
            logical_atom_local(left) && logical_atom_local(right)
        }
        Expr::Not(target) => logical_atom_local(target),
        _ => false,
    }
}

/// Returns a small relative execution-cost class for plan ordering.
///
/// Runtime and space are `O(1)`.
#[inline]
fn cost(expr: &PhysicalExpr) -> u8 {
    match expr {
        PhysicalExpr::None | PhysicalExpr::All => 0,
        PhysicalExpr::ResolvedMembership { .. } => 1,
        PhysicalExpr::And(_, _) | PhysicalExpr::Or(_, _) | PhysicalExpr::Not(_) => 2,
        PhysicalExpr::Logical(Expr::Geometric(_)) => 5,
        PhysicalExpr::Logical(_) => 3,
    }
}

/// Returns whether `column` can be lowered to dictionary symbol membership.
///
/// Explicit identifier-policy restrictions are preserved. Elements are not
/// text: a file may write iron as `FE` and a query as `Fe`, so they are
/// compared as elements rather than through the symbol dictionary.
fn symbol_column(column: Column, policy: &AnalysisPolicy) -> bool {
    if policy.identifiers == molframe_core::contract::Namespace::Explicit
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
            | Column::SegmentId
            | Column::InsertionCode
    )
}

/// Resolves membership patterns against the structure symbol dictionary.
///
/// Plain literal patterns use direct dictionary lookup and therefore avoid a
/// complete dictionary scan. Only true glob/escape patterns are compiled and
/// scanned across dictionary entries.
fn resolve_symbols(structure: &Structure, values: &[Box<str>]) -> BTreeSet<SymbolId> {
    let dictionary = &structure.data().dictionary;
    let mut symbols = BTreeSet::new();
    let mut patterns = Vec::new();

    for value in values {
        if is_plain_literal(value) {
            if let Some(symbol) = dictionary.get(value) {
                symbols.insert(symbol);
            }
        } else {
            patterns.push(Glob::new(value));
        }
    }

    if patterns.is_empty() {
        return symbols;
    }

    for (symbol, text) in dictionary.iter() {
        if patterns.iter().any(|pattern| pattern.matches(text)) {
            symbols.insert(symbol);
        }
    }

    symbols
}

/// Returns whether `pattern` has no glob or escape metacharacters.
///
/// Such patterns are semantically exact literals and can be resolved through a
/// direct dictionary lookup.
#[inline]
fn is_plain_literal(pattern: &str) -> bool {
    !pattern
        .bytes()
        .any(|byte| matches!(byte, b'*' | b'?' | b'[' | b'\\'))
}

#[cfg(test)]
#[path = "plan_tests.rs"]
mod tests;
