//! Named query definitions and their typed, closed resolution.
//!
//! A named query is a declarative definition — `pocket` means "residues within
//! five ångströms of the heme" — not the atoms that definition matched on some
//! earlier structure. Keeping the definition is what lets a later redefinition
//! of `heme` flow into `pocket` without anyone re-running `pocket` by hand, and
//! what lets one definition be evaluated against many structures.
//!
//! A query refers to a definition by name, written `$name` (or the longer
//! `group name`). Resolution replaces each reference with the referenced
//! definition's typed plan, recursively, producing a *closed* query that no
//! longer mentions any name. It works on the typed plan directly: nothing is
//! parsed again and no atom is touched, and the closed query's canonical text
//! is printed from the substituted plan rather than spliced from strings.
//!
//! Resolution is deterministic (definitions live in an ordered map), rejects an
//! unknown name, a cycle, or a chain of references deeper than
//! [`QueryAliases::MAX_DEPTH`], and costs `O(E)` in the size of the resolved
//! plan: each definition is resolved at most once per call and reused for
//! every further reference to it.
//!
//! A runtime [`crate::Groups`] map is a different thing — atoms already
//! selected by someone else, supplied at evaluation time. A query whose names
//! are meant to be runtime groups is evaluated as it is, without resolution.

use crate::Query;
use crate::ast::{Expr, GeometricExpr};
use molframe_core::diagnostic::{Code, Diagnostic};
use std::collections::{BTreeMap, BTreeSet};

/// An ordered set of named query definitions.
///
/// # Examples
///
/// ```
/// use molframe_query::{Query, QueryAliases};
///
/// let compile = |source| match Query::compile(source) {
///     Ok(query) => query,
///     Err(findings) => panic!("{findings:?}"),
/// };
/// let mut aliases = QueryAliases::new();
/// assert!(aliases.define("heme", compile("resname HEM")).is_ok());
/// let pocket = compile("byres (within 5 of $heme) and protein");
/// let Ok(closed) = aliases.resolve(&pocket) else {
///     panic!("every name is defined");
/// };
/// assert!(closed.references().is_empty());
/// assert_eq!(
///     closed.fingerprint(),
///     compile("byres (within 5 of (resname HEM)) and protein").fingerprint(),
/// );
/// ```
#[derive(Clone, Debug, Default)]
pub struct QueryAliases {
    definitions: BTreeMap<Box<str>, Query>,
}

impl QueryAliases {
    /// The deepest chain of references resolution follows before refusing.
    ///
    /// Interactive definitions nest a handful of levels deep; the bound exists
    /// so that a pathological or adversarial chain cannot exhaust the stack.
    pub const MAX_DEPTH: usize = 64;

    /// An empty set of definitions.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Whether `name` is spelled as a definable name.
    ///
    /// A name starts with an ASCII letter or underscore and continues with
    /// ASCII letters, digits and underscores, which is exactly what may follow
    /// `$` in a query.
    #[must_use]
    pub fn is_valid_name(name: &str) -> bool {
        crate::parser::is_query_name(name)
    }

    /// Defines or redefines `name`, returning the definition it replaced.
    ///
    /// The definition is stored as written: its own references are resolved
    /// when a query that uses it is resolved, so redefining a name it depends on
    /// changes what it means.
    ///
    /// # Errors
    ///
    /// Returns a syntax diagnostic when `name` is not a valid name.
    pub fn define(
        &mut self,
        name: impl Into<Box<str>>,
        query: Query,
    ) -> Result<Option<Query>, Diagnostic> {
        let name = name.into();
        if !Self::is_valid_name(&name) {
            return Err(Diagnostic::new(Code::E4001)
                .with_message("a query name is a letter or underscore followed by letters, digits or underscores")
                .with_context("name", name));
        }
        Ok(self.definitions.insert(name, query))
    }

    /// Removes `name`, returning its definition.
    pub fn remove(&mut self, name: &str) -> Option<Query> {
        self.definitions.remove(name)
    }

    /// The definition of `name`, as written.
    #[must_use]
    pub fn get(&self, name: &str) -> Option<&Query> {
        self.definitions.get(name)
    }

    /// Whether `name` is defined.
    #[must_use]
    pub fn contains(&self, name: &str) -> bool {
        self.definitions.contains_key(name)
    }

    /// Defined names in ascending order.
    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.definitions.keys().map(AsRef::as_ref)
    }

    /// Number of definitions.
    #[must_use]
    pub fn len(&self) -> usize {
        self.definitions.len()
    }

    /// Whether nothing is defined.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.definitions.is_empty()
    }

    /// Replaces every named reference in `query` with its definition.
    ///
    /// The result mentions no name, so it can be evaluated, fingerprinted or
    /// sent elsewhere without these definitions. A query with no references is
    /// returned unchanged, source text included.
    ///
    /// # Errors
    ///
    /// Returns `E4005` for a name with no definition, `E4006` for definitions
    /// that refer to each other in a cycle, and `E4007` when references nest
    /// deeper than [`Self::MAX_DEPTH`].
    pub fn resolve(&self, query: &Query) -> Result<Query, Vec<Diagnostic>> {
        if !has_reference(&query.expr) {
            return Ok(query.clone());
        }
        let mut resolver = Resolver {
            aliases: self,
            resolved: BTreeMap::new(),
            stack: Vec::new(),
            warnings: query.warnings.clone(),
        };
        let expr = resolver
            .substitute(&query.expr)
            .map_err(|finding| vec![finding])?;
        let mut warnings = resolver.warnings;
        warnings.sort_by_key(Diagnostic::code);
        let source = crate::print::print(&expr).into_boxed_str();
        Ok(Query {
            expr,
            warnings,
            source,
        })
    }
}

impl Query {
    /// Names this query refers to, in ascending order and without repeats.
    ///
    /// These are the definitions the query depends on directly; a definition's
    /// own references are its own dependencies. Runtime is `O(E)`.
    #[must_use]
    pub fn references(&self) -> Vec<&str> {
        let mut names = BTreeSet::new();
        collect_references(&self.expr, &mut names);
        names.into_iter().collect()
    }
}

struct Resolver<'a> {
    aliases: &'a QueryAliases,
    /// Each definition's closed plan, resolved at most once per call.
    resolved: BTreeMap<&'a str, Expr>,
    /// The names being resolved, outermost first, for cycle reporting.
    stack: Vec<&'a str>,
    warnings: Vec<Diagnostic>,
}

impl Resolver<'_> {
    fn substitute(&mut self, expr: &Expr) -> Result<Expr, Diagnostic> {
        Ok(match expr {
            Expr::Group(name) => self.definition(name)?,
            Expr::And(left, right) => Expr::And(
                Box::new(self.substitute(left)?),
                Box::new(self.substitute(right)?),
            ),
            Expr::Or(left, right) => Expr::Or(
                Box::new(self.substitute(left)?),
                Box::new(self.substitute(right)?),
            ),
            Expr::Not(value) => Expr::Not(Box::new(self.substitute(value)?)),
            Expr::Global(value) => Expr::Global(Box::new(self.substitute(value)?)),
            Expr::ByResidue(value) => Expr::ByResidue(Box::new(self.substitute(value)?)),
            Expr::Same { key, target } => Expr::Same {
                key: key.clone(),
                target: Box::new(self.substitute(target)?),
            },
            Expr::Bonded { depth, target } => Expr::Bonded {
                depth: *depth,
                target: Box::new(self.substitute(target)?),
            },
            Expr::Geometric(geometric) => Expr::Geometric(self.geometric(geometric)?),
            leaf => leaf.clone(),
        })
    }

    fn geometric(&mut self, expr: &GeometricExpr) -> Result<GeometricExpr, Diagnostic> {
        let mut resolved = expr.clone();
        if let Some(target) = geometric_target_mut(&mut resolved) {
            **target = self.substitute(target)?;
        }
        Ok(resolved)
    }

    fn definition(&mut self, name: &str) -> Result<Expr, Diagnostic> {
        let Some((key, query)) = self.aliases.definitions.get_key_value(name) else {
            return Err(Diagnostic::new(Code::E4005).with_context("name", name));
        };
        if let Some(expr) = self.resolved.get(key.as_ref()) {
            return Ok(expr.clone());
        }
        if let Some(start) = self.stack.iter().position(|entry| *entry == name) {
            let mut cycle = self.stack[start..].join(" -> ");
            cycle.push_str(" -> ");
            cycle.push_str(name);
            return Err(Diagnostic::new(Code::E4006).with_context("cycle", cycle));
        }
        if self.stack.len() >= QueryAliases::MAX_DEPTH {
            return Err(Diagnostic::new(Code::E4007).with_context("name", name));
        }
        self.stack.push(key);
        let expr = self.substitute(&query.expr);
        let _ = self.stack.pop();
        let expr = expr?;
        self.warnings.extend(query.warnings.iter().cloned());
        let _ = self.resolved.insert(key, expr.clone());
        Ok(expr)
    }
}

fn geometric_target_mut(expr: &mut GeometricExpr) -> Option<&mut Box<Expr>> {
    match expr {
        GeometricExpr::Within { target, .. }
        | GeometricExpr::Beyond { target, .. }
        | GeometricExpr::Around { target, .. }
        | GeometricExpr::SphereZone { target, .. }
        | GeometricExpr::SphereLayer { target, .. }
        | GeometricExpr::IsoLayer { target, .. }
        | GeometricExpr::CylinderZone { target, .. }
        | GeometricExpr::CylinderLayer { target, .. } => Some(target),
        GeometricExpr::Point { .. } => None,
    }
}

fn geometric_target(expr: &GeometricExpr) -> Option<&Expr> {
    match expr {
        GeometricExpr::Within { target, .. }
        | GeometricExpr::Beyond { target, .. }
        | GeometricExpr::Around { target, .. }
        | GeometricExpr::SphereZone { target, .. }
        | GeometricExpr::SphereLayer { target, .. }
        | GeometricExpr::IsoLayer { target, .. }
        | GeometricExpr::CylinderZone { target, .. }
        | GeometricExpr::CylinderLayer { target, .. } => Some(target),
        GeometricExpr::Point { .. } => None,
    }
}

/// Calls `visit` on each direct sub-expression of `expr`.
fn children<'e>(expr: &'e Expr, visit: &mut impl FnMut(&'e Expr)) {
    match expr {
        Expr::And(left, right) | Expr::Or(left, right) => {
            visit(left);
            visit(right);
        }
        Expr::Not(value) | Expr::Global(value) | Expr::ByResidue(value) => visit(value),
        Expr::Same { target, .. } | Expr::Bonded { target, .. } => visit(target),
        Expr::Geometric(geometric) => {
            if let Some(target) = geometric_target(geometric) {
                visit(target);
            }
        }
        Expr::All
        | Expr::None
        | Expr::Comparison { .. }
        | Expr::Membership { .. }
        | Expr::Group(_)
        | Expr::Atom { .. }
        | Expr::Macro(_)
        | Expr::Chirality(_)
        | Expr::Smarts { .. } => {}
    }
}

fn has_reference(expr: &Expr) -> bool {
    if matches!(expr, Expr::Group(_)) {
        return true;
    }
    let mut found = false;
    children(expr, &mut |child| found = found || has_reference(child));
    found
}

fn collect_references<'e>(expr: &'e Expr, names: &mut BTreeSet<&'e str>) {
    if let Expr::Group(name) = expr {
        let _ = names.insert(name);
        return;
    }
    children(expr, &mut |child| collect_references(child, names));
}

#[cfg(test)]
#[path = "aliases_tests.rs"]
mod tests;
