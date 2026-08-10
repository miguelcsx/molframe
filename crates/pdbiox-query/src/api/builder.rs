//! Compile-time checked Rust builder surface.

use crate::ast::{Column, Expr, GeometricExpr, Macro, Operator};
use std::ops::{BitAnd, BitOr, Not};

/// A typed selection expression under construction.
#[derive(Clone, PartialEq, Debug)]
pub struct Builder(pub(crate) Expr);

impl Builder {
    /// Consumes the builder into the shared selection IR.
    #[must_use]
    pub(crate) fn into_expr(self) -> Expr {
        self.0
    }
}

impl From<Expr> for Builder {
    fn from(expr: Expr) -> Self {
        Self(expr)
    }
}

impl BitAnd for Builder {
    type Output = Self;

    fn bitand(self, right: Self) -> Self::Output {
        Self(Expr::And(Box::new(self.0), Box::new(right.0)))
    }
}

impl BitOr for Builder {
    type Output = Self;

    fn bitor(self, right: Self) -> Self::Output {
        Self(Expr::Or(Box::new(self.0), Box::new(right.0)))
    }
}

impl Not for Builder {
    type Output = Self;

    fn not(self) -> Self::Output {
        Self(Expr::Not(Box::new(self.0)))
    }
}

/// A typed column awaiting a predicate.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct ColumnBuilder(Column);

impl ColumnBuilder {
    /// String or integer membership equality.
    #[must_use]
    pub fn eq(self, value: impl Into<Box<str>>) -> Builder {
        Builder(Expr::Membership {
            column: self.0,
            values: vec![value.into()],
        })
    }

    /// Numeric less-than comparison.
    #[must_use]
    pub fn lt(self, value: f64) -> Builder {
        self.compare(Operator::Less, value)
    }

    /// Numeric less-than-or-equal comparison.
    #[must_use]
    pub fn le(self, value: f64) -> Builder {
        self.compare(Operator::LessEqual, value)
    }

    /// Numeric greater-than comparison.
    #[must_use]
    pub fn gt(self, value: f64) -> Builder {
        self.compare(Operator::Greater, value)
    }

    /// Numeric greater-than-or-equal comparison.
    #[must_use]
    pub fn ge(self, value: f64) -> Builder {
        self.compare(Operator::GreaterEqual, value)
    }

    fn compare(self, operator: Operator, value: f64) -> Builder {
        Builder(Expr::Comparison {
            column: self.0,
            operator,
            value,
            absolute: false,
        })
    }
}

/// Builder functions sharing the textual parser's typed IR.
pub mod col {
    use super::{Builder, Column, ColumnBuilder, Expr, GeometricExpr, Macro};

    /// Every atom.
    #[must_use]
    pub fn all() -> Builder {
        Builder(Expr::All)
    }

    /// No atoms.
    #[must_use]
    pub fn none() -> Builder {
        Builder(Expr::None)
    }

    /// Protein atoms.
    #[must_use]
    pub fn is_protein() -> Builder {
        Builder(Expr::Macro(Macro::Protein))
    }

    /// Atom name column.
    #[must_use]
    pub const fn name() -> ColumnBuilder {
        ColumnBuilder(Column::AtomName)
    }

    /// Residue name column.
    #[must_use]
    pub const fn resname() -> ColumnBuilder {
        ColumnBuilder(Column::ResidueName)
    }

    /// Chain identifier column.
    #[must_use]
    pub const fn chain() -> ColumnBuilder {
        ColumnBuilder(Column::Chain)
    }

    /// Temperature-factor column.
    #[must_use]
    pub const fn bfactor() -> ColumnBuilder {
        ColumnBuilder(Column::BFactor)
    }

    /// Occupancy column.
    #[must_use]
    pub const fn occupancy() -> ColumnBuilder {
        ColumnBuilder(Column::Occupancy)
    }

    /// X coordinate.
    #[must_use]
    pub const fn x() -> ColumnBuilder {
        ColumnBuilder(Column::X)
    }

    /// Atoms within `radius` of `target`.
    #[must_use]
    pub fn within(radius: f32, target: Builder) -> Builder {
        Builder(Expr::Geometric(GeometricExpr::Within {
            radius,
            target: Box::new(target.0),
        }))
    }

    /// Whole residues containing `target` atoms.
    #[must_use]
    pub fn by_residue(target: Builder) -> Builder {
        Builder(Expr::ByResidue(Box::new(target.0)))
    }
}

#[cfg(test)]
#[path = "builder_tests.rs"]
mod tests;
