//! Compile-time checked Rust builder surface.

use crate::ast::{Column, Expr, GeometricExpr, Macro, Operator};
use std::ops::{BitAnd, BitOr, Not};

/// A typed selection expression under construction.
#[derive(Clone, PartialEq, Debug)]
pub struct Builder {
    expr: Expr,
    source: Box<str>,
}

impl Builder {
    /// Consumes the builder into the shared selection IR.
    #[must_use]
    pub(crate) fn into_parts(self) -> (Expr, Box<str>) {
        (self.expr, self.source)
    }

    /// Canonical textual form accepted by the query parser.
    #[must_use]
    pub fn source(&self) -> &str {
        &self.source
    }

    fn new(expr: Expr, source: impl Into<Box<str>>) -> Self {
        Self {
            expr,
            source: source.into(),
        }
    }
}

impl BitAnd for Builder {
    type Output = Self;

    fn bitand(self, right: Self) -> Self::Output {
        let source = format!("({}) and ({})", self.source, right.source);
        Self::new(Expr::And(Box::new(self.expr), Box::new(right.expr)), source)
    }
}

impl BitOr for Builder {
    type Output = Self;

    fn bitor(self, right: Self) -> Self::Output {
        let source = format!("({}) or ({})", self.source, right.source);
        Self::new(Expr::Or(Box::new(self.expr), Box::new(right.expr)), source)
    }
}

impl Not for Builder {
    type Output = Self;

    fn not(self) -> Self::Output {
        let source = format!("not ({})", self.source);
        Self::new(Expr::Not(Box::new(self.expr)), source)
    }
}

/// A typed column awaiting a predicate.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct ColumnBuilder(Column, &'static str);

impl ColumnBuilder {
    /// String or integer membership equality.
    #[must_use]
    pub fn eq(self, value: impl Into<Box<str>>) -> Builder {
        let value = value.into();
        Builder::new(
            Expr::Membership {
                column: self.0,
                values: vec![value.clone()],
            },
            format!("{} {:?}", self.1, value),
        )
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
        let symbol = match operator {
            Operator::Less => "<",
            Operator::LessEqual => "<=",
            Operator::Greater => ">",
            Operator::GreaterEqual => ">=",
            Operator::Equal => "=",
            Operator::NotEqual => "!=",
        };
        Builder::new(
            Expr::Comparison {
                column: self.0,
                operator,
                value,
                absolute: false,
            },
            format!("{} {symbol} {value}", self.1),
        )
    }
}

/// Builder functions sharing the textual parser's typed IR.
pub mod col {
    use super::{Builder, Column, ColumnBuilder, Expr, GeometricExpr, Macro};

    /// Every atom.
    #[must_use]
    pub fn all() -> Builder {
        Builder::new(Expr::All, "all")
    }

    /// No atoms.
    #[must_use]
    pub fn none() -> Builder {
        Builder::new(Expr::None, "none")
    }

    /// Protein atoms.
    #[must_use]
    pub fn is_protein() -> Builder {
        Builder::new(Expr::Macro(Macro::Protein), "protein")
    }

    macro_rules! molecular_selector {
        ($name:ident, $variant:ident, $source:literal, $doc:literal) => {
            #[doc = $doc]
            #[must_use]
            pub fn $name() -> Builder {
                Builder::new(Expr::Macro(Macro::$variant), $source)
            }
        };
    }

    molecular_selector!(backbone, Backbone, "backbone", "Protein backbone atoms.");
    molecular_selector!(
        sidechain,
        Sidechain,
        "sidechain",
        "Protein side-chain atoms."
    );
    molecular_selector!(nucleic, Nucleic, "nucleic", "Nucleic-acid atoms.");
    molecular_selector!(
        nucleic_backbone,
        NucleicBackbone,
        "nucleicbackbone",
        "Nucleic-acid backbone atoms."
    );
    molecular_selector!(
        nucleic_base,
        NucleicBase,
        "nucleicbase",
        "Nucleic-acid base atoms."
    );
    molecular_selector!(
        nucleic_sugar,
        NucleicSugar,
        "nucleicsugar",
        "Nucleic-acid sugar atoms."
    );
    molecular_selector!(water, Water, "water", "Water atoms.");
    molecular_selector!(ions, Ion, "ion", "Ion atoms.");
    molecular_selector!(lipids, Lipid, "lipid", "Lipid atoms.");
    molecular_selector!(glycans, Saccharide, "saccharide", "Saccharide atoms.");
    molecular_selector!(hetero, Hetero, "hetero", "Heterogeneous atoms.");
    molecular_selector!(hydrogen, Hydrogen, "hydrogen", "Hydrogen atoms.");
    molecular_selector!(heavy, Heavy, "heavy", "Non-hydrogen atoms.");
    molecular_selector!(polymer, Polymer, "polymer", "Polymer atoms.");
    molecular_selector!(ligands, Ligand, "ligand", "Non-polymer ligand atoms.");
    molecular_selector!(aromatic, Aromatic, "aromatic", "Aromatic atoms.");

    /// Atom name column.
    #[must_use]
    pub const fn name() -> ColumnBuilder {
        ColumnBuilder(Column::AtomName, "name")
    }

    /// Residue name column.
    #[must_use]
    pub const fn resname() -> ColumnBuilder {
        ColumnBuilder(Column::ResidueName, "resname")
    }

    /// Chain identifier column.
    #[must_use]
    pub const fn chain() -> ColumnBuilder {
        ColumnBuilder(Column::Chain, "chain")
    }

    /// Temperature-factor column.
    #[must_use]
    pub const fn bfactor() -> ColumnBuilder {
        ColumnBuilder(Column::BFactor, "bfactor")
    }

    /// Occupancy column.
    #[must_use]
    pub const fn occupancy() -> ColumnBuilder {
        ColumnBuilder(Column::Occupancy, "occupancy")
    }

    /// X coordinate.
    #[must_use]
    pub const fn x() -> ColumnBuilder {
        ColumnBuilder(Column::X, "x")
    }

    /// Atoms within `radius` of `target`.
    #[must_use]
    pub fn within(radius: f32, target: Builder) -> Builder {
        let source = format!("within {radius} of ({})", target.source);
        Builder::new(
            Expr::Geometric(GeometricExpr::Within {
                radius,
                target: Box::new(target.expr),
            }),
            source,
        )
    }

    /// Whole residues containing `target` atoms.
    #[must_use]
    pub fn by_residue(target: Builder) -> Builder {
        let source = format!("byres ({})", target.source);
        Builder::new(Expr::ByResidue(Box::new(target.expr)), source)
    }

    /// Whole residues within `radius` of `target`.
    #[must_use]
    pub fn residues_within(radius: f32, target: Builder) -> Builder {
        by_residue(within(radius, target))
    }
}

#[cfg(test)]
#[path = "builder_tests.rs"]
mod tests;
