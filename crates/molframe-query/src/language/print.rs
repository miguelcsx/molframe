//! Canonical text for a typed selection expression.
//!
//! A query built by composing other queries — by the typed builder, by the
//! boolean operators, or by substituting named definitions — must still carry
//! text that parses back to exactly the same plan, because that text is what
//! crosses process and language boundaries. This printer is the inverse of the
//! parser: every sub-expression is parenthesised, so no precedence rule is ever
//! relied upon, and every value that could be misread is quoted or escaped.
//!
//! Printing is `O(N)` in the size of the expression tree and its values.

use crate::ast::{Column, Expr, GeometricExpr, Macro, Operator, SameKey};
use std::fmt::Write as _;

/// Prints `expr` as canonical selection syntax.
pub(crate) fn print(expr: &Expr) -> String {
    let mut out = String::new();
    write_expr(&mut out, expr);
    out
}

fn write_expr(out: &mut String, expr: &Expr) {
    match expr {
        Expr::All => out.push_str("all"),
        Expr::None => out.push_str("none"),
        Expr::And(left, right) => binary(out, left, "and", right),
        Expr::Or(left, right) => binary(out, left, "or", right),
        Expr::Not(value) => prefixed(out, "not", value),
        Expr::Global(value) => prefixed(out, "global", value),
        Expr::ByResidue(value) => prefixed(out, "byres", value),
        Expr::Same { key, target } => {
            out.push_str("same ");
            out.push_str(same_key(key));
            out.push_str(" as ");
            group(out, target);
        }
        Expr::Bonded { depth, target } => {
            let _ = write!(out, "bonded {depth} ");
            group(out, target);
        }
        Expr::Geometric(geometric) => write_geometric(out, geometric),
        Expr::Comparison {
            column,
            operator,
            value,
            absolute,
        } => {
            if *absolute {
                out.push_str("prop abs ");
            }
            let _ = write!(
                out,
                "{} {} {value}",
                column_name(*column),
                operator_symbol(*operator)
            );
        }
        Expr::Membership { column, values } => {
            out.push_str(column_name(*column));
            for value in values {
                out.push(' ');
                membership_value(out, value);
            }
        }
        Expr::Group(name) => {
            if crate::parser::is_query_name(name) {
                out.push('$');
                out.push_str(name);
            } else {
                out.push_str("group ");
                quoted(out, name);
            }
        }
        Expr::Atom {
            segment,
            residue,
            name,
        } => {
            out.push_str("atom ");
            quoted(out, segment);
            let _ = write!(out, " {residue} ");
            quoted(out, name);
        }
        Expr::Macro(value) => out.push_str(macro_name(*value)),
        Expr::Chirality(value) => {
            out.push_str("chirality ");
            quoted(out, value);
        }
        Expr::Smarts { source, .. } => {
            out.push_str("smarts ");
            quoted(out, source);
        }
    }
}

fn write_geometric(out: &mut String, expr: &GeometricExpr) {
    match expr {
        GeometricExpr::Within { radius, target } => {
            let _ = write!(out, "within {radius} of ");
            group(out, target);
        }
        GeometricExpr::Beyond { radius, target } => {
            let _ = write!(out, "beyond {radius} of ");
            group(out, target);
        }
        GeometricExpr::Around { radius, target } => {
            let _ = write!(out, "around {radius} ");
            group(out, target);
        }
        GeometricExpr::SphereZone { radius, target } => {
            let _ = write!(out, "sphzone {radius} ");
            group(out, target);
        }
        GeometricExpr::SphereLayer {
            inner,
            outer,
            target,
        } => {
            let _ = write!(out, "sphlayer {inner} {outer} ");
            group(out, target);
        }
        GeometricExpr::IsoLayer {
            inner,
            outer,
            target,
        } => {
            let _ = write!(out, "isolayer {inner} {outer} ");
            group(out, target);
        }
        GeometricExpr::CylinderZone {
            radius,
            z_max,
            z_min,
            target,
        } => {
            let _ = write!(out, "cyzone {radius} {z_max} {z_min} ");
            group(out, target);
        }
        GeometricExpr::CylinderLayer {
            inner,
            outer,
            z_max,
            z_min,
            target,
        } => {
            let _ = write!(out, "cylayer {inner} {outer} {z_max} {z_min} ");
            group(out, target);
        }
        GeometricExpr::Point { point, radius } => {
            let _ = write!(out, "point {} {} {} {radius}", point[0], point[1], point[2]);
        }
    }
}

fn binary(out: &mut String, left: &Expr, keyword: &str, right: &Expr) {
    group(out, left);
    out.push(' ');
    out.push_str(keyword);
    out.push(' ');
    group(out, right);
}

fn prefixed(out: &mut String, keyword: &str, value: &Expr) {
    out.push_str(keyword);
    out.push(' ');
    group(out, value);
}

fn group(out: &mut String, expr: &Expr) {
    out.push('(');
    write_expr(out, expr);
    out.push(')');
}

/// Writes one membership pattern so the parser reads back the same pattern.
///
/// A bare word is kept bare. The two boolean keywords would end the value
/// list, so they take the parser's leading-backslash escape. Anything else that
/// the lexer would split or reinterpret is double-quoted with its quote and
/// backslash characters escaped; a pattern that itself begins with a backslash
/// gets one more, because the parser strips a single leading backslash from
/// every membership value.
fn membership_value(out: &mut String, value: &str) {
    if crate::parser::is_boundary(value) {
        out.push('\\');
        out.push_str(value);
        return;
    }
    if is_bare(value) {
        out.push_str(value);
        return;
    }
    if value.starts_with('\\') {
        out.push('"');
        out.push_str("\\\\");
        escape_into(out, value);
        out.push('"');
        return;
    }
    quoted(out, value);
}

fn quoted(out: &mut String, value: &str) {
    out.push('"');
    escape_into(out, value);
    out.push('"');
}

fn escape_into(out: &mut String, value: &str) {
    for character in value.chars() {
        if matches!(character, '"' | '\\') {
            out.push('\\');
        }
        out.push(character);
    }
}

/// Whether `value` lexes as exactly one unquoted value token and nothing else.
fn is_bare(value: &str) -> bool {
    !value.is_empty()
        && !value.bytes().any(|byte| {
            byte.is_ascii_whitespace()
                || matches!(
                    byte,
                    b'(' | b')' | b'<' | b'>' | b'=' | b'!' | b'"' | b'\'' | b'\\'
                )
        })
}

fn operator_symbol(operator: Operator) -> &'static str {
    match operator {
        Operator::Less => "<",
        Operator::LessEqual => "<=",
        Operator::Greater => ">",
        Operator::GreaterEqual => ">=",
        Operator::Equal => "==",
        Operator::NotEqual => "!=",
    }
}

fn same_key(key: &SameKey) -> &'static str {
    match key {
        SameKey::Residue => "residue",
        SameKey::Chain => "chain",
        SameKey::Model => "model",
        SameKey::Entity => "entity",
        SameKey::Fragment => "fragment",
        SameKey::Segment => "segment",
        SameKey::Column(column) => column_name(*column),
    }
}

pub(crate) const fn macro_name(value: Macro) -> &'static str {
    match value {
        Macro::Protein => "protein",
        Macro::Backbone => "backbone",
        Macro::Sidechain => "sidechain",
        Macro::Nucleic => "nucleic",
        Macro::NucleicBackbone => "nucleicbackbone",
        Macro::NucleicBase => "nucleicbase",
        Macro::NucleicSugar => "nucleicsugar",
        Macro::Water => "water",
        Macro::Ion => "ion",
        Macro::Lipid => "lipid",
        Macro::Saccharide => "saccharide",
        Macro::Hetero => "hetero",
        Macro::Hydrogen => "hydrogen",
        Macro::Heavy => "heavy",
        Macro::Polymer => "polymer",
        Macro::Ligand => "ligand",
        Macro::Aromatic => "aromatic",
        Macro::PolarHydrogen => "polar_hydrogen",
        Macro::NonpolarHydrogen => "nonpolar_hydrogen",
        Macro::Helix => "helix",
        Macro::Strand => "strand",
        Macro::Sheet => "sheet",
        Macro::AlphaHelix => "alpha_helix",
        Macro::Helix310 => "helix_310",
        Macro::PiHelix => "pi_helix",
        Macro::Polyproline => "polyproline",
        Macro::Bridge => "bridge",
        Macro::Turn => "turn",
        Macro::Bend => "bend",
        Macro::Coil => "coil",
    }
}

pub(crate) const fn column_name(column: Column) -> &'static str {
    match column {
        Column::Index => "index",
        Column::ByNumber => "bynum",
        Column::AtomSiteId => "id",
        Column::ResidueIndex => "resindex",
        Column::ChainIndex => "chainindex",
        Column::ModelIndex => "modelindex",
        Column::Chain => "chain",
        Column::LabelChain => "label_chain",
        Column::AuthChain => "auth_chain",
        Column::ResidueId => "resid",
        Column::LabelResidueId => "label_resid",
        Column::AuthResidueId => "auth_resid",
        Column::ResidueName => "resname",
        Column::LabelResidueName => "label_resname",
        Column::AuthResidueName => "auth_resname",
        Column::AtomName => "name",
        Column::LabelAtomName => "label_name",
        Column::AuthAtomName => "auth_name",
        Column::AlternateLocation => "altloc",
        Column::Model => "model",
        Column::Entity => "entity",
        Column::EntityType => "entity_type",
        Column::Element => "element",
        Column::SegmentId => "segid",
        Column::InsertionCode => "icode",
        Column::RecordType => "record_type",
        Column::X => "x",
        Column::Y => "y",
        Column::Z => "z",
        Column::Mass => "mass",
        Column::Charge => "charge",
        Column::FormalCharge => "formalcharge",
        Column::Radius => "radius",
        Column::BFactor => "bfactor",
        Column::Occupancy => "occupancy",
        Column::Plddt => "plddt",
        Column::Pae => "pae",
        Column::Assembly => "assembly",
    }
}

#[cfg(test)]
#[path = "print_tests.rs"]
mod tests;
