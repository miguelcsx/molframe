//! Precedence parser from tokens into the typed selection IR.

use crate::ast::{Column, Expr, GeometricExpr, Macro, Operator, SameKey};
use crate::lexer::{Token, TokenKind, lex};
use pdbiox_chem::SmartsPattern;
use pdbiox_core::diagnostic::{Code, Diagnostic};

#[path = "parser_helpers.rs"]
mod helpers;
use helpers::{
    canonical_name, geometric_keyword, is_boundary, mixes_boolean_precedence, syntax, unescape,
};
#[path = "parser_cursor.rs"]
mod cursor;

pub(crate) struct Parsed {
    pub(crate) expr: Expr,
    pub(crate) warnings: Vec<Diagnostic>,
}

/// Parses selection-language source into the typed expression tree.
///
/// Tokens are consumed by value through `Vec::IntoIter`, allowing owned
/// `Box<str>` token payloads to move directly into the AST instead of being
/// cloned through temporary `String`s.
pub(crate) fn parse(source: &str) -> Result<Parsed, Vec<Diagnostic>> {
    let tokens = lex(source).map_err(|finding| vec![finding])?;

    if tokens.is_empty() {
        return Err(vec![syntax(0, "selection is empty")]);
    }

    let mixed = mixes_boolean_precedence(&tokens);

    let mut parser = Parser {
        tokens: tokens.into_iter(),
        warnings: Vec::new(),
    };

    let expr = parser.parse_or().map_err(|finding| vec![finding])?;

    if !parser.tokens.as_slice().is_empty() {
        return Err(vec![parser.error_here("unexpected token after selection")]);
    }

    if mixed {
        parser.warnings.push(Diagnostic::new(Code::W4001));
    }

    Ok(Parsed {
        expr,
        warnings: parser.warnings,
    })
}

/// Stateful recursive-descent parser over an owning token iterator.
struct Parser {
    tokens: std::vec::IntoIter<Token>,
    warnings: Vec<Diagnostic>,
}

/// Geometric prefix recognized without allocating a lowercase copy.
#[derive(Clone, Copy)]
pub(super) enum GeometricKeyword {
    Within,
    Beyond,
    Around,
    SphereZone,
    SphereLayer,
    IsoLayer,
    CylinderZone,
    CylinderLayer,
    Point,
}

impl Parser {
    /// Parses the lowest-precedence `or` expression.
    fn parse_or(&mut self) -> Result<Expr, Diagnostic> {
        let mut expr = self.parse_and()?;

        while self.take_keyword("or") {
            expr = Expr::Or(Box::new(expr), Box::new(self.parse_and()?));
        }

        Ok(expr)
    }

    /// Parses left-associative `and` expressions.
    fn parse_and(&mut self) -> Result<Expr, Diagnostic> {
        let mut expr = self.parse_not()?;

        while self.take_keyword("and") {
            expr = Expr::And(Box::new(expr), Box::new(self.parse_not()?));
        }

        Ok(expr)
    }

    /// Parses unary `not` and `global` modifiers.
    fn parse_not(&mut self) -> Result<Expr, Diagnostic> {
        if self.take_keyword("not") {
            Ok(Expr::Not(Box::new(self.parse_not()?)))
        } else if self.take_keyword("global") {
            Ok(Expr::Global(Box::new(self.parse_not()?)))
        } else {
            self.parse_modifier()
        }
    }

    /// Parses hierarchy, connectivity and geometric modifiers.
    fn parse_modifier(&mut self) -> Result<Expr, Diagnostic> {
        if self.take_keyword("byres") {
            return Ok(Expr::ByResidue(Box::new(self.parse_modifier()?)));
        }

        if self.take_keyword("same") {
            let key = self.parse_same_key()?;
            self.require_keyword("as")?;

            return Ok(Expr::Same {
                key,
                target: Box::new(self.parse_modifier()?),
            });
        }

        if self.take_keyword("bonded") {
            let depth = self
                .peek_value()
                .and_then(|value| value.parse::<u32>().ok());

            let depth = match depth {
                Some(depth) => {
                    let _ = self.take_value();
                    depth
                }
                None => 1,
            };

            return Ok(Expr::Bonded {
                depth,
                target: Box::new(self.parse_modifier()?),
            });
        }

        if let Some(geometric) = self.parse_geometric()? {
            return Ok(Expr::Geometric(geometric));
        }

        self.parse_primary()
    }

    /// Parses a geometric prefix if the next value is a geometric keyword.
    ///
    /// Keyword recognition is ASCII case-insensitive and allocation-free.
    fn parse_geometric(&mut self) -> Result<Option<GeometricExpr>, Diagnostic> {
        let Some(kind) = self.peek_value().and_then(geometric_keyword) else {
            return Ok(None);
        };

        let _ = self.take_value();

        let geometric = match kind {
            GeometricKeyword::Within => {
                let radius = self.number()?;
                self.require_keyword("of")?;

                GeometricExpr::Within {
                    radius,
                    target: Box::new(self.parse_modifier()?),
                }
            }
            GeometricKeyword::Beyond => {
                let radius = self.number()?;
                self.require_keyword("of")?;

                GeometricExpr::Beyond {
                    radius,
                    target: Box::new(self.parse_modifier()?),
                }
            }
            GeometricKeyword::Around => {
                let radius = self.number()?;

                GeometricExpr::Around {
                    radius,
                    target: Box::new(self.parse_modifier()?),
                }
            }
            GeometricKeyword::SphereZone => {
                let radius = self.number()?;

                GeometricExpr::SphereZone {
                    radius,
                    target: Box::new(self.parse_modifier()?),
                }
            }
            GeometricKeyword::SphereLayer => {
                let inner = self.number()?;
                let outer = self.number()?;

                GeometricExpr::SphereLayer {
                    inner,
                    outer,
                    target: Box::new(self.parse_modifier()?),
                }
            }
            GeometricKeyword::IsoLayer => {
                let inner = self.number()?;
                let outer = self.number()?;

                GeometricExpr::IsoLayer {
                    inner,
                    outer,
                    target: Box::new(self.parse_modifier()?),
                }
            }
            GeometricKeyword::CylinderZone => {
                let radius = self.number()?;
                let z_max = self.number()?;
                let z_min = self.number()?;

                GeometricExpr::CylinderZone {
                    radius,
                    z_max,
                    z_min,
                    target: Box::new(self.parse_modifier()?),
                }
            }
            GeometricKeyword::CylinderLayer => {
                let inner = self.number()?;
                let outer = self.number()?;
                let z_max = self.number()?;
                let z_min = self.number()?;

                GeometricExpr::CylinderLayer {
                    inner,
                    outer,
                    z_max,
                    z_min,
                    target: Box::new(self.parse_modifier()?),
                }
            }
            GeometricKeyword::Point => GeometricExpr::Point {
                point: [self.number()?, self.number()?, self.number()?],
                radius: self.number()?,
            },
        };

        Ok(Some(geometric))
    }

    /// Parses a primary selection expression.
    ///
    /// Fixed keywords are recognized without lowercase allocation. Lowercase
    /// storage is created only when an uppercase identifier must be passed to
    /// macro/column lookup.
    fn parse_primary(&mut self) -> Result<Expr, Diagnostic> {
        if self.take_kind(&TokenKind::LeftParen) {
            let expr = self.parse_or()?;

            if !self.take_kind(&TokenKind::RightParen) {
                return Err(self.error_here("selection is missing a closing parenthesis"));
            }

            return Ok(expr);
        }

        let Some(value) = self.take_value() else {
            return Err(self.error_here("selection expression expected"));
        };

        if value.eq_ignore_ascii_case("all") {
            return Ok(Expr::All);
        }

        if value.eq_ignore_ascii_case("none") {
            return Ok(Expr::None);
        }

        if value.eq_ignore_ascii_case("group") {
            return Ok(Expr::Group(self.required_value("group name expected")?));
        }

        if value.eq_ignore_ascii_case("atom") {
            return self.parse_atom();
        }

        if value.eq_ignore_ascii_case("chirality") {
            return Ok(Expr::Chirality(self.required_value("chirality expected")?));
        }

        if value.eq_ignore_ascii_case("smarts") {
            let source = self.required_value("SMARTS pattern expected")?;
            let pattern = SmartsPattern::parse(&source).map_err(|error| {
                Diagnostic::new(Code::E4004).with_context("smarts", error.to_string())
            })?;
            return Ok(Expr::Smarts(pattern));
        }

        if value.eq_ignore_ascii_case("prop") {
            return self.parse_property();
        }

        let canonical = canonical_name(&value);

        if let Some(macro_name) = Macro::from_name(canonical.as_ref()) {
            return Ok(Expr::Macro(macro_name));
        }

        if let Some(column) = Column::from_name(canonical.as_ref()) {
            return self.parse_column(column);
        }

        if let Ok(number) = value.parse::<f64>() {
            return self.parse_flipped_comparison(number);
        }

        Err(Diagnostic::new(Code::E4004).with_context("keyword", value.to_string()))
    }

    /// Parses the `prop` comparison surface.
    fn parse_property(&mut self) -> Result<Expr, Diagnostic> {
        let absolute = self.take_keyword("abs");

        let Some(first) = self.take_value() else {
            return Err(self.error_here("numeric property expected"));
        };

        if let Ok(value) = first.parse::<f64>() {
            return self.parse_flipped_comparison_with(value, absolute);
        }

        let canonical = canonical_name(&first);

        let Some(column) = Column::from_name(canonical.as_ref()) else {
            return Err(Diagnostic::new(Code::E4004).with_context("property", first.to_string()));
        };

        self.comparison(column, absolute)
    }

    /// Parses a column membership or direct comparison.
    fn parse_column(&mut self, column: Column) -> Result<Expr, Diagnostic> {
        if self.peek_operator().is_some() {
            return self.comparison(column, false);
        }

        let mut values = Vec::new();

        while let Some(value) = self.peek_value() {
            if is_boundary(value) {
                break;
            }

            let Some(value) = self.take_value() else {
                break;
            };

            values.push(unescape(value));
        }

        if values.is_empty() {
            return Err(self.error_here("selection column requires at least one value"));
        }

        Ok(Expr::Membership { column, values })
    }

    /// Parses a numeric column comparison.
    fn comparison(&mut self, column: Column, absolute: bool) -> Result<Expr, Diagnostic> {
        if !column.is_numeric() {
            return Err(Diagnostic::new(Code::E4002).with_context("column", format!("{column:?}")));
        }

        let operator = self.operator()?;
        let value = f64::from(self.number()?);

        if matches!(operator, Operator::Equal | Operator::NotEqual) {
            self.warnings.push(Diagnostic::new(Code::W4002));
        }

        Ok(Expr::Comparison {
            column,
            operator,
            value,
            absolute,
        })
    }

    /// Parses a comparison whose numeric constant occurs before the property.
    fn parse_flipped_comparison(&mut self, value: f64) -> Result<Expr, Diagnostic> {
        self.parse_flipped_comparison_with(value, false)
    }

    /// Parses a flipped comparison with optional absolute-value semantics.
    fn parse_flipped_comparison_with(
        &mut self,
        value: f64,
        absolute: bool,
    ) -> Result<Expr, Diagnostic> {
        let operator = self.operator()?.flipped();
        let column_name = self.required_value("numeric property expected")?;
        let canonical = canonical_name(&column_name);

        let Some(column) = Column::from_name(canonical.as_ref()) else {
            return Err(
                Diagnostic::new(Code::E4004).with_context("property", column_name.to_string())
            );
        };

        if !column.is_numeric() {
            return Err(Diagnostic::new(Code::E4002).with_context("column", format!("{column:?}")));
        }

        if matches!(operator, Operator::Equal | Operator::NotEqual) {
            self.warnings.push(Diagnostic::new(Code::W4002));
        }

        Ok(Expr::Comparison {
            column,
            operator,
            value,
            absolute,
        })
    }

    /// Parses an exact segment/residue/atom-name selector.
    fn parse_atom(&mut self) -> Result<Expr, Diagnostic> {
        let segment = self.required_value("atom selector segment expected")?;
        let residue_text = self.required_value("atom selector residue expected")?;

        let residue = residue_text.parse::<i32>().map_err(|_| {
            Diagnostic::new(Code::E4002).with_context("value", residue_text.to_string())
        })?;

        let name = self.required_value("atom selector name expected")?;

        Ok(Expr::Atom {
            segment,
            residue,
            name,
        })
    }

    /// Parses the key following the `same` modifier.
    fn parse_same_key(&mut self) -> Result<SameKey, Diagnostic> {
        let key = self.required_value("same relation expected")?;

        if key.eq_ignore_ascii_case("residue") {
            return Ok(SameKey::Residue);
        }

        if key.eq_ignore_ascii_case("chain") {
            return Ok(SameKey::Chain);
        }

        if key.eq_ignore_ascii_case("model") {
            return Ok(SameKey::Model);
        }

        if key.eq_ignore_ascii_case("entity") {
            return Ok(SameKey::Entity);
        }

        if key.eq_ignore_ascii_case("fragment") {
            return Ok(SameKey::Fragment);
        }

        if key.eq_ignore_ascii_case("segment") {
            return Ok(SameKey::Segment);
        }

        let canonical = canonical_name(&key);

        match Column::from_name(canonical.as_ref()) {
            Some(column) => Ok(SameKey::Column(column)),
            None => Err(Diagnostic::new(Code::E4004).with_context("same key", key.to_string())),
        }
    }
}

#[cfg(test)]
#[path = "parser_tests.rs"]
mod tests;
