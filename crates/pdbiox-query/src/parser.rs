//! Precedence parser from tokens into the typed selection IR.

use crate::ast::{Column, Expr, GeometricExpr, Macro, Operator, SameKey};
use crate::lexer::{Token, TokenKind, lex};
use pdbiox_core::diagnostic::{Code, Diagnostic};

pub(super) struct Parsed {
    pub(super) expr: Expr,
    pub(super) warnings: Vec<Diagnostic>,
}

pub(super) fn parse(source: &str) -> Result<Parsed, Vec<Diagnostic>> {
    let tokens = lex(source).map_err(|finding| vec![finding])?;
    if tokens.is_empty() {
        return Err(vec![syntax(0, "selection is empty")]);
    }
    let mixed = mixes_boolean_precedence(&tokens);
    let mut parser = Parser {
        tokens,
        position: 0,
        warnings: Vec::new(),
    };
    let expr = parser.parse_or().map_err(|finding| vec![finding])?;
    if parser.position != parser.tokens.len() {
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

struct Parser {
    tokens: Vec<Token>,
    position: usize,
    warnings: Vec<Diagnostic>,
}

impl Parser {
    fn parse_or(&mut self) -> Result<Expr, Diagnostic> {
        let mut expr = self.parse_and()?;
        while self.take_keyword("or") {
            expr = Expr::Or(Box::new(expr), Box::new(self.parse_and()?));
        }
        Ok(expr)
    }

    fn parse_and(&mut self) -> Result<Expr, Diagnostic> {
        let mut expr = self.parse_not()?;
        while self.take_keyword("and") {
            expr = Expr::And(Box::new(expr), Box::new(self.parse_not()?));
        }
        Ok(expr)
    }

    fn parse_not(&mut self) -> Result<Expr, Diagnostic> {
        if self.take_keyword("not") {
            Ok(Expr::Not(Box::new(self.parse_not()?)))
        } else if self.take_keyword("global") {
            Ok(Expr::Global(Box::new(self.parse_not()?)))
        } else {
            self.parse_modifier()
        }
    }

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
                    self.position += 1;
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

    fn parse_geometric(&mut self) -> Result<Option<GeometricExpr>, Diagnostic> {
        let Some(keyword) = self.peek_value().map(str::to_ascii_lowercase) else {
            return Ok(None);
        };
        let geometric = match keyword.as_str() {
            "within" | "beyond" => {
                self.position += 1;
                let radius = self.number()?;
                self.require_keyword("of")?;
                let target = Box::new(self.parse_modifier()?);
                if keyword == "within" {
                    GeometricExpr::Within { radius, target }
                } else {
                    GeometricExpr::Beyond { radius, target }
                }
            }
            "around" | "sphzone" => {
                self.position += 1;
                let radius = self.number()?;
                let target = Box::new(self.parse_modifier()?);
                if keyword == "around" {
                    GeometricExpr::Around { radius, target }
                } else {
                    GeometricExpr::SphereZone { radius, target }
                }
            }
            "sphlayer" | "isolayer" => {
                self.position += 1;
                let inner = self.number()?;
                let outer = self.number()?;
                let target = Box::new(self.parse_modifier()?);
                if keyword == "sphlayer" {
                    GeometricExpr::SphereLayer {
                        inner,
                        outer,
                        target,
                    }
                } else {
                    GeometricExpr::IsoLayer {
                        inner,
                        outer,
                        target,
                    }
                }
            }
            "cyzone" => {
                self.position += 1;
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
            "cylayer" => {
                self.position += 1;
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
            "point" => {
                self.position += 1;
                GeometricExpr::Point {
                    point: [self.number()?, self.number()?, self.number()?],
                    radius: self.number()?,
                }
            }
            _ => return Ok(None),
        };
        Ok(Some(geometric))
    }

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
        let keyword = value.to_ascii_lowercase();
        match keyword.as_str() {
            "all" => return Ok(Expr::All),
            "none" => return Ok(Expr::None),
            "group" => {
                return Ok(Expr::Group(
                    self.required_value("group name expected")?.into(),
                ));
            }
            "atom" => return self.parse_atom(),
            "chirality" => {
                return Ok(Expr::Chirality(
                    self.required_value("chirality expected")?.into(),
                ));
            }
            "smarts" => {
                return Ok(Expr::Smarts(
                    self.required_value("SMARTS pattern expected")?.into(),
                ));
            }
            "prop" => return self.parse_property(),
            _ => {}
        }
        if let Some(macro_name) = Macro::from_name(&keyword) {
            return Ok(Expr::Macro(macro_name));
        }
        if let Some(column) = Column::from_name(&keyword) {
            return self.parse_column(column);
        }
        if let Ok(number) = value.parse::<f64>() {
            return self.parse_flipped_comparison(number);
        }
        Err(Diagnostic::new(Code::E4004).with_context("keyword", value))
    }

    fn parse_property(&mut self) -> Result<Expr, Diagnostic> {
        let absolute = self.take_keyword("abs");
        let Some(first) = self.take_value() else {
            return Err(self.error_here("numeric property expected"));
        };
        if let Ok(value) = first.parse::<f64>() {
            return self.parse_flipped_comparison_with(value, absolute);
        }
        let Some(column) = Column::from_name(&first) else {
            return Err(Diagnostic::new(Code::E4004).with_context("property", first));
        };
        self.comparison(column, absolute)
    }

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
            values.push(unescape(value).into());
        }
        if values.is_empty() {
            return Err(self.error_here("selection column requires at least one value"));
        }
        Ok(Expr::Membership { column, values })
    }

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

    fn parse_flipped_comparison(&mut self, value: f64) -> Result<Expr, Diagnostic> {
        self.parse_flipped_comparison_with(value, false)
    }

    fn parse_flipped_comparison_with(
        &mut self,
        value: f64,
        absolute: bool,
    ) -> Result<Expr, Diagnostic> {
        let operator = self.operator()?.flipped();
        let column_name = self.required_value("numeric property expected")?;
        let Some(column) = Column::from_name(&column_name) else {
            return Err(Diagnostic::new(Code::E4004).with_context("property", column_name));
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

    fn parse_atom(&mut self) -> Result<Expr, Diagnostic> {
        let segment = self.required_value("atom selector segment expected")?;
        let residue = self.required_value("atom selector residue expected")?;
        let residue = residue
            .parse::<i32>()
            .map_err(|_| Diagnostic::new(Code::E4002).with_context("value", residue))?;
        let name = self.required_value("atom selector name expected")?;
        Ok(Expr::Atom {
            segment: segment.into(),
            residue,
            name: name.into(),
        })
    }

    fn parse_same_key(&mut self) -> Result<SameKey, Diagnostic> {
        let key = self.required_value("same relation expected")?;
        Ok(match key.to_ascii_lowercase().as_str() {
            "residue" => SameKey::Residue,
            "chain" => SameKey::Chain,
            "model" => SameKey::Model,
            "entity" => SameKey::Entity,
            "fragment" => SameKey::Fragment,
            "segment" => SameKey::Segment,
            _ => match Column::from_name(&key) {
                Some(column) => SameKey::Column(column),
                None => return Err(Diagnostic::new(Code::E4004).with_context("same key", key)),
            },
        })
    }

    fn operator(&mut self) -> Result<Operator, Diagnostic> {
        let Some(operator) = self.take_operator() else {
            return Err(self.error_here("comparison operator expected"));
        };
        match operator.as_str() {
            "<" => Ok(Operator::Less),
            "<=" => Ok(Operator::LessEqual),
            ">" => Ok(Operator::Greater),
            ">=" => Ok(Operator::GreaterEqual),
            "==" => Ok(Operator::Equal),
            "!=" => Ok(Operator::NotEqual),
            _ => Err(Diagnostic::new(Code::E4001).with_context("operator", operator)),
        }
    }

    fn number(&mut self) -> Result<f32, Diagnostic> {
        let value = self.required_value("number expected")?;
        value
            .parse::<f32>()
            .map_err(|_| Diagnostic::new(Code::E4002).with_context("value", value))
    }

    fn require_keyword(&mut self, keyword: &str) -> Result<(), Diagnostic> {
        if self.take_keyword(keyword) {
            Ok(())
        } else {
            Err(self.error_here("required selection keyword is missing"))
        }
    }

    fn required_value(&mut self, message: &'static str) -> Result<String, Diagnostic> {
        self.take_value().ok_or_else(|| self.error_here(message))
    }

    fn peek_value(&self) -> Option<&str> {
        match self.tokens.get(self.position).map(|token| &token.kind) {
            Some(TokenKind::Value(value)) => Some(value),
            _ => None,
        }
    }

    fn take_value(&mut self) -> Option<String> {
        let value = self.peek_value()?.to_owned();
        self.position += 1;
        Some(value)
    }

    fn peek_operator(&self) -> Option<&str> {
        match self.tokens.get(self.position).map(|token| &token.kind) {
            Some(TokenKind::Operator(operator)) => Some(operator),
            _ => None,
        }
    }

    fn take_operator(&mut self) -> Option<String> {
        let operator = self.peek_operator()?.to_owned();
        self.position += 1;
        Some(operator)
    }

    fn take_keyword(&mut self, expected: &str) -> bool {
        if self
            .peek_value()
            .is_some_and(|value| value.eq_ignore_ascii_case(expected))
        {
            self.position += 1;
            true
        } else {
            false
        }
    }

    fn take_kind(&mut self, expected: &TokenKind) -> bool {
        if self
            .tokens
            .get(self.position)
            .is_some_and(|token| &token.kind == expected)
        {
            self.position += 1;
            true
        } else {
            false
        }
    }

    fn error_here(&self, message: &'static str) -> Diagnostic {
        let offset = self
            .tokens
            .get(self.position)
            .map_or(0, |token| token.offset);
        syntax(offset, message)
    }
}

fn mixes_boolean_precedence(tokens: &[Token]) -> bool {
    let mut depth = 0u32;
    let mut saw_and = false;
    let mut saw_or = false;
    for token in tokens {
        match &token.kind {
            TokenKind::LeftParen => depth += 1,
            TokenKind::RightParen => depth = depth.saturating_sub(1),
            TokenKind::Value(value) if depth == 0 && value.eq_ignore_ascii_case("and") => {
                saw_and = true;
            }
            TokenKind::Value(value) if depth == 0 && value.eq_ignore_ascii_case("or") => {
                saw_or = true;
            }
            _ => {}
        }
    }
    saw_and && saw_or
}

fn is_boundary(value: &str) -> bool {
    value.eq_ignore_ascii_case("and") || value.eq_ignore_ascii_case("or")
}

fn unescape(mut value: String) -> String {
    if value.starts_with('\\') {
        value.remove(0);
    }
    value
}

fn syntax(offset: usize, message: &'static str) -> Diagnostic {
    Diagnostic::new(Code::E4001)
        .with_message(message)
        .with_context("offset", offset.to_string())
}

#[cfg(test)]
#[path = "parser_tests.rs"]
mod tests;
