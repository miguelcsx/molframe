//! Token-cursor operations shared by parser productions.

use super::Parser;
use crate::ast::Operator;
use crate::lexer::{Token, TokenKind};
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::span::ByteSpan;

impl Parser {
    pub(super) fn operator(&mut self) -> Result<Operator, Diagnostic> {
        let Some(operator) = self.take_operator() else {
            return Err(self.error_here("comparison operator expected"));
        };
        match operator.as_ref() {
            "<" => Ok(Operator::Less),
            "<=" => Ok(Operator::LessEqual),
            ">" => Ok(Operator::Greater),
            ">=" => Ok(Operator::GreaterEqual),
            "==" => Ok(Operator::Equal),
            "!=" => Ok(Operator::NotEqual),
            _ => Err(Diagnostic::new(Code::E4001)
                .with_context("operator", operator.to_string())
                .at(self.last)),
        }
    }

    pub(super) fn number(&mut self) -> Result<f32, Diagnostic> {
        let value = self.required_value("number expected")?;
        value.parse::<f32>().map_err(|_| {
            Diagnostic::new(Code::E4002)
                .with_context("value", value.to_string())
                .at(self.last)
        })
    }

    pub(super) fn require_keyword(&mut self, keyword: &str) -> Result<(), Diagnostic> {
        if self.take_keyword(keyword) {
            Ok(())
        } else {
            Err(self.error_here("required selection keyword is missing"))
        }
    }

    pub(super) fn required_value(&mut self, message: &'static str) -> Result<Box<str>, Diagnostic> {
        match self.take_value() {
            Some(value) => Ok(value),
            None => Err(self.error_here(message)),
        }
    }

    pub(super) fn peek_value(&self) -> Option<&str> {
        match self.tokens.as_slice().first().map(|token| &token.kind) {
            Some(TokenKind::Value(value)) => Some(value),
            _ => None,
        }
    }

    pub(super) fn take_value(&mut self) -> Option<Box<str>> {
        if !matches!(
            self.tokens.as_slice().first().map(|token| &token.kind),
            Some(TokenKind::Value(_))
        ) {
            return None;
        }
        match self.advance() {
            Some(Token {
                kind: TokenKind::Value(value),
                ..
            }) => Some(value),
            _ => None,
        }
    }

    pub(super) fn peek_operator(&self) -> Option<&str> {
        match self.tokens.as_slice().first().map(|token| &token.kind) {
            Some(TokenKind::Operator(operator)) => Some(operator),
            _ => None,
        }
    }

    fn take_operator(&mut self) -> Option<Box<str>> {
        if !matches!(
            self.tokens.as_slice().first().map(|token| &token.kind),
            Some(TokenKind::Operator(_))
        ) {
            return None;
        }
        match self.advance() {
            Some(Token {
                kind: TokenKind::Operator(operator),
                ..
            }) => Some(operator),
            _ => None,
        }
    }

    pub(super) fn take_keyword(&mut self, expected: &str) -> bool {
        if self
            .peek_value()
            .is_some_and(|value| value.eq_ignore_ascii_case(expected))
        {
            let _ = self.advance();
            true
        } else {
            false
        }
    }

    pub(super) fn take_kind(&mut self, expected: &TokenKind) -> bool {
        if self
            .tokens
            .as_slice()
            .first()
            .is_some_and(|token| &token.kind == expected)
        {
            let _ = self.advance();
            true
        } else {
            false
        }
    }

    /// Consumes the next token, remembering where it was read from.
    ///
    /// Every production consumes through here, so a diagnostic about the value
    /// just read can point at exactly that value.
    pub(super) fn advance(&mut self) -> Option<Token> {
        let token = self.tokens.next()?;
        self.last = token.span;
        Some(token)
    }

    /// The span of the next unread token, or the end of the source.
    pub(super) fn next_span(&self) -> ByteSpan {
        match self.tokens.as_slice().first() {
            Some(token) => token.span,
            None => self.end,
        }
    }

    pub(super) fn error_here(&self, message: &'static str) -> Diagnostic {
        crate::lexer::syntax(self.next_span(), message)
    }
}
