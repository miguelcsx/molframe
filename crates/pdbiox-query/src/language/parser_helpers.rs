//! Allocation-conscious lexical helpers for the selection parser.

use super::GeometricKeyword;
use crate::lexer::{Token, TokenKind};
use pdbiox_core::diagnostic::{Code, Diagnostic};
use std::borrow::Cow;

pub(super) fn geometric_keyword(value: &str) -> Option<GeometricKeyword> {
    if value.eq_ignore_ascii_case("within") {
        Some(GeometricKeyword::Within)
    } else if value.eq_ignore_ascii_case("beyond") {
        Some(GeometricKeyword::Beyond)
    } else if value.eq_ignore_ascii_case("around") {
        Some(GeometricKeyword::Around)
    } else if value.eq_ignore_ascii_case("sphzone") {
        Some(GeometricKeyword::SphereZone)
    } else if value.eq_ignore_ascii_case("sphlayer") {
        Some(GeometricKeyword::SphereLayer)
    } else if value.eq_ignore_ascii_case("isolayer") {
        Some(GeometricKeyword::IsoLayer)
    } else if value.eq_ignore_ascii_case("cyzone") {
        Some(GeometricKeyword::CylinderZone)
    } else if value.eq_ignore_ascii_case("cylayer") {
        Some(GeometricKeyword::CylinderLayer)
    } else if value.eq_ignore_ascii_case("point") {
        Some(GeometricKeyword::Point)
    } else {
        None
    }
}

pub(super) fn canonical_name(value: &str) -> Cow<'_, str> {
    if value.bytes().any(|byte| byte.is_ascii_uppercase()) {
        Cow::Owned(value.to_ascii_lowercase())
    } else {
        Cow::Borrowed(value)
    }
}

pub(super) fn mixes_boolean_precedence(tokens: &[Token]) -> bool {
    let mut depth = 0usize;
    let mut saw_and = false;
    let mut saw_or = false;
    for token in tokens {
        match &token.kind {
            TokenKind::LeftParen => depth += 1,
            TokenKind::RightParen => {
                if let Some(parent_depth) = depth.checked_sub(1) {
                    depth = parent_depth;
                }
            }
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

#[inline]
pub(super) fn is_boundary(value: &str) -> bool {
    value.eq_ignore_ascii_case("and") || value.eq_ignore_ascii_case("or")
}

pub(super) fn unescape(value: Box<str>) -> Box<str> {
    if !value.starts_with('\\') {
        return value;
    }
    match value.get(1..) {
        Some(unescaped) => unescaped.into(),
        None => "".into(),
    }
}

pub(super) fn syntax(offset: usize, message: &'static str) -> Diagnostic {
    Diagnostic::new(Code::E4001)
        .with_message(message)
        .with_context("offset", offset.to_string())
}
