//! Allocation-conscious lexical and syntax support for the selection parser.

use super::GeometricKeyword;
use crate::ast::{Column, Expr};
use crate::lexer::{Token, TokenKind};
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
pub(crate) fn is_boundary(value: &str) -> bool {
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

/// Whether `name` can follow `$` as a named query reference.
///
/// A name starts with a letter or underscore and continues with letters,
/// digits and underscores: the same spelling rule an interactive session
/// applies to the names it defines, so every definable name is referable.
pub(crate) fn is_query_name(name: &str) -> bool {
    let mut bytes = name.bytes();
    let Some(first) = bytes.next() else {
        return false;
    };
    (first.is_ascii_alphabetic() || first == b'_')
        && bytes.all(|byte| byte.is_ascii_alphanumeric() || byte == b'_')
}

/// Builds a membership, honouring `label:` and `auth:` value prefixes.
///
/// `chain`, `resname` and `name` follow the policy's identifier namespace. A
/// value written `label:X` or `auth:X` names its namespace itself, so
/// `chain label:A` and `label_chain A` select the same atoms. Values of one
/// column may mix prefixes; the result is the union of the groups.
pub(super) fn namespaced_membership(column: Column, values: Vec<Box<str>>) -> Expr {
    let Some((label, auth)) = namespaced_columns(column) else {
        return Expr::Membership { column, values };
    };
    let mut plain = Vec::new();
    let mut labelled = Vec::new();
    let mut authored = Vec::new();
    for value in values {
        if let Some(rest) = value.strip_prefix("label:") {
            labelled.push(Box::from(rest));
        } else if let Some(rest) = value.strip_prefix("auth:") {
            authored.push(Box::from(rest));
        } else {
            plain.push(value);
        }
    }
    let mut groups = Vec::new();
    for (column, values) in [(column, plain), (label, labelled), (auth, authored)] {
        if !values.is_empty() {
            groups.push(Expr::Membership { column, values });
        }
    }
    let mut groups = groups.into_iter();
    let Some(first) = groups.next() else {
        return Expr::None;
    };
    groups.fold(first, |left, right| {
        Expr::Or(Box::new(left), Box::new(right))
    })
}

/// The label and author columns behind a policy-resolved column.
const fn namespaced_columns(column: Column) -> Option<(Column, Column)> {
    match column {
        Column::Chain => Some((Column::LabelChain, Column::AuthChain)),
        Column::ResidueName => Some((Column::LabelResidueName, Column::AuthResidueName)),
        Column::AtomName => Some((Column::LabelAtomName, Column::AuthAtomName)),
        _ => None,
    }
}
