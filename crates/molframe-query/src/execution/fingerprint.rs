//! Stable semantic identity for normalized query expressions.

use crate::ast::{Expr, GeometricExpr};
use std::fmt;
use std::hash::{Hash, Hasher};

const OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const PRIME: u64 = 0x0000_0100_0000_01b3;

/// Stable identity of a normalized molecular query.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct QueryFingerprint(u64);

impl QueryFingerprint {
    pub(super) fn of(expression: &Expr) -> Self {
        let mut encoder = Encoder(OFFSET);
        encoder.expression(expression);
        Self(encoder.finish())
    }

    /// Raw fingerprint value.
    #[must_use]
    pub const fn get(self) -> u64 {
        self.0
    }
}

impl fmt::Display for QueryFingerprint {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "query-fnv1a64:{:016x}", self.0)
    }
}

impl fmt::Debug for QueryFingerprint {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self, formatter)
    }
}

struct Encoder(u64);

impl Encoder {
    fn expression(&mut self, expression: &Expr) {
        match expression {
            Expr::All => self.tag(0),
            Expr::None => self.tag(1),
            Expr::And(left, right) => self.binary(2, left, right),
            Expr::Or(left, right) => self.binary(3, left, right),
            Expr::Not(value) => self.unary(4, value),
            Expr::Global(value) => self.unary(5, value),
            Expr::ByResidue(value) => self.unary(6, value),
            Expr::Same { key, target } => {
                self.tag(7);
                key.hash(self);
                self.expression(target);
            }
            Expr::Bonded { depth, target } => {
                self.tag(8);
                self.write_u32(*depth);
                self.expression(target);
            }
            Expr::Geometric(value) => self.geometric(value),
            Expr::Comparison {
                column,
                operator,
                value,
                absolute,
            } => {
                self.tag(10);
                column.hash(self);
                operator.hash(self);
                self.write_u64(value.to_bits());
                self.write_u8(u8::from(*absolute));
            }
            Expr::Membership { column, values } => {
                self.tag(11);
                column.hash(self);
                self.write_usize(values.len());
                for value in values {
                    self.text(value);
                }
            }
            Expr::Group(value) => self.tagged_text(12, value),
            Expr::Atom {
                segment,
                residue,
                name,
            } => {
                self.tag(13);
                self.text(segment);
                self.write_i32(*residue);
                self.text(name);
            }
            Expr::Macro(value) => {
                self.tag(14);
                value.hash(self);
            }
            Expr::Chirality(value) => self.tagged_text(15, value),
            Expr::Smarts(value) => {
                self.tag(16);
                value.hash(self);
            }
        }
    }

    fn geometric(&mut self, expression: &GeometricExpr) {
        self.tag(9);
        match expression {
            GeometricExpr::Within { radius, target } => self.radius_target(0, *radius, target),
            GeometricExpr::Beyond { radius, target } => self.radius_target(1, *radius, target),
            GeometricExpr::Around { radius, target } => self.radius_target(2, *radius, target),
            GeometricExpr::SphereZone { radius, target } => {
                self.radius_target(3, *radius, target);
            }
            GeometricExpr::SphereLayer {
                inner,
                outer,
                target,
            } => self.layer_target(4, [*inner, *outer], target),
            GeometricExpr::IsoLayer {
                inner,
                outer,
                target,
            } => self.layer_target(5, [*inner, *outer], target),
            GeometricExpr::CylinderZone {
                radius,
                z_max,
                z_min,
                target,
            } => self.layer_target(6, [*radius, *z_max, *z_min], target),
            GeometricExpr::CylinderLayer {
                inner,
                outer,
                z_max,
                z_min,
                target,
            } => self.layer_target(7, [*inner, *outer, *z_max, *z_min], target),
            GeometricExpr::Point { point, radius } => {
                self.tag(8);
                self.floats(point);
                self.write_u32(radius.to_bits());
            }
        }
    }

    fn unary(&mut self, tag: u8, value: &Expr) {
        self.tag(tag);
        self.expression(value);
    }

    fn binary(&mut self, tag: u8, left: &Expr, right: &Expr) {
        self.tag(tag);
        self.expression(left);
        self.expression(right);
    }

    fn radius_target(&mut self, tag: u8, radius: f32, target: &Expr) {
        self.tag(tag);
        self.write_u32(radius.to_bits());
        self.expression(target);
    }

    fn layer_target<const N: usize>(&mut self, tag: u8, values: [f32; N], target: &Expr) {
        self.tag(tag);
        self.floats(&values);
        self.expression(target);
    }

    fn floats(&mut self, values: &[f32]) {
        for value in values {
            self.write_u32(value.to_bits());
        }
    }

    fn tagged_text(&mut self, tag: u8, value: &str) {
        self.tag(tag);
        self.text(value);
    }

    fn text(&mut self, value: &str) {
        self.write_usize(value.len());
        self.write(value.as_bytes());
    }

    fn tag(&mut self, tag: u8) {
        self.write_u8(tag);
    }
}

impl Hasher for Encoder {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for byte in bytes {
            self.0 ^= u64::from(*byte);
            self.0 = self.0.wrapping_mul(PRIME);
        }
    }
}
