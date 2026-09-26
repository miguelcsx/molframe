//! Geometric selection productions: `within`, `around`, zones and layers.

use super::helpers::geometric_keyword;
use super::{GeometricKeyword, Parser};
use crate::ast::GeometricExpr;
use molframe_core::diagnostic::Diagnostic;

impl Parser {
    /// Parses a geometric prefix if the next value is a geometric keyword.
    ///
    /// Keyword recognition is ASCII case-insensitive and allocation-free.
    pub(super) fn parse_geometric(&mut self) -> Result<Option<GeometricExpr>, Diagnostic> {
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
}
