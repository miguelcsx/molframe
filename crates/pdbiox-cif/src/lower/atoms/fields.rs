//! Reading one coordinate row's fields, and interning what they name.
//!
//! Separated from the lowering itself so that the row loop reads as the
//! sequence of decisions it makes, rather than as those decisions interleaved
//! with the mechanics of getting a value out of a row.

use super::{AtomBuilder, AtomSiteRow, Field};
use crate::lower::diagnostics::at_source_row;
use num_traits::ToPrimitive;
use pdbiox_core::column::Presence;
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::element::Element;
use pdbiox_core::io::MissingElementPolicy;
use pdbiox_core::optional::OptionalSymbol;
use pdbiox_core::symbol::{AltId, SymbolId};

impl AtomBuilder<'_> {
    pub(super) fn position_of(&mut self, rows: &dyn AtomSiteRow) -> Option<[f32; 3]> {
        let (Some(x), Some(y), Some(z)) = (
            rows.float(Field::CartnX),
            rows.float(Field::CartnY),
            rows.float(Field::CartnZ),
        ) else {
            self.findings.push(at_source_row(
                Diagnostic::new(Code::E1202)
                    .with_message("a coordinate could not be read")
                    .in_category("atom_site"),
                rows.row(),
            ));
            return None;
        };
        let Some(position) = x
            .to_f32()
            .zip(y.to_f32())
            .zip(z.to_f32())
            .map(|((x, y), z)| [x, y, z])
        else {
            self.findings.push(at_source_row(
                Diagnostic::new(Code::E1202)
                    .with_message("a coordinate is outside the supported floating-point range")
                    .in_category("atom_site"),
                rows.row(),
            ));
            return None;
        };
        Some(position)
    }

    pub(super) fn optional_float(
        &mut self,
        rows: &dyn AtomSiteRow,
        field: Field,
        when_absent: f32,
    ) -> (f32, Presence) {
        let Some(value) = rows.float(field) else {
            return (when_absent, Presence::Unknown);
        };
        if let Some(value) = value.to_f32() {
            return (value, Presence::Present);
        }
        self.findings.push(at_source_row(
            Diagnostic::new(Code::E1202)
                .with_message("an atom value is outside the supported floating-point range")
                .in_category("atom_site")
                .with_context("item", field.item()),
            rows.row(),
        ));
        (when_absent, Presence::Unknown)
    }

    pub(super) fn element_of(&mut self, rows: &dyn AtomSiteRow, name: Option<&str>) -> Element {
        if let Some(element) = rows.text(Field::TypeSymbol).and_then(Element::from_symbol) {
            return element;
        }
        let inferred = match (self.options.missing_element_policy, name) {
            (MissingElementPolicy::InferFromAtomName, Some(name)) => Element::infer_from_name(name),
            _ => Element::UNKNOWN,
        };
        self.findings.push(at_source_row(
            Diagnostic::new(Code::W3203)
                .in_category("atom_site")
                .with_context("inferred", inferred.symbol()),
            rows.row(),
        ));
        inferred
    }

    pub(super) fn auth_name_of(
        &mut self,
        rows: &dyn AtomSiteRow,
        label: Option<&str>,
    ) -> OptionalSymbol {
        if let Some(text) = rows.identifier(Field::AuthAtomId)
            && Some(text.as_ref()) != label
        {
            return OptionalSymbol::some(self.intern(&text));
        }
        OptionalSymbol::NONE
    }

    pub(super) fn alt_of(&mut self, rows: &dyn AtomSiteRow) -> Option<AltId> {
        if let Some(text) = rows.identifier(Field::LabelAltId)
            && !text.is_empty()
        {
            return AltId::labelled(self.intern(&text));
        }
        Some(AltId::BLANK)
    }

    pub(super) fn intern(&mut self, text: &str) -> SymbolId {
        let Ok(symbol) = self.data.dictionary.intern(text) else {
            self.findings
                .push(Diagnostic::new(Code::E1901).with_message("the dictionary is full"));
            return SymbolId::from_raw(0);
        };
        symbol
    }
}
