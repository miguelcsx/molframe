//! Reading one coordinate row's fields, and interning what they name.
//!
//! Separated from the lowering itself so that the row loop reads as the
//! sequence of decisions it makes, rather than as those decisions interleaved
//! with the mechanics of getting a value out of a row.

use super::{AtomBuilder, AtomSiteRow, Field};
use crate::lower::diagnostics::at_source_row;
use molframe_core::column::Presence;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::element::Element;
use molframe_core::io::MissingElementPolicy;
use molframe_core::optional::OptionalSymbol;
use molframe_core::symbol::{AltId, SymbolId};
use num_traits::ToPrimitive;

/// Marks a cache slot that has not been filled; no real symbol takes this value.
pub(super) const UNSET_SYMBOL: SymbolId = SymbolId::from_raw(u32::MAX);

impl AtomBuilder<'_> {
    pub(super) fn position_of<R: AtomSiteRow + ?Sized>(&mut self, rows: &R) -> Option<[f32; 3]> {
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

    pub(super) fn optional_float<R: AtomSiteRow + ?Sized>(
        &mut self,
        rows: &R,
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

    pub(super) fn element_of<R: AtomSiteRow + ?Sized>(&mut self, rows: &R) -> Element {
        if let Some(element) = rows.text(Field::TypeSymbol).and_then(Element::from_symbol) {
            return element;
        }
        let name = rows.identifier(Field::LabelAtomId);
        let inferred = match (self.options.missing_element_policy, name.as_deref()) {
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

    pub(super) fn auth_name_of<R: AtomSiteRow + ?Sized>(&mut self, rows: &R) -> OptionalSymbol {
        self.symbol_of(rows, Field::AuthAtomId)
            .map_or(OptionalSymbol::NONE, OptionalSymbol::some)
    }

    pub(super) fn alt_of<R: AtomSiteRow + ?Sized>(&mut self, rows: &R) -> Option<AltId> {
        if let Some(text) = rows.identifier(Field::LabelAltId)
            && !text.is_empty()
        {
            return AltId::labelled(self.intern(&text));
        }
        Some(AltId::BLANK)
    }

    /// The symbol for a text field of this row, interning it on first sight.
    ///
    /// Dictionary-encoded readers repeat the same few strings across every
    /// row; the remembered symbol skips the two hash lookups an intern costs.
    pub(super) fn symbol_of<R: AtomSiteRow + ?Sized>(
        &mut self,
        rows: &R,
        field: Field,
    ) -> Option<SymbolId> {
        let slot = rows
            .dictionary_slot(field)
            .and_then(|slot| usize::try_from(slot).ok());
        if let Some(slot) = slot
            && let Some(symbol) = self.symbol_cache[field.position()].get(slot)
            && *symbol != UNSET_SYMBOL
        {
            return Some(*symbol);
        }
        let text = rows.identifier(field)?;
        let symbol = self.intern(&text);
        if let Some(slot) = slot {
            let cache = &mut self.symbol_cache[field.position()];
            if cache.len() <= slot {
                cache.resize(slot + 1, UNSET_SYMBOL);
            }
            cache[slot] = symbol;
        }
        Some(symbol)
    }

    /// Like [`Self::symbol_of`], with an absent value read as the empty string.
    pub(super) fn symbol_or_empty<R: AtomSiteRow + ?Sized>(
        &mut self,
        rows: &R,
        field: Field,
    ) -> SymbolId {
        match self.symbol_of(rows, field) {
            Some(symbol) => symbol,
            None => self.intern(""),
        }
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
