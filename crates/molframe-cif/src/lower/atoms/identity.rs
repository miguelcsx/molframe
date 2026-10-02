//! Residue identity: component codes, residue keys and residue opening.

use super::{AtomBuilder, AtomSiteRow, Field};
use crate::lower::diagnostics::at_source_row;
use crate::lower::keys::ResidueKey;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::index::ResidueIndex;
use molframe_core::optional::{OptionalI32, OptionalSymbol};
use molframe_core::symbol::SymbolId;
use molframe_core::topology::ResidueRecord;

impl AtomBuilder<'_> {
    /// Reports an alternate location that carries a different component code.
    ///
    /// This is a modelled point mutation: one residue with two chemical
    /// identities. Both are kept; the finding says the file did something worth
    /// knowing about rather than something wrong.
    pub(super) fn alternate_component_of<R: AtomSiteRow + ?Sized>(
        &mut self,
        rows: &R,
        residue: ResidueIndex,
    ) -> OptionalSymbol {
        let Some(existing_symbol) = self.data.topology.residues.label_comp_id(residue) else {
            return OptionalSymbol::NONE;
        };
        // Interned text is unique, so for a dictionary reader equal symbols
        // mean equal components without resolving either string.
        if rows.dictionary_slot(Field::LabelCompId).is_some()
            && self.symbol_of(rows, Field::LabelCompId) == Some(existing_symbol)
        {
            return OptionalSymbol::NONE;
        }
        let Some(comp) = rows.identifier(Field::LabelCompId) else {
            return OptionalSymbol::NONE;
        };
        let Some(existing) = self.data.dictionary.resolve(existing_symbol) else {
            return OptionalSymbol::NONE;
        };
        if existing == comp.as_ref() {
            return OptionalSymbol::NONE;
        }
        self.findings.push(at_source_row(
            Diagnostic::new(Code::W3012)
                .in_category("atom_site")
                .with_context("component", existing.to_owned())
                .with_context("alternate component", comp.to_string()),
            rows.row(),
        ));
        OptionalSymbol::some(self.intern(&comp))
    }

    pub(super) fn open_residue<R: AtomSiteRow + ?Sized>(&mut self, rows: &R, key: &ResidueKey) {
        self.close_residue();
        if self.chain != Some(key.chain) {
            self.open_chain(rows, key.chain);
        }
        self.current = Some(*key);
        self.names_in_residue.clear();

        let comp = self.symbol_or_empty(rows, Field::LabelCompId);
        let auth_comp = self
            .symbol_of(rows, Field::AuthCompId)
            .map_or(OptionalSymbol::NONE, OptionalSymbol::some);
        let ins_code = self
            .symbol_of(rows, Field::InsCode)
            .map_or(OptionalSymbol::NONE, OptionalSymbol::some);
        let het = rows.text(Field::GroupPdb) == Some("HETATM");

        if self
            .data
            .topology
            .residues
            .push(
                ResidueRecord {
                    label_comp_id: comp,
                    auth_comp_id: auth_comp,
                    label_seq_id: key.label_seq,
                    auth_seq_id: key.auth_seq,
                    ins_code,
                    het,
                },
                self.atom_position..self.atom_position,
            )
            .is_err()
        {
            self.findings.push(Diagnostic::new(Code::E3001));
            return;
        }
        self.residue_position += 1;
    }

    pub(super) fn key_of<R: AtomSiteRow + ?Sized>(&mut self, rows: &R, model: i64) -> ResidueKey {
        let chain = self
            .symbol_of(rows, Field::LabelAsymId)
            .map_or(ResidueKey::ABSENT, SymbolId::get);
        let ins_code = self
            .symbol_of(rows, Field::InsCode)
            .map_or(ResidueKey::ABSENT, SymbolId::get);
        ResidueKey {
            model,
            chain,
            label_seq: OptionalI32::from(
                rows.integer(Field::LabelSeqId)
                    .and_then(|seq| i32::try_from(seq).ok()),
            ),
            auth_seq: OptionalI32::from(
                rows.integer(Field::AuthSeqId)
                    .and_then(|seq| i32::try_from(seq).ok()),
            ),
            ins_code,
        }
    }
}
