//! Turning coordinate rows into a structure.
//!
//! This is where interpretation happens, and therefore where most findings come
//! from. Every decision the file does not force — where a residue ends, which
//! element an atom is, whether two rows are the same residue modelled twice — is
//! made here and reported when it was not forced.

use super::keys::{Boundary, ResidueKey, boundary};
use crate::document::Category;
use crate::parser::Rows;
use pdbiox_core::chunk::{AtomRecord, ChunkBuilder};
use pdbiox_core::column::Presence;
use pdbiox_core::coords::CoordinateBlock;
use pdbiox_core::diagnostic::{Code, Diagnostic, Diagnostics};
use pdbiox_core::element::Element;
use pdbiox_core::index::{EntityIndex, ResidueIndex};
use pdbiox_core::io::ReadOptions;
use pdbiox_core::optional::{OptionalI32, OptionalSymbol};
use pdbiox_core::structure::{CoordinateStore, StructureData};
use pdbiox_core::symbol::{AltId, SymbolId};
use pdbiox_core::topology::{ChainRecord, PolymerKind, ResidueRecord};

/// Builds the atoms, residues, chains and models of a structure.
pub struct AtomBuilder<'a> {
    data: &'a mut StructureData,
    findings: &'a mut Diagnostics,
    options: &'a ReadOptions,
    builder: ChunkBuilder,
    frames: Vec<CoordinateBlock>,
    /// The residue being filled, and the atom names already in it.
    current: Option<ResidueKey>,
    atom_name_epochs: Vec<u32>,
    residue_epoch: u32,
    /// The chain being filled.
    chain: Option<u32>,
    chain_first_residue: u32,
    residue_position: u32,
    atom_position: u32,
    model: i64,
    topology_locked: bool,
    reported_fallback: bool,
}

impl<'a> AtomBuilder<'a> {
    /// Starts building into `data`.
    pub fn new(
        data: &'a mut StructureData,
        findings: &'a mut Diagnostics,
        options: &'a ReadOptions,
    ) -> Self {
        Self {
            data,
            findings,
            options,
            builder: ChunkBuilder::new(),
            frames: Vec::new(),
            current: None,
            atom_name_epochs: Vec::new(),
            residue_epoch: 0,
            chain: None,
            chain_first_residue: 0,
            residue_position: 0,
            atom_position: 0,
            model: i64::MIN,
            topology_locked: false,
            reported_fallback: false,
        }
    }

    /// Reads every row of the coordinate category.
    pub fn read(mut self, category: &Category) -> CoordinateStore {
        let mut rows = Rows::new(category);

        for _ in 0..category.row_count() {
            self.row(&rows);

            if !rows.advance() {
                break;
            }
        }

        self.finish()
    }

    fn row(&mut self, rows: &Rows<'_>) {
        let model = self.model_of(rows);

        if !self.activate_model(model) {
            return;
        }

        let name_text = rows
            .identifier("label_atom_id")
            .or_else(|| rows.identifier("auth_atom_id"));

        let element = self.element_of(rows, name_text.as_deref());

        if self.options.discard_hydrogens && element.is_hydrogen() {
            return;
        }

        let atom_name = self.intern(or_empty(name_text.as_deref()));
        let key = self.key_of(rows, model);
        let residue = self.place(rows, &key, atom_name);

        self.push_atom(rows, name_text.as_deref(), element, atom_name, residue);
    }

    fn model_of(&self, rows: &Rows<'_>) -> i64 {
        // A file that does not number its models has exactly one.
        match rows.integer("pdbx_PDB_model_num") {
            Some(model) => model,
            None => self.model.max(1),
        }
    }

    fn activate_model(&mut self, model: i64) -> bool {
        if model == self.model {
            return !self.options.only_first_model || self.frames.is_empty();
        }

        if self.options.only_first_model && !self.frames.is_empty() {
            return false;
        }

        self.start_model(model);

        !self.options.only_first_model || self.frames.is_empty()
    }

    fn push_atom(
        &mut self,
        rows: &Rows<'_>,
        label: Option<&str>,
        element: Element,
        atom_name: SymbolId,
        residue: ResidueIndex,
    ) {
        let position = self.position_of(rows);
        let auth_atom_name = self.auth_name_of(rows, label);
        let alt_id = self.alt_of(rows);

        let record = AtomRecord {
            position,
            element,
            atom_name,
            auth_atom_name,
            alt_id,
            residue,
            occupancy: optional(rows.float("occupancy"), 1.0),
            b_factor: optional(rows.float("B_iso_or_equiv"), 0.0),
            formal_charge: formal_charge_of(rows),
            atom_site_id: atom_site_id_of(rows),
        };

        self.builder.push(record);
        self.atom_position += 1;
    }

    /// Where this atom's residue sits, opening a new one if the row starts one.
    fn place(&mut self, rows: &Rows<'_>, key: &ResidueKey, atom_name: SymbolId) -> ResidueIndex {
        if self.topology_locked {
            return self.locked_residue();
        }

        let repeats = self.atom_name_seen(atom_name);
        let decision = self.boundary_of(key, repeats);

        self.apply_boundary(rows, key, decision);
        self.remember_atom_name(atom_name);

        self.current_residue()
    }

    fn locked_residue(&self) -> ResidueIndex {
        match self.data.topology.residues.containing(self.atom_position) {
            Some(residue) => residue,
            None => self.current_residue(),
        }
    }

    fn boundary_of(&self, key: &ResidueKey, repeats: bool) -> Boundary {
        match self.current.as_ref() {
            Some(current) => boundary(current, key, repeats),
            None => Boundary::New,
        }
    }

    fn apply_boundary(&mut self, rows: &Rows<'_>, key: &ResidueKey, decision: Boundary) {
        match decision {
            Boundary::New => self.open_residue(rows, key),
            Boundary::NewByFileOrder => {
                self.report_file_order_fallback(rows);
                self.open_residue(rows, key);
            }
            _ => self.check_component_identity(rows),
        }
    }

    fn report_file_order_fallback(&mut self, rows: &Rows<'_>) {
        if self.reported_fallback {
            return;
        }

        self.reported_fallback = true;
        self.findings.push(
            Diagnostic::new(Code::W3011)
                .in_category("atom_site")
                .at_row(rows.row() as u32),
        );
    }

    fn atom_name_seen(&self, atom_name: SymbolId) -> bool {
        let Some(index) = atom_name_index(atom_name) else {
            return false;
        };

        self.atom_name_epochs.get(index).copied() == Some(self.residue_epoch)
    }

    fn remember_atom_name(&mut self, atom_name: SymbolId) {
        let Some(index) = atom_name_index(atom_name) else {
            return;
        };

        let Some(required_len) = index.checked_add(1) else {
            return;
        };

        if self.atom_name_epochs.len() < required_len {
            self.atom_name_epochs.resize(required_len, 0);
        }

        if let Some(epoch) = self.atom_name_epochs.get_mut(index) {
            *epoch = self.residue_epoch;
        }
    }

    fn begin_residue_epoch(&mut self) {
        match self.residue_epoch.checked_add(1) {
            Some(epoch) => {
                self.residue_epoch = epoch;
            }
            None => {
                self.atom_name_epochs.fill(0);
                self.residue_epoch = 1;
            }
        }
    }

    fn current_residue(&self) -> ResidueIndex {
        ResidueIndex::new(self.residue_position.saturating_sub(1))
    }

    /// Reports an alternate location that carries a different component code.
    ///
    /// This is a modelled point mutation: one residue with two chemical
    /// identities. Both are kept; the finding says the file did something worth
    /// knowing about rather than something wrong.
    fn check_component_identity(&mut self, rows: &Rows<'_>) {
        let Some(comp) = rows.identifier("label_comp_id") else {
            return;
        };

        let Some(existing) = self
            .data
            .topology
            .residues
            .label_comp_id(self.current_residue())
        else {
            return;
        };

        let Some(existing) = self.data.dictionary.resolve(existing) else {
            return;
        };

        if existing == comp.as_ref() {
            return;
        }

        let diagnostic = Diagnostic::new(Code::W3012)
            .in_category("atom_site")
            .at_row(rows.row() as u32)
            .with_context("component", existing)
            .with_context("alternate component", comp.as_ref());

        self.findings.push(diagnostic);
    }

    fn open_residue(&mut self, rows: &Rows<'_>, key: &ResidueKey) {
        self.close_residue();
        self.open_chain_if_needed(key.chain);

        self.current = Some(*key);
        self.begin_residue_epoch();

        let record = self.residue_record(rows, key);

        self.data
            .topology
            .residues
            .push(record, self.atom_position..self.atom_position);

        self.residue_position += 1;
    }

    fn open_chain_if_needed(&mut self, chain: u32) {
        if self.chain != Some(chain) {
            self.close_chain(chain);
        }
    }

    fn residue_record(&mut self, rows: &Rows<'_>, key: &ResidueKey) -> ResidueRecord {
        ResidueRecord {
            label_comp_id: self.required_symbol_of(rows, "label_comp_id"),
            auth_comp_id: self.optional_symbol_of(rows, "auth_comp_id"),
            label_seq_id: key.label_seq,
            auth_seq_id: key.auth_seq,
            ins_code: self.optional_symbol_of(rows, "pdbx_PDB_ins_code"),
            het: rows.text("group_PDB") == Some("HETATM"),
        }
    }

    fn required_symbol_of(&mut self, rows: &Rows<'_>, column: &str) -> SymbolId {
        match rows.identifier(column) {
            Some(text) => self.intern(text.as_ref()),
            None => self.intern(""),
        }
    }

    fn optional_symbol_of(&mut self, rows: &Rows<'_>, column: &str) -> OptionalSymbol {
        match rows.identifier(column) {
            Some(text) => OptionalSymbol::some(self.intern(text.as_ref())),
            None => OptionalSymbol::NONE,
        }
    }

    fn close_residue(&mut self) {
        let Some(position) = self.residue_position.checked_sub(1) else {
            return;
        };

        let residue = ResidueIndex::new(position);

        let Some(range) = self.data.topology.residues.atoms(residue) else {
            return;
        };

        self.data
            .topology
            .residues
            .set_atoms(residue, range.start..self.atom_position);
    }

    fn close_chain(&mut self, next: u32) {
        if let Some(label) = self.chain
            && self.residue_position > self.chain_first_residue
        {
            self.push_chain(label);
        }

        self.chain = Some(next);
        self.chain_first_residue = self.residue_position;
    }

    fn push_chain(&mut self, label: u32) {
        let entity = self.entity();

        self.data.topology.chains.push(
            ChainRecord {
                label_asym_id: SymbolId::from_raw(label),
                auth_asym_id: OptionalSymbol::some(SymbolId::from_raw(label)),
                entity,
                polymer_kind: PolymerKind::Other,
            },
            self.chain_first_residue..self.residue_position,
        );
    }

    /// The entity a chain belongs to, or the single fallback one.
    fn entity(&mut self) -> EntityIndex {
        if self.data.topology.entities.is_empty() {
            let id = self.intern("1");

            return self.data.topology.entities.push(
                id,
                pdbiox_core::topology::EntityKind::Unknown,
                OptionalSymbol::NONE,
                &[],
            );
        }

        EntityIndex::new(0)
    }

    fn start_model(&mut self, model: i64) {
        if self.model != i64::MIN {
            self.finish_model();
        }

        self.model = model;
    }

    fn finish_model(&mut self) {
        self.close_chain(u32::MAX);
        self.close_residue();

        let builder = std::mem::replace(&mut self.builder, ChunkBuilder::new());

        let (chunks, coords) = builder.finish();

        if self.data.chunks.is_empty() {
            self.data.chunks = chunks;
        }

        self.frames.push(coords);
        self.topology_locked = true;
        self.atom_position = 0;
        self.current = None;
        self.chain = None;
    }

    fn key_of(&mut self, rows: &Rows<'_>, model: i64) -> ResidueKey {
        ResidueKey {
            model,
            chain: self.chain_of(rows),
            label_seq: optional_i32(rows.integer("label_seq_id")),
            auth_seq: optional_i32(rows.integer("auth_seq_id")),
            ins_code: self.interned_raw_or_absent(rows, "pdbx_PDB_ins_code"),
        }
    }

    fn chain_of(&mut self, rows: &Rows<'_>) -> u32 {
        match rows
            .identifier("label_asym_id")
            .or_else(|| rows.identifier("auth_asym_id"))
        {
            Some(text) => self.intern(text.as_ref()).get(),
            None => ResidueKey::ABSENT,
        }
    }

    fn interned_raw_or_absent(&mut self, rows: &Rows<'_>, column: &str) -> u32 {
        match rows.identifier(column) {
            Some(text) => self.intern(text.as_ref()).get(),
            None => ResidueKey::ABSENT,
        }
    }

    fn position_of(&mut self, rows: &Rows<'_>) -> Option<[f32; 3]> {
        let (Some(x), Some(y), Some(z)) = (
            rows.float("Cartn_x"),
            rows.float("Cartn_y"),
            rows.float("Cartn_z"),
        ) else {
            self.findings.push(
                Diagnostic::new(Code::E1202)
                    .with_message("a coordinate could not be read")
                    .in_category("atom_site")
                    .at_row(rows.row() as u32),
            );

            return None;
        };

        Some([x as f32, y as f32, z as f32])
    }

    fn element_of(&mut self, rows: &Rows<'_>, name: Option<&str>) -> Element {
        if let Some(element) = rows.text("type_symbol").and_then(Element::from_symbol) {
            return element;
        }

        let inferred = match name {
            Some(name) => Element::infer_from_name(name),
            None => Element::UNKNOWN,
        };

        self.findings.push(
            Diagnostic::new(Code::W3203)
                .in_category("atom_site")
                .at_row(rows.row() as u32)
                .with_context("inferred", inferred.symbol()),
        );

        inferred
    }

    fn auth_name_of(&mut self, rows: &Rows<'_>, label: Option<&str>) -> OptionalSymbol {
        match rows.identifier("auth_atom_id") {
            Some(text) if Some(text.as_ref()) != label => {
                OptionalSymbol::some(self.intern(text.as_ref()))
            }
            _ => OptionalSymbol::NONE,
        }
    }

    fn alt_of(&mut self, rows: &Rows<'_>) -> AltId {
        match rows.identifier("label_alt_id") {
            Some(text) if !text.is_empty() => AltId::labelled(self.intern(text.as_ref())),
            _ => AltId::BLANK,
        }
    }

    fn intern(&mut self, text: &str) -> SymbolId {
        let Ok(symbol) = self.data.dictionary.intern(text) else {
            self.findings
                .push(Diagnostic::new(Code::E1901).with_message("the dictionary is full"));

            return SymbolId::from_raw(0);
        };

        symbol
    }

    fn finish(mut self) -> CoordinateStore {
        self.close_residue();
        self.close_chain(u32::MAX);

        let (chunks, coords) = self.builder.finish();

        if self.data.chunks.is_empty() {
            self.data.chunks = chunks;
        }

        self.frames.push(coords);

        self.data
            .topology
            .models
            .push(1, 0..self.data.topology.chains.len() as u32);

        coordinate_store(self.frames)
    }
}

fn atom_name_index(atom_name: SymbolId) -> Option<usize> {
    usize::try_from(atom_name.get()).ok()
}

fn formal_charge_of(rows: &Rows<'_>) -> (i8, Presence) {
    match rows.integer("pdbx_formal_charge") {
        // A charge outside a signed byte is not a formal charge; it is
        // a misread column, and is recorded as unstated rather than
        // clamped to something that looks deliberate.
        Some(charge) => match i8::try_from(charge) {
            Ok(charge) => (charge, Presence::Present),
            Err(_) => (0, Presence::Unknown),
        },
        None => (0, Presence::Inapplicable),
    }
}

fn atom_site_id_of(rows: &Rows<'_>) -> u32 {
    rows.integer("id")
        .and_then(|id| u32::try_from(id).ok())
        .map_or(0, |id| id)
}

fn optional_i32(value: Option<i64>) -> OptionalI32 {
    OptionalI32::from(value.and_then(|number| i32::try_from(number).ok()))
}

fn coordinate_store(mut frames: Vec<CoordinateBlock>) -> CoordinateStore {
    if frames.len() != 1 {
        return CoordinateStore::Dense { frames };
    }

    match frames.pop() {
        Some(block) => CoordinateStore::Single(block),
        None => CoordinateStore::Single(CoordinateBlock::new()),
    }
}

fn or_empty(value: Option<&str>) -> &str {
    match value {
        Some(value) => value,
        None => "",
    }
}

fn optional(value: Option<f64>, when_absent: f32) -> (f32, Presence) {
    match value {
        Some(value) => (value as f32, Presence::Present),
        None => (when_absent, Presence::Unknown),
    }
}
