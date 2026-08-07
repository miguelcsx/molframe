//! What is carried between lines while reading.
//!
//! A fixed-column file describes its hierarchy implicitly: a residue ends when
//! the next line names a different one. Reading it is therefore a walk with
//! memory, and this is the memory.

mod connectivity;
mod models;
mod variants;

use models::{AtomSignature, ModelRead};

use super::lines::Line;
use crate::fixed;
use crate::hybrid36;
use pdbiox_core::chunk::{AtomRecord, ChunkBuilder};
use pdbiox_core::column::Presence;
use pdbiox_core::coords::CoordinateBlock;
use pdbiox_core::diagnostic::{Code, Diagnostic, Diagnostics};
use pdbiox_core::element::Element;
use pdbiox_core::index::{AtomIndex, EntityIndex, ResidueIndex};
use pdbiox_core::io::{Format, ReadOptions, ReadResult};
use pdbiox_core::optional::{OptionalI32, OptionalSymbol};
use pdbiox_core::span::{ByteSpan, Position};
use pdbiox_core::structure::{CoordinateStore, Structure, StructureData, UnitCell};
use pdbiox_core::symbol::{AltId, SymbolId};
use pdbiox_core::topology::{ChainRecord, EntityKind, PolymerKind, ResidueRecord};
use std::collections::BTreeMap;

/// What a residue is identified by, for deciding where one ends.
#[derive(Clone, Copy, PartialEq, Eq)]
struct ResidueKey {
    chain: SymbolId,
    seq: i64,
    ins_code: OptionalSymbol,
}

/// Everything carried between lines while reading.
pub(super) struct ReadState<'a> {
    options: &'a ReadOptions,
    data: StructureData,
    builder: ChunkBuilder,
    findings: Diagnostics,
    frames: Vec<CoordinateBlock>,
    model_number: i32,
    model_numbers: Vec<i32>,
    model_read: ModelRead,
    signatures: Vec<AtomSignature>,
    chain_label: Option<SymbolId>,
    chain_first_residue: u32,
    residue_key: Option<ResidueKey>,
    residue_position: u32,
    atom_position: u32,
    entity: Option<EntityIndex>,
    saw_atoms: bool,
    topology_locked: bool,
    serial_to_atom: BTreeMap<u32, AtomIndex>,
    conect: Vec<(u32, u32, Position)>,
    variant: Format,
    partial_charges: Vec<(f64, Presence)>,
    radii: Vec<(f64, Presence)>,
    autodock_types: Vec<(SymbolId, Presence)>,
}

impl<'a> ReadState<'a> {
    pub(super) fn new(options: &'a ReadOptions, variant: Format) -> Self {
        Self {
            options,
            data: StructureData::empty(),
            builder: ChunkBuilder::new(),
            findings: Diagnostics::new(),
            frames: Vec::new(),
            model_number: 1,
            model_numbers: Vec::new(),
            model_read: ModelRead::Before,
            signatures: Vec::new(),
            chain_label: None,
            chain_first_residue: 0,
            residue_key: None,
            residue_position: 0,
            atom_position: 0,
            entity: None,
            saw_atoms: false,
            topology_locked: false,
            serial_to_atom: BTreeMap::new(),
            conect: Vec::new(),
            variant,
            partial_charges: Vec::new(),
            radii: Vec::new(),
            autodock_types: Vec::new(),
        }
    }

    pub(super) fn line(&mut self, line: &Line<'_>) {
        match fixed::record(line.text) {
            "ATOM" | "HETATM" => self.atom(line),
            "MODEL" => self.model(line),
            "ENDMDL" => self.end_model(),
            "TER" => self.close_chain(),
            "CRYST1" => self.cell(line),
            "HEADER" => self.header(line),
            "TITLE" => self.title(line),
            "CONECT" => self.conect(line),
            _ => {}
        }
    }

    fn end_model(&mut self) {
        self.close_chain();
    }

    fn cell(&mut self, line: &Line<'_>) {
        let lengths = [
            fixed::real(line.text, 7, 15),
            fixed::real(line.text, 16, 24),
            fixed::real(line.text, 25, 33),
        ];
        let angles = [
            fixed::real(line.text, 34, 40),
            fixed::real(line.text, 41, 47),
            fixed::real(line.text, 48, 54),
        ];
        let (Some(a), Some(b), Some(c)) = (lengths[0], lengths[1], lengths[2]) else {
            self.findings.push(
                Diagnostic::new(Code::E1202)
                    .with_message("cell lengths could not be read")
                    .at(ByteSpan::empty(line.at)),
            );
            return;
        };
        let (Some(alpha), Some(beta), Some(gamma)) = (angles[0], angles[1], angles[2]) else {
            return;
        };
        self.data.cell = Some(UnitCell {
            lengths: [a, b, c],
            angles: [alpha, beta, gamma],
        });
    }

    fn header(&mut self, line: &Line<'_>) {
        let id = fixed::text(line.text, 63, 66);
        if !id.is_empty() {
            self.data.entry.id = Some(id.into());
        }
    }

    fn title(&mut self, line: &Line<'_>) {
        let text = fixed::text(line.text, 11, 80);
        if text.is_empty() {
            return;
        }
        self.data.entry.title = Some(match self.data.entry.title.take() {
            Some(existing) => format!("{existing} {text}").into(),
            None => text.into(),
        });
    }

    fn atom(&mut self, line: &Line<'_>) {
        // Only the first model's atoms build topology; later models contribute
        // positions to the frames beside it.
        if matches!(self.model_read, ModelRead::Ignoring) {
            return;
        }
        if matches!(self.model_read, ModelRead::Before) {
            self.model_read = ModelRead::Reading;
            self.model_numbers.push(self.model_number);
        }

        let raw_name = fixed::raw(line.text, 13, 16);
        let element = self.element_of(line, raw_name);
        if self.options.discard_hydrogens && element.is_hydrogen() {
            return;
        }

        let chain = self.intern(fixed::text(line.text, 22, 22));
        let Some(seq) = hybrid36::decode(fixed::raw(line.text, 23, 26), 4) else {
            self.findings.push(
                Diagnostic::new(Code::E1202)
                    .with_message("residue number could not be read")
                    .at(ByteSpan::empty(line.at)),
            );
            return;
        };
        let ins_text = fixed::text(line.text, 27, 27);
        let ins_code = if ins_text.is_empty() {
            OptionalSymbol::NONE
        } else {
            OptionalSymbol::some(self.intern(ins_text))
        };
        let atom_name = self.intern(raw_name.trim());
        let component_id = self.intern(fixed::text(line.text, 18, 20));
        let alt = fixed::text(line.text, 17, 17);
        let alt_id = if alt.is_empty() {
            AltId::BLANK
        } else {
            AltId::labelled(self.intern(alt))
        };

        if !self.topology_locked {
            self.variant_atom(line);
        }

        self.observe_signature(
            line,
            AtomSignature {
                chain,
                seq,
                ins_code,
                atom_name,
                component_id,
                alt_id,
                element,
            },
        );

        let residue = if self.topology_locked {
            if let Some(residue) = self.data.topology.residues.containing(self.atom_position) {
                residue
            } else {
                self.findings.push(
                    Diagnostic::new(Code::E3401)
                        .with_message("a later model has more atoms than the first")
                        .at(ByteSpan::empty(line.at)),
                );
                return;
            }
        } else {
            self.begin_chain_if_new(chain);
            self.begin_residue_if_new(line, chain, seq, ins_code, component_id);
            ResidueIndex::new(self.residue_position.saturating_sub(1))
        };

        let position = self.position_of(line);
        let alternate_component_id = self.alternate_component_of(line, residue, component_id);
        let serial = serial_of(line);
        if !self.topology_locked && serial != 0 {
            self.serial_to_atom
                .entry(serial)
                .or_insert(AtomIndex::new(self.atom_position));
        }

        self.builder.push(AtomRecord {
            position,
            element,
            atom_name,
            auth_atom_name: OptionalSymbol::NONE,
            alternate_component_id,
            alt_id,
            residue,
            occupancy: optional_real(fixed::real(line.text, 55, 60), 1.0),
            b_factor: optional_real(fixed::real(line.text, 61, 66), 0.0),
            formal_charge: (0, Presence::Inapplicable),
            atom_site_id: serial,
        });
        self.atom_position += 1;
        self.saw_atoms = true;
    }

    /// The element a row declares, or the one its name implies.
    fn element_of(&mut self, line: &Line<'_>, raw_name: &str) -> Element {
        let declared = fixed::text(line.text, 77, 78);
        if let Some(element) = Element::from_symbol(declared) {
            return element;
        }
        let inferred = Element::infer_from_pdb_atom_name(raw_name);
        self.findings.push(
            Diagnostic::new(Code::W3203)
                .at(ByteSpan::empty(line.at))
                .with_context("atom name", raw_name.trim())
                .with_context("inferred", inferred.symbol()),
        );
        inferred
    }

    fn position_of(&mut self, line: &Line<'_>) -> Option<[f32; 3]> {
        if let (Some(x), Some(y), Some(z)) = (
            fixed::real(line.text, 31, 38),
            fixed::real(line.text, 39, 46),
            fixed::real(line.text, 47, 54),
        ) {
            return Some([x as f32, y as f32, z as f32]);
        }
        self.findings.push(
            Diagnostic::new(Code::E1202)
                .with_message("a coordinate could not be read")
                .at(ByteSpan::empty(line.at)),
        );
        None
    }

    fn alternate_component_of(
        &mut self,
        line: &Line<'_>,
        residue: ResidueIndex,
        component_id: SymbolId,
    ) -> OptionalSymbol {
        let Some(primary) = self.data.topology.residues.label_comp_id(residue) else {
            return OptionalSymbol::NONE;
        };
        if primary == component_id {
            return OptionalSymbol::NONE;
        }
        self.findings.push(
            Diagnostic::new(Code::W3012)
                .at(ByteSpan::empty(line.at))
                .with_context("component", self.component_text(primary))
                .with_context("alternate component", self.component_text(component_id)),
        );
        OptionalSymbol::some(component_id)
    }

    fn begin_chain_if_new(&mut self, chain: SymbolId) {
        if self.chain_label == Some(chain) {
            return;
        }
        self.close_chain();
        self.chain_label = Some(chain);
        self.chain_first_residue = self.residue_position;
    }

    fn begin_residue_if_new(
        &mut self,
        line: &Line<'_>,
        chain: SymbolId,
        seq: i64,
        ins_code: OptionalSymbol,
        component_id: SymbolId,
    ) {
        let key = ResidueKey {
            chain,
            seq,
            ins_code,
        };
        if self.residue_key.as_ref() == Some(&key) {
            return;
        }
        // The residue that was being filled ends where this one begins.
        self.close_residue();
        self.residue_key = Some(key);

        let het = fixed::record(line.text) == "HETATM";

        self.data.topology.residues.push(
            ResidueRecord {
                label_comp_id: component_id,
                auth_comp_id: OptionalSymbol::some(component_id),
                label_seq_id: OptionalI32::NONE,
                auth_seq_id: match i32::try_from(seq) {
                    Ok(seq) => OptionalI32::some(seq),
                    Err(_) => OptionalI32::NONE,
                },
                ins_code,
                het,
            },
            self.atom_position..self.atom_position,
        );
        self.residue_position += 1;
    }

    /// Closes the residue and chain being built, fixing their atom ranges.
    fn close_chain(&mut self) {
        self.close_residue();
        let Some(label_asym_id) = self.chain_label.take() else {
            return;
        };
        if self.residue_position == self.chain_first_residue {
            return;
        }
        let entity = self.entity_for_polymer();
        self.data.topology.chains.push(
            ChainRecord {
                label_asym_id,
                auth_asym_id: OptionalSymbol::some(label_asym_id),
                entity,
                polymer_kind: PolymerKind::Other,
            },
            self.chain_first_residue..self.residue_position,
        );
        self.residue_key = None;
    }

    /// Extends the residue being built to cover the atoms read into it.
    fn close_residue(&mut self) {
        if self.residue_position == 0 {
            return;
        }
        let residue = ResidueIndex::new(self.residue_position - 1);
        let Some(range) = self.data.topology.residues.atoms(residue) else {
            return;
        };
        self.data
            .topology
            .residues
            .set_atoms(residue, range.start..self.atom_position);
    }

    /// The single entity every chain of a fixed-column file shares.
    ///
    /// The format has no concept of a chemical species distinct from a chain, so
    /// there is nothing to distinguish here. Recording one entity rather than
    /// inventing several says exactly that.
    fn entity_for_polymer(&mut self) -> EntityIndex {
        if let Some(entity) = self.entity {
            return entity;
        }
        let id = self.intern("1");
        let entity =
            self.data
                .topology
                .entities
                .push(id, EntityKind::Unknown, OptionalSymbol::NONE, &[]);
        self.entity = Some(entity);
        entity
    }

    fn intern(&mut self, text: &str) -> SymbolId {
        if let Ok(symbol) = self.data.dictionary.intern(text) {
            symbol
        } else {
            self.findings.push(
                Diagnostic::new(Code::E1901).with_message("the identifier dictionary is full"),
            );
            SymbolId::from_raw(0)
        }
    }

    fn component_text(&self, component: SymbolId) -> Box<str> {
        match self.data.dictionary.resolve(component) {
            Some(text) => text.into(),
            None => "".into(),
        }
    }

    pub(super) fn finish(mut self) -> ReadResult {
        self.close_chain();
        self.verify_frame_len();
        if !self.saw_atoms {
            self.findings.push(
                Diagnostic::new(Code::E1001)
                    .with_message("the file contains no coordinate records"),
            );
            return Err(self.findings.finish());
        }

        self.finish_bonds();
        self.finish_variant_annotations();
        let (chunks, coords) = self.builder.finish();
        if self.data.chunks.is_empty() {
            self.data.chunks = chunks.into();
        }
        self.frames.push(coords);
        let chains = 0..self.data.topology.chains.len() as u32;
        for number in self.model_numbers {
            self.data.topology.models.push(number, chains.clone());
        }
        self.data.coords = match self.frames.len() {
            1 => match self.frames.pop() {
                Some(block) => CoordinateStore::Single(block),
                None => CoordinateStore::Single(CoordinateBlock::new()),
            },
            _ => CoordinateStore::Dense {
                frames: self.frames,
            },
        };
        let structure = Structure::new(self.data);
        self.options.finish(structure, self.findings.finish())
    }
}

/// A field that may be absent, with the value used when it is.
fn optional_real(value: Option<f64>, when_absent: f32) -> (f32, Presence) {
    match value {
        Some(value) => (value as f32, Presence::Present),
        None => (when_absent, Presence::Unknown),
    }
}

/// The file-local atom label, or the reserved unlabeled value.
fn serial_of(line: &Line<'_>) -> u32 {
    match hybrid36::decode(fixed::raw(line.text, 7, 11), 5)
        .and_then(|serial| u32::try_from(serial).ok())
    {
        Some(serial) => serial,
        None => 0,
    }
}
