//! What is carried between lines while reading.
//!
//! A fixed-column file describes its hierarchy implicitly: a residue ends when
//! the next line names a different one. Reading it is therefore a walk with
//! memory, and this is the memory.

use super::lines::Line;
use crate::fixed;
use crate::hybrid36;
use pdbiox_core::chunk::{AtomRecord, ChunkBuilder};
use pdbiox_core::column::Presence;
use pdbiox_core::coords::CoordinateBlock;
use pdbiox_core::diagnostic::{Code, Diagnostic, Diagnostics};
use pdbiox_core::element::Element;
use pdbiox_core::index::{EntityIndex, ResidueIndex};
use pdbiox_core::io::{ReadOptions, ReadResult};
use pdbiox_core::optional::{OptionalI32, OptionalSymbol};
use pdbiox_core::span::ByteSpan;
use pdbiox_core::structure::{CoordinateStore, Structure, StructureData, UnitCell};
use pdbiox_core::symbol::{AltId, SymbolId};
use pdbiox_core::topology::{ChainRecord, EntityKind, PolymerKind, ResidueRecord};

/// What a residue is identified by, for deciding where one ends.
#[derive(PartialEq, Eq)]
struct ResidueKey {
    chain: Box<str>,
    seq: i64,
    ins_code: Box<str>,
}

/// Everything carried between lines while reading.
pub(super) struct ReadState<'a> {
    options: &'a ReadOptions,
    data: StructureData,
    builder: ChunkBuilder,
    findings: Diagnostics,
    frames: Vec<CoordinateBlock>,
    model_number: i32,
    model_started: bool,
    chain_label: Option<Box<str>>,
    chain_first_residue: u32,
    residue_key: Option<ResidueKey>,
    residue_position: u32,
    atom_position: u32,
    entity: Option<EntityIndex>,
    saw_atoms: bool,
    topology_locked: bool,
}

impl<'a> ReadState<'a> {
    pub(super) fn new(options: &'a ReadOptions) -> Self {
        Self {
            options,
            data: StructureData::empty(),
            builder: ChunkBuilder::new(),
            findings: Diagnostics::new(),
            frames: Vec::new(),
            model_number: 1,
            model_started: false,
            chain_label: None,
            chain_first_residue: 0,
            residue_key: None,
            residue_position: 0,
            atom_position: 0,
            entity: None,
            saw_atoms: false,
            topology_locked: false,
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
            _ => {}
        }
    }

    fn model(&mut self, line: &Line<'_>) {
        if let Some(number) = fixed::integer(line.text, 11, 14)
            && let Ok(number) = i32::try_from(number)
        {
            self.model_number = number;
        }
        // A second model repeats the same atoms with new positions, so the
        // topology is built once and the frames accumulate beside it.
        if self.model_started {
            self.start_next_frame();
        }
        self.model_started = true;
    }

    fn end_model(&mut self) {
        self.close_chain();
    }

    /// Closes the frame just read and prepares for the next model's positions.
    ///
    /// The first model settles the topology. A later model that disagreed with
    /// it would not be a frame of the same system, and is reported rather than
    /// quietly reshaping the tables.
    fn start_next_frame(&mut self) {
        self.close_chain();
        let builder = std::mem::replace(&mut self.builder, ChunkBuilder::new());
        let (chunks, coords) = builder.finish();
        if self.data.chunks.is_empty() {
            self.data.chunks = chunks;
        }
        self.frames.push(coords);
        self.topology_locked = true;
        self.residue_key = None;
        self.chain_label = None;
        self.atom_position = 0;
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
        if self.options.only_first_model && !self.frames.is_empty() {
            return;
        }

        let raw_name = fixed::raw(line.text, 13, 16);
        let element = self.element_of(line, raw_name);
        if self.options.discard_hydrogens && element.is_hydrogen() {
            return;
        }

        let chain = fixed::text(line.text, 22, 22);
        let Some(seq) = hybrid36::decode(fixed::raw(line.text, 23, 26), 4) else {
            self.findings.push(
                Diagnostic::new(Code::E1202)
                    .with_message("residue number could not be read")
                    .at(ByteSpan::empty(line.at)),
            );
            return;
        };
        let ins_code = fixed::text(line.text, 27, 27);

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
            self.begin_residue_if_new(line, chain, seq, ins_code);
            ResidueIndex::new(self.residue_position.saturating_sub(1))
        };

        let position = if let (Some(x), Some(y), Some(z)) = (
            fixed::real(line.text, 31, 38),
            fixed::real(line.text, 39, 46),
            fixed::real(line.text, 47, 54),
        ) {
            Some([x as f32, y as f32, z as f32])
        } else {
            self.findings.push(
                Diagnostic::new(Code::E1202)
                    .with_message("a coordinate could not be read")
                    .at(ByteSpan::empty(line.at)),
            );
            None
        };

        let atom_name = self.intern(raw_name.trim());
        let alt = fixed::text(line.text, 17, 17);
        let alt_id = if alt.is_empty() {
            AltId::BLANK
        } else {
            AltId::labelled(self.intern(alt))
        };
        let serial = match hybrid36::decode(fixed::raw(line.text, 7, 11), 5)
            .and_then(|serial| u32::try_from(serial).ok())
        {
            Some(serial) => serial,
            // A serial is a file-local label, not identity; a row whose label
            // cannot be read is still a row, and loses only its label.
            None => 0,
        };

        self.builder.push(AtomRecord {
            position,
            element,
            atom_name,
            auth_atom_name: OptionalSymbol::NONE,
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

    fn begin_chain_if_new(&mut self, chain: &str) {
        if self.chain_label.as_deref() == Some(chain) {
            return;
        }
        self.close_chain();
        self.chain_label = Some(chain.into());
        self.chain_first_residue = self.residue_position;
    }

    fn begin_residue_if_new(&mut self, line: &Line<'_>, chain: &str, seq: i64, ins: &str) {
        let key = ResidueKey {
            chain: chain.into(),
            seq,
            ins_code: ins.into(),
        };
        if self.residue_key.as_ref() == Some(&key) {
            return;
        }
        // The residue that was being filled ends where this one begins.
        self.close_residue();
        self.residue_key = Some(key);

        let comp = fixed::text(line.text, 18, 20);
        let label_comp_id = self.intern(comp);
        let ins_code = if ins.is_empty() {
            OptionalSymbol::NONE
        } else {
            OptionalSymbol::some(self.intern(ins))
        };
        let het = fixed::record(line.text) == "HETATM";

        self.data.topology.residues.push(
            ResidueRecord {
                label_comp_id,
                auth_comp_id: OptionalSymbol::some(label_comp_id),
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
        let Some(label) = self.chain_label.take() else {
            return;
        };
        if self.residue_position == self.chain_first_residue {
            return;
        }
        let label_asym_id = self.intern(&label);
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

    pub(super) fn finish(mut self) -> ReadResult {
        self.close_chain();
        if !self.saw_atoms {
            self.findings.push(
                Diagnostic::new(Code::E1001)
                    .with_message("the file contains no coordinate records"),
            );
            return Err(self.findings.finish());
        }

        let (chunks, coords) = self.builder.finish();
        if self.data.chunks.is_empty() {
            self.data.chunks = chunks;
        }
        self.frames.push(coords);

        self.data
            .topology
            .models
            .push(self.model_number, 0..self.data.topology.chains.len() as u32);
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
        let findings = self.findings.finish();
        Ok((structure, findings))
    }
}

/// A field that may be absent, with the value used when it is.
fn optional_real(value: Option<f64>, when_absent: f32) -> (f32, Presence) {
    match value {
        Some(value) => (value as f32, Presence::Present),
        None => (when_absent, Presence::Unknown),
    }
}
