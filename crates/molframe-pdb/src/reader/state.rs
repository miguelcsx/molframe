//! What is carried between lines while reading.
//!
//! A fixed-column file describes its hierarchy implicitly: a residue ends when
//! the next line names a different one. Reading it is therefore a walk with
//! memory, and this is the memory.

mod connectivity;
mod finalize;
mod metadata;
mod models;
mod variants;

use models::{AtomIdentity, ModelRead};

use super::lines::Line;
use crate::fixed;
use crate::header::PdbHeaders;
use crate::hybrid36;
use molframe_core::chunk::{AtomRecord, ChunkBuilder};
use molframe_core::column::Presence;
use molframe_core::coords::CoordinateBlock;
use molframe_core::diagnostic::{Code, Diagnostic, Diagnostics};
use molframe_core::element::Element;
use molframe_core::index::{AtomIndex, EntityIndex, ResidueIndex};
use molframe_core::io::{Format, MissingElementPolicy, ReadOptions};
use molframe_core::optional::{OptionalI32, OptionalSymbol};
use molframe_core::span::{ByteSpan, Position};
use molframe_core::structure::StructureData;
use molframe_core::symbol::SymbolId;
use molframe_core::topology::{ChainRecord, EntityKind, PolymerKind, ResidueRecord};
use num_traits::ToPrimitive;
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
    frame: CoordinateBlock,
    frames: Vec<CoordinateBlock>,
    model_number: i32,
    model_numbers: Vec<i32>,
    model_read: ModelRead,
    chain_label: Option<SymbolId>,
    chain_first_residue: u32,
    residue_key: Option<ResidueKey>,
    residue_position: u32,
    atom_position: u32,
    entity: Option<EntityIndex>,
    saw_atoms: bool,
    topology_locked: bool,
    expected_chunk: usize,
    requires_ragged: bool,
    serial_to_atom: BTreeMap<u32, AtomIndex>,
    conect: Vec<(u32, u32, Position)>,
    variant: Format,
    partial_charges: Vec<(f64, Presence)>,
    radii: Vec<(f64, Presence)>,
    autodock_types: Vec<(SymbolId, Presence)>,
    headers: PdbHeaders,
}

impl<'a> ReadState<'a> {
    pub(super) fn new(options: &'a ReadOptions, variant: Format) -> Self {
        Self {
            options,
            data: StructureData::empty(),
            builder: ChunkBuilder::new(),
            findings: Diagnostics::new(),
            frame: CoordinateBlock::new(),
            frames: Vec::new(),
            model_number: 1,
            model_numbers: Vec::new(),
            model_read: ModelRead::Before,
            chain_label: None,
            chain_first_residue: 0,
            residue_key: None,
            residue_position: 0,
            atom_position: 0,
            entity: None,
            saw_atoms: false,
            topology_locked: false,
            expected_chunk: 0,
            requires_ragged: false,
            serial_to_atom: BTreeMap::new(),
            conect: Vec::new(),
            variant,
            partial_charges: Vec::new(),
            radii: Vec::new(),
            autodock_types: Vec::new(),
            headers: PdbHeaders::default(),
        }
    }

    pub(super) fn line(&mut self, line: &Line<'_>) {
        let record = fixed::record(line.text);
        self.observe_metadata(record, line);
        match record {
            "ATOM" | "HETATM" => self.atom(line),
            "MODEL" => self.model(line),
            "ENDMDL" => self.end_model(),
            "TER" => self.close_chain(),
            "CRYST1" => self.cell(line),
            "HEADER" => self.header(line),
            "TITLE" => self.title(line),
            "EXPDTA" => self.method(line),
            "CONECT" => self.conect(line),
            _ => {}
        }
    }

    fn end_model(&mut self) {
        self.close_chain();
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

        let Some(signature) = self.atom_identity(line, raw_name, element) else {
            return;
        };
        if self.topology_locked {
            self.later_atom(line, signature);
            return;
        }
        let AtomIdentity {
            chain,
            seq,
            ins_code,
            atom_name,
            component_id,
            alt_id,
            element,
            ..
        } = signature;

        self.variant_atom(line);
        self.begin_chain_if_new(chain);
        self.begin_residue_if_new(line, chain, seq, ins_code, component_id);
        let Some(position) = self.residue_position.checked_sub(1) else {
            self.findings.push(
                Diagnostic::new(Code::E1901)
                    .with_message("atom row could not be assigned to a residue")
                    .at(ByteSpan::empty(line.at)),
            );
            return;
        };
        let residue = ResidueIndex::new(position);

        let position = self.position_of(line);
        let occupancy = self.optional_real(line, 55, 60, "occupancy", 1.0);
        let b_factor = self.optional_real(line, 61, 66, "B factor", 0.0);
        let alternate_component_id = self.alternate_component_of(line, residue, component_id);
        let serial = serial_of(line);
        if serial != 0 {
            self.serial_to_atom
                .entry(serial)
                .or_insert(AtomIndex::new(self.atom_position));
        }

        self.builder.push(AtomRecord {
            position,
            element,
            atom_name,
            auth_atom_name: OptionalSymbol::some(atom_name),
            alternate_component_id,
            alt_id,
            residue,
            occupancy,
            b_factor,
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
        let inferred = match self.options.missing_element_policy {
            MissingElementPolicy::PreserveUnknown => Element::UNKNOWN,
            MissingElementPolicy::InferFromAtomName => Element::infer_from_pdb_atom_name(raw_name),
        };
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
        ) && let Some(position) = x
            .to_f32()
            .zip(y.to_f32())
            .zip(z.to_f32())
            .map(|((x, y), z)| [x, y, z])
        {
            return Some(position);
        }
        self.findings.push(
            Diagnostic::new(Code::E1202)
                .with_message("a coordinate could not be read")
                .at(ByteSpan::empty(line.at)),
        );
        None
    }

    fn optional_real(
        &mut self,
        line: &Line<'_>,
        start: usize,
        end: usize,
        field: &'static str,
        when_absent: f32,
    ) -> (f32, Presence) {
        let Some(value) = fixed::real(line.text, start, end) else {
            return (when_absent, Presence::Unknown);
        };
        if let Some(value) = value.to_f32() {
            return (value, Presence::Present);
        }
        self.findings.push(
            Diagnostic::new(Code::E1202)
                .with_message("a numeric atom field exceeds the supported floating-point range")
                .with_context("field", field)
                .at(ByteSpan::empty(line.at)),
        );
        (when_absent, Presence::Unknown)
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

        if self
            .data
            .topology
            .residues
            .push(
                ResidueRecord {
                    label_comp_id: component_id,
                    auth_comp_id: OptionalSymbol::some(component_id),
                    label_seq_id: sequence_identifier(seq),
                    auth_seq_id: sequence_identifier(seq),
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

    /// Closes the residue and chain being built, fixing their atom ranges.
    fn close_chain(&mut self) {
        self.close_residue();
        let Some(label_asym_id) = self.chain_label.take() else {
            return;
        };
        if self.residue_position == self.chain_first_residue {
            return;
        }
        let Ok(entity) = self.unclassified_entity() else {
            self.findings.push(Diagnostic::new(Code::E3001));
            return;
        };
        if self
            .data
            .topology
            .chains
            .push(
                ChainRecord {
                    label_asym_id,
                    auth_asym_id: OptionalSymbol::some(label_asym_id),
                    entity,
                    polymer_kind: PolymerKind::None,
                },
                self.chain_first_residue..self.residue_position,
            )
            .is_err()
        {
            self.findings.push(Diagnostic::new(Code::E3001));
        }
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
        if let Err(error) = self
            .data
            .topology
            .residues
            .set_atoms(residue, range.start..self.atom_position)
        {
            self.findings
                .push(Diagnostic::new(Code::E3001).with_context("cause", error.to_string()));
        }
    }

    /// The single entity every chain of a fixed-column file shares.
    ///
    /// The format has no concept of a chemical species distinct from a chain, so
    /// there is nothing to distinguish here. Recording one entity rather than
    /// inventing several says exactly that.
    fn unclassified_entity(&mut self) -> Result<EntityIndex, Diagnostic> {
        if let Some(entity) = self.entity {
            return Ok(entity);
        }
        let id = self.intern("1");
        let entity =
            self.data
                .topology
                .entities
                .push(id, EntityKind::Unknown, OptionalSymbol::NONE, &[]);
        let entity = entity.map_err(|error| {
            Diagnostic::new(Code::E3001).with_context("cause", error.to_string())
        })?;
        self.entity = Some(entity);
        Ok(entity)
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
}

fn sequence_identifier(sequence: i64) -> OptionalI32 {
    match i32::try_from(sequence) {
        Ok(sequence) => OptionalI32::some(sequence),
        Err(_) => OptionalI32::NONE,
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
