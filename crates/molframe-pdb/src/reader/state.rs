//! What is carried between lines while reading.
//!
//! A fixed-column file describes its hierarchy implicitly: a residue ends when
//! the next line names a different one. Reading it is therefore a walk with
//! memory, and this is the memory.

mod connectivity;
mod finalize;
mod hierarchy;
mod links;
mod metadata;
mod models;
mod rows;
mod secondary;
mod variants;

pub(super) use rows::{charge_in, charge_of, serial_in, serial_of, starts_with_coordinate_keyword};

use models::{AtomIdentity, ModelRead};

use super::lines::Line;
use crate::fixed;
use crate::header::PdbHeaders;
use crate::hybrid36;
use links::{LinkAtoms, LinkRequest};
use molframe_core::chunk::{AtomRecord, ChunkBuilder};
use molframe_core::column::Presence;
use molframe_core::coords::CoordinateBlock;
use molframe_core::diagnostic::{Code, Diagnostic, Diagnostics};
use molframe_core::element::Element;
use molframe_core::index::{AtomIndex, EntityIndex, ResidueIndex};
use molframe_core::io::{Format, MissingElementPolicy, ReadOptions};
use molframe_core::optional::OptionalSymbol;
use molframe_core::span::{ByteSpan, Position};
use molframe_core::structure::StructureData;
use molframe_core::symbol::SymbolId;
use num_traits::ToPrimitive;
use std::collections::BTreeMap;

/// What a residue is identified by, for deciding where one ends.
#[derive(Clone, Copy, PartialEq, Eq)]
struct ResidueKey {
    chain: SymbolId,
    seq: i64,
    ins_code: OptionalSymbol,
}

/// How the models after the first depart from it.
#[derive(Clone, Copy, Default)]
struct Divergence {
    /// The topology or a coordinate frame does not match the first model's.
    ragged: bool,
    /// A later model repeated the topology with other occupancies or B-factors.
    values: bool,
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
    /// Coordinate-looking rows whose record name could not be read.
    malformed_atom: u64,
    topology_locked: bool,
    expected_chunk: usize,
    later_models: Divergence,
    serial_to_atom: BTreeMap<u32, AtomIndex>,
    conect: Vec<(u32, u32, Position)>,
    /// `SSBOND` and `LINK` records awaiting the atoms they name.
    links: Vec<LinkRequest>,
    link_atoms: LinkAtoms,
    key_scratch: String,
    /// Ellipsoid rows awaiting their atom, joined by deposited serial.
    anisou: Vec<(u32, [f32; 6])>,
    /// The serial of the most recently read atom row, if it carried one.
    last_atom_serial: Option<u32>,
    variant: Format,
    partial_charges: Vec<(f64, Presence)>,
    radii: Vec<(f64, Presence)>,
    autodock_types: Vec<(SymbolId, Presence)>,
    headers: PdbHeaders,
    branches: usize,
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
            malformed_atom: 0,
            topology_locked: false,
            expected_chunk: 0,
            later_models: Divergence::default(),
            serial_to_atom: BTreeMap::new(),
            conect: Vec::new(),
            links: Vec::new(),
            link_atoms: LinkAtoms::new(),
            key_scratch: String::new(),
            anisou: Vec::new(),
            last_atom_serial: None,
            variant,
            partial_charges: Vec::new(),
            radii: Vec::new(),
            autodock_types: Vec::new(),
            headers: PdbHeaders::default(),
            branches: 0,
        }
    }

    pub(super) fn line(&mut self, line: &Line<'_>) {
        let record = fixed::record(line.text);
        if record == "BRANCH" && self.variant == Format::Pdbqt {
            self.branches += 1;
        }
        if !self.options.categories.keeps_record(record) {
            // A malformed coordinate row is filtered out here too, because its
            // broken record name matches no filter entry. Counting it before the
            // filter keeps the diagnostic honest about what the file held.
            if starts_with_coordinate_keyword(line.text) {
                self.malformed_atom += 1;
            }
            return;
        }
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
            "SSBOND" => self.ssbond(line),
            "LINK" => self.link(line),
            "ANISOU" => self.anisou(line),
            // A coordinate record whose serial outgrew its five columns pushes a
            // digit into column 6, so the record name reads `ATOM 1` and the line
            // matches nothing here. It is a real malformed row rather than a line
            // this reader has no interest in, and saying so is the difference
            // between a caller fixing their file and a caller being told the file
            // holds no coordinates at all.
            _ if starts_with_coordinate_keyword(line.text) => {
                self.malformed_atom += 1;
            }
            _ => {}
        }
    }

    fn end_model(&mut self) {
        self.close_chain();
    }

    /// Reads one anisotropic displacement record.
    ///
    /// The six integer fields at columns 28 to 70 are B-style displacements in
    /// units of 0.0001 ångström squared, so each divides by 10 000 to become
    /// its U component. The record belongs to the atom it names a serial for;
    /// a serial that does not match the atom just read, or a row the file
    /// introduced before any atom, describes no atom this reader has and is
    /// ignored rather than guessed onto a neighbour.
    fn anisou(&mut self, line: &Line<'_>) {
        let serial = serial_of(line);
        let Some(last) = self.last_atom_serial else {
            return;
        };
        if serial != last {
            return;
        }
        let mut components = [0.0_f32; 6];
        let columns = [(29, 35), (36, 42), (43, 49), (50, 56), (57, 63), (64, 70)];
        for (position, (from, to)) in columns.iter().enumerate() {
            let Some(value) = fixed::integer(line.text, *from, *to) else {
                return;
            };
            let Some(scaled) = value.to_f64().and_then(|value| (value / 10_000.0).to_f32()) else {
                return;
            };
            components[position] = scaled;
        }
        self.anisou.push((last, components));
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
        let (occupancy, b_factor) = if self.variant == Format::Pqr {
            ((1.0, Presence::Unknown), (0.0, Presence::Unknown))
        } else {
            (
                self.optional_real(line, line.occupancy(), "occupancy", 1.0),
                self.optional_real(line, line.b_factor(), "B factor", 0.0),
            )
        };
        let alternate_component_id = self.alternate_component_of(line, residue, component_id);
        let serial = serial_of(line);
        if serial != 0 {
            self.serial_to_atom
                .entry(serial)
                .or_insert(AtomIndex::new(self.atom_position));
            self.last_atom_serial = Some(serial);
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
            formal_charge: if matches!(self.variant, Format::Pqr | Format::Pdbqt) {
                (0, Presence::Unknown)
            } else {
                charge_of(line)
            },
            atom_site_id: serial,
        });
        self.note_link_atom(line, AtomIndex::new(self.atom_position), element);
        self.atom_position += 1;
        self.saw_atoms = true;
    }

    /// The element a row declares, or the one its name implies.
    fn element_of(&mut self, line: &Line<'_>, raw_name: &str) -> Element {
        let declared = match line.numbers {
            Some(numbers) => numbers.declared,
            None => Element::from_symbol(fixed::text(line.text, 77, 78)),
        };
        if let Some(element) = declared {
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
        let [x, y, z] = match line.numbers {
            Some(numbers) => numbers.position,
            None => [
                fixed::real(line.text, 31, 38),
                fixed::real(line.text, 39, 46),
                fixed::real(line.text, 47, 54),
            ],
        };
        if let (Some(x), Some(y), Some(z)) = (x, y, z)
            && let Some(position) = x
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
        value: Option<f64>,
        field: &'static str,
        when_absent: f32,
    ) -> (f32, Presence) {
        let Some(value) = value else {
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
}
