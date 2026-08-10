//! Turning coordinate rows into a structure.
//!
//! This is where interpretation happens, and therefore where most findings come
//! from. Every decision the file does not force — where a residue ends, which
//! element an atom is, whether two rows are the same residue modelled twice — is
//! made here and reported when it was not forced.

mod chains;
mod finish;
mod models;

use models::AtomSignature;

use super::diagnostics::at_source_row;
use super::entry::AsymEntity;
use super::keys::{Boundary, ResidueKey, boundary};
use crate::document::Category;
use crate::parser::Rows;
use num_traits::ToPrimitive;
use pdbiox_core::chunk::{AtomRecord, ChunkBuilder};
use pdbiox_core::column::Presence;
use pdbiox_core::coords::CoordinateBlock;
use pdbiox_core::diagnostic::{Code, Diagnostic, Diagnostics};
use pdbiox_core::element::Element;
use pdbiox_core::index::{EntityIndex, ResidueIndex};
use pdbiox_core::io::{AmbiguousResidueBoundaryPolicy, MissingElementPolicy, ReadOptions};
use pdbiox_core::optional::{OptionalI32, OptionalSymbol};
use pdbiox_core::structure::{CoordinateStore, StructureData};
use pdbiox_core::symbol::{AltId, SymbolId};
use pdbiox_core::topology::ResidueRecord;

/// Builds the atoms, residues, chains and models of a structure.
pub struct AtomBuilder<'a> {
    data: &'a mut StructureData,
    findings: &'a mut Diagnostics,
    options: &'a ReadOptions,
    asym_entities: &'a [AsymEntity],
    model_filter: Option<i64>,
    builder: ChunkBuilder,
    frames: Vec<CoordinateBlock>,
    model_numbers: Vec<i32>,
    signatures: Vec<AtomSignature>,
    /// The residue being filled, and the atom names already in it.
    current: Option<ResidueKey>,
    names_in_residue: Vec<(SymbolId, AltId)>,
    /// The chain being filled.
    chain: Option<u32>,
    chain_auth: OptionalSymbol,
    chain_entity: Option<EntityIndex>,
    chain_first_residue: u32,
    residue_position: u32,
    atom_position: u32,
    model: i64,
    topology_locked: bool,
    reported_boundary_inference: bool,
}

impl<'a> AtomBuilder<'a> {
    /// Starts building into `data`.
    pub fn new(
        data: &'a mut StructureData,
        findings: &'a mut Diagnostics,
        options: &'a ReadOptions,
        asym_entities: &'a [AsymEntity],
    ) -> Self {
        Self {
            data,
            findings,
            options,
            asym_entities,
            model_filter: None,
            builder: ChunkBuilder::new(),
            frames: Vec::new(),
            model_numbers: Vec::new(),
            signatures: Vec::new(),
            current: None,
            names_in_residue: Vec::new(),
            chain: None,
            chain_auth: OptionalSymbol::NONE,
            chain_entity: None,
            chain_first_residue: 0,
            residue_position: 0,
            atom_position: 0,
            model: i64::MIN,
            topology_locked: false,
            reported_boundary_inference: false,
        }
    }

    /// Restricts lowering to one deposited model.
    #[must_use]
    pub fn only_model(mut self, model: i64) -> Self {
        self.model_filter = Some(model);
        self
    }

    /// Reads every row of the coordinate category.
    pub fn read(mut self, category: &Category) -> CoordinateStore {
        let mut rows = Rows::new(category);
        if category.row_count() > 0 {
            loop {
                self.row(&rows);
                if !rows.advance() {
                    break;
                }
            }
        }
        self.finish()
    }

    fn row(&mut self, rows: &Rows<'_>) {
        // A file that does not number its models has exactly one.
        let model = match rows.integer("pdbx_PDB_model_num") {
            Some(model) => model,
            None => self.model.max(1),
        };
        if self.model_filter.is_some_and(|wanted| wanted != model) {
            return;
        }
        if model != self.model {
            if self.options.only_first_model && self.model != i64::MIN {
                return;
            }
            self.start_model(model);
        }

        let name_text = rows.identifier("label_atom_id");
        let element = self.element_of(rows, name_text.as_deref());
        if self.options.discard_hydrogens && element.is_hydrogen() {
            return;
        }

        let atom_name = self.intern(match name_text.as_deref() {
            Some(name) => name,
            None => "",
        });
        let Some(alt_id) = self.alt_of(rows) else {
            self.findings.push(
                Diagnostic::new(Code::E1901)
                    .with_message("alternate-location identifier exceeds its encoding"),
            );
            return;
        };
        let key = self.key_of(rows, model);
        let Some(residue) = self.place(rows, &key, atom_name, alt_id) else {
            return;
        };
        let alternate_component_id = self.alternate_component_of(rows, residue);

        let position = self.position_of(rows);
        let auth_atom_name = self.auth_name_of(rows, name_text.as_deref());
        let primary_component_id = match self.data.topology.residues.label_comp_id(residue) {
            Some(component) => component,
            None => SymbolId::from_raw(0),
        };
        let component_id = match alternate_component_id.get() {
            Some(component) => component,
            None => primary_component_id,
        };
        self.observe_signature(
            rows,
            AtomSignature {
                chain: key.chain,
                label_seq: key.label_seq,
                auth_seq: key.auth_seq,
                ins_code: key.ins_code,
                atom_name,
                auth_atom_name,
                component_id,
                alt_id,
                element,
            },
        );
        let occupancy = self.optional_float(rows, "occupancy", 1.0);
        let b_factor = self.optional_float(rows, "B_iso_or_equiv", 0.0);
        self.builder.push(AtomRecord {
            position,
            element,
            atom_name,
            auth_atom_name,
            alternate_component_id,
            alt_id,
            residue,
            occupancy,
            b_factor,
            formal_charge: match rows.integer("pdbx_formal_charge") {
                // A charge outside a signed byte is not a formal charge; it is
                // a misread column, and is recorded as unstated rather than
                // clamped to something that looks deliberate.
                Some(charge) => match i8::try_from(charge) {
                    Ok(charge) => (charge, Presence::Present),
                    Err(_) => (0, Presence::Unknown),
                },
                None => (0, Presence::Inapplicable),
            },
            atom_site_id: match rows.integer("id").and_then(|id| u32::try_from(id).ok()) {
                Some(id) => id,
                None => 0,
            },
        });
        self.atom_position += 1;
    }

    /// Where this atom's residue sits, opening a new one if the row starts one.
    fn place(
        &mut self,
        rows: &Rows<'_>,
        key: &ResidueKey,
        atom_name: SymbolId,
        alt_id: AltId,
    ) -> Option<ResidueIndex> {
        if self.topology_locked {
            return if let Some(residue) = self.data.topology.residues.containing(self.atom_position)
            {
                Some(residue)
            } else {
                self.findings.push(at_source_row(
                    Diagnostic::new(Code::E3401)
                        .with_message("a later model has more atoms than the first"),
                    rows.row(),
                ));
                None
            };
        }

        let repeats = self.names_in_residue.contains(&(atom_name, alt_id));
        let decision = match &self.current {
            None => Boundary::New,
            Some(current) => boundary(current, key, repeats),
        };

        if matches!(decision, Boundary::New | Boundary::NewByFileOrder) {
            if decision == Boundary::NewByFileOrder {
                self.report_ambiguous_boundary(rows);
            }
            self.open_residue(rows, key);
        }
        self.names_in_residue.push((atom_name, alt_id));
        let Some(position) = self.residue_position.checked_sub(1) else {
            self.findings.push(at_source_row(
                Diagnostic::new(Code::E1901)
                    .with_message("atom row could not be assigned to a residue"),
                rows.row(),
            ));
            return None;
        };
        Some(ResidueIndex::new(position))
    }

    fn report_ambiguous_boundary(&mut self, rows: &Rows<'_>) {
        let code = match self.options.ambiguous_residue_boundary_policy {
            AmbiguousResidueBoundaryPolicy::Reject => Code::E3015,
            AmbiguousResidueBoundaryPolicy::InferFromFileOrder => Code::W3011,
        };
        if code == Code::W3011 && self.reported_boundary_inference {
            return;
        }
        self.reported_boundary_inference = true;
        self.findings.push(at_source_row(
            Diagnostic::new(code).in_category("atom_site"),
            rows.row(),
        ));
    }

    /// Reports an alternate location that carries a different component code.
    ///
    /// This is a modelled point mutation: one residue with two chemical
    /// identities. Both are kept; the finding says the file did something worth
    /// knowing about rather than something wrong.
    fn alternate_component_of(&mut self, rows: &Rows<'_>, residue: ResidueIndex) -> OptionalSymbol {
        let Some(comp) = rows.identifier("label_comp_id") else {
            return OptionalSymbol::NONE;
        };
        let Some(existing_symbol) = self.data.topology.residues.label_comp_id(residue) else {
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

    fn open_residue(&mut self, rows: &Rows<'_>, key: &ResidueKey) {
        self.close_residue();
        if self.chain != Some(key.chain) {
            self.open_chain(rows, key.chain);
        }
        self.current = Some(*key);
        self.names_in_residue.clear();

        let comp = self.intern(match rows.identifier("label_comp_id").as_deref() {
            Some(name) => name,
            None => "",
        });
        let auth_comp = match rows.identifier("auth_comp_id") {
            Some(text) => OptionalSymbol::some(self.intern(&text)),
            None => OptionalSymbol::NONE,
        };
        let ins_code = match rows.identifier("pdbx_PDB_ins_code") {
            Some(text) => OptionalSymbol::some(self.intern(&text)),
            None => OptionalSymbol::NONE,
        };
        let het = rows.text("group_PDB") == Some("HETATM");

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

    fn key_of(&mut self, rows: &Rows<'_>, model: i64) -> ResidueKey {
        let chain = match rows.identifier("label_asym_id") {
            Some(text) => self.intern(&text).get(),
            None => ResidueKey::ABSENT,
        };
        let ins_code = match rows.identifier("pdbx_PDB_ins_code") {
            Some(text) => self.intern(&text).get(),
            None => ResidueKey::ABSENT,
        };
        ResidueKey {
            model,
            chain,
            label_seq: OptionalI32::from(
                rows.integer("label_seq_id")
                    .and_then(|seq| i32::try_from(seq).ok()),
            ),
            auth_seq: OptionalI32::from(
                rows.integer("auth_seq_id")
                    .and_then(|seq| i32::try_from(seq).ok()),
            ),
            ins_code,
        }
    }

    fn position_of(&mut self, rows: &Rows<'_>) -> Option<[f32; 3]> {
        let (Some(x), Some(y), Some(z)) = (
            rows.float("Cartn_x"),
            rows.float("Cartn_y"),
            rows.float("Cartn_z"),
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

    fn optional_float(
        &mut self,
        rows: &Rows<'_>,
        field: &'static str,
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
                .with_context("item", field),
            rows.row(),
        ));
        (when_absent, Presence::Unknown)
    }

    fn element_of(&mut self, rows: &Rows<'_>, name: Option<&str>) -> Element {
        if let Some(element) = rows.text("type_symbol").and_then(Element::from_symbol) {
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

    fn auth_name_of(&mut self, rows: &Rows<'_>, label: Option<&str>) -> OptionalSymbol {
        if let Some(text) = rows.identifier("auth_atom_id")
            && Some(text.as_ref()) != label
        {
            return OptionalSymbol::some(self.intern(&text));
        }
        OptionalSymbol::NONE
    }

    fn alt_of(&mut self, rows: &Rows<'_>) -> Option<AltId> {
        if let Some(text) = rows.identifier("label_alt_id")
            && !text.is_empty()
        {
            return AltId::labelled(self.intern(&text));
        }
        Some(AltId::BLANK)
    }

    fn intern(&mut self, text: &str) -> SymbolId {
        let Ok(symbol) = self.data.dictionary.intern(text) else {
            self.findings
                .push(Diagnostic::new(Code::E1901).with_message("the dictionary is full"));
            return SymbolId::from_raw(0);
        };
        symbol
    }
}
