//! Turning coordinate rows into a structure.
//!
//! This is where interpretation happens, and therefore where most findings come
//! from. Every decision the file does not force — where a residue ends, which
//! element an atom is, whether two rows are the same residue modelled twice — is
//! made here and reported when it was not forced.

mod chains;
mod fields;
mod finish;
mod identity;
mod models;
mod names;
mod row;

pub use row::{AtomSiteRow, Field};

use models::AtomSignature;
use names::ResidueNames;

use super::diagnostics::at_source_row;
use super::entry::AsymEntity;
use super::keys::{Boundary, ResidueKey, boundary};
use crate::document::Category;
use crate::parser::Rows;
use molframe_core::anisotropy::AnisotropicDisplacement;
use molframe_core::chunk::{AtomRecord, ChunkBuilder};
use molframe_core::column::Presence;
use molframe_core::coords::CoordinateBlock;
use molframe_core::diagnostic::{Code, Diagnostic, Diagnostics};
use molframe_core::index::{AtomIndex, EntityIndex, ResidueIndex};
use molframe_core::io::{AmbiguousResidueBoundaryPolicy, ReadOptions};
use molframe_core::optional::OptionalSymbol;
use molframe_core::structure::{CoordinateStore, StructureData};
use molframe_core::symbol::{AltId, SymbolId};
use num_traits::ToPrimitive;

/// Builds the atoms, residues, chains and models of a structure.
pub struct AtomBuilder<'a> {
    data: StructureData,
    findings: Diagnostics,
    options: &'a ReadOptions,
    asym_entities: Vec<AsymEntity>,
    model_filter: Option<i64>,
    builder: ChunkBuilder,
    frames: Vec<CoordinateBlock>,
    model_numbers: Vec<i32>,
    signatures: Vec<AtomSignature>,
    /// Ellipsoids from the inline `_atom_site.aniso_*` spelling, if any row had one.
    inline_anisotropy: Option<molframe_core::AnisotropyTableBuilder>,
    expected: Option<models::ExpectedAtoms>,
    track_identity: bool,
    /// The residue being filled, and the atom names already in it.
    current: Option<ResidueKey>,
    names_in_residue: ResidueNames,
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
    /// Symbols already interned for a dictionary slot, one table per field.
    symbol_cache: Vec<Vec<SymbolId>>,
}

impl<'a> AtomBuilder<'a> {
    /// Starts building from prepared structure state.
    pub fn new(
        data: StructureData,
        findings: Diagnostics,
        options: &'a ReadOptions,
        asym_entities: Vec<AsymEntity>,
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
            inline_anisotropy: None,
            expected: None,
            track_identity: true,
            current: None,
            names_in_residue: ResidueNames::default(),
            chain: None,
            chain_auth: OptionalSymbol::NONE,
            chain_entity: None,
            chain_first_residue: 0,
            residue_position: 0,
            atom_position: 0,
            model: i64::MIN,
            topology_locked: false,
            reported_boundary_inference: false,
            symbol_cache: vec![Vec::new(); Field::COUNT],
        }
    }

    /// Restricts lowering to one deposited model.
    #[must_use]
    pub fn only_model(mut self, model: i64) -> Self {
        self.model_filter = Some(model);
        self
    }

    /// Reads every row of the coordinate category.
    pub fn read(mut self, category: &Category) -> (StructureData, Diagnostics, CoordinateStore) {
        let mut rows = Rows::new(category);
        if category.row_count() > 0 {
            loop {
                self.feed(&rows);
                if !rows.advance() {
                    break;
                }
            }
        }
        self.finish()
    }
    /// Lowers one row, in deposition order.
    pub(crate) fn feed<R: AtomSiteRow + ?Sized>(&mut self, rows: &R) {
        // A file that does not number its models has exactly one.
        let model = match rows.integer(Field::ModelNum) {
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

        let element = self.element_of(rows);
        if self.options.discard_hydrogens && element.is_hydrogen() {
            return;
        }

        let atom_name = self.symbol_or_empty(rows, Field::LabelAtomId);
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
        let auth_atom_name = self.auth_name_of(rows);
        let primary_component_id = match self.data.topology.residues.label_comp_id(residue) {
            Some(component) => component,
            None => SymbolId::from_raw(0),
        };
        let component_id = match alternate_component_id.get() {
            Some(component) => component,
            None => primary_component_id,
        };
        let occupancy = self.optional_float(rows, Field::Occupancy, 1.0);
        let b_factor = self.optional_float(rows, Field::BFactor, 0.0);
        let record = AtomRecord {
            position,
            element,
            atom_name,
            auth_atom_name,
            alternate_component_id,
            alt_id,
            residue,
            occupancy,
            b_factor,
            formal_charge: match rows.integer(Field::FormalCharge) {
                // A charge outside a signed byte is not a formal charge; it is
                // a misread column, and is recorded as unstated rather than
                // clamped to something that looks deliberate.
                Some(charge) => match i8::try_from(charge) {
                    Ok(charge) => (charge, Presence::Present),
                    Err(_) => (0, Presence::Unknown),
                },
                None => (0, Presence::Inapplicable),
            },
            atom_site_id: identifier_or_zero(
                rows.integer(Field::Id)
                    .and_then(|id| u32::try_from(id).ok()),
            ),
        };
        self.observe_atom(
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
            &record,
        );
        self.builder.push(record);
        let atom = AtomIndex::new(self.atom_position);
        self.record_inline_tensor(rows, atom);
        self.atom_position += 1;
    }

    /// Records an inline `_atom_site.aniso_*` tensor beside its atom row.
    ///
    /// The legacy spelling kept the ellipsoid in `atom_site` itself, so one
    /// coordinate row carries both the coordinates and all six components. A
    /// row missing any component is left unattached: completing the tensor
    /// from the values that happened to be written would invent data.
    fn record_inline_tensor<R: AtomSiteRow + ?Sized>(&mut self, rows: &R, atom: AtomIndex) {
        let components = [
            Field::AnisoU11,
            Field::AnisoU22,
            Field::AnisoU33,
            Field::AnisoU12,
            Field::AnisoU13,
            Field::AnisoU23,
        ];
        let Some(u) = components
            .iter()
            .map(|field| rows.float(*field).and_then(|value| value.to_f32()))
            .collect::<Option<Vec<f32>>>()
            .map(|values| {
                [
                    values[0], values[1], values[2], values[3], values[4], values[5],
                ]
            })
        else {
            return;
        };
        self.inline_anisotropy
            .get_or_insert_with(molframe_core::AnisotropyTableBuilder::new)
            .push(AnisotropicDisplacement { atom, u });
    }

    /// Where this atom's residue sits, opening a new one if the row starts one.
    fn place<R: AtomSiteRow + ?Sized>(
        &mut self,
        rows: &R,
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

        let repeats = self.names_in_residue.contains(atom_name, alt_id);
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
        self.names_in_residue.insert(atom_name, alt_id);
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

    fn report_ambiguous_boundary<R: AtomSiteRow + ?Sized>(&mut self, rows: &R) {
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
}

fn text_or_empty(value: Option<&str>) -> &str {
    let Some(value) = value else {
        return "";
    };
    value
}

fn identifier_or_zero(value: Option<u32>) -> u32 {
    let Some(value) = value else {
        return 0;
    };
    value
}
