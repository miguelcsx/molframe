//! Writing the legacy fixed-column format, or refusing to.
//!
//! The format's fields have fixed widths, and a structure can outgrow every one
//! of them: a chain label of two characters, a serial past ninety-nine thousand,
//! a residue number past ten thousand, a coordinate too large for eight columns.
//!
//! Every one of those is checked before a single byte is written, and a
//! structure that does not fit is refused with a diagnostic naming the field and
//! offering a way out. Truncating a chain label from `AA` to `A` produces a file
//! that opens cleanly in any viewer and describes a different molecule, which is
//! the failure this library exists to prevent.

use super::capacity::{RequiredAtomFields, check_capacity, required_text, required_value};
use crate::hybrid36;
use crate::{PDB_HEADERS_EXTENSION, PdbHeaders};
use molframe_core::diagnostic::Diagnostic;
use molframe_core::io::{Select, SelectAll};
use molframe_core::structure::{AtomRef, ChainRef, ResidueRef, Structure};
use std::collections::HashMap;
use std::fmt;

#[path = "identity.rs"]
mod identity;
pub use identity::PdbIdentifierNamespace;
pub(crate) use identity::{atom_name, chain_name, component_name, residue_sequence};

/// How to write, and what to do about a structure that does not fit.
#[derive(Clone, Debug)]
pub struct PdbOptions {
    /// Rename chains on the way out, so a structure with labels the format
    /// cannot hold can still be written deliberately.
    chain_map: HashMap<Box<str>, Box<str>>,
    /// Write serials and residue numbers past the decimal range in hybrid-36.
    ///
    /// Off by default: a consumer that does not implement the scheme would read
    /// a wrong number without noticing, so it is offered rather than assumed.
    hybrid36: bool,
    /// Identifier namespace used consistently for every deposited identity.
    namespace: PdbIdentifierNamespace,
    /// Space group to name on `CRYST1` when the entry records none, for
    /// structures whose symbol lives in a typed extension rather than the entry.
    space_group: Option<Box<str>>,
    /// Whole records written after the generated header, such as the `MTRIX`
    /// lines of a non-crystallographic operator set the format crate cannot see.
    extra_records: Option<Box<str>>,
}

impl Default for PdbOptions {
    fn default() -> Self {
        Self {
            chain_map: HashMap::new(),
            hybrid36: false,
            namespace: PdbIdentifierNamespace::Label,
            space_group: None,
            extra_records: None,
        }
    }
}

impl PdbOptions {
    /// Options that refuse anything the format cannot hold.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Renames one chain on the way out.
    #[must_use]
    pub fn chain_map(mut self, from: impl Into<Box<str>>, to: impl Into<Box<str>>) -> Self {
        self.chain_map.insert(from.into(), to.into());
        self
    }

    /// Writes counts past the decimal range in hybrid-36 rather than refusing.
    #[must_use]
    pub const fn hybrid36(mut self, enabled: bool) -> Self {
        self.hybrid36 = enabled;
        self
    }

    /// Selects one identifier namespace for chains, residues and atoms.
    #[must_use]
    pub const fn namespace(mut self, namespace: PdbIdentifierNamespace) -> Self {
        self.namespace = namespace;
        self
    }

    /// Names the space group `CRYST1` carries when the entry records none.
    #[must_use]
    pub fn with_space_group(mut self, symbol: impl Into<Box<str>>) -> Self {
        self.space_group = Some(symbol.into());
        self
    }

    /// Adds whole records after the generated header, one per line.
    #[must_use]
    pub fn with_extra_records(mut self, records: impl Into<Box<str>>) -> Self {
        self.extra_records = Some(records.into());
        self
    }

    /// The extra records, when the caller supplied any.
    #[must_use]
    pub fn extra_records(&self) -> Option<&str> {
        self.extra_records.as_deref()
    }

    pub(crate) const fn identifier_namespace(&self) -> PdbIdentifierNamespace {
        self.namespace
    }

    pub(crate) const fn uses_hybrid36(&self) -> bool {
        self.hybrid36
    }

    /// The label a chain is written under.
    pub(crate) fn label_for<'a>(&'a self, label: &'a str) -> &'a str {
        match self.chain_map.get(label) {
            Some(replacement) => replacement,
            None => label,
        }
    }
}

/// Writes a structure, or explains why it cannot be written.
///
/// # Errors
///
/// Returns every capacity the structure exceeds, so a caller sees the whole
/// problem rather than fixing one field at a time.
pub fn write(structure: &Structure, options: &PdbOptions) -> Result<String, Vec<Diagnostic>> {
    write_selected(structure, options, &SelectAll)
}

/// Writes the parts of a structure a filter accepts.
///
/// # Errors
///
/// Returns every capacity the accepted parts exceed.
pub fn write_selected(
    structure: &Structure,
    options: &PdbOptions,
    select: &impl Select,
) -> Result<String, Vec<Diagnostic>> {
    let mut out = String::with_capacity(structure.atom_count() as usize * 81);
    render_selected(&mut out, structure, options, select)?;
    Ok(out)
}

pub(super) fn render_selected(
    out: &mut impl fmt::Write,
    structure: &Structure,
    options: &PdbOptions,
    select: &impl Select,
) -> Result<(), Vec<Diagnostic>> {
    let mut refusals = check_capacity(structure, options, select, RequiredAtomFields::Pdb);
    refusals.extend(super::secondary::preflight(structure, options, select));
    if !refusals.is_empty() {
        return Err(refusals);
    }

    let data = structure.data();

    let mut has_ssbond = false;
    let mut has_link = false;
    if let Some(headers) = data.extensions.get::<PdbHeaders>(PDB_HEADERS_EXTENSION) {
        has_ssbond = headers
            .records()
            .iter()
            .any(|record| record.name() == "SSBOND");
        has_link = headers
            .records()
            .iter()
            .any(|record| record.name() == "LINK");
        for record in headers.records() {
            if matches!(record.name(), "HELIX" | "SHEET")
                && structure
                    .secondary_structure()
                    .iter()
                    .any(|state| state.is_helix() || state.is_strand())
            {
                continue;
            }
            let _ = writeln!(out, "{}", record.line());
        }
    } else {
        write_generated_metadata(out, structure, options);
    }
    super::secondary::write(out, structure, options, select);
    super::connections::write_links(out, structure, options, select, has_ssbond, has_link);

    let track = super::connections::has_connect(structure);
    let mut serials: HashMap<u32, i64> = HashMap::new();
    let count = structure.model_count();
    if count <= 1 {
        write_model(
            out,
            structure,
            None,
            options,
            select,
            track.then_some(&mut serials),
        )?;
    } else {
        let mut written = 0usize;
        for number in 0..count {
            let Ok(index) = u32::try_from(number) else {
                break;
            };
            let model = molframe_core::index::ModelIndex::new(index);
            if !select.accept_model(model) {
                continue;
            }
            let Some((snapshot, local)) = structure.model_snapshot(model) else {
                continue;
            };
            // Indices in a filter are those of the first model. Ragged models
            // are independent topologies, so applying them to a later model
            // would select unrelated atoms; refuse rather than guess.
            if structure.is_ragged_ensemble() && written > 0 {
                super::models::refuse_foreign_filter(structure, select)?;
            }
            written += 1;
            // Dense models share the topology and differ only in positions;
            // ragged ones are structures of their own and are checked as such.
            let frame = if structure.is_ragged_ensemble() {
                let refusals = check_capacity(&snapshot, options, select, RequiredAtomFields::Pdb);
                if !refusals.is_empty() {
                    return Err(refusals);
                }
                None
            } else {
                let frame = structure.model_positions(local);
                if let Some(frame) = frame {
                    super::models::check_frame(&snapshot, frame, number, select)?;
                }
                frame
            };
            let _ = writeln!(out, "MODEL     {:>4}", number + 1);
            let tracked = (number == 0 && track).then_some(&mut serials);
            write_model(out, &snapshot, frame, options, select, tracked)?;
            let _ = writeln!(out, "ENDMDL");
        }
    }
    super::connections::write_connect(out, structure, options, select, &serials);
    let _ = out.write_str("END\n");
    Ok(())
}

/// One model's atoms with their `TER` records.
fn write_model(
    out: &mut impl fmt::Write,
    structure: &Structure,
    frame: Option<&[[f32; 3]]>,
    options: &PdbOptions,
    select: &impl Select,
    mut serials: Option<&mut HashMap<u32, i64>>,
) -> Result<(), Vec<Diagnostic>> {
    let mut serial = 1i64;
    for chain in structure.data().chains() {
        if !select.accept_chain(chain.index()) {
            continue;
        }
        let label = chain_label(&chain, options.namespace);
        let label = options.label_for(label);

        for residue in chain.residues() {
            if !select.accept_residue(residue.index()) {
                continue;
            }
            for atom in residue.atoms() {
                if !select.accept_atom(atom.index()) {
                    continue;
                }
                let position = match frame {
                    Some(frame) => atom
                        .position()
                        .and(frame.get(atom.index().as_usize()).copied()),
                    None => atom.position(),
                };
                if let Err(finding) = write_atom(
                    out, structure, &atom, &residue, label, serial, position, options,
                ) {
                    return Err(vec![finding]);
                }
                if let Some(serials) = serials.as_deref_mut() {
                    serials.insert(atom.index().get(), serial);
                }
                serial += 1;
            }
        }
        let _ = writeln!(out, "TER   {:>5}", serial_field(serial, options));
        serial += 1;
    }
    Ok(())
}

fn write_generated_metadata(
    out: &mut impl fmt::Write,
    structure: &Structure,
    options: &PdbOptions,
) {
    let data = structure.data();
    if let Some(id) = &data.entry.id {
        let _ = writeln!(out, "HEADER    {:<40}{:<9}   {:>4}", "", "", id);
    }
    if let Some(title) = &data.entry.title {
        let _ = writeln!(out, "TITLE     {title}");
    }
    if let Some(method) = &data.entry.method {
        let _ = writeln!(out, "EXPDTA    {method}");
    }
    if let Some(resolution) = data.entry.resolution {
        let _ = writeln!(out, "REMARK   2");
        let _ = writeln!(out, "REMARK   2 RESOLUTION. {resolution:>7.2} ANGSTROMS.");
    }
    if let Some(cell) = data.cell {
        let _ = writeln!(
            out,
            "CRYST1{:9.3}{:9.3}{:9.3}{:7.2}{:7.2}{:7.2} {:<11}",
            cell.lengths[0],
            cell.lengths[1],
            cell.lengths[2],
            cell.angles[0],
            cell.angles[1],
            cell.angles[2],
            // CRYST1 names P 1 when the entry records no space group.
            match data
                .entry
                .space_group
                .as_deref()
                .or(options.space_group.as_deref())
            {
                Some(symbol) => symbol,
                None => "P 1",
            },
        );
        super::connections::write_scale(out, &cell);
    }
    if let Some(extra) = options.extra_records() {
        for line in extra.lines() {
            let _ = writeln!(out, "{line}");
        }
    }
}

fn write_atom(
    out: &mut impl fmt::Write,
    structure: &Structure,
    atom: &AtomRef<'_>,
    residue: &ResidueRef<'_>,
    chain: &str,
    serial: i64,
    position: Option<[f32; 3]>,
    options: &PdbOptions,
) -> Result<(), Diagnostic> {
    let record = if residue.is_het() { "HETATM" } else { "ATOM  " };
    let name = required_text(
        atom_name(*atom, options.namespace),
        "atom identifier",
        *atom,
        options,
    )?;
    let element = match atom.element() {
        Some(element) => element.symbol(),
        None => "",
    };
    // An atom whose position was never recorded still has a row here, but the
    // format has no way to say "unrecorded" — so it is left out of the file
    // rather than written as an origin the experiment never observed.
    let position = required_value(position, "coordinates", *atom, options)?;
    let seq = i64::from(required_value(
        residue_sequence(*residue, options.namespace),
        "residue sequence identifier",
        *atom,
        options,
    )?);
    let component = required_text(
        component_name(*atom, *residue, options.namespace),
        "component identifier",
        *atom,
        options,
    )?;
    let occupancy = required_value(atom.occupancy(), "occupancy", *atom, options)?;
    let b_factor = required_value(atom.b_factor(), "B factor", *atom, options)?;

    let _ = writeln!(
        out,
        "{record}{serial:>5} {name:<4}{alt:1}{comp:>3} {chain:>1}{seq}{ins:1}   \
         {x:8.3}{y:8.3}{z:8.3}{occ:6.2}{b:6.2}          {element:>2}",
        serial = serial_field(serial, options),
        name = atom_name_field(name),
        alt = alt_field(structure, atom),
        comp = component,
        seq = residue_field(seq, options),
        ins = insertion_code(*residue),
        x = f64::from(position[0]),
        y = f64::from(position[1]),
        z = f64::from(position[2]),
        occ = f64::from(occupancy),
        b = f64::from(b_factor),
    );
    // A tensor recorded for this atom follows under the same serial and the
    // same identity columns, as the record pairs by serial.
    if let Some(u) = structure.data().anisotropy.for_atom(atom.index()) {
        super::aniso::Record {
            atom,
            structure,
            residue,
            name,
            component,
            chain,
            seq: &residue_field(seq, options),
            serial: &serial_field(serial, options),
            u: &u,
        }
        .write(out);
    }
    Ok(())
}

/// The four-column name field.
///
/// A name shorter than four characters starts in the second column, which is the
/// convention that lets a reader tell a one-letter element from a two-letter one.
pub(crate) fn atom_name_field(name: &str) -> String {
    if name.len() >= 4 {
        name.to_owned()
    } else {
        format!(" {name:<3}")
    }
}

/// The one-column alternate-location field.
///
/// A label longer than one character has already been refused by the capacity
/// check, so anything reaching here fits.
pub(crate) fn alt_field<'a>(structure: &'a Structure, atom: &AtomRef<'a>) -> &'a str {
    let Some(alt) = atom.alt_id() else { return " " };
    let Some(symbol) = alt.symbol() else {
        return " ";
    };
    match structure.resolve(symbol) {
        Some(label) => label,
        None => " ",
    }
}

pub(crate) fn serial_field(serial: i64, options: &PdbOptions) -> String {
    match options
        .hybrid36
        .then(|| hybrid36::encode(serial, 5))
        .flatten()
    {
        Some(field) => field,
        None => format!("{serial:>5}"),
    }
}

pub(crate) fn residue_field(seq: i64, options: &PdbOptions) -> String {
    match options.hybrid36.then(|| hybrid36::encode(seq, 4)).flatten() {
        Some(field) => field,
        None => format!("{seq:>4}"),
    }
}

pub(crate) fn chain_label<'a>(
    chain: &'a ChainRef<'a>,
    namespace: PdbIdentifierNamespace,
) -> &'a str {
    let Some(name) = chain_name(*chain, namespace) else {
        return "";
    };
    name
}

pub(crate) fn insertion_code(residue: ResidueRef<'_>) -> &str {
    let Some(code) = residue.ins_code() else {
        return " ";
    };
    code
}

#[cfg(test)]
#[path = "pdb_tests.rs"]
mod tests;
