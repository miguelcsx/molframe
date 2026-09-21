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

use crate::hybrid36;
use crate::{PDB_HEADERS_EXTENSION, PdbHeaders};
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::io::{Select, SelectAll};
use molframe_core::structure::{AtomRef, ChainRef, ResidueRef, Structure};
use std::collections::HashMap;
use std::fmt;

#[path = "identity.rs"]
mod identity;
pub use identity::PdbIdentifierNamespace;
pub(crate) use identity::{atom_name, chain_name, component_name, residue_sequence};

/// Bounds whose rounding would overflow an eight-column, three-decimal field.
/// The negative side loses one digit to the sign, so the limits are asymmetric.
const MIN_COORDINATE: f64 = -999.999_5;
const MAX_COORDINATE: f64 = 9_999.999_5;

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
}

impl Default for PdbOptions {
    fn default() -> Self {
        Self {
            chain_map: HashMap::new(),
            hybrid36: false,
            namespace: PdbIdentifierNamespace::Label,
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

    pub(crate) const fn identifier_namespace(&self) -> PdbIdentifierNamespace {
        self.namespace
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
    let refusals = check_capacity(structure, options, select, RequiredAtomFields::Pdb);
    if !refusals.is_empty() {
        return Err(refusals);
    }

    let data = structure.data();

    if let Some(headers) = data.extensions.get::<PdbHeaders>(PDB_HEADERS_EXTENSION) {
        for record in headers.records() {
            let _ = writeln!(out, "{}", record.line());
        }
    } else {
        write_generated_metadata(out, structure);
    }

    let mut serial = 1i64;
    for chain in data.chains() {
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
                if let Err(finding) =
                    write_atom(out, structure, &atom, &residue, label, serial, options)
                {
                    return Err(vec![finding]);
                }
                serial += 1;
            }
        }
        let _ = writeln!(out, "TER   {:>5}", serial_field(serial, options));
        serial += 1;
    }
    let _ = out.write_str("END\n");
    Ok(())
}

fn write_generated_metadata(out: &mut impl fmt::Write, structure: &Structure) {
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
    if let Some(cell) = data.cell {
        let _ = writeln!(
            out,
            "CRYST1{:9.3}{:9.3}{:9.3}{:7.2}{:7.2}{:7.2}",
            cell.lengths[0],
            cell.lengths[1],
            cell.lengths[2],
            cell.angles[0],
            cell.angles[1],
            cell.angles[2],
        );
    }
}

fn write_atom(
    out: &mut impl fmt::Write,
    structure: &Structure,
    atom: &AtomRef<'_>,
    residue: &ResidueRef<'_>,
    chain: &str,
    serial: i64,
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
    let position = required_value(atom.position(), "coordinates", *atom, options)?;
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

#[derive(Clone, Copy)]
pub(crate) enum RequiredAtomFields {
    Pdb,
    Pqr,
    Pdbqt,
}

/// Everything about the structure that the format cannot hold.
pub(crate) fn check_capacity(
    structure: &Structure,
    options: &PdbOptions,
    select: &impl Select,
    required: RequiredAtomFields,
) -> Vec<Diagnostic> {
    let data = structure.data();
    let mut refusals = Vec::new();

    for chain in data.chains() {
        if !select.accept_chain(chain.index()) {
            continue;
        }
        let original = chain_label(&chain, options.namespace);
        let written = options.label_for(original);
        if written.chars().count() > 1 {
            refusals.push(
                Diagnostic::new(Code::E4102)
                    .with_context("chain", original.to_string())
                    .with_context("written as", written.to_string()),
            );
        }
        for residue in chain.residues() {
            if !select.accept_residue(residue.index()) {
                continue;
            }
            check_residue(&residue, options, select, required, &mut refusals);
        }
    }

    let atoms = i64::from(structure.atom_count());
    if !options.hybrid36 && hybrid36::needs_hybrid36(atoms, 5) {
        refusals.push(
            Diagnostic::new(Code::E4101)
                .with_context("atoms", atoms.to_string())
                .with_context("field holds", "99999"),
        );
    }
    refusals
}

fn check_residue(
    residue: &ResidueRef<'_>,
    options: &PdbOptions,
    select: &impl Select,
    required: RequiredAtomFields,
    refusals: &mut Vec<Diagnostic>,
) {
    if let Some(seq) = residue_sequence(*residue, options.namespace)
        && !options.hybrid36
        && hybrid36::needs_hybrid36(i64::from(seq), 4)
    {
        refusals.push(
            Diagnostic::new(Code::E4103)
                .with_context("residue", seq.to_string())
                .with_context("field holds", "9999"),
        );
    }
    for atom in residue.atoms() {
        if !select.accept_atom(atom.index()) {
            continue;
        }
        check_required_atom_fields(atom, *residue, options, required, refusals);
        let Some(position) = atom.position() else {
            continue;
        };
        if position.iter().any(|value| {
            let value = f64::from(*value);
            !value.is_finite() || !(MIN_COORDINATE..MAX_COORDINATE).contains(&value)
        }) {
            refusals.push(
                Diagnostic::new(Code::E4104)
                    .with_context("atom", atom.index().to_string())
                    .with_context("position", format!("{position:?}")),
            );
        }
    }
}

fn check_required_atom_fields(
    atom: AtomRef<'_>,
    residue: ResidueRef<'_>,
    options: &PdbOptions,
    required: RequiredAtomFields,
    refusals: &mut Vec<Diagnostic>,
) {
    let common = [
        (atom.position().is_some(), "coordinates"),
        (
            atom_name(atom, options.namespace).is_some(),
            "atom identifier",
        ),
        (
            component_name(atom, residue, options.namespace).is_some(),
            "component identifier",
        ),
        (
            residue_sequence(residue, options.namespace).is_some(),
            "residue sequence identifier",
        ),
    ];
    for (present, field) in common {
        if !present {
            refusals.push(missing_field(field, atom, options));
        }
    }
    if matches!(
        required,
        RequiredAtomFields::Pdb | RequiredAtomFields::Pdbqt
    ) {
        for (present, field) in [
            (atom.occupancy().is_some(), "occupancy"),
            (atom.b_factor().is_some(), "B factor"),
        ] {
            if !present {
                refusals.push(missing_field(field, atom, options));
            }
        }
    }
}

fn required_text<'a>(
    value: Option<&'a str>,
    field: &'static str,
    atom: AtomRef<'_>,
    options: &PdbOptions,
) -> Result<&'a str, Diagnostic> {
    value.ok_or_else(|| missing_field(field, atom, options))
}

fn required_value<T>(
    value: Option<T>,
    field: &'static str,
    atom: AtomRef<'_>,
    options: &PdbOptions,
) -> Result<T, Diagnostic> {
    value.ok_or_else(|| missing_field(field, atom, options))
}

fn missing_field(field: &'static str, atom: AtomRef<'_>, options: &PdbOptions) -> Diagnostic {
    Diagnostic::new(Code::E4105)
        .with_context("required field", field)
        .with_context("atom", atom.index().to_string())
        .with_context("namespace", options.namespace.as_str())
}

#[cfg(test)]
#[path = "pdb_tests.rs"]
mod tests;
