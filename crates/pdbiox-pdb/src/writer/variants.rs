//! Writing the PQR and PDBQT extensions without duplicating topology logic.

use crate::writer::{
    PdbOptions, RequiredAtomFields, alt_field, atom_name, atom_name_field, chain_label,
    check_capacity, component_name, residue_field, residue_sequence, serial_field,
};
use pdbiox_core::annotation::{
    ATOM_RADIUS_ANNOTATION, AUTODOCK_TYPE_ANNOTATION, AnnotationColumn, AtomAnnotation,
    PARTIAL_CHARGE_ANNOTATION,
};
use pdbiox_core::column::Presence;
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::io::{SelectAll, TextOutput};
use pdbiox_core::structure::{AtomRef, ResidueRef, Structure};
use pdbiox_core::symbol::SymbolId;
use std::fmt;
use std::io::Write as IoWrite;

#[derive(Clone, Copy)]
enum Variant {
    Pqr,
    Pdbqt,
}

struct VariantColumns<'a> {
    charges: &'a AnnotationColumn<f64>,
    radii: Option<&'a AnnotationColumn<f64>>,
    atom_types: Option<&'a AnnotationColumn<SymbolId>>,
}

struct LineContext<'a> {
    chain: &'a str,
    options: &'a PdbOptions,
    variant: Variant,
    columns: &'a VariantColumns<'a>,
}

/// Writes PQR coordinates, requiring charge and radius annotations.
///
/// # Errors
///
/// Returns a conversion diagnostic when the structure exceeds fixed-column
/// capacity or lacks a present charge/radius value for any atom.
pub fn write_pqr(structure: &Structure, options: &PdbOptions) -> Result<String, Vec<Diagnostic>> {
    write_variant(structure, options, Variant::Pqr)
}

/// Writes a rigid PDBQT model, requiring charge and `AutoDock` type annotations.
///
/// # Errors
///
/// Returns a conversion diagnostic when the structure exceeds fixed-column
/// capacity or lacks a present charge/type value for any atom.
pub fn write_pdbqt(structure: &Structure, options: &PdbOptions) -> Result<String, Vec<Diagnostic>> {
    write_variant(structure, options, Variant::Pdbqt)
}

/// Streams PQR coordinates without retaining the complete text.
///
/// # Errors
///
/// Returns fixed-column, annotation, or destination diagnostics.
pub fn write_pqr_to<W: IoWrite>(
    structure: &Structure,
    options: &PdbOptions,
    output: &mut W,
) -> Result<(), Vec<Diagnostic>> {
    write_variant_to(structure, options, Variant::Pqr, output)
}

/// Streams a rigid PDBQT model without retaining the complete text.
///
/// # Errors
///
/// Returns fixed-column, annotation, or destination diagnostics.
pub fn write_pdbqt_to<W: IoWrite>(
    structure: &Structure,
    options: &PdbOptions,
    output: &mut W,
) -> Result<(), Vec<Diagnostic>> {
    write_variant_to(structure, options, Variant::Pdbqt, output)
}

fn write_variant(
    structure: &Structure,
    options: &PdbOptions,
    variant: Variant,
) -> Result<String, Vec<Diagnostic>> {
    let mut out = String::with_capacity(structure.atom_count() as usize * 82);
    render_variant(&mut out, structure, options, variant)?;
    Ok(out)
}

fn write_variant_to<W: IoWrite>(
    structure: &Structure,
    options: &PdbOptions,
    variant: Variant,
    output: &mut W,
) -> Result<(), Vec<Diagnostic>> {
    let mut output = TextOutput::new(output);
    render_variant(&mut output, structure, options, variant)?;
    output.finish().map_err(|error| {
        vec![Diagnostic::new(Code::E7901).with_context("reason", error.to_string())]
    })
}

fn render_variant(
    out: &mut impl fmt::Write,
    structure: &Structure,
    options: &PdbOptions,
    variant: Variant,
) -> Result<(), Vec<Diagnostic>> {
    let required = match variant {
        Variant::Pqr => RequiredAtomFields::Pqr,
        Variant::Pdbqt => RequiredAtomFields::Pdbqt,
    };
    let refusals = check_capacity(structure, options, &SelectAll, required);
    if !refusals.is_empty() {
        return Err(refusals);
    }
    let charges = real_column(structure, PARTIAL_CHARGE_ANNOTATION)?;
    let radii = match variant {
        Variant::Pqr => Some(real_column(structure, ATOM_RADIUS_ANNOTATION)?),
        Variant::Pdbqt => None,
    };
    let atom_types = match variant {
        Variant::Pdbqt => Some(symbol_column(structure, AUTODOCK_TYPE_ANNOTATION)?),
        Variant::Pqr => None,
    };
    let columns = VariantColumns {
        charges,
        radii,
        atom_types,
    };
    require_present(structure, &columns)?;

    if matches!(variant, Variant::Pdbqt) {
        let _ = out.write_str("ROOT\n");
    }
    let mut serial = 1_i64;
    for chain in structure.data().chains() {
        let original_label = chain_label(&chain, options.identifier_namespace());
        let label = options.label_for(original_label);
        let context = LineContext {
            chain: label,
            options,
            variant,
            columns: &columns,
        };
        for residue in chain.residues() {
            for atom in residue.atoms() {
                if let Err(finding) = write_atom(out, structure, &atom, &residue, serial, &context)
                {
                    return Err(vec![finding]);
                }
                serial += 1;
            }
        }
    }
    match variant {
        Variant::Pqr => {
            let _ = out.write_str("END\n");
        }
        Variant::Pdbqt => {
            let _ = out.write_str("ENDROOT\nTORSDOF 0\n");
        }
    }
    Ok(())
}

fn write_atom(
    out: &mut impl fmt::Write,
    structure: &Structure,
    atom: &AtomRef<'_>,
    residue: &ResidueRef<'_>,
    serial: i64,
    context: &LineContext<'_>,
) -> Result<(), Diagnostic> {
    let position = atom
        .position()
        .ok_or_else(|| missing_field("coordinates", *atom))?;
    let record = if residue.is_het() { "HETATM" } else { "ATOM  " };
    let name = atom_name(*atom, context.options.identifier_namespace())
        .ok_or_else(|| missing_field("atom identifier", *atom))?;
    let name = atom_name_field(name);
    let component = component_name(*atom, *residue, context.options.identifier_namespace())
        .ok_or_else(|| missing_field("component identifier", *atom))?;
    let sequence = i64::from(
        residue_sequence(*residue, context.options.identifier_namespace())
            .ok_or_else(|| missing_field("residue sequence identifier", *atom))?,
    );
    let charge = present_real(context.columns.charges, atom.index().get())
        .ok_or_else(|| missing_field(PARTIAL_CHARGE_ANNOTATION, *atom))?;
    let _ = write!(
        out,
        "{record}{serial:>5} {name:<4}{alt:1}{component:>3} {chain:>1}{sequence}{ins:1}   \
         {x:8.3}{y:8.3}{z:8.3}",
        chain = context.chain,
        serial = serial_field(serial, context.options),
        alt = alt_field(structure, atom),
        sequence = residue_field(sequence, context.options),
        ins = match residue.ins_code() {
            Some(code) => code,
            None => " ",
        },
        x = f64::from(position[0]),
        y = f64::from(position[1]),
        z = f64::from(position[2]),
    );
    match context.variant {
        Variant::Pqr => {
            let radius = context
                .columns
                .radii
                .and_then(|column| present_real(column, atom.index().get()));
            let _ = writeln!(
                out,
                "{charge:8.4}{radius:7.4}",
                charge = charge,
                radius = radius.ok_or_else(|| missing_field(ATOM_RADIUS_ANNOTATION, *atom))?
            );
        }
        Variant::Pdbqt => {
            let atom_type = context
                .columns
                .atom_types
                .and_then(|column| present_symbol(column, atom.index().get()))
                .and_then(|symbol| structure.resolve(symbol))
                .ok_or_else(|| missing_field(AUTODOCK_TYPE_ANNOTATION, *atom))?;
            let occupancy = atom
                .occupancy()
                .ok_or_else(|| missing_field("occupancy", *atom))?;
            let b_factor = atom
                .b_factor()
                .ok_or_else(|| missing_field("B factor", *atom))?;
            let _ = writeln!(
                out,
                "{occ:6.2}{b:6.2}      {charge:6.3} {atom_type:<2}",
                occ = f64::from(occupancy),
                b = f64::from(b_factor),
                charge = charge,
                atom_type = atom_type,
            );
        }
    }
    Ok(())
}

fn real_column<'a>(
    structure: &'a Structure,
    name: &str,
) -> Result<&'a AnnotationColumn<f64>, Vec<Diagnostic>> {
    match structure.annotations().get(name) {
        Some(AtomAnnotation::Real(column)) => Ok(column),
        _ => Err(missing(name)),
    }
}

fn symbol_column<'a>(
    structure: &'a Structure,
    name: &str,
) -> Result<&'a AnnotationColumn<SymbolId>, Vec<Diagnostic>> {
    match structure.annotations().get(name) {
        Some(AtomAnnotation::Symbol(column)) => Ok(column),
        _ => Err(missing(name)),
    }
}

fn require_present(
    structure: &Structure,
    columns: &VariantColumns<'_>,
) -> Result<(), Vec<Diagnostic>> {
    for atom in 0..structure.atom_count() {
        if present_real(columns.charges, atom).is_none()
            || columns
                .radii
                .is_some_and(|column| present_real(column, atom).is_none())
            || columns
                .atom_types
                .is_some_and(|column| present_symbol(column, atom).is_none())
        {
            return Err(missing("a required per-atom value"));
        }
    }
    Ok(())
}

fn present_real(column: &AnnotationColumn<f64>, atom: u32) -> Option<f64> {
    column
        .get(atom)
        .and_then(|(value, presence)| (presence == Presence::Present).then_some(value))
}

fn present_symbol(column: &AnnotationColumn<SymbolId>, atom: u32) -> Option<SymbolId> {
    column
        .get(atom)
        .and_then(|(value, presence)| (presence == Presence::Present).then_some(value))
}

fn missing(annotation: &str) -> Vec<Diagnostic> {
    vec![
        Diagnostic::new(Code::E4105)
            .with_context("required annotation", annotation)
            .with_context("target", "PDB-derived charged coordinate format"),
    ]
}

fn missing_field(field: &str, atom: AtomRef<'_>) -> Diagnostic {
    Diagnostic::new(Code::E4105)
        .with_context("required field", field)
        .with_context("atom", atom.index().to_string())
}

#[cfg(test)]
#[path = "variants_tests.rs"]
mod tests;
