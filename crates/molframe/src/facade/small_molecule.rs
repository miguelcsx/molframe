//! SDF, MOL2 and core-CIF small molecules through the format facade.
//!
//! Each is one molecule in one residue, as `molframe_chem` lowers them. A
//! multi-record file reads its first record and says so.

use molframe_chem::{
    MolAtomMetadata, MolBondMetadata, MolRecord, MolVersion, mol_record_to_structure,
    mol2_record_to_structure, parse_mol_record, parse_mol2_record, parse_sdf_records,
    structure_to_molecule, write_sdf,
};
use molframe_core::io::InputBuffer;
use molframe_core::{Code, Diagnostic, Structure as CoreStructure};
use std::io::Write;

use super::Findings;

type Read = Result<(CoreStructure, Vec<Diagnostic>), Findings>;

fn text_of(input: &InputBuffer) -> Result<&str, Findings> {
    std::str::from_utf8(input.as_bytes()).map_err(|error| {
        Findings::from(
            Diagnostic::new(Code::E1102)
                .with_message("the file is not valid UTF-8 text")
                .with_context("decoder", error.to_string()),
        )
    })
}

fn malformed(format: &str, error: &dyn std::fmt::Display) -> Findings {
    Findings::from(
        Diagnostic::new(Code::E1102)
            .with_message(format!("{format} could not be decoded"))
            .with_context("decoder", error.to_string()),
    )
}

/// Reads the first SDF (or MOL) record.
pub(super) fn read_sdf(input: &InputBuffer) -> Read {
    let text = text_of(input)?;
    let records = if text.contains("$$$$") {
        parse_sdf_records(text).map_err(|error| malformed("SDF", &error))?
    } else {
        vec![parse_mol_record(text).map_err(|error| malformed("MOL", &error))?]
    };
    let Some(first) = records.first() else {
        return Err(malformed("SDF", &"the file holds no records"));
    };
    let mut findings = Vec::new();
    if records.len() > 1 {
        findings.push(
            Diagnostic::new(Code::W1001)
                .with_message("the SDF holds several records; only the first is read")
                .with_context("records", records.len().to_string()),
        );
    }
    let structure = mol_record_to_structure(first).map_err(Findings::from)?;
    Ok((structure, findings))
}

/// Reads the first MOL2 molecule.
pub(super) fn read_mol2(input: &InputBuffer) -> Read {
    let text = text_of(input)?;
    let mut findings = Vec::new();
    let molecules = text.matches("@<TRIPOS>MOLECULE").count();
    let first = match text.match_indices("@<TRIPOS>MOLECULE").nth(1) {
        Some((second, _)) => &text[..second],
        None => text,
    };
    if molecules > 1 {
        findings.push(
            Diagnostic::new(Code::W1001)
                .with_message("the MOL2 file holds several molecules; only the first is read")
                .with_context("records", molecules.to_string()),
        );
    }
    let record = parse_mol2_record(first).map_err(|error| malformed("MOL2", &error))?;
    let structure = mol2_record_to_structure(&record).map_err(Findings::from)?;
    Ok((structure, findings))
}
/// Reads the first core-CIF block as a small-molecule structure.
#[cfg(feature = "mmcif")]
pub(super) fn read_small_cif(input: &InputBuffer) -> Read {
    let (document, findings) = molframe_cif::parse(input).map_err(Findings::from)?;
    // Legacy underscore tags are only read when the file uses them.
    let legacy = input
        .as_bytes()
        .windows(b"_atom_site_fract_x".len())
        .any(|window| window == b"_atom_site_fract_x");
    let options = if legacy {
        molframe_cif::SmallCifOptions::ddl1()
    } else {
        molframe_cif::SmallCifOptions::ddl2()
    };
    let model = molframe_cif::lower_small_cif_with_options(&document, options)
        .map_err(|error| malformed("small-molecule CIF", &error))?;
    let structure = molframe_chem::small_cif_to_structure(&model).map_err(Findings::from)?;
    Ok((structure, findings))
}

/// Writes the structure as one SDF record, with the formal charges it carries.
pub(super) fn write_sdf_to<W: Write>(
    structure: &CoreStructure,
    output: &mut W,
) -> Result<(), Findings> {
    let molecule = structure_to_molecule(structure).map_err(Findings::from)?;
    let name = structure
        .data()
        .residues()
        .next()
        .and_then(molframe_core::structure::ResidueRef::name)
        .map_or_else(String::new, str::to_owned);
    let atom_metadata = structure
        .data()
        .atoms()
        .map(|atom| MolAtomMetadata {
            formal_charge: atom.formal_charge().filter(|charge| *charge != 0),
            ..MolAtomMetadata::default()
        })
        .collect();
    let record = MolRecord {
        name: name.into(),
        program: "molframe".into(),
        comment: "".into(),
        version: MolVersion::V2000,
        bond_metadata: vec![MolBondMetadata::default(); molecule.bonds.len()],
        molecule,
        atom_metadata,
        properties: Vec::new(),
    };
    let text = write_sdf(&[record]).map_err(|error| {
        Findings::from(
            Diagnostic::new(Code::E1102)
                .with_message("the molecule cannot be written as SDF")
                .with_context("cause", error.to_string()),
        )
    })?;
    output.write_all(text.as_bytes()).map_err(|error| {
        Findings::from(
            Diagnostic::new(Code::E1102)
                .with_message("the SDF could not be written")
                .with_context("cause", error.to_string()),
        )
    })
}
