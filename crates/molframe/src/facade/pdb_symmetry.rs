//! Symmetry and assemblies for legacy PDB files.
//!
//! The format names its space group only in `CRYST1` and its biological
//! assemblies only in `REMARK 350`; both resolve here into the extensions a CIF
//! read attaches, so a legacy file is as usable crystallographically.

use molframe_core::{Diagnostic, Structure};

/// Attaches the symmetry set and the assemblies a legacy file declares.
///
/// An unknown or absent space-group symbol leaves the symmetry unattached: a
/// name this catalogue does not know is not evidence of any symmetry. An
/// assembly the file declares inconsistently is reported, not guessed.
#[cfg(feature = "crystal")]
pub(super) fn attach_pdb_symmetry(
    structure: Structure,
    findings: Vec<Diagnostic>,
) -> (Structure, Vec<Diagnostic>) {
    let (mut structure, mut findings) = attach_space_group(structure, findings);
    match biomolecule_assemblies(&structure) {
        Ok(Some(assemblies)) => {
            structure = structure.with_extension(molframe_xtal::ASSEMBLIES_EXTENSION, assemblies);
        }
        Ok(None) => {}
        Err(finding) => findings.push(finding),
    }
    (structure, findings)
}

#[cfg(feature = "crystal")]
fn attach_space_group(
    structure: Structure,
    findings: Vec<Diagnostic>,
) -> (Structure, Vec<Diagnostic>) {
    let Some(symbol) = structure.data().entry.space_group.clone() else {
        return (structure, findings);
    };
    match molframe_xtal::space_group_by_hermann_mauguin(&symbol) {
        Ok(setting) => (
            structure.with_extension(molframe_xtal::SYMMETRY_EXTENSION, setting.symmetry_set()),
            findings,
        ),
        Err(_) => (structure, findings),
    }
}

/// The assemblies `REMARK 350` defines, as the set an mmCIF read would build.
///
/// Operator identifiers are `<assembly>.<group>.<serial>`, because BIOMT
/// serials restart for every assembly and every chain group.
#[cfg(feature = "crystal")]
fn biomolecule_assemblies(
    structure: &Structure,
) -> Result<Option<molframe_xtal::AssemblySet>, Diagnostic> {
    use molframe_pdb::PdbHeadersExt;
    use molframe_xtal::{AssemblyDef, AssemblySet, Generator, OperExpression, Operator};

    let Some(headers) = structure.pdb_headers() else {
        return Ok(None);
    };
    let biomolecules = headers.biomolecules();
    if biomolecules.is_empty() {
        return Ok(None);
    }
    let mut definitions = Vec::new();
    let mut operators = Vec::new();
    for biomolecule in &biomolecules {
        let mut generators = Vec::new();
        for (group_number, group) in biomolecule.groups.iter().enumerate() {
            let ids: Vec<String> = group
                .operations
                .iter()
                .map(|operation| {
                    format!(
                        "{}.{}.{}",
                        biomolecule.id,
                        group_number + 1,
                        operation.serial
                    )
                })
                .collect();
            for (id, operation) in ids.iter().zip(&group.operations) {
                operators.push(Operator::from_matrix(
                    id.as_str(),
                    operation.rotation,
                    operation.translation,
                ));
            }
            if ids.is_empty() {
                continue;
            }
            generators.push(Generator {
                oper_expression: OperExpression::parse(&ids.join(","))?,
                asym_ids: group.chains.clone().into(),
            });
        }
        definitions.push(AssemblyDef {
            id: biomolecule.id.to_string().into(),
            details: None,
            method: Some("REMARK 350".into()),
            oligomeric: None,
            generators,
        });
    }
    AssemblySet::from_definitions(definitions, operators).map(Some)
}

#[cfg(not(feature = "crystal"))]
pub(super) fn attach_pdb_symmetry(
    structure: Structure,
    findings: Vec<Diagnostic>,
) -> (Structure, Vec<Diagnostic>) {
    (structure, findings)
}
