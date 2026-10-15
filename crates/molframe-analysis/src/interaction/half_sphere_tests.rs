use super::half_sphere_exposure;
use molframe_chem::PolymerAtomRole;
use molframe_core::annotation::{AnnotationColumn, AtomAnnotation};
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;
use molframe_core::{ExecutionContext, Presence};
use molframe_spatial::SpatialBackend;

// Residue 1's side chain points along +x. Residue 2 sits on the +x side (upper),
// residue 3 on the −x side (lower).
const SOURCE: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C CA ALA A 1 0 0 0\n\
ATOM 2 C CB ALA A 1 1 0 0\n\
ATOM 3 C CA ALA A 2 5 0 0\n\
ATOM 4 C CB ALA A 2 6 0 0\n\
ATOM 5 C CA ALA A 3 -5 0 0\n\
ATOM 6 C CB ALA A 3 -6 0 0\n";

fn structure() -> Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => with_roles(&structure),
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

fn with_roles(structure: &Structure) -> Structure {
    let entries: Vec<_> = structure
        .data()
        .atoms()
        .map(|atom| {
            let role = match atom.name() {
                Some("CA") => PolymerAtomRole::PROTEIN_ALPHA_CARBON,
                Some("CB") => PolymerAtomRole::PROTEIN_BETA_CARBON,
                _ => PolymerAtomRole::UNKNOWN,
            };
            (role.code(), Presence::Present)
        })
        .collect();
    let mut data = structure.data().clone();
    let _ = data.annotations.insert(
        molframe_core::POLYMER_ATOM_ROLE_ANNOTATION,
        AtomAnnotation::Integer(AnnotationColumn::from_entries(entries).expect("small column")),
    );
    Structure::new(data)
}

#[test]
fn a_residue_counts_neighbours_on_each_side_of_its_side_chain() {
    let Ok(exposure) = half_sphere_exposure(
        &structure(),
        13.0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    assert_eq!(exposure.len(), 3);
    // Residue 1 (index 0): residue 2 is up (+x), residue 3 is down (−x).
    assert_eq!(exposure[0].residue.get(), 0);
    assert_eq!((exposure[0].upper, exposure[0].lower), (1, 1));
}

#[test]
fn residues_without_a_beta_carbon_are_skipped() {
    // Glycine has no CB, so it never appears in the result.
    let source = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C CA GLY A 1 0 0 0\n";
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let (structure, _) = match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok(result) => result,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    let Ok(exposure) = half_sphere_exposure(
        &with_roles(&structure),
        13.0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    assert!(exposure.is_empty());
}

/// The governed path must report a residue it skipped as missing coverage rather
/// than returning an empty table.
///
/// A caller whose profile supplies only backbone roles gets no row for any
/// residue, and an empty table with `Complete` quality would read as "no residue
/// is crowded" where the truth is "nothing was measured". Reporting the skipped
/// residues puts the question to `missing_atoms`, which is where every other
/// analysis leaves it, and the descriptor requires a side-chain direction so an
/// input with no direction at all has no answer instead of an empty one.
#[test]
fn the_governed_kernel_reports_skipped_residues_as_missing_coverage() {
    use crate::policy_execution::{analyse_structure, half_sphere_exposure_kernel};
    use molframe_core::contract::{AnalysisPolicy, Status};

    // Three residues: two with both roles, one glycine without a CB.
    let source = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C CA ALA A 1 0 0 0\n\
ATOM 2 C CB ALA A 1 1 0 0\n\
ATOM 3 C CA ALA A 2 5 0 0\n\
ATOM 4 C CB ALA A 2 6 0 0\n\
ATOM 5 C CA GLY A 3 9 0 0\n";
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let (structure, _) = match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok(result) => result,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    let annotated = with_roles(&structure);
    // Two of the three residues carry a direction; the glycine does not.
    let kernel = half_sphere_exposure_kernel(13.0, SpatialBackend::BruteForce);
    let policy = AnalysisPolicy::default();
    let result = analyse_structure(&annotated, &policy, &kernel, &ExecutionContext::default())
        .unwrap_or_else(|error| panic!("valid: {error}"));

    let rows = result.value().cloned().unwrap_or_default();
    assert_eq!(
        rows.len(),
        2,
        "only the two residues with a CB are measured"
    );
    assert_eq!(
        result.coverage.used, 2,
        "the two measured residues are covered"
    );
    assert_eq!(
        result.coverage.missing, 1,
        "the glycine is reported as missing rather than silently dropped"
    );
    assert_eq!(result.coverage.intended, 3);
    assert_eq!(result.status(), Status::Partial);
}

/// An input with no side-chain direction at all has no answer, not an empty one.
#[test]
fn a_structure_without_any_side_chain_direction_is_indeterminate() {
    use crate::policy_execution::{analyse_structure, half_sphere_exposure_kernel};
    use molframe_core::contract::AnalysisPolicy;

    // Alpha carbons only: no residue has an outward direction.
    let source = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C CA ALA A 1 0 0 0\n\
ATOM 2 C CA ALA A 2 5 0 0\n";
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let (structure, _) = match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok(result) => result,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    let annotated = with_roles(&structure);
    assert!(
        !crate::policy_execution::Requirement::PolymerSideChainDirection.satisfied_by(&annotated),
        "alpha carbons alone are not a side-chain direction"
    );

    let kernel = half_sphere_exposure_kernel(13.0, SpatialBackend::BruteForce);
    let policy = AnalysisPolicy::default();
    let result = analyse_structure(&annotated, &policy, &kernel, &ExecutionContext::default());
    let Ok(analysis) = result else {
        // Refusing outright is also correct: there is nothing to measure.
        return;
    };
    // If it runs, it must not claim a complete empty answer.
    assert!(
        !analysis.value().is_some_and(Vec::is_empty),
        "an empty table must never be reported as a complete answer"
    );
}

/// Solvent and other non-protein residues are outside the analysis's domain.
///
/// Half-sphere exposure describes how buried a protein residue is, so a bound
/// water is not a protein residue that failed to be measured. Counting it would
/// report a run that covered every measurable residue as half incomplete, and
/// `missing_atoms` would then reject a perfectly good answer.
#[test]
fn solvent_residues_are_outside_the_coverage_domain() {
    use crate::policy_execution::{analyse_structure, half_sphere_exposure_kernel};
    use molframe_core::contract::{AnalysisPolicy, Status};

    // Two protein residues carrying both roles, plus a water and an ion.
    let source = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C CA ALA A 1 0 0 0\n\
ATOM 2 C CB ALA A 1 1 0 0\n\
ATOM 3 C CA ALA A 2 5 0 0\n\
ATOM 4 C CB ALA A 2 6 0 0\n\
HETATM 5 O HOH B 3 20 0 0 0\n\
HETATM 6 ZN ZN C 4 30 0 0 0\n";
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let (structure, _) = match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok(result) => result,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    // The profile names CA and CB, so the two alanines are the protein domain. The
    // water and the zinc carry no protein role and are outside it.
    let annotated = with_roles(&structure);

    let kernel = half_sphere_exposure_kernel(13.0, SpatialBackend::BruteForce);
    let result = analyse_structure(
        &annotated,
        &AnalysisPolicy::default(),
        &kernel,
        &ExecutionContext::default(),
    )
    .unwrap_or_else(|error| panic!("valid: {error}"));

    assert_eq!(
        result.coverage.intended, 2,
        "only the two protein residues are in the domain; water and the ion are not"
    );
    assert_eq!(result.coverage.used, 2);
    assert_eq!(result.coverage.missing, 0);
    assert_eq!(
        result.status(),
        Status::Complete,
        "a run that measured every residue in its domain is complete"
    );
}
