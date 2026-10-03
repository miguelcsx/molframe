use super::*;
use crate::validation_inputs::ValidationInputError;
use crate::{AnalysisPolicy, ReadOptions, read_bytes};

const PDB: &str = "\
ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N
ATOM      2  CA  ALA A   1      12.560   6.195  -6.504  1.00  0.00           C
ATOM      3  C   ALA A   1      13.100   7.600  -6.504  1.00  0.00           C
ATOM      4  N   GLY B   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      5  CA  GLY B   1       2.000   2.500   3.200  1.00  0.00           C
END
";

fn structure() -> Structure {
    match read_bytes(PDB.into(), Some("t.pdb"), &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(error) => panic!("fixture failed: {error}"),
    }
}

fn write(text: &str) -> tempfile::NamedTempFile {
    let file = match tempfile::Builder::new().suffix(".toml").tempfile() {
        Ok(file) => file,
        Err(error) => panic!("temporary file failed: {error}"),
    };
    if let Err(error) = std::fs::write(file.path(), text) {
        panic!("write failed: {error}")
    }
    file
}

const ZERO: &str = "[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]";

fn group(id: &str, selection: &str) -> String {
    format!(
        "[[group]]\nid = '{id}'\nselection = '{selection}'\norigin = [1.0, 2.0, 3.0]\nt = [[0.1, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.3]]\nl = {ZERO}\ns = {ZERO}\n"
    )
}

#[test]
fn planes_resolve_their_selections_in_the_structure() {
    let file = write(
        "[[plane]]\nid = 'a'\nselection = 'chain A'\n[[plane]]\nid = 'b'\nselection = 'chain B'\n",
    );
    let planes = match read_plane_restraints(file.path(), &structure(), &AnalysisPolicy::default())
    {
        Ok(planes) => planes,
        Err(error) => panic!("planes failed: {error}"),
    };
    assert_eq!(planes.len(), 2);
    assert_eq!(planes[0].id, "a");
    assert_eq!(planes[0].atoms.len(), 3);
    assert_eq!(planes[1].atoms.len(), 2);
}

#[test]
fn a_selection_matching_nothing_is_refused_rather_than_covering_nothing() {
    let file = write("[[plane]]\nid = 'a'\nselection = 'chain Z'\n");
    assert!(matches!(
        read_plane_restraints(file.path(), &structure(), &AnalysisPolicy::default()),
        Err(ValidationInputError::InvalidValue { .. })
    ));
}

#[test]
fn a_selection_that_does_not_parse_names_its_plane() {
    let file = write("[[plane]]\nid = 'broken'\nselection = 'chain ((('\n");
    match read_plane_restraints(file.path(), &structure(), &AnalysisPolicy::default()) {
        Err(ValidationInputError::Selection { id, .. }) => assert_eq!(id, "broken"),
        other => panic!("expected a selection error, got {other:?}"),
    }
}

#[test]
fn a_repeated_plane_id_is_refused() {
    let file = write(
        "[[plane]]\nid = 'a'\nselection = 'chain A'\n[[plane]]\nid = 'a'\nselection = 'chain B'\n",
    );
    assert!(matches!(
        read_plane_restraints(file.path(), &structure(), &AnalysisPolicy::default()),
        Err(ValidationInputError::InvalidValue { .. })
    ));
}

#[test]
fn tls_groups_carry_the_declared_model() {
    let file = write(&format!(
        "{}{}",
        group("A", "chain A"),
        group("B", "chain B")
    ));
    let groups = match read_tls_groups(file.path(), &structure(), &AnalysisPolicy::default()) {
        Ok(groups) => groups,
        Err(error) => panic!("groups failed: {error}"),
    };
    assert_eq!(groups.len(), 2);
    assert_eq!(groups[0].id, "A");
    assert_eq!(groups[0].atoms.len(), 3);
    for (read, expected) in groups[0].model.origin.iter().zip([1.0, 2.0, 3.0]) {
        assert!((read - expected).abs() < f64::EPSILON);
    }
    assert!((groups[0].model.translation[1][1] - 0.2).abs() < f64::EPSILON);
}

#[test]
fn a_tensor_that_is_not_three_by_three_is_a_schema_error() {
    let file = write(&group("A", "chain A").replace(
        "t = [[0.1, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.3]]",
        "t = [[0.1, 0.0], [0.0, 0.2]]",
    ));
    assert!(matches!(
        read_tls_groups(file.path(), &structure(), &AnalysisPolicy::default()),
        Err(ValidationInputError::Toml(_))
    ));
}

#[test]
fn a_nonfinite_tensor_entry_is_refused() {
    let file = write(&group("A", "chain A").replace("0.1,", "nan,"));
    assert!(matches!(
        read_tls_groups(file.path(), &structure(), &AnalysisPolicy::default()),
        Err(ValidationInputError::InvalidValue { .. })
    ));
}

#[test]
fn read_planes_drive_the_planarity_kernel_end_to_end() {
    use crate::ValidationExt;
    use molframe_validate::PlanarityOptions;
    const FLAT_AND_BENT: &str = "\
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  C2  LIG A   1       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  C3  LIG A   1       0.000   1.000   0.000  1.00  0.00           C
ATOM      4  C4  LIG A   1       1.000   1.000   0.000  1.00  0.00           C
ATOM      5  C5  LIG B   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      6  C6  LIG B   1       1.000   0.000   0.000  1.00  0.00           C
ATOM      7  C7  LIG B   1       0.000   1.000   0.000  1.00  0.00           C
ATOM      8  C8  LIG B   1       1.000   1.000   2.000  1.00  0.00           C
END
";
    let structure = match read_bytes(FLAT_AND_BENT.into(), Some("t.pdb"), &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(error) => panic!("fixture failed: {error}"),
    };
    let file = write(
        "[[plane]]\nid = 'flat'\nselection = 'chain A'\n[[plane]]\nid = 'bent'\nselection = 'chain B'\n",
    );
    let restraints =
        match read_plane_restraints(file.path(), &structure, &AnalysisPolicy::default()) {
            Ok(restraints) => restraints,
            Err(error) => panic!("planes failed: {error}"),
        };
    let report = match structure.plane_restraint_outliers(
        &restraints,
        PlanarityOptions {
            maximum_deviation: 0.1,
            plane_fit: molframe_geom::EigenOptions::standard(),
        },
    ) {
        Ok(report) => report,
        Err(error) => panic!("planarity failed: {error}"),
    };
    assert_eq!(report.assessed, 2);
    assert_eq!(report.flags.len(), 1);
    assert_eq!(report.flags[0].id, "bent");
}
