use super::*;

#[test]
fn every_declared_blosum_and_pam_profile_loads() {
    for level in (30..=90).step_by(5) {
        assert!(
            load_matrix(MatrixProfile::Blosum(level)).is_ok(),
            "BLOSUM{level}"
        );
    }
    for distance in (10..=500).step_by(10) {
        assert!(
            load_matrix(MatrixProfile::Pam(distance)).is_ok(),
            "PAM{distance}"
        );
    }
}

#[test]
fn nucleotide_and_identity_scores_match_the_bundled_tables() {
    let Ok(nucleotide) = load_matrix(MatrixProfile::Nuc44) else {
        panic!("valid NUC.4.4");
    };
    assert_eq!(nucleotide.get(b'A', b'A'), 5);
    assert_eq!(nucleotide.get(b'A', b'T'), -4);
    let Ok(identity) = load_matrix(MatrixProfile::Identity) else {
        panic!("valid identity matrix");
    };
    assert_eq!(identity.get(b'W', b'W'), 1);
    assert_eq!(identity.get(b'W', b'A'), -10_000);
}

#[test]
fn profile_provenance_is_exact_and_inspectable() {
    let Ok(matrix) = load_matrix(MatrixProfile::Blosum(30)) else {
        panic!("valid BLOSUM30");
    };
    assert_eq!(matrix.identity().name(), "BLOSUM30");
    assert_eq!(matrix.identity().version(), "ncbi-blast-1997-08-26");
    assert_eq!(
        matrix.identity().source(),
        "https://ftp.ncbi.nih.gov/blast/matrices/BLOSUM30"
    );
    assert_eq!(matrix.identity().retrieved(), "2026-08-08");
}

#[test]
fn unsupported_family_parameters_are_errors() {
    assert!(matches!(
        load_matrix(MatrixProfile::Blosum(31)),
        Err(MatrixError::UnknownProfile)
    ));
    assert!(matches!(
        load_matrix(MatrixProfile::Pam(25)),
        Err(MatrixError::UnknownProfile)
    ));
}
