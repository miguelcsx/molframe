use super::{PsfError, parse_psf, write_psf};

const SOURCE: &str = "PSF EXT XPLOR\n\n       1 !NTITLE\n REMARKS test topology\n\n       4 !NATOM\n\
       1 SEG 1 GLY N  NH1 -0.30 14.007 0\n\
       2 SEG 1 GLY CA CT1  0.10 12.011 0\n\
       3 SEG 1 GLY C  C    0.50 12.011 0\n\
       4 SEG 1 GLY O  O   -0.30 15.999 0\n\n\
       3 !NBOND: bonds\n       1 2  2 3  3 4\n\n\
       2 !NTHETA: angles\n       1 2 3  2 3 4\n\n\
       1 !NPHI: dihedrals\n       1 2 3 4\n";

#[test]
fn ext_xplor_atoms_and_bonded_sections_round_trip() {
    let expected = parse_psf(SOURCE).unwrap_or_else(|error| panic!("PSF failed: {error}"));
    assert_eq!(expected.atoms[0].atom_type.as_ref(), "NH1");
    assert_eq!(expected.bonds, [[0, 1], [1, 2], [2, 3]]);
    assert_eq!(expected.angles, [[0, 1, 2], [1, 2, 3]]);
    let observed = parse_psf(&write_psf(&expected))
        .unwrap_or_else(|error| panic!("written PSF failed: {error}"));
    assert_eq!(observed, expected);
}

#[test]
fn connectivity_outside_natom_is_rejected() {
    let source = SOURCE.replace("3 4\n", "3 5\n");
    assert_eq!(parse_psf(&source), Err(PsfError::IndexOutOfRange));
}
