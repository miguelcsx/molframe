use super::one_letter;

#[test]
fn the_twenty_amino_acids_have_their_letters() {
    let pairs = [
        ("ALA", 'A'),
        ("ARG", 'R'),
        ("ASN", 'N'),
        ("ASP", 'D'),
        ("CYS", 'C'),
        ("GLN", 'Q'),
        ("GLU", 'E'),
        ("GLY", 'G'),
        ("HIS", 'H'),
        ("ILE", 'I'),
        ("LEU", 'L'),
        ("LYS", 'K'),
        ("MET", 'M'),
        ("PHE", 'F'),
        ("PRO", 'P'),
        ("SER", 'S'),
        ("THR", 'T'),
        ("TRP", 'W'),
        ("TYR", 'Y'),
        ("VAL", 'V'),
    ];
    for (name, letter) in pairs {
        assert_eq!(one_letter(name), letter, "{name}");
    }
}

#[test]
fn modified_and_unknown_residues_are_not_guessed() {
    assert_eq!(one_letter("MSE"), 'M');
    assert_eq!(one_letter("XYZ"), 'X');
    assert_eq!(one_letter("hoh"), 'X');
    assert_eq!(one_letter("dg"), 'G');
}
