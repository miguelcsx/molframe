use super::*;

#[test]
fn symbols_parse_case_insensitively_because_some_formats_upper_case_them() {
    assert_eq!(Element::from_symbol("FE"), Some(Element::IRON));
    assert_eq!(Element::from_symbol("Fe"), Some(Element::IRON));
    assert_eq!(Element::from_symbol(" C "), Some(Element::CARBON));
    assert_eq!(Element::from_symbol(""), None);
    assert_eq!(Element::from_symbol("Xx"), None);
    assert_eq!(Element::from_symbol("Carbon"), None);
}

#[test]
fn hydrogen_isotopes_parse_as_hydrogen_rather_than_as_unknown() {
    assert_eq!(Element::from_symbol("D"), Some(Element::HYDROGEN));
    assert_eq!(Element::from_symbol("T"), Some(Element::HYDROGEN));
}

#[test]
fn a_leading_space_distinguishes_c_alpha_from_calcium() {
    assert_eq!(Element::infer_from_pdb_atom_name(" CA "), Element::CARBON);
    assert_eq!(Element::infer_from_pdb_atom_name("CA  "), Element::CALCIUM);
    assert_eq!(Element::infer_from_pdb_atom_name(" N  "), Element::NITROGEN);
    assert_eq!(Element::infer_from_pdb_atom_name("ZN  "), Element::ZINC);
    assert_eq!(Element::infer_from_pdb_atom_name(" SE "), Element::SULFUR);
}

#[test]
fn names_beginning_with_a_digit_fall_through_to_their_first_letter() {
    assert_eq!(Element::infer_from_pdb_atom_name("1HB "), Element::HYDROGEN);
    assert_eq!(Element::infer_from_pdb_atom_name("2HD1"), Element::HYDROGEN);
}

#[test]
fn a_name_without_column_context_reads_two_letters_first() {
    assert_eq!(Element::infer_from_name("CA"), Element::CALCIUM);
    assert_eq!(Element::infer_from_name("C1"), Element::CARBON);
    assert_eq!(Element::infer_from_name("SE"), Element::SELENIUM);
    assert_eq!(Element::infer_from_name("1H"), Element::HYDROGEN);
    assert_eq!(Element::infer_from_name(""), Element::UNKNOWN);
}

#[test]
fn atomic_numbers_round_trip_and_out_of_range_values_are_unknown() {
    for z in 0..=Element::MAX_ATOMIC_NUMBER {
        assert_eq!(Element::from_atomic_number(z).atomic_number(), z);
    }
    assert!(Element::from_atomic_number(200).is_unknown());
    assert!(Element::UNKNOWN.symbol().is_empty());
}

#[test]
fn keys_are_derived_from_the_symbol_table_at_compile_time() {
    assert_eq!(KEYS[0], [0, 0]);
    assert_eq!(KEYS[6], [b'C', 0]);
    assert_eq!(KEYS[26], [b'F', b'e']);
}
