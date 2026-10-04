use super::*;

#[test]
fn every_format_name_round_trips_and_suffixes_are_accepted() {
    for name in Format::NAMES {
        assert_eq!(name.parse::<Format>().map(Format::name), Ok(*name));
    }
    assert_eq!("CIF".parse(), Ok(Format::Mmcif));
    assert_eq!("ent".parse(), Ok(Format::Pdb));
    assert_eq!("mol".parse(), Ok(Format::Sdf));
    let refused = "xyz".parse::<Format>().expect_err("not a format");
    assert_eq!(refused.field, "format");
}

#[test]
fn read_decisions_are_spelled_in_the_shared_vocabulary() {
    assert_eq!("recover".parse(), Ok(ParseMode::Recover));
    assert_eq!(
        "infer_from_atom_name".parse(),
        Ok(MissingElementPolicy::InferFromAtomName)
    );
    assert_eq!(
        AmbiguousResidueBoundaryPolicy::InferFromFileOrder.to_string(),
        "infer-from-file-order"
    );
    assert!("lenient".parse::<ParseMode>().is_err());
}
