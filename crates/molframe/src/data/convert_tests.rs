use super::*;

#[test]
fn ids_are_resolved_and_unknown_ones_are_refused() {
    assert_eq!(Converter::named("tsv"), Ok(Converter::Tsv));
    assert_eq!(
        Converter::named("stereo_chemical_props"),
        Ok(Converter::StereoChemicalProps)
    );
    assert!(Converter::named("csv").is_err());
}

#[test]
fn tsv_keeps_rows_and_drops_blank_lines_and_comments() {
    let rows = Converter::Tsv
        .rows("# note\nC\t1.7\n\n  N\t1.55  \n")
        .expect("rows");
    assert_eq!(rows, ["C\t1.7", "N\t1.55"]);
}
