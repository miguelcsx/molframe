use super::*;

#[test]
fn builder_composition_produces_the_same_ir_shape_as_text() {
    let built =
        col::is_protein() & col::bfactor().gt(30.0) & col::within(5.0, col::resname().eq("ATP"));
    let parsed = match crate::parser::parse("protein and bfactor > 30 and within 5 of resname ATP")
    {
        Ok(parsed) => parsed.expr,
        Err(findings) => panic!("parse failed: {findings:?}"),
    };
    assert_eq!(built.clone().into_parts().0, parsed);
    assert!(crate::parser::parse(built.source()).is_ok());
}

#[test]
fn molecular_helpers_share_the_textual_query_ir() {
    let built = col::residues_within(5.0, col::ligands())
        & col::heavy()
        & !col::water()
        & (col::nucleic() | col::glycans());
    let parsed = match crate::parser::parse(
        "byres within 5 of ligand and heavy and not water and (nucleic or saccharide)",
    ) {
        Ok(parsed) => parsed.expr,
        Err(findings) => panic!("parse failed: {findings:?}"),
    };
    assert_eq!(built.clone().into_parts().0, parsed);
    assert!(crate::parser::parse(built.source()).is_ok());
}
