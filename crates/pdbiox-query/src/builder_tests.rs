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
    assert_eq!(built.into_expr(), parsed);
}
