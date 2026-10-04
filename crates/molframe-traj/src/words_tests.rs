use super::*;
use molframe_core::Diagnostic;

#[test]
fn each_choice_is_read_from_its_word() {
    assert_eq!("single".parse::<Linkage>(), Ok(Linkage::Single));
    assert_eq!("average".parse::<Linkage>(), Ok(Linkage::Average));
    assert_eq!(
        "reject".parse::<RemainderPolicy>(),
        Ok(RemainderPolicy::Reject)
    );
    assert_eq!(
        "fitted_rmsd".parse::<PathFrameMetric>(),
        Ok(PathFrameMetric::FittedRmsd)
    );
    assert_eq!(
        "cartesian-rmsd".parse::<PathFrameMetric>(),
        Ok(PathFrameMetric::CartesianRmsd)
    );
}

#[test]
fn an_unknown_word_names_the_choices() {
    let Err(error) = "ward".parse::<Linkage>() else {
        panic!("ward is not a linkage");
    };
    assert!(error.to_string().contains("single, complete and average"));
    assert_eq!(Diagnostic::from(error).code(), Code::E5101);
}
