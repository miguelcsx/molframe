use super::*;

const HEADER: &str = "name: lab\nversion: 1\nsource: survey\nlicence: CC0-1.0\n";

fn table(body: &str) -> Result<Provided<IonicRadii>, TableTextError> {
    IonicRadii::decode(
        format!("{HEADER}{body}").as_bytes(),
        TableIdentity::supplied("x", "x", "x", "x"),
    )
}

#[test]
fn entries_keep_charge_coordination_spin_and_reliability() {
    let Ok(provided) =
        table("Na 1 VI - 1.02 1.16 yes\nFe 2 VI HS 0.78 0.92 yes\nFe 2 VI LS 0.61 0.75 -\n")
    else {
        panic!("a well-formed table parses");
    };
    let iron = provided.table.entries(Element::IRON);
    assert_eq!(iron.len(), 2);
    assert_eq!(iron[0].spin, IonicSpin::High);
    assert_eq!(iron[1].spin, IonicSpin::Low);
    assert_eq!(iron[1].most_reliable, None);
    assert_eq!(&*iron[0].coordination, "VI");
    assert!(provided.table.entries(Element::CARBON).is_empty());
}

#[test]
fn every_defect_is_refused_with_its_line() {
    for (body, line) in [
        ("Na 1 VI - 1.02 1.16\n", 5),
        ("Xx 1 VI - 1.02 1.16 yes\n", 5),
        ("Na one VI - 1.02 1.16 yes\n", 5),
        ("Na 1 VI up 1.02 1.16 yes\n", 5),
        ("Na 1 VI - big 1.16 yes\n", 5),
        ("Na 1 VI - 1.02 1.16 maybe\n", 5),
        ("Fe 2 VI - 1.0 1.0 yes\nNa 1 VI - 1.0 1.0 yes\n", 6),
    ] {
        let Err(error) = table(body) else {
            panic!("{body:?} must be refused");
        };
        assert!(
            error.to_string().contains(&format!("line {line}")),
            "{body:?}: {error}"
        );
    }
    assert!(table("").is_err());
}
