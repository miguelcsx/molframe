use super::*;

const HEADER: &str = "name: lab\nversion: 1\nsource: survey\nlicence: CC0-1.0\n";

fn table(body: &str) -> Result<Provided<BondReach>, TableTextError> {
    BondReach::decode(
        format!("{HEADER}{body}").as_bytes(),
        TableIdentity::supplied("x", "x", "x", "x"),
    )
}

#[test]
fn a_table_gives_reach_and_the_metal_flag() {
    let Ok(provided) = table("C 1.75\nZn 2.85 metal\n") else {
        panic!("a well-formed table parses");
    };
    assert_eq!(provided.table.reach(Element::CARBON), Some(1.75));
    assert!(!provided.table.is_metal(Element::CARBON));
    assert!(provided.table.is_metal(Element::from_atomic_number(30)));
    assert_eq!(provided.table.reach(Element::OXYGEN), None);
}

#[test]
fn every_defect_is_refused_with_its_line() {
    for body in [
        "C\n",
        "Xx 1.0\n",
        "C big\n",
        "C 0\n",
        "C 99\n",
        "C 1.7 bent\n",
        "C 1.7\nC 1.8\n",
    ] {
        let Err(error) = table(body) else {
            panic!("{body:?} must be refused");
        };
        assert!(error.to_string().contains("line"), "{body:?}: {error}");
    }
    assert!(table("").is_err());
}
