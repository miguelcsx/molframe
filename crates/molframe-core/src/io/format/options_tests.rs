use super::*;

#[test]
fn names_compare_ignoring_ascii_case() {
    let filter = CategoryFilter::except(["Struct_Conf"]);
    assert!(!filter.admits("struct_conf"));
    assert!(filter.admits("struct_conn"));
    assert!(CategoryFilter::only(["CONECT"]).admits("conect"));
}

#[test]
fn required_names_survive_every_filter() {
    let none = CategoryFilter::only(Vec::<Box<str>>::new());
    for name in REQUIRED_CIF_CATEGORIES {
        assert!(!none.admits(name));
        assert!(none.keeps_category(name));
        assert!(CategoryFilter::except([name]).keeps_category(name));
    }
    for name in REQUIRED_PDB_RECORDS {
        assert!(none.keeps_record(&name.to_ascii_lowercase()));
    }
    assert!(!none.keeps_category("struct_conf"));
    assert!(CategoryFilter::All.keeps_record("CONECT"));
}
