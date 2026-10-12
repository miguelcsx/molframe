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

#[test]
fn a_filter_that_keeps_connections_keeps_the_operators_that_explain_them() {
    // An operator table dropped by the filter would leave a nonidentity symmetry
    // code with nothing to check against, and the reader would have to assume the
    // local identity. Keeping `struct_conn` therefore keeps its dependencies.
    let connections = CategoryFilter::only(["struct_conn"]);
    for name in SYMMETRY_DEPENDENCIES {
        assert!(!connections.admits(name));
        assert!(connections.keeps_category(name), "{name}");
    }
    // Without the connections, the tables are the caller's business again.
    let unrelated = CategoryFilter::only(["entry"]);
    assert!(!unrelated.keeps_category("struct_conn"));
    for name in SYMMETRY_DEPENDENCIES {
        assert!(!unrelated.keeps_category(name), "{name}");
    }
}
