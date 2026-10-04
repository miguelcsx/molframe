use super::{ReflectionColumnType, ReflectionTable, UnknownColumnType};

#[test]
fn every_column_type_is_read_back_from_its_own_name() {
    for column_type in [
        ReflectionColumnType::MillerIndex,
        ReflectionColumnType::Amplitude,
        ReflectionColumnType::Intensity,
        ReflectionColumnType::StandardDeviation,
        ReflectionColumnType::Phase,
        ReflectionColumnType::Flag,
        ReflectionColumnType::Real,
        ReflectionColumnType::Text,
    ] {
        assert_eq!(column_type.name().parse(), Ok(column_type));
    }
}

#[test]
fn an_unknown_column_type_is_named_in_the_refusal() {
    assert_eq!(
        "weights".parse::<ReflectionColumnType>(),
        Err(UnknownColumnType("weights".into()))
    );
}

#[test]
fn a_table_that_names_its_space_group_gets_that_groups_operations() {
    let mut table = ReflectionTable {
        space_group_name: Some("P 21 21 21".into()),
        ..ReflectionTable::default()
    };
    assert!(table.resolve_symmetry_operations().is_ok());
    assert_eq!(table.symmetry_operations.len(), 4);
    assert!(
        table
            .symmetry_operations
            .iter()
            .any(|operation| &**operation == "x,y,z")
    );
    let mut by_number = ReflectionTable {
        space_group_number: Some(1),
        ..ReflectionTable::default()
    };
    assert!(by_number.resolve_symmetry_operations().is_ok());
    assert_eq!(by_number.symmetry_operations.len(), 1);
}

#[test]
fn a_table_that_lists_operations_or_names_no_group_is_left_alone() {
    let mut listed = ReflectionTable {
        space_group_name: Some("P 1".into()),
        symmetry_operations: vec!["x,y,z".into(), "-x,-y,-z".into()],
        ..ReflectionTable::default()
    };
    assert!(listed.resolve_symmetry_operations().is_ok());
    assert_eq!(listed.symmetry_operations.len(), 2);
    let mut none = ReflectionTable::default();
    assert!(none.resolve_symmetry_operations().is_ok());
    assert!(none.symmetry_operations.is_empty());
    let mut unknown = ReflectionTable {
        space_group_name: Some("not a group".into()),
        ..ReflectionTable::default()
    };
    assert!(unknown.resolve_symmetry_operations().is_err());
}
