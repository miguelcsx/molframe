use super::*;

#[test]
fn contact_table_keeps_each_native_contact_in_its_column() {
    let contacts = vec![
        pdbiox::analysis::Contact {
            first: pdbiox::AtomIndex::new(2),
            second: pdbiox::AtomIndex::new(7),
            distance: 3.25,
        },
        pdbiox::analysis::Contact {
            first: pdbiox::AtomIndex::new(9),
            second: pdbiox::AtomIndex::new(11),
            distance: 4.5,
        },
    ];

    let table = PyContactTable::from(contacts);

    assert_eq!(table.first, vec![2, 9]);
    assert_eq!(table.second, vec![7, 11]);
    assert_eq!(table.distance, vec![3.25, 4.5]);
}
