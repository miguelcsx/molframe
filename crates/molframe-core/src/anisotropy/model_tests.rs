use super::*;

/// Asserts two six-component tensors agree componentwise.
fn assert_tensor(actual: [f32; 6], expected: [f32; 6]) {
    for (actual, expected) in actual.iter().zip(expected) {
        assert!((actual - expected).abs() <= f32::EPSILON);
    }
}

#[test]
fn tensor_components_keep_the_u11_u22_u33_u12_u13_u23_order() {
    let mut builder = AnisotropyTableBuilder::new();
    builder.push(AnisotropicDisplacement {
        atom: AtomIndex::new(0),
        u: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
    });
    let table = builder.finish();
    let Some(stored) = table.for_atom(AtomIndex::new(0)) else {
        panic!("the stored tensor was lost");
    };
    assert_tensor(stored, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    let Some(record) = table.get(AnisotropyIndex::new(0)) else {
        panic!("the only row was lost")
    };
    assert_tensor(record.u, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
}

#[test]
fn rows_are_sorted_by_atom_and_the_first_duplicate_wins() {
    let mut builder = AnisotropyTableBuilder::new();
    builder.push(AnisotropicDisplacement {
        atom: AtomIndex::new(3),
        u: [3.0; 6],
    });
    builder.push(AnisotropicDisplacement {
        atom: AtomIndex::new(1),
        u: [1.0; 6],
    });
    builder.push(AnisotropicDisplacement {
        atom: AtomIndex::new(3),
        u: [9.0; 6],
    });
    let table = builder.finish();
    let records: Vec<_> = table.iter().collect();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].atom, AtomIndex::new(1));
    assert_eq!(records[1].atom, AtomIndex::new(3));
    assert_tensor(records[1].u, [3.0; 6]);
    assert_eq!(table.len(), 2);
    assert!(!table.is_empty());
}

#[test]
fn for_atom_locates_rows_by_binary_search() {
    let mut builder = AnisotropyTableBuilder::new();
    for (atom, value) in [(0_u32, 0.0_f32), (7, 7.0), (9, 9.0)] {
        builder.push(AnisotropicDisplacement {
            atom: AtomIndex::new(atom),
            u: [value; 6],
        });
    }
    let table = builder.finish();
    let Some(empty) = table.for_atom(AtomIndex::new(0)) else {
        panic!("the first tensor was lost");
    };
    assert_tensor(empty, [0.0; 6]);
    let Some(seven) = table.for_atom(AtomIndex::new(7)) else {
        panic!("the second tensor was lost");
    };
    assert_tensor(seven, [7.0; 6]);
    assert!(table.for_atom(AtomIndex::new(9)).is_some());
    // A gap in the sparse coverage resolves to nothing, never to a neighbour.
    assert!(table.for_atom(AtomIndex::new(1)).is_none());
    assert!(table.for_atom(AtomIndex::new(10)).is_none());
}

#[test]
fn availability_survives_construction_when_asked() {
    let mut builder = AnisotropyTableBuilder::new();
    builder.push(AnisotropicDisplacement {
        atom: AtomIndex::new(0),
        u: [1.0; 6],
    });
    let table = builder.finish_with_availability(false);
    assert_eq!(table.len(), 1);
    assert!(!table.is_available());
    let mut builder = AnisotropyTableBuilder::new();
    builder.push(AnisotropicDisplacement {
        atom: AtomIndex::new(0),
        u: [1.0; 6],
    });
    let table = builder.finish();
    assert!(table.is_available());
    assert!(!AnisotropyTable::default().is_available());
    assert!(AnisotropyTable::default().is_empty());
}
