use super::ContactArrowTable;
use molframe_analysis::{Contact, ContactTable};
use molframe_core::AtomIndex;
use std::sync::Arc;

#[test]
fn contact_columns_remain_aligned_in_arrow() {
    let contacts = [
        Contact {
            first: AtomIndex::new(1),
            second: AtomIndex::new(2),
            distance: 3.0,
        },
        Contact {
            first: AtomIndex::new(4),
            second: AtomIndex::new(9),
            distance: 1.5,
        },
    ]
    .into_iter()
    .collect::<ContactTable>();
    let batches = ContactArrowTable::new(Arc::new(contacts))
        .record_batches()
        .expect("valid contact table");
    assert_eq!(batches.len(), 1);
    assert_eq!(batches[0].num_rows(), 2);
    assert_eq!(batches[0].num_columns(), 3);
}
