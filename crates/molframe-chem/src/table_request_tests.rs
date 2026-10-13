use super::*;
use crate::{TableIdentity, TableTextError};

#[derive(Debug, PartialEq)]
struct Words(Vec<String>);

impl TableKind for Words {
    const KIND: &'static str = "words";
    const EXTENSION: &'static str = "txt";

    fn decode(bytes: &[u8], fallback: TableIdentity) -> Result<Provided<Self>, TableTextError> {
        let text = crate::text_of("words table", bytes)?;
        Ok(Provided {
            identity: fallback,
            table: Self(text.split_whitespace().map(str::to_owned).collect()),
        })
    }
}

fn empty_store() -> TableStore {
    let dir = std::env::temp_dir().join(format!("molframe-request-{}", std::process::id()));
    TableStore::new(dir.join("nothing-here"))
}

#[test]
fn a_given_table_resolves_without_touching_the_store() {
    let request = TableRequest::given(Provided {
        identity: TableIdentity::supplied("w", "1", "test", "CC0-1.0"),
        table: Words(vec!["a".into()]),
    });
    let Ok(table) = request.resolve(&empty_store(), "--words") else {
        panic!("a given table resolves");
    };
    assert_eq!(table.table, Words(vec!["a".into()]));
}

#[test]
fn the_default_request_defers_to_the_profile_and_fails_typed_without_one() {
    let request = TableRequest::<Words>::default();
    assert!(matches!(
        request.resolve(&empty_store(), "--words"),
        Err(TableError::Missing { .. })
    ));
}

#[test]
fn a_named_request_that_is_not_in_the_cache_is_missing() {
    let request = TableRequest::<Words>::named("ghost");
    assert!(matches!(
        request.resolve(&empty_store(), "--words"),
        Err(TableError::Missing { .. })
    ));
}
