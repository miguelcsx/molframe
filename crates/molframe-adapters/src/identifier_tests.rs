use super::{IdentifierFetchError, fetch_identifier};

#[test]
fn an_empty_identifier_is_rejected_before_provider_resolution() {
    assert_eq!(
        fetch_identifier("  "),
        Err(IdentifierFetchError::EmptyIdentifier)
    );
}

#[test]
fn an_identifier_reports_provider_unavailability_instead_of_guessing_a_url() {
    assert_eq!(
        fetch_identifier("1HHO"),
        Err(IdentifierFetchError::ProviderUnavailable {
            identifier: "1HHO".to_owned(),
        })
    );
}
