use super::*;

#[test]
fn checksum_vocabulary_is_strict_and_case_normalised() {
    let uppercase = "A".repeat(64);
    assert!(matches!(normalise_checksum(&uppercase), Ok(value) if value == "a".repeat(64)));
    assert!(matches!(
        normalise_checksum("abc"),
        Err(DownloadError::InvalidChecksum)
    ));
    assert!(matches!(
        normalise_checksum(&"x".repeat(64)),
        Err(DownloadError::InvalidChecksum)
    ));
}

#[test]
fn digest_is_the_standard_sha256_hex_representation() {
    assert_eq!(
        digest(b"abc"),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    );
}

#[test]
fn request_limits_are_explicit_and_validated_before_network_access() {
    let result = fetch_verified(
        "",
        &"0".repeat(64),
        DownloadOptions {
            max_bytes: 1,
            timeout: std::time::Duration::from_secs(1),
            redirect_limit: 0,
        },
    );
    assert!(matches!(result, Err(DownloadError::InvalidRequest)));
}
