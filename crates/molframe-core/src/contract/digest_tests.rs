use super::ContentDigest;

#[test]
fn the_digest_is_the_published_sha256_of_the_bytes() {
    // FIPS 180-4 example vectors.
    assert_eq!(
        ContentDigest::of(b"abc").to_string(),
        "sha256:ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    );
    assert_eq!(
        ContentDigest::of(b"").to_string(),
        "sha256:e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    );
    assert_eq!(
        ContentDigest::of(b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq").to_string(),
        "sha256:248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"
    );
}

#[test]
fn a_digest_reads_back_from_its_text_and_refuses_anything_else() {
    let digest = ContentDigest::of(b"1crn");
    assert_eq!(ContentDigest::from_hex(&digest.to_string()), Some(digest));
    let bare = digest.to_string().replace("sha256:", "");
    assert_eq!(ContentDigest::from_hex(&bare), Some(digest));
    assert_eq!(ContentDigest::from_hex("sha256:abc"), None);
    assert_eq!(ContentDigest::from_hex(&"z".repeat(64)), None);
}

#[test]
fn one_changed_byte_changes_the_digest() {
    assert_ne!(ContentDigest::of(b"ATOM  1"), ContentDigest::of(b"ATOM  2"));
}
