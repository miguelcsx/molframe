//! Decompression bombs, forged sizes and damaged streams: every one ends in a
//! registered diagnostic, with output never larger than the configured cap.

use super::super::InputBuffer;
use super::*;
use std::io::Write;

fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    encoder.write_all(bytes).unwrap();
    encoder.finish().unwrap()
}

fn zeros_gzip(length: usize) -> Vec<u8> {
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    let block = vec![0_u8; 1 << 20];
    for _ in 0..length / block.len() {
        encoder.write_all(&block).unwrap();
    }
    encoder.finish().unwrap()
}

fn refused(result: Result<InputBuffer, Diagnostic>) -> Code {
    result.expect_err("input must be refused").code()
}

struct SplitMix(u64);

impl SplitMix {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn below(&mut self, bound: usize) -> usize {
        usize::try_from(self.next() % u64::try_from(bound).unwrap()).unwrap()
    }
}

#[test]
fn a_real_bomb_is_refused_by_the_ratio_ceiling() {
    let bomb = zeros_gzip(32 << 20);
    assert!(bomb.len() < 256 << 10, "the fixture must be a real bomb");
    let limits = Limits {
        compression_ratio: 50,
        ..Limits::default()
    };
    assert_eq!(
        refused(InputBuffer::from_encoded_bytes(bomb.clone(), limits)),
        Code::E1901
    );
    assert_eq!(
        refused(InputBuffer::from_reader(bomb.as_slice(), limits)),
        Code::E1901
    );
}

#[test]
fn a_real_bomb_is_refused_by_the_absolute_ceiling_alone() {
    let bomb = zeros_gzip(32 << 20);
    let limits = Limits {
        decompressed_bytes: 1 << 20,
        compression_ratio: u64::MAX,
        ..Limits::default()
    };
    assert_eq!(
        refused(InputBuffer::from_encoded_bytes(bomb.clone(), limits)),
        Code::E1901
    );
    assert_eq!(
        refused(InputBuffer::from_reader(bomb.as_slice(), limits)),
        Code::E1901
    );
}

#[test]
fn collection_never_hands_back_more_than_the_cap() {
    let payload = vec![7_u8; 100_000];
    let compressed = gzip(&payload);
    for cap in [0_u64, 1, 4095, 65_535, 65_536, 99_999] {
        let limits = Limits {
            decompressed_bytes: cap,
            compression_ratio: u64::MAX,
            ..Limits::default()
        };
        let result = expand_gzip_reader(compressed.as_slice(), 1, limits, None);
        assert_eq!(
            result.err().map(|d| d.code()),
            Some(Code::E1901),
            "cap {cap}"
        );
    }
    let limits = Limits {
        decompressed_bytes: 100_000,
        compression_ratio: u64::MAX,
        ..Limits::default()
    };
    let exact = expand_gzip_reader(compressed.as_slice(), 1, limits, None).unwrap();
    assert_eq!(exact, payload);
}

#[test]
fn a_forged_size_trailer_cannot_force_an_allocation() {
    let mut forged = gzip(b"tiny");
    let end = forged.len();
    forged[end - 4..].copy_from_slice(&u32::MAX.to_le_bytes());
    let code = refused(InputBuffer::from_encoded_bytes(forged, Limits::default()));
    assert!(matches!(code, Code::E7101 | Code::E1901));
}

#[test]
fn every_truncation_of_a_valid_stream_is_refused() {
    let compressed = gzip(&b"data_example\n".repeat(500));
    for cut in 0..compressed.len() {
        let result = InputBuffer::from_encoded_bytes(compressed[..cut].to_vec(), Limits::default());
        if cut >= 4 {
            let code = refused(result);
            assert!(matches!(code, Code::E7101 | Code::E1901), "cut {cut}");
        }
    }
}

#[test]
fn mutated_and_random_gzip_input_never_panics_or_exceeds_the_cap() {
    let cap = 1 << 16;
    let limits = Limits {
        decompressed_bytes: cap,
        ..Limits::default()
    };
    let seed = gzip(&b"data_example\n".repeat(200));
    let mut rng = SplitMix(0x5EED);
    for round in 0..3000 {
        let mut bytes = if round % 2 == 0 {
            seed.clone()
        } else {
            let mut raw = vec![0x1f, 0x8b, 0x08];
            raw.extend((0..rng.below(300)).map(|_| u8::try_from(rng.next() & 0xff).unwrap()));
            raw
        };
        for _ in 0..=rng.below(4) {
            let at = rng.below(bytes.len());
            bytes[at] ^= 1 << rng.below(8);
        }
        match InputBuffer::from_encoded_bytes(bytes, limits) {
            Ok(input) => assert!(input.as_bytes().len() as u64 <= cap),
            Err(diagnostic) => {
                assert!(matches!(diagnostic.code(), Code::E7101 | Code::E1901));
            }
        }
    }
}
