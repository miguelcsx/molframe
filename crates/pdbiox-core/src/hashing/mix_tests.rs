use super::*;
use std::hash::{BuildHasher, Hash};

fn hash_of<T: Hash>(value: &T) -> u64 {
    IdentityBuildHasher::default().hash_one(value)
}

#[test]
fn distinct_integers_hash_to_distinct_values() {
    let first = hash_of(&1u32);
    let second = hash_of(&2u32);

    assert_ne!(first, second);
    assert_ne!(first, hash_of(&0u32));
}

#[test]
fn a_pair_does_not_collide_with_its_own_reverse() {
    assert_ne!(hash_of(&(3u32, 7u32)), hash_of(&(7u32, 3u32)));
}

#[test]
fn hashing_is_a_pure_function_of_the_key() {
    assert_eq!(hash_of(&(11u32, 13u32)), hash_of(&(11u32, 13u32)));
    assert_eq!(hash_of(&99u64), hash_of(&99u64));
}

#[test]
fn the_multiplier_is_fixed_rather_than_process_seeded() {
    // Pinning the value proves the hash is reproducible across runs and
    // machines, which is what makes map iteration order deterministic.
    let mut hasher = IdentityHasher::default();
    hasher.write_u32(1);

    assert_eq!(hasher.finish(), {
        let mixed = 1u64.wrapping_mul(0x9E37_79B9_7F4A_7C15);
        mixed ^ (mixed >> 29)
    });
}

#[test]
fn low_bits_differ_for_keys_that_differ_only_in_low_bits() {
    // A hash table selects its bucket from the low bits, so consecutive atom
    // indices must not all land in one bucket.
    let buckets: std::collections::BTreeSet<u64> =
        (0..64u32).map(|index| hash_of(&index) & 0x3F).collect();

    assert!(
        buckets.len() > 40,
        "only {} distinct buckets",
        buckets.len()
    );
}

#[test]
fn byte_writes_absorb_every_chunk_and_the_remainder() {
    let mut full = IdentityHasher::default();
    full.write(&[1, 2, 3, 4, 5, 6, 7, 8, 9]);

    let mut short = IdentityHasher::default();
    short.write(&[1, 2, 3, 4, 5, 6, 7, 8]);

    assert_ne!(full.finish(), short.finish());
}

#[test]
fn the_maps_behave_as_ordinary_collections() {
    let mut map: IdentityHashMap<(u32, u32), f32> = IdentityHashMap::default();
    map.insert((1, 2), 0.5);
    map.insert((2, 1), 1.5);

    assert_eq!(map.get(&(1, 2)), Some(&0.5));
    assert_eq!(map.get(&(2, 1)), Some(&1.5));
    assert_eq!(map.len(), 2);

    let mut set: IdentityHashSet<u32> = IdentityHashSet::default();
    assert!(set.insert(7));
    assert!(!set.insert(7));
}
