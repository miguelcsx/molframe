use super::*;
use crate::{Sample, structure};

/// Bit patterns of a position, so equality is exact rather than approximate.
fn bits(position: [f32; 3]) -> [u32; 3] {
    [
        position[0].to_bits(),
        position[1].to_bits(),
        position[2].to_bits(),
    ]
}

fn tiny_tile() -> Tile {
    Tile::from_structure(&structure(Sample::Tiny), Seed::new(1))
}

#[test]
fn a_tile_carries_every_atom_of_its_source() {
    let source = structure(Sample::Tiny);
    let tile = Tile::from_structure(&source, Seed::new(1));
    let Ok(atoms) = u32::try_from(tile.atoms().len()) else {
        panic!("the tiny sample does not fit a u32 atom count")
    };
    assert_eq!(atoms, source.atom_count());
}

#[test]
fn the_name_table_holds_far_fewer_entries_than_the_atom_count() {
    let tile = tiny_tile();
    assert!(
        tile.names().len() * 4 < tile.atoms().len(),
        "name table {} is not small against {} atoms",
        tile.names().len(),
        tile.atoms().len()
    );
}

#[test]
fn a_name_position_past_the_table_reads_as_empty() {
    let tile = tiny_tile();
    assert_eq!(tile.name(u32::MAX), "");
}

#[test]
fn copies_for_rounds_up_to_cover_the_requested_atoms() {
    let tile = tiny_tile();
    let Ok(per_copy) = u64::try_from(tile.atoms().len()) else {
        panic!("the tiny sample does not fit a u64 atom count")
    };
    assert_eq!(tile.copies_for(per_copy), 1);
    assert_eq!(tile.copies_for(per_copy + 1), 2);
    assert!(tile.copies_for(per_copy) < tile.copies_for(per_copy + 1));
    assert!(tile.atoms_in(tile.copies_for(per_copy + 1)) > per_copy);
}

#[test]
fn the_same_seed_places_a_copy_the_same_way() {
    let tile = tiny_tile();
    let Some(atom) = tile.atoms().first().copied() else {
        panic!("the tiny sample has no atoms")
    };

    let mut left_seed = tile.jitter_seed(3);
    let mut right_seed = tile.jitter_seed(3);
    let left = tile.placement(3, 8).apply(atom.position, &mut left_seed);
    let right = tile.placement(3, 8).apply(atom.position, &mut right_seed);
    assert_eq!(bits(left), bits(right));
}

#[test]
fn distinct_copies_are_placed_apart() {
    let tile = tiny_tile();
    let Some(atom) = tile.atoms().first().copied() else {
        panic!("the tiny sample has no atoms")
    };

    let mut first_seed = tile.jitter_seed(0);
    let mut second_seed = tile.jitter_seed(1);
    let first = tile.placement(0, 8).apply(atom.position, &mut first_seed);
    let second = tile.placement(1, 8).apply(atom.position, &mut second_seed);

    let separation = f64::from(first[0] - second[0])
        .hypot(f64::from(first[1] - second[1]))
        .hypot(f64::from(first[2] - second[2]));
    assert!(separation > 1.0, "copies coincide: {separation}");
}

#[test]
fn a_placed_position_stays_finite() {
    let tile = tiny_tile();
    let mut seed = tile.jitter_seed(11);
    let placement = tile.placement(11, 64);
    for atom in tile.atoms() {
        let placed = placement.apply(atom.position, &mut seed);
        assert!(
            placed.iter().all(|value| value.is_finite()),
            "placed position is not finite: {placed:?}"
        );
    }
}

#[test]
fn jitter_stays_within_the_configured_half_width() {
    let tile = tiny_tile().with_jitter(0.0);
    let Some(atom) = tile.atoms().first().copied() else {
        panic!("the tiny sample has no atoms")
    };
    let mut left_seed = tile.jitter_seed(0);
    let mut right_seed = tile.jitter_seed(0);
    let placement = tile.placement(0, 1);
    assert_eq!(
        bits(placement.apply(atom.position, &mut left_seed)),
        bits(placement.apply(atom.position, &mut right_seed))
    );
}
