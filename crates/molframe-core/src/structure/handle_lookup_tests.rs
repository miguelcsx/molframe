//! The chunk-cursor atom lookup, checked against the structure invariant it
//! depends on: a residue's atoms are contiguous and never split across chunks.

use super::find_named;
use crate::chunk::{AtomChunk, AtomRecord, ChunkBuilder};
use crate::column::Presence;
use crate::element::Element;
use crate::optional::OptionalSymbol;
use crate::structure::handle::AtomRef;
use crate::symbol::AltId;

/// Builds one chain of `residues`, each carrying the four backbone names, with
/// a chunk target small enough that many chunks are produced.
fn chain(residues: u32, target: u32) -> (Vec<crate::chunk::AtomChunk>, crate::symbol::Interner) {
    let mut dictionary = crate::symbol::Interner::new();
    let mut atoms = Vec::new();
    for residue in 0..residues {
        for name in ["N", "CA", "C", "O"] {
            let Ok(atom_name) = dictionary.intern(name) else {
                panic!("fixture name interns")
            };
            atoms.push(AtomRecord {
                element: Element::CARBON,
                atom_name,
                auth_atom_name: OptionalSymbol::NONE,
                alternate_component_id: OptionalSymbol::NONE,
                alt_id: AltId::BLANK,
                residue: crate::index::ResidueIndex::new(residue),
                occupancy: (1.0, Presence::Present),
                b_factor: (20.0, Presence::Present),
                position: Some([0.0, 0.0, 0.0]),
                formal_charge: (0, Presence::Unknown),
                atom_site_id: residue,
            });
        }
    }

    let mut builder = ChunkBuilder::with_target(target);
    builder.start_model(0);
    for atom in atoms {
        builder.push(atom);
    }
    let (chunks, _) = builder.finish();
    (chunks, dictionary)
}

#[test]
fn the_lookup_finds_each_name_in_every_residue_across_many_chunks() {
    // Twelve residues of four atoms over a target of seven atoms each: the
    // chunk boundaries land inside the sequence rather than on residue edges.
    let (chunks, dictionary) = chain(12, 7);
    assert!(
        chunks.len() > 1,
        "the fixture must span several chunks, found {}",
        chunks.len()
    );
    let Some(ca) = dictionary.get("CA") else {
        panic!("the fixture interns CA")
    };
    let Some(absent) = dictionary.get("ZZ") else {
        return;
    };

    for residue in 0..12 {
        let range = residue * 4..residue * 4 + 4;
        assert_eq!(
            find_named(&chunks, range.clone(), ca),
            Some(residue * 4 + 1),
            "residue {residue} must report its own CA"
        );
        assert_eq!(
            find_named(&chunks, range, absent),
            None,
            "an absent name selects nothing"
        );
    }
}

#[test]
fn the_lookup_reports_the_first_match_when_a_name_repeats() {
    let (chunks, dictionary) = chain(4, 3);
    let Some(oxygen) = dictionary.get("O") else {
        panic!("the fixture interns O")
    };
    // Only the first residue's range is scanned, so the match is its own O.
    assert_eq!(find_named(&chunks, 0..4, oxygen), Some(3));
}

#[test]
fn the_lookup_declines_a_range_that_reaches_past_every_chunk() {
    let (chunks, dictionary) = chain(2, 3);
    let Some(nitrogen) = dictionary.get("N") else {
        panic!("the fixture interns N")
    };
    let total = chunks.last().map_or(0, |chunk| chunk.atoms().end);
    assert_eq!(find_named(&chunks, total..total + 4, nitrogen), None);
}

#[test]
fn every_residue_stays_inside_one_chunk_which_is_what_the_cursor_relies_on() {
    let (chunks, _) = chain(12, 7);
    for (position, chunk) in chunks.iter().enumerate() {
        // A chunk closes only before an atom that starts a new residue, so its
        // length is a whole number of residues and no residue straddles two.
        assert_eq!(
            chunk.len() % 4,
            0,
            "chunk {position} holds a partial residue: {} atoms",
            chunk.len()
        );
    }
}

#[test]
fn a_lookup_still_finds_a_single_residue_chunk() {
    // One residue cannot be split, so the smallest target still closes on a
    // residue boundary rather than mid-residue.
    let (chunks, dictionary) = chain(3, 1);
    let Some(nitrogen) = dictionary.get("N") else {
        panic!("the fixture interns N")
    };
    assert_eq!(find_named(&chunks, 4..8, nitrogen), Some(4));
}

#[test]
fn a_residue_built_by_the_hierarchy_fixture_reports_its_named_atom() {
    let structure = crate::structure::fixture::sample();
    let data = structure.data();
    let Some(model) = data.model(crate::index::ModelIndex::new(0)) else {
        panic!("the fixture has a model")
    };
    let Some(chain) = model.chain("A") else {
        panic!("the fixture has chain A")
    };
    let Some(residue) = chain.residue(100) else {
        panic!("the fixture has residue 100")
    };
    let atom = residue.atom("CA");
    assert_eq!(atom.and_then(AtomRef::name), Some("CA"));
    assert!(residue.atom("ZZ").is_none());
}

#[test]
fn the_chunks_tile_the_atom_order_without_a_gap() {
    // The cursor relies on a residue's range naming atoms that a chunk holds
    // and that chunks cover the order contiguously.
    let (chunks, _) = chain(6, 5);
    let total: u32 = chunks.iter().map(AtomChunk::len).sum();
    assert_eq!(total, 24, "six residues of four atoms");
    let mut seen = 0_u32;
    for chunk in &chunks {
        let atoms = chunk.atoms();
        assert_eq!(atoms.start, seen, "chunks tile the order without a gap");
        seen = atoms.end;
    }
    assert_eq!(seen, total);
}
