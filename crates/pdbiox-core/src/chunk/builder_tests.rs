use super::*;
use crate::element::Element;
use crate::index::ResidueIndex;
use crate::optional::OptionalSymbol;
use crate::symbol::AltId;

fn atom(residue: u32, element: Element, serial: u32) -> AtomRecord {
    AtomRecord {
        position: Some([serial as f32, 0.0, 0.0]),
        element,
        atom_name: SymbolId::from_raw(0),
        auth_atom_name: OptionalSymbol::NONE,
        alt_id: AltId::BLANK,
        residue: ResidueIndex::new(residue),
        occupancy: (1.0, Presence::Present),
        b_factor: (10.0 + serial as f32, Presence::Present),
        formal_charge: (0, Presence::Inapplicable),
        atom_site_id: serial,
    }
}

#[test]
fn a_chunk_never_ends_in_the_middle_of_a_residue() {
    let mut builder = ChunkBuilder::with_target(10);
    // Residues of seven atoms each: the target falls inside residue 1.
    for residue in 0..6u32 {
        for offset in 0..7u32 {
            builder.push(atom(residue, Element::CARBON, residue * 7 + offset));
        }
    }
    let (chunks, _) = builder.finish();
    assert!(chunks.len() > 1, "the target should have closed a chunk");
    for chunk in &chunks {
        let first = chunk.atoms().start;
        let last = chunk.atoms().end - 1;
        assert_eq!(first % 7, 0, "chunk starts at a residue boundary");
        assert_eq!(last % 7, 6, "chunk ends at a residue boundary");
    }
}

#[test]
fn chunks_cover_every_atom_once_and_in_order() {
    let mut builder = ChunkBuilder::with_target(8);
    for serial in 0..100u32 {
        builder.push(atom(serial / 3, Element::CARBON, serial));
    }
    let (chunks, coords) = builder.finish();
    assert_eq!(coords.len(), 100);

    let mut expected = 0u32;
    for chunk in &chunks {
        assert_eq!(chunk.atoms().start, expected);
        expected = chunk.atoms().end;
    }
    assert_eq!(expected, 100);
}

#[test]
fn a_chunks_summary_reports_the_elements_it_actually_contains() {
    let mut builder = ChunkBuilder::new();
    builder.push(atom(0, Element::CARBON, 0));
    builder.push(atom(0, Element::ZINC, 1));
    let (chunks, _) = builder.finish();

    let stats = chunks.first().map(AtomChunk::stats);
    assert_eq!(
        stats.map(|stats| stats.elements.contains(Element::ZINC)),
        Some(true)
    );
    assert_eq!(
        stats.map(|stats| stats.elements.contains(Element::IRON)),
        Some(false)
    );
    assert_eq!(stats.map(|stats| stats.has_hydrogen), Some(false));
}

#[test]
fn an_atom_with_no_recorded_position_keeps_its_row_and_is_marked_absent() {
    let mut builder = ChunkBuilder::new();
    let mut without = atom(0, Element::OXYGEN, 0);
    without.position = None;
    builder.push(without);
    builder.push(atom(0, Element::OXYGEN, 1));

    let (chunks, coords) = builder.finish();
    assert_eq!(coords.len(), 2, "the row exists even without a position");
    let Some(chunk) = chunks.first() else {
        panic!("expected one chunk")
    };
    assert!(!chunk.has_position(0));
    assert!(chunk.has_position(1));
    assert!(chunk.stats().has_missing_coords);
}

#[test]
#[allow(clippy::float_cmp, reason = "reads back stored values")]
fn a_missing_position_does_not_poison_the_bounding_box() {
    let mut builder = ChunkBuilder::new();
    let mut without = atom(0, Element::OXYGEN, 0);
    without.position = None;
    builder.push(without);
    builder.push(atom(0, Element::OXYGEN, 5));

    let (chunks, _) = builder.finish();
    let Some(bounds) = chunks.first().map(|chunk| chunk.stats().bounds) else {
        panic!("expected one chunk")
    };
    assert_eq!(bounds.min, [5.0, 0.0, 0.0]);
    assert_eq!(bounds.max, [5.0, 0.0, 0.0]);
}

#[test]
fn a_new_model_starts_a_new_chunk_so_the_model_column_can_vanish() {
    let mut builder = ChunkBuilder::new();
    builder.start_model(0);
    builder.push(atom(0, Element::CARBON, 0));
    builder.start_model(1);
    builder.push(atom(1, Element::CARBON, 1));

    let (chunks, _) = builder.finish();
    assert_eq!(chunks.len(), 2);
    assert_eq!(chunks.first().map(AtomChunk::model), Some(0));
    assert_eq!(chunks.get(1).map(AtomChunk::model), Some(1));
}

#[test]
fn a_chunk_reads_back_every_value_it_was_given() {
    let mut builder = ChunkBuilder::new();
    for serial in 0..20u32 {
        builder.push(atom(
            serial / 4,
            Element::from_atomic_number(6 + serial as u8),
            serial,
        ));
    }
    let (chunks, coords) = builder.finish();
    let Some(chunk) = chunks.first() else {
        panic!("expected one chunk")
    };

    for local in 0..20u32 {
        assert_eq!(
            chunk.element(local),
            Some(Element::from_atomic_number(6 + local as u8))
        );
        assert_eq!(chunk.atom_site_id(local), Some(local));
        assert_eq!(
            chunk.b_factor(local).map(|value| value.0),
            Some(10.0 + local as f32)
        );
        assert_eq!(chunk.alt_id(local), Some(AltId::BLANK));
    }
    assert_eq!(chunk.positions(&coords).map(<[[f32; 3]]>::len), Some(20));
}

#[test]
fn the_depositor_atom_name_column_is_absent_when_no_file_carried_one() {
    let mut builder = ChunkBuilder::new();
    builder.push(atom(0, Element::CARBON, 0));
    let (chunks, _) = builder.finish();
    assert_eq!(
        chunks.first().and_then(|chunk| chunk.auth_atom_name(0)),
        None
    );
}
