use crate::AssemblyExt;
use crate::view::tests::{ENTRY, attached};
use molframe_core::annotation::AtomAnnotation;
use molframe_core::structure::ChainRef;
use molframe_core::{INSTANCE_ID_ANNOTATION, Structure};

fn many() -> Structure {
    attached(&ENTRY.replace("'(T)(R)'", "'(T,R,T,R,T)'"))
}

fn labels(structure: &Structure) -> Vec<String> {
    structure
        .data()
        .chains()
        .filter_map(ChainRef::label)
        .map(str::to_owned)
        .collect()
}

fn instance_ids(structure: &Structure) -> Vec<i64> {
    let Some(AtomAnnotation::Integer(column)) = structure.annotations().get(INSTANCE_ID_ANNOTATION)
    else {
        panic!("instance identity annotation absent")
    };
    column.values().to_vec()
}

#[test]
fn chunks_concatenate_to_the_whole_expansion() {
    let source = many();
    let view = source.assembly("1").expect("assembly");
    let whole = view.materialize().expect("whole");
    for bound in [1, 2, 3, 4, 100] {
        let pieces = view
            .chunks(bound)
            .collect::<Result<Vec<_>, _>>()
            .expect("pieces");
        let atoms: usize = pieces.iter().map(|piece| piece.atom_count() as usize).sum();
        assert_eq!(atoms, whole.atom_count() as usize);
        let joined: Vec<[f32; 3]> = pieces
            .iter()
            .flat_map(|piece| piece.positions().iter().copied())
            .collect();
        assert_eq!(joined, whole.positions().to_vec(), "bound {bound}");
        let names: Vec<String> = pieces.iter().flat_map(labels).collect();
        assert_eq!(names, labels(&whole), "bound {bound}");
        let ids: Vec<i64> = pieces.iter().flat_map(instance_ids).collect();
        assert_eq!(ids, instance_ids(&whole), "bound {bound}");
        let sites: Vec<Option<u32>> = pieces
            .iter()
            .flat_map(|piece| piece.data().atoms().map(|atom| atom.atom_site_id()))
            .collect();
        let expected: Vec<Option<u32>> = whole
            .data()
            .atoms()
            .map(|atom| atom.atom_site_id())
            .collect();
        assert_eq!(sites, expected, "bound {bound}");
    }
}

#[test]
fn a_placement_is_never_split_and_pieces_respect_the_bound() {
    let source = many();
    let view = source.assembly("1").expect("assembly");
    let sizes: Vec<u32> = view
        .chunks(3)
        .map(|piece| piece.expect("piece").atom_count())
        .collect();
    // Each placement holds two atoms; a piece closes at four, never at three.
    assert_eq!(sizes, [4, 4, 2]);
    assert_eq!(view.chunks(0).count(), 5);
}
