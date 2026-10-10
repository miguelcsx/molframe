use molframe_core::annotation::{AnnotationColumn, AtomAnnotation};
use molframe_core::column::Presence;
use molframe_core::structure::Structure;

pub(crate) fn charges(structure: &Structure, charges: &[(u32, i64)]) -> Structure {
    let entries = (0..structure.atom_count()).map(|atom| {
        charges
            .iter()
            .find(|(index, _)| *index == atom)
            .map_or((0, Presence::Inapplicable), |(_, charge)| {
                (*charge, Presence::Present)
            })
    });
    let Ok(column) = AnnotationColumn::from_entries(entries) else {
        panic!("test annotation fits the column index domain");
    };
    with_annotation(
        structure,
        molframe_core::FORMAL_CHARGE_ANNOTATION,
        AtomAnnotation::Integer(column),
    )
}

fn with_annotation(structure: &Structure, name: &str, annotation: AtomAnnotation) -> Structure {
    let mut data = structure.data().clone();
    data.annotations.insert(name, annotation);
    Structure::new(data)
}
