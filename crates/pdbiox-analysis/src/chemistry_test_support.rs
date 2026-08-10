use std::collections::BTreeSet;

use pdbiox_core::annotation::{AnnotationColumn, AtomAnnotation};
use pdbiox_core::column::Presence;
use pdbiox_core::structure::Structure;

pub(crate) fn aromatic(structure: &Structure, atoms: impl IntoIterator<Item = u32>) -> Structure {
    let selected: BTreeSet<_> = atoms.into_iter().collect();
    let entries = (0..structure.atom_count()).map(|atom| {
        if selected.contains(&atom) {
            (true, Presence::Present)
        } else {
            (false, Presence::Inapplicable)
        }
    });
    let Ok(column) = AnnotationColumn::from_entries(entries) else {
        panic!("test annotation fits the column index domain");
    };
    with_annotation(
        structure,
        pdbiox_core::AROMATIC_ATOM_ANNOTATION,
        AtomAnnotation::Boolean(column),
    )
}

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
        pdbiox_core::FORMAL_CHARGE_ANNOTATION,
        AtomAnnotation::Integer(column),
    )
}

fn with_annotation(structure: &Structure, name: &str, annotation: AtomAnnotation) -> Structure {
    let mut data = structure.data().clone();
    data.annotations.insert(name, annotation);
    Structure::new(data)
}
