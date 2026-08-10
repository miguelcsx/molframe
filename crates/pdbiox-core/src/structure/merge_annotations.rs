//! Concatenation of typed custom annotation columns.

use super::Structure;
use super::merge::{remap_symbol, single};
use crate::annotation::{AnnotationColumn, AtomAnnotation, AtomAnnotations};
use crate::column::Presence;
use crate::diagnostic::{Code, Diagnostic};
use crate::symbol::{Interner, SymbolId};
use std::collections::BTreeSet;

pub(super) fn merge_annotations(
    sources: &[Structure],
    dictionary: &mut Interner,
) -> Result<AtomAnnotations, Vec<Diagnostic>> {
    let names: BTreeSet<&str> = sources
        .iter()
        .flat_map(|source| source.annotations().iter().map(|(name, _)| name))
        .collect();
    let mut merged = AtomAnnotations::default();
    for name in names {
        let prototype = sources
            .iter()
            .find_map(|source| source.annotations().get(name));
        let Some(prototype) = prototype else {
            continue;
        };
        let column = merge_annotation(name, prototype, sources, dictionary)?;
        let _ = merged.insert(name, column);
    }
    Ok(merged)
}

fn merge_annotation(
    name: &str,
    prototype: &AtomAnnotation,
    sources: &[Structure],
    dictionary: &mut Interner,
) -> Result<AtomAnnotation, Vec<Diagnostic>> {
    match prototype {
        AtomAnnotation::Boolean(_) => merge_plain(name, sources, false, AtomAnnotation::Boolean),
        AtomAnnotation::Integer(_) => merge_plain(name, sources, 0_i64, AtomAnnotation::Integer),
        AtomAnnotation::Real(_) => merge_plain(name, sources, 0.0_f64, AtomAnnotation::Real),
        AtomAnnotation::Symbol(_) => merge_symbols(name, sources, dictionary),
    }
}

trait AnnotationType: Copy {
    fn column(annotation: &AtomAnnotation) -> Option<&AnnotationColumn<Self>>;
}

impl AnnotationType for bool {
    fn column(annotation: &AtomAnnotation) -> Option<&AnnotationColumn<Self>> {
        match annotation {
            AtomAnnotation::Boolean(column) => Some(column),
            _ => None,
        }
    }
}

impl AnnotationType for i64 {
    fn column(annotation: &AtomAnnotation) -> Option<&AnnotationColumn<Self>> {
        match annotation {
            AtomAnnotation::Integer(column) => Some(column),
            _ => None,
        }
    }
}

impl AnnotationType for f64 {
    fn column(annotation: &AtomAnnotation) -> Option<&AnnotationColumn<Self>> {
        match annotation {
            AtomAnnotation::Real(column) => Some(column),
            _ => None,
        }
    }
}

fn merge_plain<T: AnnotationType>(
    name: &str,
    sources: &[Structure],
    absent: T,
    wrap: fn(AnnotationColumn<T>) -> AtomAnnotation,
) -> Result<AtomAnnotation, Vec<Diagnostic>> {
    let mut entries = Vec::new();
    for source in sources {
        match source.annotations().get(name) {
            Some(annotation) => {
                let Some(column) = T::column(annotation) else {
                    return Err(annotation_type_error(name));
                };
                entries.extend((0..column.len()).filter_map(|atom| column.get(atom)));
            }
            None => {
                entries.extend((0..source.atom_count()).map(|_| (absent, Presence::Inapplicable)));
            }
        }
    }
    AnnotationColumn::from_entries(entries)
        .map(wrap)
        .map_err(|_| single(Diagnostic::new(Code::E6009)))
}

fn merge_symbols(
    name: &str,
    sources: &[Structure],
    dictionary: &mut Interner,
) -> Result<AtomAnnotation, Vec<Diagnostic>> {
    let mut entries = Vec::new();
    for source in sources {
        match source.annotations().get(name) {
            Some(AtomAnnotation::Symbol(column)) => {
                for atom in 0..column.len() {
                    let Some((symbol, presence)) = column.get(atom) else {
                        continue;
                    };
                    let remapped = remap_symbol(source, symbol, dictionary).map_err(single)?;
                    entries.push((remapped, presence));
                }
            }
            Some(_) => return Err(annotation_type_error(name)),
            None => entries.extend(
                (0..source.atom_count()).map(|_| (SymbolId::from_raw(0), Presence::Inapplicable)),
            ),
        }
    }
    AnnotationColumn::from_entries(entries)
        .map(AtomAnnotation::Symbol)
        .map_err(|_| single(Diagnostic::new(Code::E6009)))
}

fn annotation_type_error(name: &str) -> Vec<Diagnostic> {
    single(Diagnostic::new(Code::E3013).with_context("annotation", name))
}
