//! Per-atom CCD annotations and stereochemical observation.

use crate::{Component, StereoConfiguration};
use molframe_core::annotation::{AnnotationColumn, AtomAnnotation};
use molframe_core::column::Presence;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::structure::{ResidueRef, StructureData};

#[derive(Clone, Copy)]
pub(super) struct AtomChemistry {
    pub(super) kind: crate::ComponentKind,
    pub(super) aromatic: bool,
    pub(super) charge: i8,
    pub(super) donor: bool,
    pub(super) acceptor: bool,
    pub(super) stereo: Option<StereoConfiguration>,
}

pub(super) fn observed_stereo(
    component: &Component,
    residue: ResidueRef<'_>,
    centre_name: &str,
    declared: Option<StereoConfiguration>,
) -> Option<StereoConfiguration> {
    let declared = declared?;
    if declared == StereoConfiguration::Mixed {
        return Some(declared);
    }
    let reference = component
        .ideal_coordinates
        .as_deref()
        .or(component.model_coordinates.as_deref())?;
    let centre = component
        .atoms
        .iter()
        .position(|atom| atom.name.as_ref() == centre_name)?;
    let neighbours = first_three_neighbours(component, centre)?;
    let reference_volume = signed_volume(
        *reference.get(centre)?,
        [
            *reference.get(neighbours[0])?,
            *reference.get(neighbours[1])?,
            *reference.get(neighbours[2])?,
        ],
    );
    let observed_centre = residue.atom(centre_name)?.position()?;
    let observed_volume = signed_volume(
        observed_centre,
        [
            residue
                .atom(component.atoms.get(neighbours[0])?.name.as_ref())?
                .position()?,
            residue
                .atom(component.atoms.get(neighbours[1])?.name.as_ref())?
                .position()?,
            residue
                .atom(component.atoms.get(neighbours[2])?.name.as_ref())?
                .position()?,
        ],
    );
    if reference_volume.abs() <= f64::EPSILON || observed_volume.abs() <= f64::EPSILON {
        return None;
    }
    if reference_volume.is_sign_positive() == observed_volume.is_sign_positive() {
        Some(declared)
    } else {
        Some(match declared {
            StereoConfiguration::R => StereoConfiguration::S,
            StereoConfiguration::S => StereoConfiguration::R,
            StereoConfiguration::Mixed => StereoConfiguration::Mixed,
        })
    }
}

fn first_three_neighbours(component: &Component, centre: usize) -> Option<[usize; 3]> {
    let centre_name = component.atoms.get(centre)?.name.as_ref();
    let mut neighbours = [usize::MAX; 3];
    let mut count = 0;
    for bond in component.bonds.iter() {
        let other = if bond.atom_a.as_ref() == centre_name {
            component
                .atoms
                .iter()
                .position(|atom| atom.name == bond.atom_b)?
        } else if bond.atom_b.as_ref() == centre_name {
            component
                .atoms
                .iter()
                .position(|atom| atom.name == bond.atom_a)?
        } else {
            continue;
        };
        if neighbours[..count].contains(&other) {
            continue;
        }
        neighbours[count] = other;
        count += 1;
        if count == neighbours.len() {
            return Some(neighbours);
        }
    }
    None
}

fn signed_volume(centre: [f32; 3], neighbours: [[f32; 3]; 3]) -> f64 {
    let vectors = neighbours.map(|point| {
        [
            f64::from(point[0] - centre[0]),
            f64::from(point[1] - centre[1]),
            f64::from(point[2] - centre[2]),
        ]
    });
    let cross = [
        vectors[1][1] * vectors[2][2] - vectors[1][2] * vectors[2][1],
        vectors[1][2] * vectors[2][0] - vectors[1][0] * vectors[2][2],
        vectors[1][0] * vectors[2][1] - vectors[1][1] * vectors[2][0],
    ];
    vectors[0][0] * cross[0] + vectors[0][1] * cross[1] + vectors[0][2] * cross[2]
}

pub(super) fn attach_annotations(
    data: &mut StructureData,
    chemistry: &[Option<AtomChemistry>],
) -> Result<(), Diagnostic> {
    let unknown = data.dictionary.intern("unknown").map_err(dictionary_full)?;
    let mut kinds = Vec::with_capacity(chemistry.len());
    let mut aromatic = Vec::with_capacity(chemistry.len());
    let mut charges = Vec::with_capacity(chemistry.len());
    let mut donors = Vec::with_capacity(chemistry.len());
    let mut acceptors = Vec::with_capacity(chemistry.len());
    let mut stereo = Vec::with_capacity(chemistry.len());
    for annotation in chemistry {
        if let Some(annotation) = annotation {
            kinds.push((annotation.kind.code(), Presence::Present));
            aromatic.push((annotation.aromatic, Presence::Present));
            charges.push((i64::from(annotation.charge), Presence::Present));
            donors.push((annotation.donor, Presence::Present));
            acceptors.push((annotation.acceptor, Presence::Present));
            match annotation.stereo {
                Some(configuration) => {
                    let symbol = data
                        .dictionary
                        .intern(stereo_name(configuration))
                        .map_err(dictionary_full)?;
                    stereo.push((symbol, Presence::Present));
                }
                None => stereo.push((unknown, Presence::Inapplicable)),
            }
        } else {
            kinds.push((crate::ComponentKind::Unknown.code(), Presence::Unknown));
            aromatic.push((false, Presence::Unknown));
            charges.push((0, Presence::Unknown));
            donors.push((false, Presence::Unknown));
            acceptors.push((false, Presence::Unknown));
            stereo.push((unknown, Presence::Unknown));
        }
    }
    let annotations =
        [
            (
                molframe_core::COMPONENT_KIND_ANNOTATION,
                AtomAnnotation::Integer(AnnotationColumn::from_entries(kinds).map_err(
                    |error| annotation_capacity().with_context("cause", error.to_string()),
                )?),
            ),
            (
                molframe_core::AROMATIC_ATOM_ANNOTATION,
                AtomAnnotation::Boolean(AnnotationColumn::from_entries(aromatic).map_err(
                    |error| annotation_capacity().with_context("cause", error.to_string()),
                )?),
            ),
            (
                molframe_core::FORMAL_CHARGE_ANNOTATION,
                AtomAnnotation::Integer(AnnotationColumn::from_entries(charges).map_err(
                    |error| annotation_capacity().with_context("cause", error.to_string()),
                )?),
            ),
            (
                molframe_core::HBOND_DONOR_ANNOTATION,
                AtomAnnotation::Boolean(AnnotationColumn::from_entries(donors).map_err(
                    |error| annotation_capacity().with_context("cause", error.to_string()),
                )?),
            ),
            (
                molframe_core::HBOND_ACCEPTOR_ANNOTATION,
                AtomAnnotation::Boolean(AnnotationColumn::from_entries(acceptors).map_err(
                    |error| annotation_capacity().with_context("cause", error.to_string()),
                )?),
            ),
            (
                molframe_core::STEREO_CONFIGURATION_ANNOTATION,
                AtomAnnotation::Symbol(AnnotationColumn::from_entries(stereo).map_err(
                    |error| annotation_capacity().with_context("cause", error.to_string()),
                )?),
            ),
        ];
    for (name, annotation) in annotations {
        let _ = data.annotations.insert(name, annotation);
    }
    Ok(())
}

fn annotation_capacity() -> Diagnostic {
    Diagnostic::new(Code::E6009).with_message("annotation exceeds the supported atom range")
}

fn dictionary_full(_: molframe_core::DictionaryFull) -> Diagnostic {
    Diagnostic::new(Code::E1901).with_context("limit", "dictionary entries")
}

const fn stereo_name(stereo: StereoConfiguration) -> &'static str {
    match stereo {
        StereoConfiguration::R => "R",
        StereoConfiguration::S => "S",
        StereoConfiguration::Mixed => "mixed",
    }
}
