use super::{
    RamachandranBasin, RamachandranError, RamachandranOptions, RamachandranRegion, classify,
    ramachandran, ramachandran_outliers,
};
use crate::{ReferenceDistribution, ReferenceLibrary};
use pdbiox_chem::PolymerAtomRole;
use pdbiox_core::BondOrder;
use pdbiox_core::bond::{BondProvenance, BondRecord, BondTableBuilder};
use pdbiox_core::index::AtomIndex;
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_core::structure::Structure;
use pdbiox_core::{AnnotationColumn, AtomAnnotation};

fn library() -> ReferenceLibrary {
    let alpha = ReferenceDistribution::grid(
        "right_alpha",
        vec![-180.0, 0.0, 180.0],
        vec![-180.0, 0.0, 180.0],
        vec![90.0, 10.0, 0.0, 0.0],
    )
    .unwrap_or_else(|error| panic!("alpha grid failed: {error}"));
    let beta = ReferenceDistribution::grid(
        "beta",
        vec![-180.0, 0.0, 180.0],
        vec![-180.0, 0.0, 180.0],
        vec![1.0, 80.0, 0.0, 0.0],
    )
    .unwrap_or_else(|error| panic!("beta grid failed: {error}"));
    ReferenceLibrary::new("rama", "2026.08", [alpha, beta])
        .unwrap_or_else(|error| panic!("library failed: {error}"))
}

fn options(library: &ReferenceLibrary) -> RamachandranOptions<'_> {
    let alpha = RamachandranBasin::new(RamachandranRegion::AlphaHelixRight, "right_alpha")
        .unwrap_or_else(|error| panic!("alpha basin failed: {error}"));
    let beta = RamachandranBasin::new(RamachandranRegion::BetaSheet, "beta")
        .unwrap_or_else(|error| panic!("beta basin failed: {error}"));
    RamachandranOptions::new(library, [alpha, beta], 0.5)
        .unwrap_or_else(|error| panic!("options failed: {error}"))
}

#[test]
fn classification_uses_selected_versioned_grids() {
    let library = library();
    let options = options(&library);
    let (alpha, assessment) = classify(-63.0, -43.0, &options)
        .unwrap_or_else(|error| panic!("classification failed: {error}"));
    let (beta, _) = classify(-135.0, 135.0, &options)
        .unwrap_or_else(|error| panic!("classification failed: {error}"));
    assert_eq!(alpha, RamachandranRegion::AlphaHelixRight);
    assert_eq!(beta, RamachandranRegion::BetaSheet);
    assert_eq!(assessment.set.as_ref(), "rama");
    assert_eq!(assessment.version.as_ref(), "2026.08");
    assert_eq!(assessment.distribution.as_ref(), "right_alpha");
}

#[test]
fn probability_threshold_is_explicit_and_controls_outliers() {
    let library = library();
    let alpha = RamachandranBasin::new(RamachandranRegion::AlphaHelixRight, "right_alpha")
        .unwrap_or_else(|error| panic!("basin failed: {error}"));
    let options = RamachandranOptions::new(&library, [alpha], 1.0)
        .unwrap_or_else(|error| panic!("options failed: {error}"));
    let (region, assessment) = classify(-60.0, -40.0, &options)
        .unwrap_or_else(|error| panic!("classification failed: {error}"));
    assert_eq!(region, RamachandranRegion::Outlier);
    assert!(assessment.probability < 1.0);
}

#[test]
fn invalid_or_ambiguous_policies_are_rejected() {
    let library = library();
    assert_eq!(
        RamachandranOptions::new(&library, [], 0.1),
        Err(RamachandranError::NoBasins)
    );
    let one = RamachandranBasin::new(RamachandranRegion::BetaSheet, "beta")
        .unwrap_or_else(|error| panic!("basin failed: {error}"));
    assert_eq!(
        RamachandranOptions::new(&library, [one.clone(), one], 0.1),
        Err(RamachandranError::DuplicateBasin)
    );
    assert_eq!(
        RamachandranBasin::new(RamachandranRegion::Outlier, "outlier"),
        Err(RamachandranError::OutlierBasin)
    );
    let alpha = RamachandranBasin::new(RamachandranRegion::AlphaHelixRight, "beta")
        .unwrap_or_else(|error| panic!("alpha basin failed: {error}"));
    let beta = RamachandranBasin::new(RamachandranRegion::BetaSheet, "beta")
        .unwrap_or_else(|error| panic!("beta basin failed: {error}"));
    assert_eq!(
        RamachandranOptions::new(&library, [alpha, beta], 0.1),
        Err(RamachandranError::DuplicateDistribution)
    );
    let missing = RamachandranBasin::new(RamachandranRegion::BetaSheet, "absent")
        .unwrap_or_else(|error| panic!("missing basin failed: {error}"));
    assert_eq!(
        RamachandranOptions::new(&library, [missing], 0.1),
        Err(RamachandranError::Reference(crate::ReferenceError::Missing))
    );
}

const SOURCE: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 N N ALA A 1 0 0 0\n\
ATOM 2 C CA ALA A 1 1.46 0 0\n\
ATOM 3 C C ALA A 1 2.0 1.2 0\n\
ATOM 4 N N ALA A 2 3.3 1.4 0\n\
ATOM 5 C CA ALA A 2 4.0 2.5 0\n\
ATOM 6 C C ALA A 2 5.4 2.5 0\n\
ATOM 7 N N ALA A 3 6.0 3.6 0\n\
ATOM 8 C CA ALA A 3 7.0 3.6 0\n\
ATOM 9 C C ALA A 3 8.0 4.0 0\n";

#[test]
fn structure_records_retain_reference_provenance() {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    let (structure, _): (Structure, _) = match pdbiox_cif::read(&input, &ReadOptions::new()) {
        Ok(result) => result,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    let mut data = structure.data().clone();
    let mut bonds = BondTableBuilder::new();
    for (carbon, nitrogen) in [(2, 3), (5, 6)] {
        bonds.push(BondRecord {
            atom_a: AtomIndex::new(carbon),
            atom_b: AtomIndex::new(nitrogen),
            order: BondOrder::Single,
            provenance: BondProvenance::User,
        });
    }
    data.bonds = bonds.finish();
    let roles = [
        PolymerAtomRole::PROTEIN_NITROGEN,
        PolymerAtomRole::PROTEIN_ALPHA_CARBON,
        PolymerAtomRole::PROTEIN_CARBONYL_CARBON,
    ];
    data.annotations.insert(
        pdbiox_core::POLYMER_ATOM_ROLE_ANNOTATION,
        AtomAnnotation::Integer(
            AnnotationColumn::from_values(
                (0..9)
                    .map(|index| roles[index % roles.len()].code())
                    .collect(),
            )
            .expect("small annotation column"),
        ),
    );
    let structure = Structure::new(data);
    let library = library();
    let options = options(&library);
    let records = ramachandran(&structure, &options)
        .unwrap_or_else(|error| panic!("Ramachandran failed: {error}"));
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].residue.get(), 1);
    assert_eq!(records[0].reference.version.as_ref(), "2026.08");

    let outliers = ramachandran_outliers(&structure, &options)
        .unwrap_or_else(|error| panic!("outlier view failed: {error}"));
    assert!(
        outliers
            .iter()
            .all(|record| record.region == RamachandranRegion::Outlier)
    );
    assert!(outliers.len() <= records.len());
}
