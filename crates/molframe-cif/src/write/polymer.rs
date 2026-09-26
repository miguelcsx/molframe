//! The `_entity_poly.type` a writer declares for each polymer entity.
//!
//! The reader takes a chain's polymer kind from its entity's declared type, so
//! a writer that drops the declaration hands back a structure whose proteins
//! are only "some polymer". Declaring it is what keeps a round trip faithful.

use molframe_core::index::EntityIndex;
use molframe_core::structure::Structure;
use molframe_core::topology::PolymerKind;

/// Each polymer entity with a kind this format can declare, in entity order.
///
/// The kind is taken from the entity's first chain that has one. An entity
/// whose chains are only an unspecified polymer declares nothing, because
/// `other` would state more than the structure knows.
#[doc(hidden)]
#[must_use]
pub fn declared_polymer_types(structure: &Structure) -> Vec<(EntityIndex, &'static str)> {
    let topology = &structure.data().topology;
    let mut declared = Vec::new();
    for entity in topology.entities.iter() {
        let kind = topology.chains.iter().find_map(|chain| {
            (topology.chains.entity(chain) == Some(entity))
                .then(|| topology.chains.polymer_kind(chain))
                .flatten()
                .and_then(declared_type)
        });
        if let Some(kind) = kind {
            declared.push((entity, kind));
        }
    }
    declared
}

/// The PDBx/mmCIF `_entity_poly.type` value for a polymer kind, where it has one.
///
/// `Protein` is written as the L form, which nearly every deposition declares.
fn declared_type(kind: PolymerKind) -> Option<&'static str> {
    match kind {
        PolymerKind::Protein => Some("polypeptide(L)"),
        PolymerKind::Dna => Some("polydeoxyribonucleotide"),
        PolymerKind::Rna => Some("polyribonucleotide"),
        PolymerKind::NucleicHybrid => Some("polydeoxyribonucleotide/polyribonucleotide hybrid"),
        PolymerKind::Saccharide => Some("polysaccharide(D)"),
        _ => None,
    }
}
