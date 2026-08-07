use super::*;
use crate::chunk::{AtomRecord, ChunkBuilder};
use crate::column::Presence;
use crate::contract::AltlocPolicy;
use crate::element::Element;
use crate::index::ResidueIndex;
use crate::optional::{OptionalI32, OptionalSymbol};
use crate::topology::{ChainRecord, EntityKind, PolymerKind, ResidueRecord};

fn disordered() -> Structure {
    let mut data = super::super::StructureData::empty();
    let entity_id = data.dictionary.intern("1").expect("small dictionary");
    let component = data.dictionary.intern("ALA").expect("small dictionary");
    let chain = data.dictionary.intern("A").expect("small dictionary");
    let alt_a = AltId::labelled(data.dictionary.intern("A").expect("small dictionary"));
    let alt_b = AltId::labelled(data.dictionary.intern("B").expect("small dictionary"));
    let entity = data.topology.entities.push(
        entity_id,
        EntityKind::Polymer,
        OptionalSymbol::NONE,
        &[component],
    );
    let residue = data.topology.residues.push(
        ResidueRecord {
            label_comp_id: component,
            auth_comp_id: OptionalSymbol::NONE,
            label_seq_id: OptionalI32::some(1),
            auth_seq_id: OptionalI32::some(1),
            ins_code: OptionalSymbol::NONE,
            het: false,
        },
        0..5,
    );
    data.topology.chains.push(
        ChainRecord {
            label_asym_id: chain,
            auth_asym_id: OptionalSymbol::some(chain),
            entity,
            polymer_kind: PolymerKind::Protein,
        },
        0..1,
    );
    data.topology.models.push(1, 0..1);

    let mut builder = ChunkBuilder::new();
    builder.start_model(0);
    for (index, (name, alt, occupancy)) in [
        ("N", AltId::BLANK, 1.0),
        ("CA", alt_a, 0.9),
        ("CA", alt_b, 0.1),
        ("CB", alt_a, 0.1),
        ("CB", alt_b, 0.9),
    ]
    .into_iter()
    .enumerate()
    {
        let atom_name = data.dictionary.intern(name).expect("small dictionary");
        builder.push(AtomRecord {
            position: Some([index as f32, 0.0, 0.0]),
            element: Element::CARBON,
            atom_name,
            auth_atom_name: OptionalSymbol::NONE,
            alternate_component_id: OptionalSymbol::NONE,
            alt_id: alt,
            residue: ResidueIndex::new(residue.get()),
            occupancy: (occupancy, Presence::Present),
            b_factor: (0.0, Presence::Present),
            formal_charge: (0, Presence::Inapplicable),
            atom_site_id: index as u32 + 1,
        });
    }
    let (chunks, coords) = builder.finish();
    data.chunks = chunks.into();
    data.coords = super::super::CoordinateStore::Single(coords);
    Structure::new(data)
}

fn selected(structure: &Structure, altloc: AltlocPolicy) -> (Vec<u32>, Vec<Code>, Status) {
    let policy = AnalysisPolicy::default().with_altloc(altloc);
    let result = structure.resolve_altlocs(&policy);
    (
        result.value.iter().collect(),
        result.warnings.iter().map(Diagnostic::code).collect(),
        result.status,
    )
}

#[test]
fn every_altloc_policy_is_explicit_and_deterministic() {
    let structure = disordered();
    assert_eq!(
        selected(&structure, AltlocPolicy::KeepAll).0,
        [0, 1, 2, 3, 4]
    );
    assert_eq!(
        selected(&structure, AltlocPolicy::ConformerConsistent).0,
        [0, 1, 3]
    );
    assert_eq!(selected(&structure, AltlocPolicy::First).0, [0, 1, 3]);
    assert_eq!(
        selected(&structure, AltlocPolicy::HighestOccupancyPerResidue).0,
        [0, 1, 3]
    );
    assert_eq!(
        selected(&structure, AltlocPolicy::Label("B".into())).0,
        [0, 2, 4]
    );
}

#[test]
fn hazardous_per_atom_selection_warns_and_can_mix_labels() {
    let structure = disordered();
    let (atoms, warnings, _) = selected(&structure, AltlocPolicy::HighestOccupancyPerAtom);
    assert_eq!(atoms, [0, 1, 4]);
    assert_eq!(warnings, [Code::W6001]);
}

#[test]
fn unresolved_and_missing_labels_are_visible_in_the_result() {
    let structure = disordered();
    assert_eq!(
        selected(&structure, AltlocPolicy::KeepAll).2,
        Status::Ambiguous
    );
    let (atoms, warnings, _) = selected(&structure, AltlocPolicy::Label("Z".into()));
    assert_eq!(atoms, [0]);
    assert_eq!(warnings, [Code::W4003]);
}
