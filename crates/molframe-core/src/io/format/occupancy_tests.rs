use super::super::options::{OccupancyPolicy, ParseMode, ReadOptions};
use crate::chunk::{AtomRecord, ChunkBuilder};
use crate::column::Presence;
use crate::diagnostic::Code;
use crate::element::Element;
use crate::index::ResidueIndex;
use crate::optional::{OptionalI32, OptionalSymbol};
use crate::structure::{CoordinateStore, Structure, StructureData};
use crate::symbol::{AltId, SymbolId};
use crate::topology::ResidueRecord;

fn structure(occupancies: &[f32]) -> Structure {
    let mut data = StructureData::empty();
    let mut builder = ChunkBuilder::new();
    let count = u32::try_from(occupancies.len()).expect("small");
    data.topology
        .residues
        .push(
            ResidueRecord {
                label_comp_id: SymbolId::from_raw(0),
                auth_comp_id: OptionalSymbol::NONE,
                label_seq_id: OptionalI32::NONE,
                auth_seq_id: OptionalI32::NONE,
                ins_code: OptionalSymbol::NONE,
                het: false,
            },
            0..count,
        )
        .expect("small residue table");
    for (serial, occupancy) in (1..).zip(occupancies) {
        builder.push(AtomRecord {
            position: Some([0.0; 3]),
            element: Element::CARBON,
            atom_name: SymbolId::from_raw(0),
            auth_atom_name: OptionalSymbol::NONE,
            alternate_component_id: OptionalSymbol::NONE,
            alt_id: AltId::BLANK,
            residue: ResidueIndex::new(0),
            occupancy: (*occupancy, Presence::Present),
            b_factor: (0.0, Presence::Present),
            formal_charge: (0, Presence::Inapplicable),
            atom_site_id: serial,
        });
    }
    let (chunks, coords) = builder.finish();
    data.chunks = chunks.into();
    data.coords = CoordinateStore::Single(coords);
    Structure::new(data)
}

fn occupancies(structure: &Structure) -> Vec<f32> {
    structure
        .data()
        .chunks
        .iter()
        .flat_map(|chunk| (0..chunk.len()).filter_map(|i| chunk.occupancy(i).map(|o| o.0)))
        .collect()
}

#[test]
fn keep_preserves_values_and_reports_once_without_failing_a_permissive_read() {
    let options = ReadOptions::new();
    let (read, findings) = options
        .finish(structure(&[1.0, 1.5, 2.0]), [])
        .expect("permissive read keeps");
    assert_eq!(occupancies(&read), [1.0, 1.5, 2.0]);
    let hits: Vec<_> = findings
        .iter()
        .filter(|f| f.code() == Code::E3009)
        .collect();
    assert_eq!(hits.len(), 1);
}

#[test]
fn keep_fails_a_strict_read() {
    let options = ReadOptions::new().mode(ParseMode::Strict);
    assert!(options.finish(structure(&[1.5]), []).is_err());
}

#[test]
fn clamp_rewrites_out_of_range_values_and_reports_a_warning() {
    let options = ReadOptions::new()
        .mode(ParseMode::Strict)
        .occupancy_policy(OccupancyPolicy::Clamp);
    let (read, findings) = options
        .finish(structure(&[0.5, 1.5, -0.25]), [])
        .expect("clamping repairs the file");
    assert_eq!(occupancies(&read), [0.5, 1.0, 0.0]);
    assert!(findings.iter().any(|f| f.code() == Code::W3014));
}

#[test]
fn reject_fails_even_a_permissive_read() {
    let options = ReadOptions::new().occupancy_policy(OccupancyPolicy::Reject);
    assert!(options.finish(structure(&[1.5]), []).is_err());
}

#[test]
fn in_range_occupancy_reports_nothing_under_every_policy() {
    for policy in [
        OccupancyPolicy::Keep,
        OccupancyPolicy::Clamp,
        OccupancyPolicy::Reject,
    ] {
        let options = ReadOptions::new().occupancy_policy(policy);
        let (_, findings) = options.finish(structure(&[0.0, 1.0]), []).expect("clean");
        assert!(findings.is_empty());
    }
}
