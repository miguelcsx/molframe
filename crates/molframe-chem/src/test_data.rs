//! Test helpers: decode the tables shipped in the repository's `data/` directory.

use crate::{
    BondReach, Component, ComponentAtom, ComponentBond, ComponentKind, IonicRadii, MetalSites,
    MetalTables, Provided, RadiusSource, TableIdentity, TableKind, TableTextError, VdwRadii,
};
use molframe_core::{BondOrder, Element};
use std::sync::{Arc, OnceLock};

/// The bytes of `data/<name>`, panicking with the reason when unreadable.
pub(crate) fn shipped_bytes(name: &str) -> Vec<u8> {
    let path = format!("{}/../../data/{name}", env!("CARGO_MANIFEST_DIR"));
    match std::fs::read(&path) {
        Ok(bytes) => bytes,
        Err(error) => panic!("{path}: {error}"),
    }
}

/// Decodes `data/<name>`, panicking with the reason when it is refused.
pub(crate) fn shipped<T: TableKind>(name: &str) -> Provided<T> {
    let path = format!("{}/../../data/{name}", env!("CARGO_MANIFEST_DIR"));
    let bytes = shipped_bytes(name);
    match T::decode(&bytes, TableIdentity::supplied("x", "x", "x", "x")) {
        Ok(provided) => provided,
        Err(error) => panic!("{path}: {error}"),
    }
}

/// The shipped Bondi table as a radius source.
pub(crate) fn bondi() -> RadiusSource {
    RadiusSource::from(shipped::<VdwRadii>("bondi.radii"))
}

/// The shipped bonding reach table, decoded once.
pub(crate) fn reach() -> &'static BondReach {
    static REACH: OnceLock<BondReach> = OnceLock::new();
    REACH.get_or_init(|| shipped::<BondReach>("macromolecular-reach.reach").table)
}

/// The shipped metal and Shannon tables.
pub(crate) fn metal_tables() -> MetalTables {
    MetalTables {
        sites: Arc::new(shipped::<MetalSites>("macromolecular-metals.metals")),
        ionic: Some(Arc::new(shipped::<IonicRadii>("shannon-ionic.ionic"))),
    }
}

/// A neutral, non-aromatic, non-leaving atom with no stereo configuration.
pub(crate) fn plain_atom(name: &str, element: Element) -> ComponentAtom {
    ComponentAtom {
        name: name.into(),
        alternate_name: None,
        element,
        charge: 0,
        aromatic: false,
        leaving: false,
        stereo: None,
    }
}

/// A bond with no stereo configuration, aromatic exactly when its order is.
pub(crate) fn plain_bond(atom_a: &str, atom_b: &str, order: BondOrder) -> ComponentBond {
    ComponentBond {
        atom_a: atom_a.into(),
        atom_b: atom_b.into(),
        order,
        aromatic: order == BondOrder::Aromatic,
        stereo: None,
    }
}

/// A component built from plain atoms and bonds, without coordinates.
pub(crate) fn plain_component(
    id: &str,
    kind: ComponentKind,
    atoms: &[(&str, Element)],
    bonds: &[(&str, &str, BondOrder)],
) -> Component {
    Component {
        id: id.into(),
        name: id.into(),
        kind,
        parent: None,
        one_letter_code: None,
        formula: None,
        atoms: atoms
            .iter()
            .map(|&(name, element)| plain_atom(name, element))
            .collect::<Vec<_>>()
            .into(),
        bonds: bonds
            .iter()
            .map(|&(atom_a, atom_b, order)| plain_bond(atom_a, atom_b, order))
            .collect::<Vec<_>>()
            .into(),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}

/// A whitespace-separated word table that refuses text containing `bad`.
#[derive(Debug, PartialEq)]
pub(crate) struct Words(pub(crate) Vec<String>);

impl TableKind for Words {
    const KIND: &'static str = "words";
    const EXTENSION: &'static str = "txt";

    fn decode(bytes: &[u8], fallback: TableIdentity) -> Result<Provided<Self>, TableTextError> {
        let text = crate::text_of("words table", bytes)?;
        if text.contains("bad") {
            return Err(TableTextError::at("words table", 1, "bad word"));
        }
        Ok(Provided {
            identity: fallback,
            table: Self(text.split_whitespace().map(str::to_owned).collect()),
        })
    }
}
