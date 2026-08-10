//! Explicit connectivity lowering from `struct_conn`.

use crate::document::{Category, CifValue, DataBlock};
use crate::lower::diagnostics::at_source_row;
use pdbiox_core::bond::{BondOrder, BondProvenance, BondRecord, BondTableBuilder};
use pdbiox_core::diagnostic::{Code, Diagnostic, Diagnostics};
use pdbiox_core::index::AtomIndex;
use pdbiox_core::structure::{AtomRef, ChainRef, ResidueRef, StructureData};
use std::collections::BTreeMap;

type Key<'a> = (&'a str, Option<i32>, &'a str, &'a str, Option<&'a str>);
type AtomMap<'a> = BTreeMap<Key<'a>, Option<AtomIndex>>;

pub(super) fn read(block: &DataBlock, data: &mut StructureData, findings: &mut Diagnostics) {
    let Some(category) = block.category("struct_conn") else {
        return;
    };
    let (label, auth) = atom_maps(data);
    let mut bonds = BondTableBuilder::new();
    for row in 0..category.row_count() {
        if !is_connectivity(category.text("conn_type_id", row)) {
            continue;
        }
        let atom_a = endpoint(category, row, 1, &label, &auth);
        let atom_b = endpoint(category, row, 2, &label, &auth);
        let (Some(atom_a), Some(atom_b)) = (atom_a, atom_b) else {
            findings.push(at_source_row(
                Diagnostic::new(Code::E3006).in_category("struct_conn"),
                row,
            ));
            continue;
        };
        bonds.push(BondRecord {
            atom_a,
            atom_b,
            order: order(category.text("pdbx_value_order", row)),
            provenance: provenance(category.text("details", row)),
        });
    }
    data.bonds = bonds.finish();
}

fn atom_maps(data: &StructureData) -> (AtomMap<'_>, AtomMap<'_>) {
    let mut label = BTreeMap::new();
    let mut auth = BTreeMap::new();
    for chain in data.chains() {
        for residue in chain.residues() {
            for atom in residue.atoms() {
                if let Some(key) = label_key(chain, residue, atom) {
                    insert(&mut label, key, atom.index());
                }
                if let Some(key) = auth_key(chain, residue, atom) {
                    insert(&mut auth, key, atom.index());
                }
            }
        }
    }
    (label, auth)
}

fn label_key<'a>(
    chain: ChainRef<'a>,
    residue: ResidueRef<'a>,
    atom: AtomRef<'a>,
) -> Option<Key<'a>> {
    Some((
        chain.label()?,
        residue.label_seq_id(),
        atom.component_name()?,
        atom.name()?,
        alt_text(atom),
    ))
}

fn auth_key<'a>(
    chain: ChainRef<'a>,
    residue: ResidueRef<'a>,
    atom: AtomRef<'a>,
) -> Option<Key<'a>> {
    Some((
        chain.auth_label()?,
        residue.auth_seq_id(),
        match residue.auth_name() {
            Some(name) => name,
            None => atom.component_name()?,
        },
        match atom.auth_name() {
            Some(name) => name,
            None => atom.name()?,
        },
        alt_text(atom),
    ))
}

fn alt_text(atom: AtomRef<'_>) -> Option<&str> {
    atom.alt_label()
}

fn insert<'a>(map: &mut AtomMap<'a>, key: Key<'a>, atom: AtomIndex) {
    match map.entry(key) {
        std::collections::btree_map::Entry::Vacant(entry) => {
            entry.insert(Some(atom));
        }
        std::collections::btree_map::Entry::Occupied(mut entry) => {
            entry.insert(None);
        }
    }
}

fn endpoint(
    category: &Category,
    row: usize,
    partner: u8,
    label: &AtomMap<'_>,
    auth: &AtomMap<'_>,
) -> Option<AtomIndex> {
    let label_key = row_key(category, row, partner, "label");
    if let Some(key) = label_key.as_ref()
        && let Some(atom) = label.get(&as_borrowed(key)).copied().flatten()
    {
        return Some(atom);
    }
    let auth_key = row_key(category, row, partner, "auth")?;
    auth.get(&as_borrowed(&auth_key)).copied().flatten()
}

type OwnedKey = (Box<str>, Option<i32>, Box<str>, Box<str>, Option<Box<str>>);

fn row_key(category: &Category, row: usize, partner: u8, namespace: &str) -> Option<OwnedKey> {
    let field = |name: &str| {
        category
            .identifier(name, row)
            .map(|value| value.into_owned().into())
    };
    let prefix = format!("ptnr{partner}_{namespace}_");
    Some((
        field(&format!("{prefix}asym_id"))?,
        category
            .value(&format!("{prefix}seq_id"), row)
            .and_then(CifValue::as_integer)
            .and_then(|value| i32::try_from(value).ok()),
        field(&format!("{prefix}comp_id"))?,
        field(&format!("{prefix}atom_id"))?,
        category
            .identifier(&format!("pdbx_ptnr{partner}_label_alt_id"), row)
            .map(|value| value.into_owned().into()),
    ))
}

fn as_borrowed(key: &OwnedKey) -> Key<'_> {
    (&key.0, key.1, &key.2, &key.3, key.4.as_deref())
}

fn is_connectivity(kind: Option<&str>) -> bool {
    match kind {
        Some(kind) => {
            let lower = kind.to_ascii_lowercase();
            lower.starts_with("covale")
                || lower.starts_with("disulf")
                || lower.starts_with("metalc")
                || lower.starts_with("modres")
        }
        None => true,
    }
}

fn order(value: Option<&str>) -> BondOrder {
    match value.map(str::to_ascii_uppercase).as_deref() {
        Some("SING" | "SINGLE") => BondOrder::Single,
        Some("DOUB" | "DOUBLE") => BondOrder::Double,
        Some("TRIP" | "TRIPLE") => BondOrder::Triple,
        Some("QUAD" | "QUADRUPLE") => BondOrder::Quadruple,
        Some("AROM" | "AROMATIC") => BondOrder::Aromatic,
        _ => BondOrder::Unknown,
    }
}

fn provenance(value: Option<&str>) -> BondProvenance {
    match value {
        Some("pdbiox:ccd") => BondProvenance::ChemicalComponentDictionary,
        Some("pdbiox:inferred-distance") => BondProvenance::InferredDistance,
        Some("pdbiox:user") => BondProvenance::User,
        _ => BondProvenance::File,
    }
}

#[cfg(test)]
#[path = "bonds_tests.rs"]
mod tests;
