//! Explicit connectivity lowering from `struct_conn`.

use crate::document::{Category, CifValue, DataBlock};
use crate::lower::diagnostics::at_source_row;
use pdbiox_core::bond::{BondOrder, BondProvenance, BondRecord, BondTableBuilder};
use pdbiox_core::diagnostic::{Code, Diagnostic, Diagnostics};
use pdbiox_core::index::AtomIndex;
use pdbiox_core::structure::{AtomRef, ChainRef, ResidueRef, StructureData};
use pdbiox_core::symbol::SymbolId;
use std::collections::HashMap;

/// A connectivity endpoint expressed in structure-local integer identifiers.
///
/// Only keys referenced by `struct_conn` enter a map. A giant structure with a
/// handful of explicit bonds therefore uses memory proportional to the bonds,
/// rather than allocating two hash-table entries for every atom.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct AtomKey {
    asym: SymbolId,
    sequence: Option<i32>,
    component: SymbolId,
    atom: SymbolId,
    alternate: Option<SymbolId>,
}

#[derive(Clone, Copy, Default)]
enum Match {
    #[default]
    Missing,
    Unique(AtomIndex),
    Ambiguous,
}

impl Match {
    fn observe(&mut self, atom: AtomIndex) {
        *self = match self {
            Self::Missing => Self::Unique(atom),
            Self::Unique(_) | Self::Ambiguous => Self::Ambiguous,
        };
    }

    const fn unique(self) -> Option<AtomIndex> {
        match self {
            Self::Unique(atom) => Some(atom),
            Self::Missing | Self::Ambiguous => None,
        }
    }
}

type AtomMap = HashMap<AtomKey, Match>;

#[derive(Clone, Copy)]
struct Endpoint {
    label: Option<AtomKey>,
    auth: Option<AtomKey>,
}

impl Endpoint {
    fn resolve(self, label: &AtomMap, auth: &AtomMap) -> Option<AtomIndex> {
        self.label
            .and_then(|key| label.get(&key).copied().and_then(Match::unique))
            .or_else(|| {
                self.auth
                    .and_then(|key| auth.get(&key).copied().and_then(Match::unique))
            })
    }
}

struct Connection {
    row: usize,
    atom_a: Endpoint,
    atom_b: Endpoint,
    order: BondOrder,
    provenance: BondProvenance,
}

pub(super) fn read(block: &DataBlock, data: &mut StructureData, findings: &mut Diagnostics) {
    let Some(category) = block.category("struct_conn") else {
        return;
    };
    let connections = connections(category, data);
    let (mut label, mut auth) = requested_maps(&connections);
    match_atoms(data, &mut label, &mut auth);

    let mut bonds = BondTableBuilder::new();
    for connection in connections {
        let atom_a = connection.atom_a.resolve(&label, &auth);
        let atom_b = connection.atom_b.resolve(&label, &auth);
        let (Some(atom_a), Some(atom_b)) = (atom_a, atom_b) else {
            findings.push(at_source_row(
                Diagnostic::new(Code::E3006).in_category("struct_conn"),
                connection.row,
            ));
            continue;
        };
        bonds.push(BondRecord {
            atom_a,
            atom_b,
            order: connection.order,
            provenance: connection.provenance,
        });
    }
    data.bonds = bonds.finish();
}

fn connections(category: &Category, data: &StructureData) -> Vec<Connection> {
    let mut output = Vec::with_capacity(category.row_count());
    for row in 0..category.row_count() {
        if !is_connectivity(category.text("conn_type_id", row)) {
            continue;
        }
        output.push(Connection {
            row,
            atom_a: endpoint(category, row, 1, data),
            atom_b: endpoint(category, row, 2, data),
            order: order(category.text("pdbx_value_order", row)),
            provenance: provenance(category.text("details", row)),
        });
    }
    output
}

fn endpoint(category: &Category, row: usize, partner: u8, data: &StructureData) -> Endpoint {
    Endpoint {
        label: row_key(category, row, partner, Namespace::Label, data),
        auth: row_key(category, row, partner, Namespace::Author, data),
    }
}

#[derive(Clone, Copy)]
enum Namespace {
    Label,
    Author,
}

fn row_key(
    category: &Category,
    row: usize,
    partner: u8,
    namespace: Namespace,
    data: &StructureData,
) -> Option<AtomKey> {
    let items = endpoint_items(partner, namespace)?;
    Some(AtomKey {
        asym: symbol(category, items.asym, row, data)?,
        sequence: sequence(category, items.sequence, items.fallback_sequence, row),
        component: symbol(category, items.component, row, data)?,
        atom: symbol(category, items.atom, row, data)?,
        alternate: category
            .identifier(items.alternate, row)
            .and_then(|value| data.dictionary.get(&value)),
    })
}

fn sequence(category: &Category, primary: &str, fallback: Option<&str>, row: usize) -> Option<i32> {
    integer(category, primary, row)
        .or_else(|| fallback.and_then(|item| integer(category, item, row)))
}

fn integer(category: &Category, item: &str, row: usize) -> Option<i32> {
    category
        .value(item, row)
        .and_then(CifValue::as_integer)
        .and_then(|value| i32::try_from(value).ok())
}

fn symbol(category: &Category, item: &str, row: usize, data: &StructureData) -> Option<SymbolId> {
    category
        .identifier(item, row)
        .and_then(|value| data.dictionary.get(&value))
}

struct EndpointItems {
    asym: &'static str,
    sequence: &'static str,
    fallback_sequence: Option<&'static str>,
    component: &'static str,
    atom: &'static str,
    alternate: &'static str,
}

fn endpoint_items(partner: u8, namespace: Namespace) -> Option<EndpointItems> {
    match (partner, namespace) {
        (1, Namespace::Label) => Some(EndpointItems {
            asym: "ptnr1_label_asym_id",
            sequence: "ptnr1_label_seq_id",
            fallback_sequence: Some("ptnr1_auth_seq_id"),
            component: "ptnr1_label_comp_id",
            atom: "ptnr1_label_atom_id",
            alternate: "pdbx_ptnr1_label_alt_id",
        }),
        (2, Namespace::Label) => Some(EndpointItems {
            asym: "ptnr2_label_asym_id",
            sequence: "ptnr2_label_seq_id",
            fallback_sequence: Some("ptnr2_auth_seq_id"),
            component: "ptnr2_label_comp_id",
            atom: "ptnr2_label_atom_id",
            alternate: "pdbx_ptnr2_label_alt_id",
        }),
        (1, Namespace::Author) => Some(EndpointItems {
            asym: "ptnr1_auth_asym_id",
            sequence: "ptnr1_auth_seq_id",
            fallback_sequence: None,
            component: "ptnr1_auth_comp_id",
            atom: "ptnr1_auth_atom_id",
            alternate: "pdbx_ptnr1_label_alt_id",
        }),
        (2, Namespace::Author) => Some(EndpointItems {
            asym: "ptnr2_auth_asym_id",
            sequence: "ptnr2_auth_seq_id",
            fallback_sequence: None,
            component: "ptnr2_auth_comp_id",
            atom: "ptnr2_auth_atom_id",
            alternate: "pdbx_ptnr2_label_alt_id",
        }),
        _ => None,
    }
}

fn requested_maps(connections: &[Connection]) -> (AtomMap, AtomMap) {
    let capacity = connections.len().saturating_mul(2);
    let mut label = HashMap::with_capacity(capacity);
    let mut auth = HashMap::with_capacity(capacity);
    for connection in connections {
        request(&mut label, connection.atom_a.label);
        request(&mut label, connection.atom_b.label);
        request(&mut auth, connection.atom_a.auth);
        request(&mut auth, connection.atom_b.auth);
    }
    (label, auth)
}

fn request(map: &mut AtomMap, key: Option<AtomKey>) {
    if let Some(key) = key {
        map.entry(key).or_default();
    }
}

fn match_atoms(data: &StructureData, label: &mut AtomMap, auth: &mut AtomMap) {
    for chain in data.chains() {
        for residue in chain.residues() {
            for atom in residue.atoms() {
                observe(label, label_key(chain, residue, atom), atom.index());
                observe(auth, auth_key(chain, residue, atom), atom.index());
            }
        }
    }
}

fn observe(map: &mut AtomMap, key: Option<AtomKey>, atom: AtomIndex) {
    if let Some(key) = key
        && let Some(found) = map.get_mut(&key)
    {
        found.observe(atom);
    }
}

fn label_key(chain: ChainRef<'_>, residue: ResidueRef<'_>, atom: AtomRef<'_>) -> Option<AtomKey> {
    Some(AtomKey {
        asym: chain.label_asym_id()?,
        sequence: residue.label_seq_id().or_else(|| residue.auth_seq_id()),
        component: atom.component_id()?,
        atom: atom.name_symbol()?,
        alternate: atom.alt_id()?.symbol(),
    })
}

fn auth_key(chain: ChainRef<'_>, residue: ResidueRef<'_>, atom: AtomRef<'_>) -> Option<AtomKey> {
    Some(AtomKey {
        asym: chain.auth_asym_id()?,
        sequence: residue.auth_seq_id(),
        component: match residue.auth_comp_id() {
            Some(component) => component,
            None => atom.component_id()?,
        },
        atom: match atom.auth_name_symbol() {
            Some(name) => name,
            None => atom.name_symbol()?,
        },
        alternate: atom.alt_id()?.symbol(),
    })
}

fn is_connectivity(kind: Option<&str>) -> bool {
    match kind {
        Some(kind) => {
            starts_with_ignore_ascii_case(kind, "covale")
                || starts_with_ignore_ascii_case(kind, "disulf")
                || starts_with_ignore_ascii_case(kind, "metalc")
                || starts_with_ignore_ascii_case(kind, "modres")
        }
        None => true,
    }
}

fn starts_with_ignore_ascii_case(value: &str, prefix: &str) -> bool {
    value
        .get(..prefix.len())
        .is_some_and(|start| start.eq_ignore_ascii_case(prefix))
}

fn order(value: Option<&str>) -> BondOrder {
    let Some(value) = value else {
        return BondOrder::Unknown;
    };
    if value.eq_ignore_ascii_case("SING") || value.eq_ignore_ascii_case("SINGLE") {
        BondOrder::Single
    } else if value.eq_ignore_ascii_case("DOUB") || value.eq_ignore_ascii_case("DOUBLE") {
        BondOrder::Double
    } else if value.eq_ignore_ascii_case("TRIP") || value.eq_ignore_ascii_case("TRIPLE") {
        BondOrder::Triple
    } else if value.eq_ignore_ascii_case("QUAD") || value.eq_ignore_ascii_case("QUADRUPLE") {
        BondOrder::Quadruple
    } else if value.eq_ignore_ascii_case("AROM") || value.eq_ignore_ascii_case("AROMATIC") {
        BondOrder::Aromatic
    } else {
        BondOrder::Unknown
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
