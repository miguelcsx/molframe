//! Direct canonical connectivity encoding.

use super::column::{IntegerColumnBuilder, TextColumnBuilder};
use crate::container::EncodedCategory;
use pdbiox_cif::{CanonicalProjection, CanonicalValue};
use pdbiox_core::BondOrder;
use pdbiox_core::diagnostic::Diagnostic;
use pdbiox_core::structure::{AtomRef, ResidueRef, Structure};

pub(super) fn encode(
    projection: CanonicalProjection<'_>,
) -> Result<Option<EncodedCategory>, Diagnostic> {
    let structure = projection.structure();
    let bonds = &structure.data().bonds;
    if bonds.is_empty() {
        return Ok(None);
    }
    let Some(connection_type) = projection.connection_type_id() else {
        return Err(super::structure::projection_diagnostic(
            &pdbiox_cif::CifWriteError::MissingConnectionTypeId,
        ));
    };
    let rows = bonds.len();
    let mut ids = IntegerColumnBuilder::new("id", rows);
    let mut types = TextColumnBuilder::new("conn_type_id", rows);
    let mut a = EndpointColumns::new(1, rows);
    let mut b = EndpointColumns::new(2, rows);
    let mut orders = TextColumnBuilder::new("pdbx_value_order", rows);
    let mut details = TextColumnBuilder::new("details", rows);
    for (position, bond) in bonds.iter().enumerate() {
        let number =
            i64::try_from(position + 1).map_or(CanonicalValue::Unknown, CanonicalValue::Present);
        ids.push(number);
        types.push(CanonicalValue::Present(connection_type));
        let endpoint_a = endpoint(structure, bond.atom_a.get(), position, 1)
            .map_err(|error| super::structure::projection_diagnostic(&error))?;
        let endpoint_b = endpoint(structure, bond.atom_b.get(), position, 2)
            .map_err(|error| super::structure::projection_diagnostic(&error))?;
        a.push(&endpoint_a);
        b.push(&endpoint_b);
        orders.push(order(bond.order));
        details.push(CanonicalValue::Unknown);
    }
    let mut columns = vec![ids.finish()?, types.finish()?];
    columns.extend(a.finish()?);
    columns.extend(b.finish()?);
    columns.push(orders.finish()?);
    columns.push(details.finish()?);
    Ok(Some(EncodedCategory {
        name: "_struct_conn".to_owned(),
        row_count: rows,
        columns,
    }))
}

struct EndpointColumns {
    chain: TextColumnBuilder,
    sequence: IntegerColumnBuilder,
    component: TextColumnBuilder,
    atom: TextColumnBuilder,
    alternate: TextColumnBuilder,
}

impl EndpointColumns {
    fn new(side: u8, rows: usize) -> Self {
        let names = if side == 1 {
            (
                "ptnr1_label_asym_id",
                "ptnr1_label_seq_id",
                "ptnr1_label_comp_id",
                "ptnr1_label_atom_id",
                "pdbx_ptnr1_label_alt_id",
            )
        } else {
            (
                "ptnr2_label_asym_id",
                "ptnr2_label_seq_id",
                "ptnr2_label_comp_id",
                "ptnr2_label_atom_id",
                "pdbx_ptnr2_label_alt_id",
            )
        };
        Self {
            chain: TextColumnBuilder::new(names.0, rows),
            sequence: IntegerColumnBuilder::new(names.1, rows),
            component: TextColumnBuilder::new(names.2, rows),
            atom: TextColumnBuilder::new(names.3, rows),
            alternate: TextColumnBuilder::new(names.4, rows),
        }
    }

    fn push(&mut self, endpoint: &Endpoint<'_>) {
        self.chain.push(CanonicalValue::Present(endpoint.chain));
        self.sequence.push(match endpoint.residue.label_seq_id() {
            Some(value) => CanonicalValue::Present(i64::from(value)),
            None => CanonicalValue::Inapplicable,
        });
        self.component
            .push(CanonicalValue::Present(endpoint.component));
        self.atom.push(CanonicalValue::Present(endpoint.atom));
        self.alternate.push(match endpoint.source.alt_label() {
            Some(value) if !value.is_empty() => CanonicalValue::Present(value),
            Some(_) | None => CanonicalValue::Inapplicable,
        });
    }

    fn finish(self) -> Result<Vec<crate::container::EncodedColumn>, Diagnostic> {
        Ok(vec![
            self.chain.finish()?,
            self.sequence.finish()?,
            self.component.finish()?,
            self.atom.finish()?,
            self.alternate.finish()?,
        ])
    }
}

struct Endpoint<'a> {
    source: AtomRef<'a>,
    residue: ResidueRef<'a>,
    chain: &'a str,
    component: &'a str,
    atom: &'a str,
}

fn endpoint(
    structure: &Structure,
    index: u32,
    bond: usize,
    side: u8,
) -> Result<Endpoint<'_>, pdbiox_cif::CifWriteError> {
    let Some(source) = structure.atom(pdbiox_core::AtomIndex::new(index)) else {
        return Err(pdbiox_cif::CifWriteError::InvalidBondAtom { bond, atom: index });
    };
    let Some(residue) = source.residue() else {
        return Err(missing(bond, side, "residue"));
    };
    let chain = {
        structure
            .data()
            .topology
            .chains
            .containing(residue.index().get())
            .and_then(|index| structure.chain(index))
    };
    let Some(chain) = chain
        .and_then(pdbiox_core::structure::ChainRef::label)
        .filter(|value| !value.is_empty())
    else {
        return Err(missing(bond, side, "label_asym_id"));
    };
    let Some(component) = source.component_name().filter(|value| !value.is_empty()) else {
        return Err(missing(bond, side, "label_comp_id"));
    };
    let Some(atom) = source.name().filter(|value| !value.is_empty()) else {
        return Err(missing(bond, side, "label_atom_id"));
    };
    Ok(Endpoint {
        source,
        residue,
        chain,
        component,
        atom,
    })
}

const fn missing(bond: usize, endpoint: u8, field: &'static str) -> pdbiox_cif::CifWriteError {
    pdbiox_cif::CifWriteError::MissingBondField {
        bond,
        endpoint,
        field,
    }
}

fn order(value: BondOrder) -> CanonicalValue<&'static str> {
    match value {
        BondOrder::Single => CanonicalValue::Present("SING"),
        BondOrder::Double => CanonicalValue::Present("DOUB"),
        BondOrder::Triple => CanonicalValue::Present("TRIP"),
        BondOrder::Quadruple => CanonicalValue::Present("QUAD"),
        BondOrder::Aromatic => CanonicalValue::Present("AROM"),
        BondOrder::Polymeric | BondOrder::Unknown => CanonicalValue::Unknown,
    }
}
