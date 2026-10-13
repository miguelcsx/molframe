//! NCS ensembles, domains and the operators pairing them.
//!
//! The operator table says how coordinates relate; these records say which
//! residues, in which chains, the operators relate. Together they bind an
//! operator to the molecules it acts on.

use crate::category_transform::malformed;
use molframe_cif::{Category, Document};
use molframe_core::{Code, Diagnostic};

/// One end of a domain range, as the residue it names.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct NcsResidue {
    /// Chain identifier in the label namespace.
    pub asym_id: Box<str>,
    /// Residue name.
    pub comp_id: Option<Box<str>>,
    /// Residue number in the label namespace.
    pub seq_id: Option<i32>,
    /// Chain identifier in the author namespace.
    pub auth_asym_id: Option<Box<str>>,
    /// Residue number in the author namespace.
    pub auth_seq_id: Option<i32>,
}

/// One contiguous residue range of a domain.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NcsDomainLimit {
    /// Component number within the domain, when the source numbers them.
    pub component: Option<Box<str>>,
    /// First residue of the range.
    pub begin: NcsResidue,
    /// Last residue of the range.
    pub end: NcsResidue,
}

/// A set of residue ranges that one operator maps onto another domain.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NcsDomain {
    /// Source domain identifier.
    pub id: Box<str>,
    /// Ensemble the domain belongs to.
    pub ensemble: Box<str>,
    /// Optional source description.
    pub details: Option<Box<str>>,
    /// Residue ranges in source order.
    pub limits: Vec<NcsDomainLimit>,
}

/// A group of domains related by NCS operators.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NcsEnsemble {
    /// Source ensemble identifier.
    pub id: Box<str>,
    /// Optional source description.
    pub details: Option<Box<str>>,
    /// Optional point-group name.
    pub point_group: Option<Box<str>>,
}

/// An operator declared to carry one domain onto another.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NcsPairing {
    /// Ensemble the pairing belongs to.
    pub ensemble: Box<str>,
    /// Domain the operator is applied to.
    pub domain_1: Box<str>,
    /// Domain the operator produces.
    pub domain_2: Box<str>,
    /// Operator identifier.
    pub operator: Box<str>,
}

/// Domain-level NCS records of one data block.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct NcsDomains {
    /// Ensembles in source order.
    pub ensembles: Vec<NcsEnsemble>,
    /// Domains in source order, with their limits attached.
    pub domains: Vec<NcsDomain>,
    /// Operator pairings in source order.
    pub pairings: Vec<NcsPairing>,
}

impl NcsDomains {
    /// Whether the block recorded no domain-level NCS information.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.ensembles.is_empty() && self.domains.is_empty() && self.pairings.is_empty()
    }

    /// One domain by identifier.
    #[must_use]
    pub fn domain(&self, id: &str) -> Option<&NcsDomain> {
        self.domains.iter().find(|domain| &*domain.id == id)
    }
}

pub(crate) fn lower_domains(document: &Document) -> Result<NcsDomains, Vec<Diagnostic>> {
    let Some(block) = document.first_block() else {
        return Ok(NcsDomains::default());
    };
    let mut findings = Vec::new();
    let mut out = NcsDomains::default();
    if let Some(category) = block.category("struct_ncs_ens") {
        for row in 0..category.row_count() {
            match text(category, "id", row) {
                Some(id) => out.ensembles.push(NcsEnsemble {
                    id,
                    details: text(category, "details", row),
                    point_group: text(category, "point_group", row),
                }),
                None => findings.push(malformed(category, row, "id", Code::E6014)),
            }
        }
    }
    if let Some(category) = block.category("struct_ncs_dom") {
        for row in 0..category.row_count() {
            match text(category, "id", row) {
                Some(id) => out.domains.push(NcsDomain {
                    id,
                    ensemble: text(category, "pdbx_ens_id", row).unwrap_or_default(),
                    details: text(category, "details", row),
                    limits: Vec::new(),
                }),
                None => findings.push(malformed(category, row, "id", Code::E6014)),
            }
        }
    }
    if let Some(category) = block.category("struct_ncs_dom_lim") {
        lower_limits(category, &mut out, &mut findings);
    }
    if let Some(category) = block.category("struct_ncs_ens_gen") {
        for row in 0..category.row_count() {
            let fields = (
                text(category, "dom_id_1", row),
                text(category, "dom_id_2", row),
                text(category, "oper_id", row),
            );
            let (Some(domain_1), Some(domain_2), Some(operator)) = fields else {
                findings.push(malformed(category, row, "oper_id", Code::E6014));
                continue;
            };
            out.pairings.push(NcsPairing {
                ensemble: text(category, "ens_id", row).unwrap_or_default(),
                domain_1,
                domain_2,
                operator,
            });
        }
    }
    if findings.is_empty() {
        Ok(out)
    } else {
        Err(findings)
    }
}

fn lower_limits(category: &Category, out: &mut NcsDomains, findings: &mut Vec<Diagnostic>) {
    for row in 0..category.row_count() {
        let (Some(dom_id), Some(begin), Some(end)) = (
            text(category, "dom_id", row),
            residue(category, "beg", row),
            residue(category, "end", row),
        ) else {
            findings.push(malformed(category, row, "dom_id", Code::E6014));
            continue;
        };
        let limit = NcsDomainLimit {
            component: text(category, "pdbx_component_id", row),
            begin,
            end,
        };
        let ensemble = text(category, "pdbx_ens_id", row);
        let position = out.domains.iter().position(|domain| {
            domain.id == dom_id && ensemble.as_ref().is_none_or(|e| *e == domain.ensemble)
        });
        match position {
            Some(position) => out.domains[position].limits.push(limit),
            None => findings.push(malformed(category, row, "dom_id", Code::E6014)),
        }
    }
}

fn residue(category: &Category, side: &str, row: usize) -> Option<NcsResidue> {
    Some(NcsResidue {
        asym_id: text(category, &format!("{side}_label_asym_id"), row)?,
        comp_id: text(category, &format!("{side}_label_comp_id"), row),
        seq_id: integer(category, &format!("{side}_label_seq_id"), row),
        auth_asym_id: text(category, &format!("{side}_auth_asym_id"), row),
        auth_seq_id: integer(category, &format!("{side}_auth_seq_id"), row),
    })
}

fn text(category: &Category, item: &str, row: usize) -> Option<Box<str>> {
    category
        .identifier(item, row)
        .filter(|value| !value.trim().is_empty())
        .map(|value| value.into_owned().into_boxed_str())
}

fn integer(category: &Category, item: &str, row: usize) -> Option<i32> {
    category
        .value(item, row)
        .and_then(molframe_cif::CifValue::as_integer)
        .and_then(|value| i32::try_from(value).ok())
}

#[cfg(test)]
#[path = "ncs_domain_tests.rs"]
mod tests;
