//! Text of the domain-level NCS categories.

use crate::cif_text::optional;
use crate::ncs_domain::{NcsDomains, NcsResidue};
use molframe_cif::quote_text;
use std::fmt::Write as _;

/// Appends `_struct_ncs_ens`, `_struct_ncs_dom`, `_struct_ncs_dom_lim` and
/// `_struct_ncs_ens_gen` for whichever of them hold rows.
pub(crate) fn domains_cif(out: &mut String, domains: &NcsDomains) {
    if !domains.ensembles.is_empty() {
        out.push_str(
            "loop_\n_struct_ncs_ens.id\n_struct_ncs_ens.details\n_struct_ncs_ens.point_group\n",
        );
        for ensemble in &domains.ensembles {
            let _ = writeln!(
                out,
                "{} {} {}",
                quote_text(&ensemble.id),
                optional(ensemble.details.as_deref()),
                optional(ensemble.point_group.as_deref()),
            );
        }
        out.push_str("#\n");
    }
    if !domains.domains.is_empty() {
        out.push_str(
            "loop_\n_struct_ncs_dom.id\n_struct_ncs_dom.pdbx_ens_id\n_struct_ncs_dom.details\n",
        );
        for domain in &domains.domains {
            let _ = writeln!(
                out,
                "{} {} {}",
                quote_text(&domain.id),
                optional(Some(&domain.ensemble)),
                optional(domain.details.as_deref()),
            );
        }
        out.push_str("#\n");
    }
    limits_cif(out, domains);
    if !domains.pairings.is_empty() {
        out.push_str(
            "loop_\n_struct_ncs_ens_gen.ens_id\n_struct_ncs_ens_gen.dom_id_1\n\
             _struct_ncs_ens_gen.dom_id_2\n_struct_ncs_ens_gen.oper_id\n",
        );
        for pairing in &domains.pairings {
            let _ = writeln!(
                out,
                "{} {} {} {}",
                optional(Some(&pairing.ensemble)),
                quote_text(&pairing.domain_1),
                quote_text(&pairing.domain_2),
                quote_text(&pairing.operator),
            );
        }
        out.push_str("#\n");
    }
}

const LIMIT_ITEMS: [&str; 2] = ["beg", "end"];
const RESIDUE_ITEMS: [&str; 5] = [
    "label_asym_id",
    "label_comp_id",
    "label_seq_id",
    "auth_asym_id",
    "auth_seq_id",
];

fn limits_cif(out: &mut String, domains: &NcsDomains) {
    if domains
        .domains
        .iter()
        .all(|domain| domain.limits.is_empty())
    {
        return;
    }
    out.push_str(
        "loop_\n_struct_ncs_dom_lim.pdbx_ens_id\n_struct_ncs_dom_lim.dom_id\n\
         _struct_ncs_dom_lim.pdbx_component_id\n",
    );
    for side in LIMIT_ITEMS {
        for item in RESIDUE_ITEMS {
            let _ = writeln!(out, "_struct_ncs_dom_lim.{side}_{item}");
        }
    }
    for domain in &domains.domains {
        for limit in &domain.limits {
            let _ = write!(
                out,
                "{} {} {} ",
                optional(Some(&domain.ensemble)),
                quote_text(&domain.id),
                optional(limit.component.as_deref()),
            );
            residue_cif(out, &limit.begin);
            out.push(' ');
            residue_cif(out, &limit.end);
            out.push('\n');
        }
    }
    out.push_str("#\n");
}

fn residue_cif(out: &mut String, residue: &NcsResidue) {
    let number = |value: Option<i32>| value.map_or_else(|| "?".to_owned(), |n| n.to_string());
    let _ = write!(
        out,
        "{} {} {} {} {}",
        quote_text(&residue.asym_id),
        optional(residue.comp_id.as_deref()),
        number(residue.seq_id),
        optional(residue.auth_asym_id.as_deref()),
        number(residue.auth_seq_id),
    );
}
