//! Canonical secondary ranges from the authoritative per-residue column.
//! Preflight uses O(residues log residues) time and one O(residues) identity buffer;
//! output visits residues linearly without retaining ranges or allocating row strings.

use super::options::CifWriteError;
use super::value::quoted;
use molframe_core::SecondaryStructure as Ss;
use molframe_core::structure::{ChainRef, ResidueRef, Structure};
use std::fmt::{self, Display, Formatter};

pub(super) fn preflight(structure: &Structure) -> Result<(), CifWriteError> {
    if let Some(models) = structure.ragged_models() {
        if models.iter().any(|model| {
            model
                .secondary_structure()
                .iter()
                .copied()
                .any(representable)
        }) {
            return Err(CifWriteError::ModelSpecificSecondary);
        }
        return Ok(());
    }
    if !structure
        .secondary_structure()
        .iter()
        .copied()
        .any(representable)
    {
        return Ok(());
    }
    let mut identities = Vec::with_capacity(structure.data().topology.residues.len());
    for chain in structure.data().chains() {
        for residue in chain.residues() {
            if let (Some(chain), Some(sequence)) = (chain.label(), residue.label_seq_id()) {
                identities.push((chain, sequence, residue.index().get()));
            }
            let kind = state(structure, residue);
            if !representable(kind) {
                continue;
            }
            for (field, present) in [
                ("label_asym_id", chain.label().is_some()),
                ("label_comp_id", residue.name().is_some()),
                ("label_seq_id", residue.label_seq_id().is_some()),
            ] {
                if !present {
                    return Err(CifWriteError::MissingSecondaryField {
                        residue: residue.index().get(),
                        field,
                    });
                }
            }
        }
    }
    identities.sort_unstable();
    if let Some(pair) = identities
        .windows(2)
        .find(|pair| pair[0].0 == pair[1].0 && pair[0].1 == pair[1].1)
    {
        return Err(CifWriteError::AmbiguousSecondaryIdentity { residue: pair[1].2 });
    }
    Ok(())
}

pub(super) fn write(out: &mut impl fmt::Write, structure: &Structure) {
    let mut conf_id = 0;
    let mut sheet_id = 0;
    for sheet in [false, true] {
        ranges(structure, |chain, first, last, kind| {
            if (kind == Ss::Strand) != sheet {
                return;
            }
            let id = if sheet { &mut sheet_id } else { &mut conf_id };
            if *id == 0 {
                let _ = out.write_str(if sheet {
                    "loop_\n_struct_sheet_range.sheet_id\n_struct_sheet_range.id\n\
                     _struct_sheet_range.beg_label_comp_id\n_struct_sheet_range.beg_label_asym_id\n\
                     _struct_sheet_range.beg_label_seq_id\n_struct_sheet_range.end_label_comp_id\n\
                     _struct_sheet_range.end_label_asym_id\n_struct_sheet_range.end_label_seq_id\n"
                } else {
                    "loop_\n_struct_conf.conf_type_id\n_struct_conf.id\n\
                     _struct_conf.pdbx_PDB_helix_class\n_struct_conf.beg_label_comp_id\n\
                     _struct_conf.beg_label_asym_id\n_struct_conf.beg_label_seq_id\n\
                     _struct_conf.end_label_comp_id\n_struct_conf.end_label_asym_id\n\
                     _struct_conf.end_label_seq_id\n"
                });
            }
            *id += 1;
            if sheet {
                let _ = writeln!(
                    out,
                    "S{id} 1 {} {}",
                    Endpoint(first, chain),
                    Endpoint(last, chain)
                );
            } else {
                let (conf_type, class) = match kind {
                    Ss::AlphaHelix => ("HELX_P", "1"),
                    Ss::ThreeTenHelix => ("HELX_P", "5"),
                    Ss::PiHelix => ("HELX_P", "3"),
                    Ss::PolyProline => ("HELX_P", "10"),
                    Ss::OtherHelix => ("HELX_P", "2"),
                    Ss::Turn => ("TURN_P", "."),
                    _ => return,
                };
                let _ = writeln!(
                    out,
                    "{conf_type} C{id} {class} {} {}",
                    Endpoint(first, chain),
                    Endpoint(last, chain)
                );
            }
        });
        if (sheet && sheet_id > 0) || (!sheet && conf_id > 0) {
            let _ = out.write_str("#\n");
        }
    }
}

const fn representable(kind: Ss) -> bool {
    kind.is_helix() || matches!(kind, Ss::Strand | Ss::Turn)
}

fn state(structure: &Structure, residue: ResidueRef<'_>) -> Ss {
    match structure
        .secondary_structure()
        .get(residue.index().as_usize())
    {
        Some(kind) => *kind,
        None => Ss::Unknown,
    }
}

fn ranges<'a>(
    structure: &'a Structure,
    mut visit: impl FnMut(ChainRef<'a>, ResidueRef<'a>, ResidueRef<'a>, Ss),
) {
    for chain in structure.data().chains() {
        let mut pending: Option<(ResidueRef<'_>, ResidueRef<'_>, Ss)> = None;
        for residue in chain.residues() {
            let kind = state(structure, residue);
            if let Some((first, last, old)) = pending {
                if kind == old {
                    pending = Some((first, residue, kind));
                    continue;
                }
                visit(chain, first, last, old);
            }
            pending = representable(kind).then_some((residue, residue, kind));
        }
        if let Some((first, last, kind)) = pending {
            visit(chain, first, last, kind);
        }
    }
}

struct Endpoint<'a>(ResidueRef<'a>, ChainRef<'a>);

impl Display for Endpoint<'_> {
    fn fmt(&self, out: &mut Formatter<'_>) -> fmt::Result {
        let residue = self.0;
        let (Some(component), Some(chain), Some(sequence)) =
            (residue.name(), self.1.label(), residue.label_seq_id())
        else {
            return Err(fmt::Error);
        };
        write!(out, "{} {} {sequence}", quoted(component), quoted(chain))
    }
}

#[cfg(test)]
#[path = "secondary_tests.rs"]
mod tests;
