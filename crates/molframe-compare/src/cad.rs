//! Contact Area Difference (CAD) scoring over mapped residue contact areas.

use molframe_core::ExecutionContext;
use std::collections::BTreeMap;

/// Area assigned to one unordered residue contact.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactArea {
    /// First residue identifier in the comparison mapping.
    pub first: u32,
    /// Second residue identifier in the comparison mapping.
    pub second: u32,
    /// Contact area in square ångström.
    pub area: f64,
}

/// Per-contact contribution to a CAD score.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CadContact {
    /// First mapped residue.
    pub first: u32,
    /// Second mapped residue.
    pub second: u32,
    /// Reference contact area.
    pub reference_area: f64,
    /// Model contact area, or zero when absent.
    pub model_area: f64,
    /// Absolute area difference capped by the reference area.
    pub lost_area: f64,
}

/// Residue-local CAD score and its normalization terms.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LocalCad {
    /// Mapped residue identifier.
    pub residue: u32,
    /// Sum of its reference contact areas.
    pub reference_area: f64,
    /// Sum of capped absolute differences.
    pub lost_area: f64,
    /// `1 - lost_area / reference_area`.
    pub score: f64,
}

/// Global and per-item CAD result.
#[derive(Clone, Debug, PartialEq)]
pub struct CadScore {
    /// Global CAD score in `[0, 1]`.
    pub score: f64,
    /// Reference area used as the global denominator.
    pub reference_area: f64,
    /// Capped area difference used as the global numerator.
    pub lost_area: f64,
    /// Every reference contact and its contribution, sorted by residue pair.
    pub contacts: Vec<CadContact>,
    /// Per-residue scores sorted by residue identifier.
    pub local: Vec<LocalCad>,
}

/// Invalid contact-area input.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum CadError {
    /// An area is negative or non-finite, or a contact names one residue twice.
    #[error("CAD contact areas must be finite, non-negative and join two residues")]
    InvalidContact,
}

/// Invalid surface-to-residue CAD construction input.
#[derive(Debug, thiserror::Error)]
pub enum CadConstructionError {
    /// Atom arrays do not share one length.
    #[error("positions, radii and residue identifiers must have equal length")]
    LengthMismatch,
    /// Surface construction refused geometry or sampling controls.
    #[error(transparent)]
    Surface(#[from] molframe_surface::SasaError),
}

/// Constructs residue contact areas from solvent-excluded atom patches.
///
/// Atom contacts within one residue are omitted. All other atom-pair patches
/// are coalesced into stable unordered residue pairs. Chemistry remains explicit
/// in the caller-supplied radii, probe and surface density.
///
/// # Errors
///
/// Returns length or surface-construction errors.
pub fn cad_contact_areas(
    positions: &[[f32; 3]],
    radii: &[f32],
    residues: &[u32],
    probe: f32,
    density: f32,
    context: &ExecutionContext,
) -> Result<Vec<ContactArea>, CadConstructionError> {
    if positions.len() != residues.len() || positions.len() != radii.len() {
        return Err(CadConstructionError::LengthMismatch);
    }
    let mut areas = BTreeMap::new();
    for contact in molframe_surface::atom_contact_areas(positions, radii, probe, density, context)?
    {
        let first = residues[contact.first];
        let second = residues[contact.second];
        if first == second || contact.area == 0.0 {
            continue;
        }
        let pair = if first < second {
            (first, second)
        } else {
            (second, first)
        };
        *areas.entry(pair).or_insert(0.0) += contact.area;
    }
    Ok(areas
        .into_iter()
        .map(|((first, second), area)| ContactArea {
            first,
            second,
            area,
        })
        .collect())
}

/// Computes global and local CAD scores from corresponding contact-area maps.
///
/// Duplicate or reversed pairs are coalesced by summing their areas. For each
/// reference contact, the absolute model difference is capped at the reference
/// area, which keeps every contribution and the final score in `[0, 1]`.
/// Model-only pairs do not create a new denominator term; their geometric effect
/// belongs in the contact-area construction that precedes this metric.
///
/// # Errors
///
/// Refuses self-contacts and negative or non-finite areas.
pub fn cad_score(reference: &[ContactArea], model: &[ContactArea]) -> Result<CadScore, CadError> {
    let reference = normalize(reference)?;
    let model = normalize(model)?;
    let mut contacts = Vec::with_capacity(reference.len());
    let mut local: BTreeMap<u32, (f64, f64)> = BTreeMap::new();
    let mut reference_area = 0.0;
    let mut lost_area = 0.0;
    for ((first, second), target_area) in reference {
        let model_area = area_or_zero(model.get(&(first, second)).copied());
        let difference = (target_area - model_area).abs().min(target_area);
        reference_area += target_area;
        lost_area += difference;
        contacts.push(CadContact {
            first,
            second,
            reference_area: target_area,
            model_area,
            lost_area: difference,
        });
        accumulate_local(&mut local, first, target_area, difference);
        accumulate_local(&mut local, second, target_area, difference);
    }
    let local = local
        .into_iter()
        .map(|(residue, (area, lost))| LocalCad {
            residue,
            reference_area: area,
            lost_area: lost,
            score: normalized(area, lost),
        })
        .collect();
    Ok(CadScore {
        score: normalized(reference_area, lost_area),
        reference_area,
        lost_area,
        contacts,
        local,
    })
}

fn area_or_zero(area: Option<f64>) -> f64 {
    let Some(area) = area else {
        return 0.0;
    };
    area
}

fn normalize(areas: &[ContactArea]) -> Result<BTreeMap<(u32, u32), f64>, CadError> {
    let mut normalized = BTreeMap::new();
    for contact in areas {
        if contact.first == contact.second || !contact.area.is_finite() || contact.area < 0.0 {
            return Err(CadError::InvalidContact);
        }
        let pair = if contact.first < contact.second {
            (contact.first, contact.second)
        } else {
            (contact.second, contact.first)
        };
        *normalized.entry(pair).or_insert(0.0) += contact.area;
    }
    Ok(normalized)
}

fn accumulate_local(local: &mut BTreeMap<u32, (f64, f64)>, residue: u32, area: f64, lost: f64) {
    let totals = local.entry(residue).or_insert((0.0, 0.0));
    totals.0 += area;
    totals.1 += lost;
}

fn normalized(reference: f64, lost: f64) -> f64 {
    if reference <= 0.0 {
        1.0
    } else {
        (1.0 - lost / reference).clamp(0.0, 1.0)
    }
}

#[cfg(test)]
#[path = "cad_tests.rs"]
mod tests;
