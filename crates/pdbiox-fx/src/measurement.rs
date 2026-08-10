use crate::{AlignedMotif, Constraint, Motif};
use pdbiox_core::index::AtomIndex;
use pdbiox_core::structure::Structure;
use std::collections::BTreeMap;

use crate::numeric::usize_to_f64;

/// Typed observable produced by one constraint.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum MeasurementValue {
    /// Distance, angle, dihedral or plane deviation.
    Scalar(f64),
    /// Signed tetrahedral handedness: `-1` or `1`.
    Chirality(i8),
    /// Number of partners satisfying a coordination cutoff.
    Count(usize),
}

impl MeasurementValue {
    /// Numeric representation used by verdict profiles.
    #[must_use]
    pub fn numeric(self) -> f64 {
        match self {
            Self::Scalar(value) => value,
            Self::Chirality(sign) => f64::from(sign),
            Self::Count(count) => usize_to_f64(count),
        }
    }
}

/// Why a constraint could not produce a defensible value.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IndeterminateReason {
    /// At least one mapped atom has no recorded coordinate.
    MissingCoordinate,
    /// The selected points do not define the requested geometry.
    DegenerateGeometry,
}

/// Fully decomposed outcome of one named constraint.
#[derive(Clone, Debug, PartialEq)]
pub struct ConstraintMeasurement {
    /// Constraint name from the motif.
    pub name: Box<str>,
    /// Observed value, absent only for an indeterminate constraint.
    pub value: Option<MeasurementValue>,
    /// Distance from the acceptable region; zero means satisfied.
    pub deviation: Option<f64>,
    /// Intrinsic constraint result, separate from a profile verdict.
    pub satisfied: Option<bool>,
    /// Exact atom alternative selected, in constraint order.
    pub atoms: Vec<AtomIndex>,
    /// Cause when no alternative was measurable.
    pub indeterminate: Option<IndeterminateReason>,
}

/// Constraint results in declaration order.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MeasurementSet {
    /// One result per declared constraint.
    pub constraints: Vec<ConstraintMeasurement>,
}

/// Explicit bounds and numerical controls for constraint measurement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MeasurementOptions {
    /// Maximum equivalent-atom combinations evaluated for one constraint.
    pub maximum_alternatives: usize,
    /// Numerical controls for plane-fitting constraints.
    pub plane_fit: pdbiox_geom::EigenOptions,
}

impl MeasurementSet {
    /// Raw observables keyed by constraint name for profile evaluation.
    #[must_use]
    pub fn metrics(&self) -> BTreeMap<Box<str>, f64> {
        self.constraints
            .iter()
            .filter_map(|result| {
                result
                    .value
                    .map(|value| (result.name.clone(), value.numeric()))
            })
            .collect()
    }
}

/// Failure during bounded atom-equivalence evaluation.
#[derive(Clone, Copy, Debug, thiserror::Error, PartialEq, Eq)]
#[non_exhaustive]
pub enum MeasurementError {
    /// A constraint's atom alternatives exceeded the caller-selected bound.
    #[error("constraint measurement exceeded the limit of {limit} alternatives")]
    LimitExceeded {
        /// Maximum combinations allowed per constraint.
        limit: usize,
    },
    /// A constraint's geometric decomposition failed.
    #[error("constraint geometry failed: {0:?}")]
    Geometry(pdbiox_geom::EigenError),
}

/// Measures every constraint against one aligned mapping.
///
/// Equivalent atoms are enumerated in atom-index order. The measurable
/// alternative with the smallest deviation is retained, with atom indices as a
/// stable tie-break. No score is collapsed across constraints.
///
/// # Errors
///
/// Returns [`MeasurementError::LimitExceeded`] before evaluating more than the
/// configured number of alternatives for any one constraint, or a geometry
/// error when a numerical decomposition fails.
pub fn measure_constraints(
    structure: &Structure,
    motif: &Motif,
    aligned: &AlignedMotif,
    options: MeasurementOptions,
) -> Result<MeasurementSet, MeasurementError> {
    let mut constraints = Vec::with_capacity(motif.constraints().len());
    for named in motif.constraints() {
        let sites = named.constraint.sites();
        let alternatives: Vec<_> = sites
            .iter()
            .map(|site| {
                aligned
                    .mapping
                    .atoms
                    .get(*site)
                    .map_or_else(Vec::new, Clone::clone)
            })
            .collect();
        combination_count(&alternatives, options.maximum_alternatives)?;
        let mut state = ChoiceState::new();
        enumerate_choices(
            &sites,
            &alternatives,
            0,
            &mut Vec::new(),
            &mut state,
            |atoms| {
                observe(
                    structure,
                    &named.constraint,
                    aligned,
                    atoms,
                    options.plane_fit,
                )
            },
        )?;
        constraints.push(state.finish(named.name.clone()));
    }
    Ok(MeasurementSet { constraints })
}

fn combination_count(
    alternatives: &[Vec<AtomIndex>],
    limit: usize,
) -> Result<usize, MeasurementError> {
    let mut count = 1usize;
    for choices in alternatives {
        count = count
            .checked_mul(choices.len())
            .filter(|count| *count <= limit)
            .ok_or(MeasurementError::LimitExceeded { limit })?;
    }
    Ok(count)
}

fn enumerate_choices(
    sites: &[&crate::AtomSite],
    alternatives: &[Vec<AtomIndex>],
    depth: usize,
    chosen: &mut Vec<AtomIndex>,
    state: &mut ChoiceState,
    mut measure: impl FnMut(&[AtomIndex]) -> Result<Observation, MeasurementError> + Copy,
) -> Result<(), MeasurementError> {
    if depth == alternatives.len() {
        state.consider(measure(chosen)?, chosen);
        return Ok(());
    }
    for &atom in &alternatives[depth] {
        if chosen
            .iter()
            .enumerate()
            .any(|(index, &selected)| selected == atom && sites[index] != sites[depth])
        {
            continue;
        }
        chosen.push(atom);
        enumerate_choices(sites, alternatives, depth + 1, chosen, state, measure)?;
        chosen.pop();
    }
    Ok(())
}

#[derive(Clone, Copy)]
struct Observation {
    value: Option<MeasurementValue>,
    deviation: f64,
    satisfied: bool,
    reason: Option<IndeterminateReason>,
}

fn observe(
    structure: &Structure,
    constraint: &Constraint,
    aligned: &AlignedMotif,
    atoms: &[AtomIndex],
    plane_fit: pdbiox_geom::EigenOptions,
) -> Result<Observation, MeasurementError> {
    let positions: Option<Vec<_>> = atoms
        .iter()
        .map(|&atom| {
            structure
                .data()
                .atom(atom)
                .and_then(pdbiox_core::structure::AtomRef::position)
                .map(|position| aligned.transform.apply(position))
        })
        .collect();
    let Some(points) = positions else {
        return Ok(Observation::missing());
    };
    let observation = match constraint {
        Constraint::Distance {
            target, tolerance, ..
        } => scalar(
            pdbiox_geom::distance(points[0], points[1]),
            *target,
            *tolerance,
            false,
        ),
        Constraint::Angle {
            target, tolerance, ..
        } => pdbiox_geom::angle(points[0], points[1], points[2])
            .map(pdbiox_geom::degrees)
            .map_or_else(Observation::degenerate, |value| {
                scalar(value, *target, *tolerance, false)
            }),
        Constraint::Dihedral {
            target, tolerance, ..
        } => pdbiox_geom::dihedral(points[0], points[1], points[2], points[3])
            .map(pdbiox_geom::degrees)
            .map_or_else(Observation::degenerate, |value| {
                scalar(value, *target, *tolerance, true)
            }),
        Constraint::Chirality { positive, .. } => chirality(&points, *positive),
        Constraint::Planarity { tolerance, .. } => {
            pdbiox_geom::plane_deviation_with_options(&points, plane_fit)
                .map_err(MeasurementError::Geometry)?
                .map_or_else(Observation::degenerate, |value| {
                    upper_bound(value, *tolerance)
                })
        }
        Constraint::Coordination {
            count,
            max_distance,
            ..
        } => {
            let observed = points[1..]
                .iter()
                .filter(|&&partner| pdbiox_geom::distance(points[0], partner) <= *max_distance)
                .count();
            let deviation = usize_to_f64(observed.abs_diff(*count));
            Observation {
                value: Some(MeasurementValue::Count(observed)),
                deviation,
                satisfied: observed == *count,
                reason: None,
            }
        }
        Constraint::StericExclusion { min_distance, .. } => {
            let value = pdbiox_geom::distance(points[0], points[1]);
            Observation {
                value: Some(MeasurementValue::Scalar(value)),
                deviation: (*min_distance - value).max(0.0),
                satisfied: value >= *min_distance,
                reason: None,
            }
        }
    };
    Ok(observation)
}

impl Observation {
    fn missing() -> Self {
        Self {
            value: None,
            deviation: f64::INFINITY,
            satisfied: false,
            reason: Some(IndeterminateReason::MissingCoordinate),
        }
    }

    fn degenerate() -> Self {
        Self {
            value: None,
            deviation: f64::INFINITY,
            satisfied: false,
            reason: Some(IndeterminateReason::DegenerateGeometry),
        }
    }
}

fn scalar(value: f64, target: f64, tolerance: f64, circular: bool) -> Observation {
    let difference = if circular {
        ((value - target + 180.0).rem_euclid(360.0) - 180.0).abs()
    } else {
        (value - target).abs()
    };
    Observation {
        value: Some(MeasurementValue::Scalar(value)),
        deviation: (difference - tolerance).max(0.0),
        satisfied: difference <= tolerance,
        reason: None,
    }
}

fn upper_bound(value: f64, limit: f64) -> Observation {
    Observation {
        value: Some(MeasurementValue::Scalar(value)),
        deviation: (value - limit).max(0.0),
        satisfied: value <= limit,
        reason: None,
    }
}

fn chirality(points: &[[f32; 3]], positive: bool) -> Observation {
    let first = pdbiox_geom::displacement(points[0], points[1]);
    let second = pdbiox_geom::displacement(points[0], points[2]);
    let third = pdbiox_geom::displacement(points[0], points[3]);
    let triple = pdbiox_geom::dot(pdbiox_geom::cross(first, second), third);
    if triple.abs() <= f64::EPSILON {
        return Observation::degenerate();
    }
    let sign = if triple.is_sign_positive() { 1 } else { -1 };
    let satisfied = (sign > 0) == positive;
    Observation {
        value: Some(MeasurementValue::Chirality(sign)),
        deviation: if satisfied { 0.0 } else { 2.0 },
        satisfied,
        reason: None,
    }
}

struct ChoiceState {
    best: Option<(Observation, Vec<AtomIndex>)>,
    reason: Option<IndeterminateReason>,
}

impl ChoiceState {
    fn new() -> Self {
        Self {
            best: None,
            reason: None,
        }
    }

    fn consider(&mut self, observation: Observation, atoms: &[AtomIndex]) {
        if !observation.deviation.is_finite() {
            if self.reason != Some(IndeterminateReason::MissingCoordinate) {
                self.reason = observation.reason;
            }
            return;
        }
        let replace = self.best.as_ref().is_none_or(|(best, best_atoms)| {
            observation
                .deviation
                .total_cmp(&best.deviation)
                .then_with(|| atoms.cmp(best_atoms))
                .is_lt()
        });
        if replace {
            self.best = Some((observation, atoms.to_vec()));
        }
    }

    fn finish(self, name: Box<str>) -> ConstraintMeasurement {
        if let Some((observation, atoms)) = self.best {
            return ConstraintMeasurement {
                name,
                value: observation.value,
                deviation: Some(observation.deviation),
                satisfied: Some(observation.satisfied),
                atoms,
                indeterminate: None,
            };
        }
        let reason = match self.reason {
            Some(reason) => reason,
            None => IndeterminateReason::DegenerateGeometry,
        };
        ConstraintMeasurement {
            name,
            value: None,
            deviation: None,
            satisfied: None,
            atoms: Vec::new(),
            indeterminate: Some(reason),
        }
    }
}

#[cfg(test)]
#[path = "measurement_tests.rs"]
mod tests;
