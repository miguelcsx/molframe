//! Governed streamed SASA: radii from a topology under the policy's rule.
//!
//! [`SasaStreamPlan`] is the part of a streamed surface analysis that needs the
//! topology: which atoms are measured and with which radii. The trajectory is then
//! read through exactly those atoms, so a frame never carries an atom that has no
//! radius, and the totals equal the per-atom areas of the batch calculation summed
//! over the same atoms.

use crate::sasa_atoms::measured_radii;
use crate::stream_surface::{SasaStreamError, sasa_stream};
use molframe_chem::RadiusSelection;
use molframe_core::contract::{Analysis, AnalysisPolicy, Coverage, PolicyField};
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::{ExecutionContext, Structure, diagnostic_from};
use molframe_surface::{SasaMethod, SasaOptions};
use molframe_traj::{TrajectoryError, TrajectoryReader};

/// Total accessible area per streamed frame, in square ångström.
#[derive(Clone, Debug, PartialEq)]
pub struct SasaSeries {
    /// Source frame ordinal of every value.
    pub frames: Vec<u64>,
    /// Simulation time in picoseconds; `NaN` where unrecorded.
    pub times: Vec<f64>,
    /// Total solvent-accessible area of every frame.
    pub areas: Vec<f64>,
}

/// The atoms and radii a streamed SASA analysis measures.
#[derive(Debug)]
pub struct SasaStreamPlan {
    atoms: Vec<usize>,
    radii: Vec<f32>,
    topology_atoms: usize,
    missing: usize,
    probe: f32,
    samples: u16,
    assumption: molframe_core::contract::Assumption,
    table: String,
    policy: AnalysisPolicy,
}

diagnostic_from!(SasaStreamError, |error| match error {
    SasaStreamError::Trajectory(inner) => Diagnostic::from(inner).code(),
    SasaStreamError::Surface(inner) => Diagnostic::from(inner).code(),
    SasaStreamError::Radii(inner) => Diagnostic::from(inner).code(),
    SasaStreamError::Method => Code::E5101,
    SasaStreamError::NoAtoms => Code::E5103,
    SasaStreamError::TopologyAtoms { .. } => Code::E5102,
});

impl SasaStreamPlan {
    /// Chooses the measured atoms of `topology` and their radii under `policy`.
    ///
    /// `options` and `selection` mean what they do for the batch SASA kernel: the table
    /// is checked against the policy (`E7103` when they conflict) and recorded as an
    /// assumption. Hydrogens, other conformers and atoms the heteroatom options leave
    /// out are not measured.
    ///
    /// # Errors
    ///
    /// [`SasaStreamError::Method`] for Lee-Richards, [`SasaStreamError::Radii`] for a
    /// table conflict, [`SasaStreamError::NoAtoms`] when nothing is left to measure.
    pub fn new(
        topology: &Structure,
        options: &SasaOptions,
        selection: &RadiusSelection,
        policy: &AnalysisPolicy,
    ) -> Result<Self, SasaStreamError> {
        let SasaMethod::ShrakeRupley { points } = options.method else {
            return Err(SasaStreamError::Method);
        };
        let assumption = crate::policy_execution::radii_assumption(options, selection, policy)
            .map_err(SasaStreamError::Radii)?;
        let measured = measured_radii(topology, options);
        let atoms: Vec<usize> = (0..measured.radii.len())
            .filter(|&atom| measured.radii[atom].is_finite())
            .collect();
        if atoms.is_empty() {
            return Err(SasaStreamError::NoAtoms);
        }
        let radii = atoms.iter().map(|&atom| measured.radii[atom]).collect();
        Ok(Self {
            atoms,
            radii,
            topology_atoms: measured.radii.len(),
            missing: measured.missing,
            probe: options.probe,
            samples: points,
            assumption,
            table: options.radii.identity().name.to_string(),
            policy: policy.clone(),
        })
    }

    /// Topology atom indices to read from the trajectory, ascending.
    #[must_use]
    pub fn atoms(&self) -> &[usize] {
        &self.atoms
    }

    /// Streams the total area of every frame of `reader`, which must read exactly
    /// [`Self::atoms`] of a trajectory with `trajectory_atoms` atoms.
    ///
    /// # Errors
    ///
    /// [`SasaStreamError::TopologyAtoms`] when the trajectory and the topology differ
    /// in size; otherwise as [`sasa_stream`].
    pub fn run<R: TrajectoryReader + ?Sized>(
        &self,
        reader: &mut R,
        trajectory_atoms: usize,
        context: &ExecutionContext,
        frame_workspace_bytes: usize,
    ) -> Result<Analysis<SasaSeries>, SasaStreamError> {
        if trajectory_atoms != self.topology_atoms {
            return Err(SasaStreamError::TopologyAtoms {
                topology: self.topology_atoms,
                trajectory: trajectory_atoms,
            });
        }
        let mut series = SasaSeries {
            frames: Vec::new(),
            times: Vec::new(),
            areas: Vec::new(),
        };
        sasa_stream(
            reader,
            &self.radii,
            self.probe,
            self.samples,
            context,
            frame_workspace_bytes,
            |frame, time, area| {
                let frame = u64::try_from(frame).map_err(|_| TrajectoryError::IdentityOverflow)?;
                series.frames.push(frame);
                series.times.push(match time {
                    Some(time) => time,
                    None => f64::NAN,
                });
                series.areas.push(area);
                Ok(())
            },
        )?;
        Ok(self.envelope(series))
    }

    fn envelope(&self, series: SasaSeries) -> Analysis<SasaSeries> {
        let used = count(self.atoms.len());
        let analysis = if self.missing == 0 {
            Analysis::complete(series, Coverage::complete(used), &self.policy)
        } else {
            let coverage = Coverage {
                intended: count(self.atoms.len() + self.missing),
                used,
                missing: count(self.missing),
                ambiguous: 0,
            };
            Analysis::partial(series, coverage, &self.policy).with_warning(
                Diagnostic::new(Code::W6004)
                    .about_field(PolicyField::VdwRadii.name())
                    .with_message(format!(
                        "the {} radius set has no radius for {} included atoms; they are not measured and do not occlude (choose a set that covers them, such as alvarez)",
                        self.table, self.missing
                    )),
            )
        };
        analysis.with_assumption(self.assumption.clone())
    }
}

fn count(value: usize) -> u32 {
    match u32::try_from(value) {
        Ok(count) => count,
        Err(_) => u32::MAX,
    }
}

#[cfg(test)]
#[path = "stream_plan_tests.rs"]
mod tests;
