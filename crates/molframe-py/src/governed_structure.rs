//! Governed analyses of one structure: contact maps, interaction geometry,
//! interfaces, exposure and nucleic-acid torsions.
//!
//! Every scientific threshold is a required keyword: there is no universally
//! right cutoff, and the result records the parameters it was computed with.

use crate::analysis_result::PyAnalysis;
use crate::bindings::PyStructure;
use crate::execution::PyExecutionContext;
use crate::governed::run;
use crate::policy::{PyAnalysisPolicy, policy_of};
use crate::table::TableBuilder;
use molframe::analysis::{
    CationPiOptions, HydrogenBondOptions, PiStackingOptions, WaterBridgeOptions,
};
use molframe::geometry::EigenOptions;
use pyo3::prelude::*;
use pyo3::types::PyDict;

fn residues(column: &[molframe::ResidueIndex]) -> Vec<u32> {
    column.iter().map(|index| index.get()).collect()
}

fn atoms(column: &[molframe::AtomIndex]) -> Vec<u32> {
    column.iter().map(|index| index.get()).collect()
}

/// Residue pairs within `cutoff` ångström of each other, with their closest
/// atom distance.
///
/// Pairs closer than `min_separation` residues in the structure's own order are
/// dropped; cross-chain pairs are always kept.
#[pyfunction]
#[pyo3(signature = (structure, *, cutoff, min_separation, backend="auto", policy=None, context=None))]
pub(crate) fn contact_map(
    py: Python<'_>,
    structure: &PyStructure,
    cutoff: f32,
    min_separation: u32,
    backend: &str,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let kernel = molframe::analysis::contact_map_kernel(
        cutoff,
        min_separation,
        crate::backend::parse(backend)?,
    );
    run(
        py,
        structure,
        &policy_of(policy),
        &kernel,
        |py, map| {
            let contacts = map.contacts();
            let table = TableBuilder::new(py, contacts.len())
                .indices("first_residue", &residues(contacts.first()))
                .indices("second_residue", &residues(contacts.second()))
                .single("min_distance", contacts.min_distance())
                .finish();
            Ok(Py::new(py, table)?.into_any())
        },
        context,
    )
}

/// Ring-ring stacking between aromatic residues of CCD-annotated chemistry.
///
/// `kind` is `0` for parallel planes and `1` for T-shaped ones. Plane fits use
/// the named standard small-matrix eigensolver profile.
#[pyfunction]
#[pyo3(signature = (
    structure,
    *,
    max_centre_distance,
    max_parallel_angle,
    min_perpendicular_angle,
    policy=None,
    context=None,
))]
pub(crate) fn pi_stacking(
    py: Python<'_>,
    structure: &PyStructure,
    max_centre_distance: f32,
    max_parallel_angle: f64,
    min_perpendicular_angle: f64,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let kernel = molframe::analysis::pi_stacking_kernel(PiStackingOptions {
        maximum_centre_distance: max_centre_distance,
        maximum_parallel_angle: max_parallel_angle,
        minimum_perpendicular_angle: min_perpendicular_angle,
        plane_fit: EigenOptions::standard(),
    });
    run(
        py,
        structure,
        &policy_of(policy),
        &kernel,
        |py, table| {
            let kinds: Vec<u32> = table
                .kind()
                .iter()
                .map(|kind| match kind {
                    molframe::analysis::StackingKind::Parallel => 0,
                    molframe::analysis::StackingKind::TShaped => 1,
                })
                .collect();
            let table = TableBuilder::new(py, table.len())
                .indices("first_residue", &residues(table.first()))
                .indices("second_residue", &residues(table.second()))
                .single("centre_distance", table.centre_distance())
                .double("angle_degrees", table.angle())
                .indices("kind", &kinds)
                .finish();
            Ok(Py::new(py, table)?.into_any())
        },
        context,
    )
}

/// Cations over the face of an aromatic ring, from CCD formal charges.
#[pyfunction]
#[pyo3(signature = (structure, *, max_distance, max_face_angle, policy=None, context=None))]
pub(crate) fn cation_pi(
    py: Python<'_>,
    structure: &PyStructure,
    max_distance: f32,
    max_face_angle: f64,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let kernel = molframe::analysis::cation_pi_kernel(CationPiOptions {
        maximum_distance: max_distance,
        maximum_face_angle: max_face_angle,
        plane_fit: EigenOptions::standard(),
    });
    run(
        py,
        structure,
        &policy_of(policy),
        &kernel,
        |py, table| {
            let table = TableBuilder::new(py, table.len())
                .indices("cation_residue", &residues(table.cation_residue()))
                .indices("ring_residue", &residues(table.ring_residue()))
                .single("distance", table.distance())
                .finish();
            Ok(Py::new(py, table)?.into_any())
        },
        context,
    )
}

/// Pairs of non-solvent atoms bridged by one water, from hydrogen bonds.
#[pyfunction]
#[pyo3(signature = (
    structure,
    *,
    max_distance,
    min_angle,
    backend="auto",
    policy=None,
    context=None,
))]
pub(crate) fn water_bridges(
    py: Python<'_>,
    structure: &PyStructure,
    max_distance: f32,
    min_angle: f64,
    backend: &str,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let kernel = molframe::analysis::water_bridges_kernel(WaterBridgeOptions {
        hydrogen_bonds: HydrogenBondOptions {
            maximum_donor_acceptor_distance: max_distance,
            minimum_angle_degrees: min_angle,
            backend: crate::backend::parse(backend)?,
            periodic: false,
        },
    });
    run(
        py,
        structure,
        &policy_of(policy),
        &kernel,
        |py, table| {
            let table = TableBuilder::new(py, table.len())
                .indices("water", &atoms(table.water()))
                .indices("first", &atoms(table.first()))
                .indices("second", &atoms(table.second()))
                .finish();
            Ok(Py::new(py, table)?.into_any())
        },
        context,
    )
}

/// Residues of two chains that lie within `cutoff` ångström of each other.
#[pyfunction]
#[pyo3(signature = (
    structure,
    *,
    first_chain,
    second_chain,
    cutoff,
    backend="auto",
    policy=None,
    context=None,
))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn chain_interface(
    py: Python<'_>,
    structure: &PyStructure,
    first_chain: &str,
    second_chain: &str,
    cutoff: f32,
    backend: &str,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let kernel = molframe::analysis::chain_interface_kernel(
        first_chain,
        second_chain,
        cutoff,
        crate::backend::parse(backend)?,
    );
    run(
        py,
        structure,
        &policy_of(policy),
        &kernel,
        |py, found| {
            let table = TableBuilder::new(py, found.len())
                .indices("residue", &residues(&found))
                .finish();
            Ok(Py::new(py, table)?.into_any())
        },
        context,
    )
}

/// Half-sphere exposure: α-carbons within `radius` ångström on each side of
/// every residue's side-chain direction.
#[pyfunction]
#[pyo3(signature = (structure, *, radius, backend="auto", policy=None, context=None))]
pub(crate) fn half_sphere_exposure(
    py: Python<'_>,
    structure: &PyStructure,
    radius: f32,
    backend: &str,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let kernel =
        molframe::analysis::half_sphere_exposure_kernel(radius, crate::backend::parse(backend)?);
    run(
        py,
        structure,
        &policy_of(policy),
        &kernel,
        |py, rows| {
            let residue: Vec<u32> = rows.iter().map(|row| row.residue.get()).collect();
            let upper: Vec<u32> = rows.iter().map(|row| row.upper).collect();
            let lower: Vec<u32> = rows.iter().map(|row| row.lower).collect();
            let table = TableBuilder::new(py, rows.len())
                .indices("residue", &residue)
                .indices("upper", &upper)
                .indices("lower", &lower)
                .finish();
            Ok(Py::new(py, table)?.into_any())
        },
        context,
    )
}

/// The seven backbone and glycosidic torsions of every nucleotide, in degrees.
///
/// A torsion that cannot be defined (a missing neighbour or atom) is `NaN`,
/// never a placeholder angle. Requires polymer atom-role annotations.
#[pyfunction]
#[pyo3(signature = (structure, *, policy=None, context=None))]
pub(crate) fn nucleic_torsions(
    py: Python<'_>,
    structure: &PyStructure,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let kernel = molframe::analysis::nucleic_torsions_kernel();
    run(
        py,
        structure,
        &policy_of(policy),
        &kernel,
        |py, rows| {
            let residue: Vec<u32> = rows.iter().map(|row| row.residue.get()).collect();
            let column =
                |pick: fn(&molframe::analysis::NucleicTorsions) -> Option<f64>| -> Vec<f64> {
                    rows.iter()
                        .map(|row| match pick(row) {
                            Some(angle) => angle,
                            None => f64::NAN,
                        })
                        .collect()
                };
            let table = TableBuilder::new(py, rows.len())
                .indices("residue", &residue)
                .double("alpha", &column(|row| row.alpha))
                .double("beta", &column(|row| row.beta))
                .double("gamma", &column(|row| row.gamma))
                .double("delta", &column(|row| row.delta))
                .double("epsilon", &column(|row| row.epsilon))
                .double("zeta", &column(|row| row.zeta))
                .double("chi", &column(|row| row.chi))
                .finish();
            Ok(Py::new(py, table)?.into_any())
        },
        context,
    )
}

/// The fraction of the reference's contacts that `target` keeps (Q).
///
/// A reference contact is a pair within `cutoff` ångström; it is kept when the
/// same pair is within `tolerance * cutoff` in `target`. The two structures must
/// share atom numbering. The value is a dictionary with `native`, `kept` and
/// `fraction`.
#[pyfunction]
#[pyo3(signature = (
    reference,
    target,
    *,
    cutoff,
    tolerance,
    backend="auto",
    policy=None,
    context=None,
))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn native_contacts(
    py: Python<'_>,
    reference: &PyStructure,
    target: &PyStructure,
    cutoff: f32,
    tolerance: f32,
    backend: &str,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let target = target.inner.clone();
    let kernel = molframe::analysis::native_contact_fraction_kernel(
        target.engine(),
        cutoff,
        tolerance,
        crate::backend::parse(backend)?,
    );
    run(
        py,
        reference,
        &policy_of(policy),
        &kernel,
        |py, found| {
            let value = PyDict::new(py);
            value.set_item("native", found.native)?;
            value.set_item("kept", found.kept)?;
            value.set_item("fraction", found.fraction)?;
            Ok(value.into_any().unbind())
        },
        context,
    )
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(contact_map, module)?)?;
    module.add_function(wrap_pyfunction!(pi_stacking, module)?)?;
    module.add_function(wrap_pyfunction!(cation_pi, module)?)?;
    module.add_function(wrap_pyfunction!(water_bridges, module)?)?;
    module.add_function(wrap_pyfunction!(chain_interface, module)?)?;
    module.add_function(wrap_pyfunction!(half_sphere_exposure, module)?)?;
    module.add_function(wrap_pyfunction!(nucleic_torsions, module)?)?;
    module.add_function(wrap_pyfunction!(native_contacts, module)?)
}
