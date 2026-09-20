//! CLI projection over intrinsic geometry and ensemble library workflows.

use crate::commands::open;
use crate::exit::Exit;
use crate::report::{Context, Json};
use molframe::chemistry::{RadiusSet, vdw_radius};
use molframe::geometry::PeriodicAngle;
use molframe::surface::{edge_geodesic_distances, surface_curvatures, surface_patch};
use molframe::trajectory::{
    CartesianFit, cartesian_pca, diffusion_map, dihedral_pca, pairwise_fitted_rmsd,
};
use molframe::{ModelIndex, Structure};
use std::fmt::Write as _;
use std::path::Path;

pub(super) fn torsions(path: &Path, ccd: &Path, ccd_version: &str, context: Context) -> Exit {
    let structure = match open(path, context) {
        Ok(value) => value,
        Err(exit) => return exit,
    };
    let structure = match crate::chemistry::annotate(&structure, ccd, ccd_version, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let values = match molframe::structure_backbone_torsions(structure.engine()) {
        Ok(values) => values,
        Err(findings) => {
            context.findings(&findings, &path.display().to_string());
            return Exit::Consistency;
        }
    };
    let complete = values
        .iter()
        .filter(|record| {
            record.torsions.phi.is_some()
                && record.torsions.psi.is_some()
                && record.torsions.omega.is_some()
        })
        .count();
    emit_summary(
        context,
        &[
            ("residues", values.len().to_string()),
            ("complete", complete.to_string()),
            ("missing", (values.len() - complete).to_string()),
            ("torsion_set", "phi,psi,omega".to_owned()),
        ],
    );
    Exit::Success
}

pub(super) fn helix(
    path: &Path,
    ccd: &Path,
    ccd_version: &str,
    eigen: molframe::geometry::EigenOptions,
    context: Context,
) -> Exit {
    let structure = match open(path, context) {
        Ok(value) => value,
        Err(exit) => return exit,
    };
    let structure = match crate::chemistry::annotate(&structure, ccd, ccd_version, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let mut alpha_count = 0usize;
    let mut frame_count = 0usize;
    let mut helix_count = 0usize;
    let mut chain_count = 0usize;
    let traces = match molframe::structure_protein_alpha_traces(structure.engine()) {
        Ok(traces) => traces,
        Err(findings) => {
            context.findings(&findings, &path.display().to_string());
            return Exit::Consistency;
        }
    };
    for trace in traces {
        let chain_positions = trace.positions;
        alpha_count += chain_positions.iter().flatten().count();
        frame_count += molframe::geometry::backbone_frames(&chain_positions)
            .iter()
            .flatten()
            .count();
        match molframe::geometry::helix_geometry_with_options(&chain_positions, eigen) {
            Ok(Some(_)) => helix_count += 1,
            Ok(None) => {}
            Err(error) => {
                eprintln!("helix eigensolver failed: {error:?}");
                return Exit::Consistency;
            }
        }
        chain_count += 1;
    }
    emit_summary(
        context,
        &[
            ("ca_atoms", alpha_count.to_string()),
            ("frenet_frames", frame_count.to_string()),
            ("helix_chains", helix_count.to_string()),
            ("coverage", format!("{helix_count}/{chain_count}")),
        ],
    );
    Exit::Success
}

#[derive(Clone, Copy)]
pub(super) struct SurfaceOptions<'a> {
    pub(super) probe: f32,
    pub(super) resolution: f32,
    pub(super) max_cells: usize,
    pub(super) radii: RadiusSet,
    pub(super) source_vertex: Option<u32>,
    pub(super) patch_radius: Option<f64>,
    pub(super) output: Option<&'a Path>,
}

pub(super) fn surface(path: &Path, options: SurfaceOptions<'_>, context: Context) -> Exit {
    let patch_request = match (options.source_vertex, options.patch_radius) {
        (Some(source), Some(radius)) => Some((source, radius)),
        (None, None) => None,
        _ => {
            eprintln!("--source-vertex and --patch-radius must be supplied together");
            return Exit::Usage;
        }
    };
    let structure = match open(path, context) {
        Ok(value) => value,
        Err(exit) => return exit,
    };
    let mut positions = Vec::new();
    let mut radii = Vec::new();
    for index in 0..structure.atom_count() {
        let Some(atom) = structure.atom_at(index as usize) else {
            continue;
        };
        let (Some(position), Some(element)) = (atom.position(), atom.element()) else {
            continue;
        };
        let Some(radius) = vdw_radius(element, options.radii) else {
            eprintln!("selected radius set has no value for atom {index}");
            return Exit::Indeterminate;
        };
        positions.push(position);
        radii.push(radius);
    }
    let Ok(surface) = molframe::surface::solvent_excluded_surface_with_options(
        &positions,
        &radii,
        options.probe,
        molframe::surface::SurfaceGridOptions {
            resolution: options.resolution,
            max_cells: options.max_cells,
            max_workspace_bytes: molframe::surface::SurfaceGridOptions::STANDARD_WORKSPACE_BYTES,
        },
    ) else {
        eprintln!("surface construction failed for the supplied probe or resolution");
        return Exit::Consistency;
    };
    let mesh = surface.indexed_mesh();
    let Ok(curvature) = surface_curvatures(&mesh) else {
        eprintln!("surface mesh is not manifold enough for curvature");
        return Exit::Consistency;
    };
    if let Some(output) = options.output {
        if output.extension().and_then(std::ffi::OsStr::to_str) != Some("obj") {
            eprintln!("surface output currently requires an .obj destination");
            return Exit::Usage;
        }
        if let Err(error) = molframe::surface::write_obj(output, &mesh) {
            eprintln!("surface output failed: {error}");
            return Exit::Failure;
        }
    }
    let mut rows = vec![
        ("atoms", positions.len().to_string()),
        ("vertices", mesh.vertices.len().to_string()),
        ("faces", mesh.faces.len().to_string()),
        ("components", mesh.report.components.to_string()),
        ("curvature_values", curvature.len().to_string()),
        ("radii_set", radius_set_name(options.radii).to_owned()),
        ("probe", options.probe.to_string()),
        ("resolution", options.resolution.to_string()),
    ];
    if let Some((source, radius)) = patch_request {
        let Ok(geodesics) = edge_geodesic_distances(&mesh, source) else {
            eprintln!("surface source vertex is invalid");
            return Exit::Usage;
        };
        let Ok(patch_vertices) = surface_patch(&mesh, source, radius) else {
            eprintln!("surface patch radius is invalid");
            return Exit::Usage;
        };
        rows.push(("patch_vertices", patch_vertices.len().to_string()));
        rows.push((
            "unreachable_vertices",
            geodesics
                .distances
                .iter()
                .filter(|distance| !distance.is_finite())
                .count()
                .to_string(),
        ));
        rows.push(("geodesic", "mesh-edge-dijkstra".to_owned()));
    }
    if let Some(output) = options.output {
        rows.push(("output", output.display().to_string()));
    }
    emit_summary(context, &rows);
    Exit::Success
}

pub(super) fn pca(
    path: &Path,
    components: usize,
    memory_limit: usize,
    no_fit: bool,
    fit_tolerance: Option<f64>,
    fit_iterations: Option<usize>,
    context: Context,
) -> Exit {
    let structure = match open(path, context) {
        Ok(value) => value,
        Err(exit) => return exit,
    };
    let Some(frames) = dense_frames(&structure) else {
        eprintln!("ensemble input must be a dense multi-model structure");
        return Exit::Consistency;
    };
    let fit = match (no_fit, fit_tolerance, fit_iterations) {
        (true, None, None) => CartesianFit::None,
        (false, Some(tolerance), Some(max_iterations)) => CartesianFit::IterativeMean {
            tolerance,
            max_iterations,
        },
        (true, _, _) => {
            eprintln!("--fit-tolerance and --fit-iterations cannot be used with --no-fit");
            return Exit::Usage;
        }
        (false, _, _) => {
            eprintln!("fitted PCA requires --fit-tolerance and --fit-iterations");
            return Exit::Usage;
        }
    };
    let Ok(result) = cartesian_pca(&frames, fit, components, memory_limit) else {
        eprintln!("PCA failed: check model correspondence, components and memory limit");
        return Exit::Consistency;
    };
    emit_ensemble(
        context,
        "cartesian-pca",
        frames.len(),
        memory_limit,
        if no_fit { "none" } else { "iterative-mean" },
        &result.eigenvalues,
    );
    Exit::Success
}

pub(super) fn torsion_pca(
    path: &Path,
    components: usize,
    memory_limit: usize,
    ccd: &Path,
    ccd_version: &str,
    context: Context,
) -> Exit {
    let structure = match open(path, context) {
        Ok(value) => value,
        Err(exit) => return exit,
    };
    let structure = match crate::chemistry::annotate(&structure, ccd, ccd_version, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let mut observations = Vec::new();
    for model in 0..structure.model_count() {
        let Ok(model) = u32::try_from(model) else {
            eprintln!("model count exceeds the supported index range");
            return Exit::Resource;
        };
        let mut row = Vec::new();
        let records = match molframe::structure_backbone_torsions_model(
            structure.engine(),
            ModelIndex::new(model),
        ) {
            Ok(records) => records,
            Err(findings) => {
                context.findings(&findings, &path.display().to_string());
                return Exit::Consistency;
            }
        };
        for record in records {
            for angle in [
                record.torsions.phi,
                record.torsions.psi,
                record.torsions.omega,
            ] {
                let Some(angle) = angle else {
                    eprintln!(
                        "torsion PCA requires the same complete torsion mapping in every model"
                    );
                    return Exit::Consistency;
                };
                let Ok(periodic) = PeriodicAngle::from_radians(angle) else {
                    return Exit::Consistency;
                };
                row.push(periodic);
            }
        }
        observations.push(row);
    }
    let Ok(result) = dihedral_pca(&observations, components, memory_limit) else {
        return Exit::Consistency;
    };
    emit_ensemble(
        context,
        "torsion-pca",
        observations.len(),
        memory_limit,
        "cos-sin",
        &result.eigenvalues,
    );
    Exit::Success
}

const fn radius_set_name(radius_set: RadiusSet) -> &'static str {
    match radius_set {
        RadiusSet::Bondi => "bondi",
        RadiusSet::AmberUnited => "amber-united",
        RadiusSet::Charmm => "charmm",
        RadiusSet::Alvarez => "alvarez",
    }
}

pub(super) fn diffusion(
    path: &Path,
    epsilon: f64,
    time: u32,
    dimensions: usize,
    memory_limit: usize,
    context: Context,
) -> Exit {
    let structure = match open(path, context) {
        Ok(value) => value,
        Err(exit) => return exit,
    };
    let Some(frames) = dense_frames(&structure) else {
        return Exit::Consistency;
    };
    let Ok(distances) = pairwise_fitted_rmsd(&frames, memory_limit) else {
        return Exit::Consistency;
    };
    let Ok(result) = diffusion_map(&distances, epsilon, time, dimensions) else {
        return Exit::Consistency;
    };
    emit_summary(
        context,
        &[
            ("algorithm", "diffusion-map".to_owned()),
            ("frames", frames.len().to_string()),
            ("metric", "fitted-rmsd".to_owned()),
            ("alignment", "qcp".to_owned()),
            ("epsilon", epsilon.to_string()),
            ("time", time.to_string()),
            ("dimensions", result.eigenvalues.len().to_string()),
            (
                "graph_components",
                result
                    .graph_components
                    .iter()
                    .max()
                    .map_or(0, |value| value + 1)
                    .to_string(),
            ),
            ("memory_limit", memory_limit.to_string()),
        ],
    );
    Exit::Success
}

fn dense_frames(structure: &Structure) -> Option<Vec<Vec<[f32; 3]>>> {
    (0..structure.model_count())
        .map(|model| {
            let Ok(model) = u32::try_from(model) else {
                return None;
            };
            structure
                .model_coordinates(ModelIndex::new(model))
                .map(<[_]>::to_vec)
        })
        .collect()
}

fn emit_ensemble(
    context: Context,
    algorithm: &str,
    frames: usize,
    memory_limit: usize,
    alignment: &str,
    eigenvalues: &[f64],
) {
    emit_summary(
        context,
        &[
            ("algorithm", algorithm.to_owned()),
            ("frames", frames.to_string()),
            ("alignment", alignment.to_owned()),
            ("selection", "all-atoms".to_owned()),
            ("memory_limit", memory_limit.to_string()),
            (
                "eigenvalues",
                eigenvalues
                    .iter()
                    .map(|value| format!("{value:.8}"))
                    .collect::<Vec<_>>()
                    .join(","),
            ),
        ],
    );
}

fn emit_summary(context: Context, fields: &[(&str, String)]) {
    if context.is_json() {
        let mut object = Json::new();
        for (key, value) in fields {
            object.text(key, value);
        }
        context.result(&object.finish());
    } else if let Some(delimiter) = context.delimiter() {
        let header = fields
            .iter()
            .map(|(key, _)| *key)
            .collect::<Vec<_>>()
            .join(&delimiter.to_string());
        let values = fields
            .iter()
            .map(|(_, value)| value.as_str())
            .collect::<Vec<_>>()
            .join(&delimiter.to_string());
        context.result(&format!("{header}\n{values}"));
    } else {
        let mut output = String::new();
        for (index, (key, value)) in fields.iter().enumerate() {
            if index != 0 {
                output.push('\n');
            }
            let _ = write!(output, "{key:<18} {value}");
        }
        context.result(&output);
    }
}
