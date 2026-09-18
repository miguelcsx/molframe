//! Routing for geometry, ensemble, trajectory and system command groups.

use crate::exit::Exit;
use crate::report::Context;
use crate::{EnsembleCommand, GeometryCommand, ProfileCommand, SystemCommand, TrajectoryCommand};

pub(super) fn profile(command: ProfileCommand, context: Context) -> Exit {
    match command {
        ProfileCommand::List => {
            let profile = molframe::ProfileId::DEFAULT.as_str();
            if context.is_json() {
                let mut json = crate::report::Json::new();
                json.raw("profiles", &format!("[\"{profile}\"]"));
                context.result(&json.finish());
            } else {
                context.result(profile);
            }
            Exit::Success
        }
    }
}

pub(super) fn geometry(command: GeometryCommand, context: Context) -> Exit {
    match command {
        GeometryCommand::Torsions { input, chemistry } => crate::chemistry::with_ccd(
            chemistry.ccd.as_deref(),
            chemistry.ccd_version.as_deref(),
            context,
            |ccd, version, context| {
                crate::intrinsic_commands::torsions(&input, ccd, version, context)
            },
        ),
        GeometryCommand::Helix {
            input,
            chemistry,
            eigen_relative_tolerance,
            eigen_maximum_sweeps,
        } => crate::chemistry::with_ccd(
            chemistry.ccd.as_deref(),
            chemistry.ccd_version.as_deref(),
            context,
            |ccd, version, context| {
                crate::intrinsic_commands::helix(
                    &input,
                    ccd,
                    version,
                    molframe::EigenOptions {
                        relative_tolerance: eigen_relative_tolerance,
                        maximum_sweeps: eigen_maximum_sweeps,
                    },
                    context,
                )
            },
        ),
        GeometryCommand::Surface { args } => surface(&args, context),
    }
}

pub(super) fn surface(args: &crate::SurfaceArguments, context: Context) -> Exit {
    crate::intrinsic_commands::surface(
        &args.input,
        crate::intrinsic_commands::SurfaceOptions {
            probe: args.probe,
            resolution: args.resolution,
            max_cells: args.max_cells,
            radii: args.radii.into(),
            source_vertex: args.source_vertex,
            patch_radius: args.patch_radius,
            output: args.output.as_deref(),
        },
        context,
    )
}

pub(super) fn ensemble(command: EnsembleCommand, context: Context) -> Exit {
    match command {
        EnsembleCommand::Pca {
            input,
            components,
            memory_limit,
            no_fit,
            fit_tolerance,
            fit_iterations,
        } => crate::intrinsic_commands::pca(
            &input,
            components,
            memory_limit,
            no_fit,
            fit_tolerance,
            fit_iterations,
            context,
        ),
        EnsembleCommand::TorsionPca {
            input,
            components,
            memory_limit,
            chemistry,
        } => crate::chemistry::with_ccd(
            chemistry.ccd.as_deref(),
            chemistry.ccd_version.as_deref(),
            context,
            |ccd, version, context| {
                crate::intrinsic_commands::torsion_pca(
                    &input,
                    components,
                    memory_limit,
                    ccd,
                    version,
                    context,
                )
            },
        ),
        EnsembleCommand::Diffusion {
            input,
            epsilon,
            time,
            dimensions,
            memory_limit,
        } => crate::intrinsic_commands::diffusion(
            &input,
            epsilon,
            time,
            dimensions,
            memory_limit,
            context,
        ),
    }
}

pub(super) fn trajectory(command: TrajectoryCommand, context: Context) -> Exit {
    match command {
        TrajectoryCommand::Info {
            input,
            topology,
            gsd_length_scale,
        } => crate::trajectory::info(&input, topology.as_deref(), gsd_length_scale, context),
        TrajectoryCommand::Convert {
            input,
            output,
            stride,
            gsd_length_scale,
            trz_title,
        } => crate::trajectory::convert(
            &input,
            &output,
            stride,
            gsd_length_scale,
            trz_title.as_deref(),
            context,
        ),
        TrajectoryCommand::Extract {
            input,
            output,
            frames,
            gsd_length_scale,
            trz_title,
        } => crate::trajectory::extract(
            &input,
            &output,
            &frames,
            gsd_length_scale,
            trz_title.as_deref(),
            context,
        ),
        TrajectoryCommand::Contacts { input, cutoff } => {
            crate::trajectory::contacts(&input, cutoff, context)
        }
        TrajectoryCommand::Sasa {
            input,
            topology,
            probe,
            samples,
            radii,
        } => crate::trajectory::sasa(&input, &topology, probe, samples, radii.into(), context),
        TrajectoryCommand::Rmsd {
            input,
            reference,
            no_fit,
            gsd_length_scale,
        } => crate::trajectory::rmsd(&input, reference, no_fit, gsd_length_scale, context),
    }
}

pub(super) fn system(command: SystemCommand, context: Context) -> Exit {
    match command {
        SystemCommand::Info { input } => crate::system::info(&input, context),
        SystemCommand::Copy { input, output } => crate::system::copy(&input, &output, context),
    }
}
