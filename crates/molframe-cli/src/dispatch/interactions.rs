//! Routing for the interaction analyses: arguments become library options.

use crate::InteractionCommand;
use crate::args::{HydrogenBondArguments, InteractionInput, PlaneFitArguments};
use crate::exit::Exit;
use crate::interaction_commands::{self, Source};
use crate::report::Context;
use molframe::analysis::{CationPiOptions, HydrogenBondOptions, PiStackingOptions};
use molframe::geometry::EigenOptions;
use molframe::spatial::SpatialBackend;

pub(super) fn interactions(command: InteractionCommand, context: Context) -> Exit {
    match command {
        InteractionCommand::Hbonds { args } => {
            with_source(&args.source, context, |source, context| {
                interaction_commands::hbonds(source, hydrogen_bonds(&args), context)
            })
        }
        InteractionCommand::WaterBridges { args } => {
            with_source(&args.source, context, |source, context| {
                interaction_commands::water_bridges(source, hydrogen_bonds(&args), context)
            })
        }
        InteractionCommand::SaltBridges { args } => {
            with_source(&args.source, context, |source, context| {
                interaction_commands::salt_bridges(source, args.max_distance, context)
            })
        }
        InteractionCommand::PiStacking { args } => {
            with_source(&args.source, context, |source, context| {
                interaction_commands::pi_stacking(
                    source,
                    PiStackingOptions {
                        maximum_centre_distance: args.max_centre_distance,
                        maximum_parallel_angle: args.max_parallel_angle,
                        minimum_perpendicular_angle: args.min_perpendicular_angle,
                        plane_fit: plane_fit(&args.plane_fit),
                    },
                    context,
                )
            })
        }
        InteractionCommand::CationPi { args } => {
            with_source(&args.source, context, |source, context| {
                interaction_commands::cation_pi(
                    source,
                    CationPiOptions {
                        maximum_distance: args.max_distance,
                        maximum_face_angle: args.max_face_angle,
                        plane_fit: plane_fit(&args.plane_fit),
                    },
                    context,
                )
            })
        }
    }
}

fn with_source(
    input: &InteractionInput,
    context: Context,
    run: impl FnOnce(Source<'_>, Context) -> Exit,
) -> Exit {
    crate::chemistry::with_ccd(
        input.chemistry.ccd.as_deref(),
        input.chemistry.ccd_version.as_deref(),
        context,
        |ccd, version, context| {
            run(
                Source {
                    input: &input.input,
                    ccd,
                    ccd_version: version,
                },
                context,
            )
        },
    )
}

fn hydrogen_bonds(args: &HydrogenBondArguments) -> HydrogenBondOptions {
    HydrogenBondOptions {
        maximum_donor_acceptor_distance: args.max_donor_acceptor_distance,
        minimum_angle_degrees: args.min_angle_degrees,
        backend: SpatialBackend::Auto,
        periodic: args.periodic,
    }
}

const fn plane_fit(args: &PlaneFitArguments) -> EigenOptions {
    EigenOptions {
        relative_tolerance: args.eigen_relative_tolerance,
        maximum_sweeps: args.eigen_maximum_sweeps,
    }
}
