//! Exhaustive routing from parsed commands to library-backed projections.

use crate::exit::Exit;
use crate::report::Context;
use crate::{CcdCommand, Command};
use clap::CommandFactory as _;

mod domains;
use domains::{ensemble, geometry, profile, surface, system, trajectory};

pub(crate) fn execute(command: Command, context: Context) -> Exit {
    if context.is_table_file()
        && !matches!(
            command,
            Command::Select {
                output: Some(_),
                count: false,
                ..
            }
        )
    {
        eprintln!("Arrow and Parquet output require `select -o PATH`");
        return Exit::Usage;
    }
    match command {
        Command::Info { input, detail } => crate::commands::info(&input, detail, context),
        Command::Convert {
            input,
            output,
            chain_map,
            hybrid36,
            preserve,
            cif_block_id,
            generate_cif_connection_ids,
            cif_connection_type,
        } => {
            let mut cif_options = pdbiox::CifWriteOptions::new();
            if let Some(block_id) = cif_block_id {
                cif_options = cif_options.with_block_id(block_id);
            }
            if generate_cif_connection_ids {
                cif_options = cif_options.with_generated_connection_ids();
            }
            if let Some(connection_type) = cif_connection_type {
                cif_options = cif_options.with_connection_type_id(connection_type);
            }
            crate::commands::convert(
                &input,
                &output,
                &chain_map,
                hybrid36,
                preserve,
                &cif_options,
                context,
            )
        }
        Command::Validate {
            input,
            checks,
            chemistry,
            bond_tolerance,
            planarity_tolerance,
            plane_relative_tolerance,
            plane_maximum_sweeps,
            clash_tolerance,
            radii,
            altloc_expected_sum,
            altloc_tolerance,
            b_factor_z_score,
        } => crate::validation_commands::validate(
            &input,
            crate::validation_commands::ValidationOptions {
                checks: &checks,
                ccd: chemistry.ccd.as_deref(),
                ccd_version: chemistry.ccd_version.as_deref(),
                bond_tolerance,
                planarity_tolerance,
                plane_relative_tolerance,
                plane_maximum_sweeps,
                clash_tolerance,
                radii: radii.map(Into::into),
                altloc_expected_sum,
                altloc_tolerance,
                b_factor_z_score,
            },
            context,
        ),
        Command::Measure { input } => crate::commands::measure(&input, context),
        Command::Rmsd {
            mobile,
            reference,
            no_fit,
            on,
        } => crate::comparison_commands::rmsd(&mobile, &reference, no_fit, on.as_deref(), context),
        Command::Policy => crate::commands::policy(context),
        Command::Profile { command } => profile(command, context),
        Command::Geometry { command } => geometry(command, context),
        Command::Surface { args } => surface(&args, context),
        Command::Ensemble { command } => ensemble(command, context),
        Command::Trajectory { command } => trajectory(command, context),
        Command::System { command } => system(command, context),
        command => execute_data(command, context),
    }
}

fn execute_data(command: Command, context: Context) -> Exit {
    match command {
        Command::Categories { input } => crate::inspect::categories(&input, context),
        Command::Sequence { command } => crate::sequence_commands::execute(command, context),
        Command::Assemblies { input } => crate::inspect::assemblies(&input, context),
        Command::Show {
            input,
            category,
            limit,
        } => crate::inspect::show(&input, &category, limit, context),
        Command::Contacts {
            input,
            between,
            cutoff,
        } => crate::analysis_commands::contacts(&input, between.as_deref(), cutoff, context),
        Command::Neighbors {
            input,
            query,
            cutoff,
        } => crate::analysis_commands::neighbors(&input, &query, cutoff, context),
        Command::Sse {
            input,
            chemistry,
            electrostatic_prefactor,
            hydrogen_bond_energy,
            amide_hydrogen_distance,
            minimum_sequence_separation,
            helix_offset,
            turn_offsets,
        } => crate::chemistry::with_ccd(
            chemistry.ccd.as_deref(),
            chemistry.ccd_version.as_deref(),
            context,
            |ccd, version, context| {
                crate::analysis_commands::sse(
                    &input,
                    crate::analysis_commands::SseOptions {
                        ccd,
                        ccd_version: version,
                        electrostatic_prefactor,
                        hydrogen_bond_energy,
                        amide_hydrogen_distance,
                        minimum_sequence_separation,
                        helix_offset,
                        turn_offsets: &turn_offsets,
                    },
                    context,
                )
            },
        ),
        Command::Torsions { input, chemistry } => crate::chemistry::with_ccd(
            chemistry.ccd.as_deref(),
            chemistry.ccd_version.as_deref(),
            context,
            |ccd, version, context| {
                crate::intrinsic_commands::torsions(&input, ccd, version, context)
            },
        ),
        Command::Interfaces {
            input,
            between,
            cutoff,
        } => crate::analysis_commands::interfaces(&input, &between, cutoff, context),
        Command::Sasa {
            input,
            probe,
            points,
            radii,
        } => crate::analysis_commands::sasa(&input, probe, points, radii.into(), context),
        command => execute_comparison(command, context),
    }
}

fn execute_comparison(command: Command, context: Context) -> Exit {
    match command {
        Command::Compare {
            model,
            reference,
            metrics,
            lddt_radius,
            lddt_minimum_distance,
            lddt_tolerances,
            lddt_empty,
            receptor,
            ligand,
            contact_distance,
            ligand_scale,
            interface_scale,
        } => crate::comparison_commands::compare(
            &model,
            &reference,
            crate::comparison_commands::ComparisonOptions {
                metrics: &metrics,
                lddt_radius,
                lddt_minimum_distance,
                lddt_tolerances: &lddt_tolerances,
                lddt_empty: lddt_empty.map(Into::into),
                receptor: receptor.as_deref(),
                ligand: ligand.as_deref(),
                contact_distance,
                ligand_scale,
                interface_scale,
            },
            context,
        ),
        Command::Diff {
            left,
            right,
            coordinate_tolerance,
        } => crate::diff_commands::diff(&left, &right, coordinate_tolerance, context),
        Command::Select {
            input,
            query,
            output,
            count,
        } => crate::selection_commands::select(&input, &query, output.as_deref(), count, context),
        Command::Superpose {
            mobile,
            reference,
            output,
            on,
            collinear_relative_tolerance,
            eigen_relative_tolerance,
            eigen_maximum_sweeps,
        } => crate::comparison_commands::superpose(
            &mobile,
            &reference,
            &output,
            &on,
            pdbiox::SuperposeOptions {
                collinear_relative_tolerance,
                eigen: pdbiox::EigenOptions {
                    relative_tolerance: eigen_relative_tolerance,
                    maximum_sweeps: eigen_maximum_sweeps,
                },
            },
            context,
        ),
        Command::MapChains {
            model,
            reference,
            chemistry,
            min_identity,
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
        } => crate::chemistry::with_ccd(
            chemistry.ccd.as_deref(),
            chemistry.ccd_version.as_deref(),
            context,
            |ccd, version, context| {
                crate::comparison_commands::map_chains(
                    &model,
                    &reference,
                    ccd,
                    version,
                    min_identity,
                    pdbiox::seq::Scoring {
                        match_score,
                        mismatch_score,
                        gap_open,
                        gap_extend,
                    },
                    context,
                )
            },
        ),
        command => execute_utilities(command, context),
    }
}

fn execute_utilities(command: Command, context: Context) -> Exit {
    match command {
        Command::Completions { shell } => crate::man_commands::completions(shell),
        Command::Ccd { command } => ccd(command, context),
        Command::Fetch {
            id,
            url_template,
            sha256,
            max_bytes,
            timeout_seconds,
            redirect_limit,
            output,
        } => crate::network_commands::fetch(
            &id,
            &url_template,
            &sha256,
            crate::network_commands::NetworkOptions {
                max_bytes,
                timeout_seconds,
                redirect_limit,
            },
            &output,
            context,
        ),
        Command::Audit { args } => crate::audit_commands::audit_selection(&args, context),
        Command::Fx { command } => crate::fx_commands::execute(command, context),
        Command::Batch { command } => crate::batch_commands::execute(command, context),
        Command::Man { outdir } => crate::man_commands::generate(&crate::Cli::command(), &outdir),
        _ => {
            eprintln!("internal command routing error");
            Exit::Failure
        }
    }
}

fn ccd(command: CcdCommand, context: Context) -> Exit {
    match command {
        CcdCommand::Get {
            component,
            chemistry,
        } => crate::chemistry::with_ccd(
            chemistry.ccd.as_deref(),
            chemistry.ccd_version.as_deref(),
            context,
            |ccd, version, context| crate::chemistry::component(ccd, version, &component, context),
        ),
        CcdCommand::Update {
            url,
            sha256,
            max_bytes,
            timeout_seconds,
            redirect_limit,
            version,
            output,
            replace,
        } => crate::network_commands::update_ccd(
            &url,
            &sha256,
            crate::network_commands::NetworkOptions {
                max_bytes,
                timeout_seconds,
                redirect_limit,
            },
            &version,
            &output,
            replace,
            context,
        ),
    }
}
