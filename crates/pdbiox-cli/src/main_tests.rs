use super::{Cli, Command, TrajectoryCommand};
use clap::Parser;

#[test]
fn traj_alias_has_the_same_declarative_command() {
    let cli = Cli::try_parse_from(["pdbiox", "traj", "info", "run.xtc"])
        .expect("trajectory alias should parse");
    assert!(matches!(
        cli.command,
        Command::Trajectory {
            command: TrajectoryCommand::Info { .. }
        }
    ));
}

#[test]
fn zero_stride_is_rejected_during_argument_parsing() {
    let result = Cli::try_parse_from([
        "pdbiox",
        "trajectory",
        "convert",
        "in.xtc",
        "out.dcd",
        "--stride",
        "0",
    ]);
    assert!(result.is_err());
}

#[test]
fn output_format_is_optional_for_context_sensitive_defaulting() {
    let cli = Cli::try_parse_from(["pdbiox", "policy"]).expect("policy should parse");
    assert_eq!(cli.global.format, None);
}

#[test]
fn diff_requires_an_explicit_coordinate_tolerance() {
    assert!(Cli::try_parse_from(["pdbiox", "diff", "left.cif", "right.cif"]).is_err());
    let parsed = Cli::try_parse_from([
        "pdbiox",
        "diff",
        "left.cif",
        "right.cif",
        "--coordinate-tolerance",
        "0.01",
    ]);
    assert!(matches!(
        parsed,
        Ok(Cli {
            command: Command::Diff { .. },
            ..
        })
    ));
}

#[test]
fn select_accepts_declarative_query_and_binary_table_output() {
    let parsed = Cli::try_parse_from([
        "pdbiox",
        "select",
        "entry.cif",
        "protein and chain A",
        "--format",
        "parquet",
        "-o",
        "protein.parquet",
    ]);
    assert!(matches!(
        parsed,
        Ok(Cli {
            command: Command::Select { .. },
            ..
        })
    ));
}

#[test]
fn rmsd_accepts_one_declarative_selection_for_both_structures() {
    let parsed = Cli::try_parse_from([
        "pdbiox",
        "rmsd",
        "model.cif",
        "reference.cif",
        "--on",
        "protein and backbone",
    ]);
    assert!(matches!(
        parsed,
        Ok(Cli {
            command: Command::Rmsd { on: Some(_), .. },
            ..
        })
    ));
}

#[test]
fn dockq_accepts_only_explicit_scientific_controls() {
    let parsed = Cli::try_parse_from([
        "pdbiox",
        "compare",
        "model.cif",
        "native.cif",
        "--metrics",
        "dockq",
        "--receptor",
        "A",
        "--ligand",
        "B",
        "--contact-distance",
        "5.0",
        "--ligand-scale",
        "8.5",
        "--interface-scale",
        "1.5",
    ]);
    assert!(matches!(
        parsed,
        Ok(Cli {
            command: Command::Compare { .. },
            ..
        })
    ));
}

#[test]
fn contacts_accepts_two_declarative_selections() {
    let parsed = Cli::try_parse_from([
        "pdbiox",
        "contacts",
        "entry.cif",
        "--between",
        "chain A",
        "chain B",
        "--cutoff",
        "4.5",
    ]);
    assert!(matches!(
        parsed,
        Ok(Cli {
            command: Command::Contacts {
                between: Some(_),
                ..
            },
            ..
        })
    ));
}

#[test]
fn surface_has_a_short_top_level_form_with_explicit_physics() {
    let parsed = Cli::try_parse_from([
        "pdbiox",
        "surface",
        "entry.cif",
        "--probe",
        "1.4",
        "--resolution",
        "0.5",
        "--max-cells",
        "1000000",
        "--radii",
        "bondi",
        "-o",
        "surface.obj",
    ]);
    assert!(matches!(
        parsed,
        Ok(Cli {
            command: Command::Surface { .. },
            ..
        })
    ));
}

#[test]
fn sse_refuses_an_implicit_scientific_definition() {
    assert!(
        Cli::try_parse_from([
            "pdbiox",
            "sse",
            "entry.cif",
            "--ccd",
            "components.cif",
            "--ccd-version",
            "2026-08-01",
        ])
        .is_err()
    );
}

#[test]
fn advanced_workflows_have_declarative_entrypoints() {
    assert!(
        Cli::try_parse_from([
            "pdbiox",
            "audit",
            "entry.cif",
            "protein",
            "--model-values",
            "first,all",
            "--max-runs",
            "16",
        ])
        .is_ok()
    );
    assert!(
        Cli::try_parse_from([
            "pdbiox",
            "fx",
            "evaluate",
            "entry.cif",
            "--motif",
            "motif.toml",
            "--ccd",
            "components.cif",
            "--ccd-version",
            "v1",
            "--mapping-limit",
            "64",
            "--measurement-limit",
            "256",
            "--plane-relative-tolerance",
            "1e-14",
            "--plane-maximum-sweeps",
            "24",
        ])
        .is_ok()
    );
    assert!(Cli::try_parse_from(["pdbiox", "sequence", "extract", "entry.cif"]).is_ok());
    assert!(
        Cli::try_parse_from([
            "pdbiox",
            "sequence",
            "extract",
            "entry.cif",
            "--ccd",
            "components.cif"
        ])
        .is_err()
    );
    assert!(
        Cli::try_parse_from([
            "pdbiox", "sequence", "convert", "-", "--from", "fasta", "--to", "clustal"
        ])
        .is_ok()
    );
    assert!(
        Cli::try_parse_from([
            "pdbiox",
            "sequence",
            "kmer",
            "sequences.fa",
            "--k",
            "5",
            "--operation",
            "minimizers"
        ])
        .is_err()
    );
    assert!(Cli::try_parse_from(["pdbiox", "batch", "info", "*.cif", "--threads", "4",]).is_ok());
    assert!(Cli::try_parse_from(["pdbiox", "man", "--outdir", "man"]).is_ok());
}

#[test]
fn global_result_and_provenance_destinations_parse() {
    let parsed = Cli::try_parse_from([
        "pdbiox",
        "info",
        "entry.cif",
        "--result-output",
        "result.json",
        "--provenance",
        "result.provenance.json",
        "--format",
        "json",
        "--infer-elements",
        "--infer-residue-boundaries",
    ]);
    assert!(parsed.is_ok());
}

#[test]
fn global_policy_overrides_are_typed_and_strict() {
    assert!(
        Cli::try_parse_from([
            "pdbiox",
            "info",
            "entry.cif",
            "--assembly",
            "crystal:8.5",
            "--model",
            "index:2",
            "--altloc",
            "label:A",
            "--identifiers",
            "auth",
        ])
        .is_ok()
    );
    assert!(
        Cli::try_parse_from(["pdbiox", "info", "entry.cif", "--assembly", "crystal:nan"]).is_err()
    );
    assert!(
        Cli::try_parse_from(["pdbiox", "info", "entry.cif", "--identifiers", "automatic"]).is_err()
    );
}
