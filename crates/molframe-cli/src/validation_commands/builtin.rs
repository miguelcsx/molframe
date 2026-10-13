//! Checks backed by reference tables.
//!
//! Each check resolves its tables from the `--*-table` file, else from the
//! profile's choice in the cache (a missing table is a coded error that says how
//! to supply it), runs the Rust kernel and reports, besides its findings, which table or policy
//! decided them (name, version, source and licence).

use super::{Row, Subject};
use crate::exit::Exit;
use molframe::chemistry::{TableIdentity, TableRequest};
use molframe::validation::{
    BuiltinError, ClashScoreOptions, Layout, ReferenceSet, clash_score, geometry_scores,
    nucleic_table, percent, protein_table,
};

fn failure(what: &str, error: &BuiltinError) -> Exit {
    eprintln!("{what} failed: {error}");
    Exit::Policy
}

fn named<T>(path: Option<&std::path::PathBuf>) -> TableRequest<T> {
    match path {
        Some(path) => TableRequest::named(path.to_string_lossy().into_owned()),
        None => TableRequest::profile(),
    }
}

fn reference(subject: &Subject<'_>) -> ReferenceSet {
    let args = &subject.options.references;
    let mut set = ReferenceSet::default();
    set.rama = named(args.rama_table.as_ref());
    set.rotamers = named(args.rotamer_table.as_ref());
    set.geometry = named(args.geometry_table.as_ref());
    set.chi = named(args.chi_table.as_ref());
    set.nucleic = named(args.nucleic_table.as_ref());
    set
}

fn table_rows(check: &str, tables: &[TableIdentity]) -> impl Iterator<Item = Row> {
    tables.iter().map(move |info| {
        Row::info(
            check,
            format!("table:{}", info.name),
            format!(
                "version={} source={} licence={}",
                info.version, info.source, info.licence
            ),
        )
    })
}

fn rate(check: &str, item: &str, part: usize, whole: usize) -> Row {
    Row::info(
        check,
        item,
        format!("{part}/{whole} ({:.2}%)", percent(part, whole)),
    )
}

pub(super) fn protein_rows(subject: &Subject<'_>) -> Result<Vec<Row>, Exit> {
    let set = reference(subject);
    let layout = Layout::new(subject.structure.engine());
    let table = protein_table(&layout, &set, subject.context.execution)
        .map_err(|error| failure("protein validation", &error))?;
    let summary = table.summary();
    let mut rows = vec![
        rate(
            "ramachandran",
            "favoured",
            summary.rama_favoured,
            summary.rama_evaluated,
        ),
        rate(
            "ramachandran",
            "allowed",
            summary.rama_allowed,
            summary.rama_evaluated,
        ),
        rate(
            "ramachandran",
            "outlier",
            summary.rama_outliers,
            summary.rama_evaluated,
        ),
        rate(
            "rotamer",
            "outlier",
            summary.rotamer_outliers,
            summary.rotamer_evaluated,
        ),
        rate(
            "omega",
            "cis",
            summary.cis_peptides,
            summary.omega_evaluated,
        ),
        rate(
            "omega",
            "twisted",
            summary.twisted_peptides,
            summary.omega_evaluated,
        ),
    ];
    for residue in 0..layout.len() {
        let label = format!("{}#{residue}", layout.name(residue));
        if table.rama_class[residue] == 0 {
            let value = format!(
                "phi={:.1} psi={:.1}",
                table.phi[residue], table.psi[residue]
            );
            rows.push(Row::new("ramachandran-outlier", label.clone(), value));
        }
        if table.rotamer_class[residue] == 0 {
            rows.push(Row::new("rotamer-outlier", label.clone(), "outlier"));
        }
        if table.omega_class[residue] == 1 || table.omega_class[residue] == 2 {
            let kind = if table.omega_class[residue] == 1 {
                "cis"
            } else {
                "twisted"
            };
            let value = format!("{kind} omega={:.1}", table.omega[residue]);
            rows.push(Row::new("omega", label, value));
        }
    }
    rows.extend(table_rows("protein", &table.tables));
    Ok(rows)
}

pub(super) fn geometry_z_rows(subject: &Subject<'_>) -> Result<Vec<Row>, Exit> {
    let set = reference(subject);
    let layout = Layout::new(subject.structure.engine());
    let cutoff = subject.options.references.geometry_z_cutoff;
    let scores = geometry_scores(&layout, &set, cutoff, subject.context.execution)
        .map_err(|error| failure("geometry z-score", &error))?;
    let mut rows = vec![
        Row::info(
            "geometry-z",
            "bonds",
            format!("n={} rmsz={:.3}", scores.bond_count, scores.bond_rmsz),
        ),
        Row::info(
            "geometry-z",
            "angles",
            format!("n={} rmsz={:.3}", scores.angle_count, scores.angle_rmsz),
        ),
    ];
    for outlier in &scores.outliers {
        let kind = if outlier.is_angle() { "angle" } else { "bond" };
        rows.push(Row::new(
            format!("geometry-z-{kind}"),
            format!(
                "{}#{}",
                layout.name(outlier.residue as usize),
                outlier.residue
            ),
            format!(
                "z={:.2} observed={:.3} expected={:.3}",
                outlier.z, outlier.observed, outlier.expected
            ),
        ));
    }
    rows.extend(table_rows("geometry-z", &scores.tables));
    Ok(rows)
}

pub(super) fn clashscore_rows(subject: &Subject<'_>) -> Result<Vec<Row>, Exit> {
    let mut options = ClashScoreOptions::default();
    if let Some(tolerance) = subject.options.clash_tolerance {
        options.tolerance = tolerance;
    }
    let score = clash_score(
        subject.structure.engine(),
        options,
        subject.context.execution,
    )
    .map_err(|error| failure("clashscore", &error))?;
    let applied = score.options;
    Ok(vec![
        Row::info(
            "clashscore",
            "score",
            format!(
                "{:.2} per 1000 atoms ({} clashes, {} atoms)",
                score.score, score.clashes, score.atoms
            ),
        ),
        Row::info(
            "clashscore",
            "worst-overlap",
            format!("{:.2}", score.worst_overlap),
        ),
        Row::info(
            "clashscore",
            "policy",
            format!(
                "hydrogen_aware={} tolerance={} radii={:?} polar_margin={}",
                score.hydrogen_aware, applied.tolerance, applied.radius_set, applied.polar_margin
            ),
        ),
    ])
}

pub(super) fn nucleic_rows(subject: &Subject<'_>) -> Result<Vec<Row>, Exit> {
    let set = reference(subject);
    let layout = Layout::new(subject.structure.engine());
    let table = nucleic_table(&layout, &set, subject.context.execution)
        .map_err(|error| failure("nucleic conformation", &error))?;
    let summary = table.summary();
    let mut rows = vec![
        rate(
            "pucker",
            "c3-endo",
            summary.c3_endo,
            summary.pucker_evaluated,
        ),
        rate(
            "pucker",
            "c2-endo",
            summary.c2_endo,
            summary.pucker_evaluated,
        ),
        rate(
            "pucker",
            "other",
            summary.other_pucker,
            summary.pucker_evaluated,
        ),
        rate(
            "suite",
            "flagged",
            summary.suite_flagged,
            summary.suite_evaluated,
        ),
    ];
    for residue in 0..layout.len() {
        if table.suite_class[residue] == 1 {
            rows.push(Row::new(
                "suite-flag",
                format!("{}#{residue}", layout.name(residue)),
                format!("flags={:#05b}", table.suite_flags[residue]),
            ));
        }
    }
    rows.extend(table_rows("nucleic", &table.tables));
    Ok(rows)
}
