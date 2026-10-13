//! Feature-level round-trip audit of real files: alternate locations,
//! insertion codes, models, chain and sequence naming, B-factors, charges,
//! hetero flags, bonds and header metadata must survive cif, bcif and pdb
//! in every direction, and whatever a format cannot hold must be reported.

use molframe::formats::pdb::{PdbIdentifierNamespace, PdbOptions};
use molframe::{ReadOptions, Structure};
use std::collections::BTreeSet;
use std::fmt::Write as _;
use std::path::PathBuf;

const CORPUS: [&str; 7] = [
    "1crn.cif",
    "1ubq.cif",
    "4hhb.cif",
    "1aon.bcif",
    "1crn.pdb",
    "4hhb.pdb",
    "2m7c.pdb",
];
const OPTIONAL: [&str; 9] = [
    "2m7c.cif", "1igt.cif", "1a02.cif", "1bna.cif", "2hhb.cif", "3nir.cif", "1ejg.cif", "1d3z.cif",
    "5xnl.cif",
];

fn read(bytes: Vec<u8>, name: &str) -> Structure {
    match molframe::read_bytes(bytes, Some(name), &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(failure) => panic!("{name}: {failure:?}"),
    }
}

fn corpus() -> Vec<(String, Structure)> {
    let bench = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../molframe-bench/data");
    let extra = PathBuf::from("/tmp/rt-data/rcsb");
    let mut found = Vec::new();
    for name in CORPUS {
        found.push((name.to_owned(), bench.join(name)));
    }
    for name in OPTIONAL {
        found.push((name.to_owned(), extra.join(name)));
    }
    found
        .into_iter()
        .filter_map(|(name, path)| {
            let bytes = std::fs::read(path).ok()?;
            let structure = read(bytes, &name);
            Some((name, structure))
        })
        .collect()
}

fn auth() -> PdbOptions {
    PdbOptions::new().namespace(PdbIdentifierNamespace::Auth)
}

fn to_pdb(name: &str, structure: &Structure) -> (Structure, Vec<molframe::Diagnostic>) {
    let (text, findings) = match molframe::write_pdb_with_findings(structure, &auth()) {
        Ok(written) => written,
        Err(failure) => panic!("{name}: pdb write {failure:?}"),
    };
    (read(text.into_bytes(), "again.pdb"), findings.to_vec())
}

fn via_cif(name: &str, structure: &Structure) -> Structure {
    match molframe::write_mmcif(structure) {
        Ok(text) => read(text.into_bytes(), "again.cif"),
        Err(failure) => panic!("{name}: cif write {failure:?}"),
    }
}

fn via_bcif(name: &str, structure: &Structure) -> Structure {
    match molframe::write_bcif(structure) {
        Ok(bytes) => read(bytes, "again.bcif"),
        Err(failure) => panic!("{name}: bcif write {failure:?}"),
    }
}

/// What one atom row carries. `naming` adds the label namespace that PDB
/// cannot hold.
fn atom_rows(structure: &Structure, naming: bool) -> Vec<String> {
    let mut rows = Vec::new();
    for model in structure.engine().data().models() {
        for chain in model.chains() {
            for residue in chain.residues() {
                for atom in residue.atoms() {
                    let label = if naming {
                        format!("{:?}|{:?}", chain.label(), residue.label_seq_id())
                    } else {
                        String::new()
                    };
                    rows.push(format!(
                        "{label}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}",
                        chain.auth_label(),
                        residue.auth_seq_id(),
                        residue.ins_code().filter(|c| !c.is_empty()),
                        residue.auth_name(),
                        residue.is_het(),
                        atom.auth_name(),
                        atom.alt_label().filter(|c| !c.is_empty()),
                        atom.occupancy(),
                        atom.b_factor(),
                        atom.formal_charge().unwrap_or(0),
                        atom.element(),
                    ));
                }
            }
        }
    }
    rows
}

/// Positions of every model, in atom order.
fn model_positions(structure: &Structure) -> Vec<Vec<[f32; 3]>> {
    let data = structure.engine().data();
    if let Some(ragged) = structure.ragged_models() {
        return ragged.iter().map(single_positions).collect();
    }
    (0..structure.model_count())
        .map(|m| {
            data.coords
                .block(molframe::ModelIndex::new(u32::try_from(m).unwrap_or(0)))
                .map(|block| block.as_slice().to_vec())
                .unwrap_or_default()
        })
        .collect()
}

fn single_positions(structure: &Structure) -> Vec<[f32; 3]> {
    (0..structure.atom_count())
        .filter_map(|i| {
            structure
                .engine()
                .atom(molframe::AtomIndex::new(i))?
                .position()
        })
        .collect()
}

fn close(left: &[Vec<[f32; 3]>], right: &[Vec<[f32; 3]>], tolerance: f32) -> bool {
    left.len() == right.len()
        && left.iter().zip(right).all(|(a, b)| {
            a.len() == b.len()
                && a.iter()
                    .zip(b)
                    .all(|(p, q)| (0..3).all(|k| (p[k] - q[k]).abs() <= tolerance))
        })
}

/// Bonds as unordered pairs of atom identities.
fn bond_set(structure: &Structure) -> BTreeSet<String> {
    let data = structure.engine().data();
    let mut names = std::collections::HashMap::new();
    for chain in data.chains() {
        for residue in chain.residues() {
            for atom in residue.atoms() {
                names.insert(
                    atom.index(),
                    format!(
                        "{:?}.{:?}{:?}.{:?}.{:?}",
                        chain.auth_label(),
                        residue.auth_seq_id(),
                        residue.ins_code().filter(|c| !c.is_empty()),
                        atom.auth_name(),
                        atom.alt_label().filter(|c| !c.is_empty()),
                    ),
                );
            }
        }
    }
    data.bonds
        .iter()
        .filter_map(|bond| {
            let a = names.get(&bond.atom_a)?.clone();
            let b = names.get(&bond.atom_b)?.clone();
            let (low, high) = if a <= b { (a, b) } else { (b, a) };
            Some(format!("{low}~{high}"))
        })
        .collect()
}

/// The space group as the file spells it, wherever the reader keeps it.
fn space_group(structure: &Structure) -> Option<String> {
    let engine = structure.engine();
    engine
        .data()
        .entry
        .space_group
        .as_deref()
        .map(str::to_owned)
        .or_else(|| {
            engine
                .extensions()
                .get::<molframe::crystal::SymmetrySet>(molframe::crystal::SYMMETRY_EXTENSION)?
                .hermann_mauguin
                .as_deref()
                .map(str::to_owned)
        })
}

fn header(structure: &Structure) -> String {
    let data = structure.engine().data();
    format!(
        "{:?}|{:?}|{:?}|{:?}|{:?}",
        data.entry.title,
        data.entry.resolution,
        space_group(structure),
        data.cell,
        data.entry.id
    )
}

fn count_where(structure: &Structure, test: impl Fn(&molframe::AtomRef<'_>) -> bool) -> usize {
    (0..structure.atom_count())
        .filter_map(|i| structure.engine().atom(molframe::AtomIndex::new(i)))
        .filter(|a| test(a))
        .count()
}

#[test]
fn the_corpus_exercises_every_audited_feature() {
    // The downloads carry the alternate locations and insertion codes the
    // committed files lack; without them there is nothing to guard.
    if !PathBuf::from("/tmp/rt-data/rcsb/3nir.cif").exists() {
        return;
    }
    let files = corpus();
    let any = |test: &dyn Fn(&Structure) -> bool| files.iter().any(|(_, s)| test(s));
    assert!(any(&|s| count_where(s, |a| a
        .alt_label()
        .is_some_and(|l| !l.is_empty()))
        > 0));
    assert!(any(&|s| s
        .engine()
        .data()
        .residues()
        .any(|r| r.ins_code().is_some_and(|c| !c.is_empty()))));
    assert!(any(&|s| s.model_count() > 1));
    assert!(any(&|s| !s.engine().data().bonds.is_empty()));
    assert!(any(&|s| count_where(s, |a| a
        .formal_charge()
        .is_some_and(|c| c != 0))
        > 0));
}

#[test]
fn text_and_binary_mmcif_agree_on_every_feature() {
    for (name, source) in corpus() {
        // A source read from PDB may repeat a chain label; the first write
        // numbers the repeats, so later hops must then be exact.
        let cif = via_cif(&name, &source);
        same_rows(&name, &source, &cif, false);
        let bcif = via_bcif(&name, &cif);
        same_rows(&name, &cif, &bcif, true);
        assert!(
            close(&model_positions(&cif), &model_positions(&bcif), 0.0),
            "{name}"
        );
        assert_eq!(bond_set(&source), bond_set(&bcif), "{name}");
        assert_eq!(header(&source), header(&bcif), "{name}");
    }
}

#[test]
fn mmcif_through_pdb_keeps_what_pdb_can_hold() {
    for (name, source) in corpus() {
        let (back, _) = to_pdb(&name, &source);
        let again = via_cif(&name, &back);
        same_rows(&name, &source, &again, false);
        assert!(
            close(&model_positions(&source), &model_positions(&again), 5e-4),
            "{name}"
        );
        assert_eq!(source.model_count(), again.model_count(), "{name}");
    }
}

#[test]
fn pdb_through_mmcif_and_back_is_stable() {
    for (name, source) in corpus() {
        let (once, _) = to_pdb(&name, &source);
        let cif = via_cif(&name, &once);
        let (twice, _) = to_pdb(&name, &cif);
        same_rows(&name, &once, &twice, false);
        assert_eq!(bond_set(&once), bond_set(&twice), "{name}");
        assert_eq!(header(&once), header(&twice), "{name}");
    }
}

#[test]
fn bonds_survive_pdb_when_a_record_exists_and_are_reported_otherwise() {
    for (name, source) in corpus() {
        let (back, findings) = to_pdb(&name, &source);
        let lost: BTreeSet<_> = bond_set(&source)
            .difference(&bond_set(&back))
            .cloned()
            .collect();
        if !lost.is_empty() {
            let reported = findings.iter().any(|f| f.message().contains("bond"));
            assert!(reported, "{name}: {} bonds lost silently", lost.len());
        }
    }
}

#[test]
fn header_metadata_survives_pdb_or_is_reported() {
    for (name, source) in corpus() {
        let (back, findings) = to_pdb(&name, &source);
        let (left, right) = (source.engine().data(), back.engine().data());
        let mut lost = Vec::new();
        if left.entry.title != right.entry.title {
            lost.push("title");
        }
        if left.entry.resolution != right.entry.resolution {
            lost.push("resolution");
        }
        if space_group(&source) != space_group(&back) {
            lost.push("space group");
        }
        if left.cell != right.cell {
            lost.push("cell");
        }
        let short =
            |id: &Option<Box<str>>| id.as_deref().map(|i| i.chars().take(4).collect::<String>());
        if short(&left.entry.id) != short(&right.entry.id) {
            lost.push("entry id");
        }
        assert!(
            lost.is_empty() || findings.iter().any(|f| f.message().contains("metadata")),
            "{name}: {lost:?} lost without a finding"
        );
    }
}

/// Fails with the first differing atom row rather than the whole table.
fn same_rows(name: &str, left: &Structure, right: &Structure, naming: bool) {
    let (a, b) = (atom_rows(left, naming), atom_rows(right, naming));
    assert_eq!(a.len(), b.len(), "{name}: atom count");
    if let Some(at) = (0..a.len()).find(|&i| a[i] != b[i]) {
        panic!("{name}: row {at}\n  {}\n  {}", a[at], b[at]);
    }
}

const SYNTHETIC_RESIDUES: [(&str, &str, i32, &str); 5] = [
    ("CYS", "SG", 1, " "),
    ("CYS", "SG", 2, " "),
    ("ALA", "CA", 3, " "),
    ("GLY", "CA", 3, "A"),
    ("ALA", "CA", 4, " "),
];

fn atom_line(
    serial: u32,
    record: &str,
    atom: (&str, &str, &str, i32, &str),
    at: [f64; 3],
) -> String {
    let (name, alt, residue, seq, ins) = atom;
    let element = name.trim_start_matches(|c: char| c.is_ascii_digit());
    let element: String = element.chars().take(1).collect();
    format!(
        "{record:<6}{serial:>5} {name:<4}{alt}{residue:>3} A{seq:>4}{ins}   {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}          {element:>2}  ",
        at[0], at[1], at[2], 1.0, 20.5
    )
}

/// A small structure with two models, an alternate location, an insertion
/// code, a charged metal, a disulfide and a metal link.
fn synthetic_pdb() -> String {
    let mut text =
        String::from("HEADER    TEST                                    01-JAN-00   9XYZ\n");
    text.push_str("TITLE     SYNTHETIC ROUND TRIP FIXTURE\n");
    text.push_str("CRYST1   50.000   60.000   70.000  90.00  90.00  90.00 P 21 21 21    4\n");
    let _ = writeln!(
        text,
        "HELIX  {:>3} {:>3} ALA A {:>4}  ALA A {:>4}  {:>2}{:>35}{:>5}",
        1, 1, 3, 4, 1, "", 3
    );
    text.push_str("SSBOND   1 CYS A    1    CYS A    2\n");
    text.push_str("LINK        ZN    ZN A 100                 SG  CYS A   1\n");
    for model in 1..=2 {
        let _ = writeln!(text, "MODEL     {model:>4}");
        let shift = f64::from(model) * 0.5;
        let mut serial = 1;
        for (index, (residue, name, seq, ins)) in SYNTHETIC_RESIDUES.iter().enumerate() {
            let x = 3.8 * f64::from(u32::try_from(index).unwrap_or(0)) + shift;
            let alts: &[&str] = if *seq == 3 && *ins == " " {
                &["A", "B"]
            } else {
                &[" "]
            };
            for alt in alts {
                let mut line = atom_line(
                    serial,
                    "ATOM",
                    (name, alt, residue, *seq, ins),
                    [x, 1.0, 2.0],
                );
                if *alt == "A" {
                    line.replace_range(54..60, "  0.60");
                }
                if *alt == "B" {
                    line.replace_range(54..60, "  0.40");
                }
                text.push_str(&line);
                text.push('\n');
                serial += 1;
            }
        }
        text.push_str("TER\n");
        let mut zinc = atom_line(
            serial,
            "HETATM",
            ("ZN", " ", "ZN", 100, " "),
            [4.0 + shift, 1.0, 3.0],
        );
        zinc.replace_range(76..78, "ZN");
        zinc.replace_range(78..80, "2+");
        text.push_str(&zinc);
        text.push_str("\nENDMDL\n");
    }
    text.push_str("CONECT    2    7\nEND\n");
    text
}

#[test]
fn a_synthetic_structure_keeps_every_feature_in_every_direction() {
    let source = read(synthetic_pdb().into_bytes(), "synthetic.pdb");
    assert_eq!(source.model_count(), 2);
    let (pdb, _) = to_pdb("synthetic", &source);
    let cif = via_cif("synthetic", &source);
    let bcif = via_bcif("synthetic", &cif);
    let again = via_cif("synthetic", &pdb);
    for other in [&pdb, &cif, &bcif, &again] {
        same_rows("synthetic", &source, other, false);
        assert!(close(
            &model_positions(&source),
            &model_positions(other),
            5e-4
        ));
        assert_eq!(source.model_count(), other.model_count());
        assert_eq!(bond_set(&source), bond_set(other));
        assert_eq!(header(&source), header(other));
    }
    let rows = atom_rows(&source, false);
    assert!(rows.iter().any(|r| r.contains("Some(\"A\")|Some(0.6)")));
    assert!(
        rows.iter()
            .any(|r| r.contains("Some(\"A\")|Some(3)|Some(\"A\")"))
    );
    assert!(
        rows.iter()
            .any(|r| r.contains("Some(\"ZN\")") && r.contains("|2|"))
    );
    assert!(bond_set(&source).len() >= 3);
    // The helix crosses the insertion code, which a numbering that ignores
    // codes could not address.
    let helix = |s: &Structure| {
        s.engine()
            .secondary_structure()
            .iter()
            .filter(|k| k.is_helix())
            .count()
    };
    assert!(helix(&source) > 0);
    assert_eq!(helix(&source), helix(&again));
}
