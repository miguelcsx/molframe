//! A default write keeps what the input recorded, or refuses; it never drops
//! it and exits cleanly.

use molframe::ReadOptions;
use molframe::crystal::AssemblyExt;

const UBIQUITIN_PDB: &[u8] = include_bytes!("../../molframe-bench/data/1ubq.pdb");
const FOUR_HHB_BCIF: &[u8] = include_bytes!("../../molframe-bench/data/4hhb.bcif");
const NMR_ENSEMBLE: &[u8] = include_bytes!("../../molframe-bench/data/2m7c.pdb");

fn read(bytes: &[u8], name: &str) -> molframe::Structure {
    match molframe::read_bytes(bytes.to_vec(), Some(name), &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("{name} must read: {findings:?}"),
    }
}

fn assemblies(structure: &molframe::Structure) -> usize {
    structure
        .assembly_set()
        .map_or(0, molframe::crystal::AssemblySet::len)
}

#[test]
fn default_mmcif_keeps_space_group_resolution_assemblies_and_ncs() {
    let text = molframe::write_mmcif(&read(FOUR_HHB_BCIF, "4hhb.bcif")).expect("default write");
    assert!(
        text.contains("_symmetry.space_group_name_H-M"),
        "space group"
    );
    assert!(text.contains("_refine.ls_d_res_high"), "resolution");
    assert!(
        text.contains("_refine.pdbx_refine_id"),
        "resolution needs a refine id"
    );
    assert!(
        text.contains("_pdbx_struct_assembly.oligomeric_details"),
        "assemblies"
    );
    assert!(text.contains("_struct_ncs_oper"), "NCS operators");
    assert!(
        text.contains("_struct_conf.beg_auth_seq_id"),
        "secondary structure"
    );
    assert!(text.contains("_struct_conn.id"), "connections");
    let source = read(FOUR_HHB_BCIF, "4hhb.bcif");
    let again = read(text.as_bytes(), "again.cif");
    assert_eq!(assemblies(&again), assemblies(&source));
    assert!(assemblies(&again) > 0);
}

#[test]
fn default_bcif_keeps_the_same_metadata() {
    let bytes = molframe::write_bcif(&read(FOUR_HHB_BCIF, "4hhb.bcif")).expect("default write");
    let again = read(&bytes, "again.bcif");
    assert_eq!(
        assemblies(&again),
        assemblies(&read(FOUR_HHB_BCIF, "4hhb.bcif"))
    );
    assert!(assemblies(&again) > 0);
}

#[test]
fn a_pdb_read_structure_writes_mmcif_without_an_entry_id() {
    let text = molframe::write_mmcif(&read(UBIQUITIN_PDB, "1ubq.pdb")).expect("default write");
    assert!(text.starts_with("data_"));
}

#[test]
fn cif_to_pdb_keeps_cell_space_group_links_and_operators() {
    let structure = read(FOUR_HHB_BCIF, "4hhb.bcif");
    let text = molframe::write_pdb(
        &structure,
        &molframe::formats::pdb::PdbOptions::new()
            .namespace(molframe::formats::pdb::PdbIdentifierNamespace::Auth),
    )
    .expect("auth write");
    let cryst = text
        .lines()
        .find(|line| line.starts_with("CRYST1"))
        .expect("CRYST1");
    assert!(cryst.contains("P 1 21 1"), "{cryst}");
    assert_eq!(
        text.lines().filter(|line| line.starts_with("LINK")).count(),
        4
    );
    assert!(text.lines().any(|line| line.starts_with("SCALE1")));
    assert!(text.lines().any(|line| line.starts_with("MTRIX1")));
    assert!(
        text.lines()
            .any(|line| line.contains("RESOLUTION.    1.74"))
    );
}

#[test]
fn every_model_of_an_ensemble_is_written() {
    let structure = read(NMR_ENSEMBLE, "2m7c.pdb");
    let models = structure.model_count();
    assert!(models > 1);
    let text = molframe::write_pdb(&structure, &molframe::formats::pdb::PdbOptions::new())
        .or_else(|_| {
            molframe::write_pdb(
                &structure,
                &molframe::formats::pdb::PdbOptions::new()
                    .namespace(molframe::formats::pdb::PdbIdentifierNamespace::Auth),
            )
        })
        .expect("ensemble write");
    assert_eq!(
        text.lines().filter(|l| l.starts_with("MODEL")).count(),
        models
    );
    assert_eq!(
        text.lines().filter(|l| l.starts_with("ENDMDL")).count(),
        models
    );
    let again = read(text.as_bytes(), "again.pdb");
    assert_eq!(again.model_count(), models);
}

#[test]
fn a_refused_write_leaves_no_empty_file() {
    // 4HHB's haem has no label sequence number, so the label namespace refuses.
    let structure = read(FOUR_HHB_BCIF, "4hhb.bcif");
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("out.pdb");
    let options = molframe::WriteOptions::canonical().with_pdb(
        molframe::formats::pdb::PdbOptions::new()
            .namespace(molframe::formats::pdb::PdbIdentifierNamespace::Label),
    );
    let result = molframe::write_with_options(&path, &structure, &options);
    assert!(result.is_err(), "the label namespace cannot hold haem");
    assert!(
        !path.exists(),
        "a failed write must not leave a file behind"
    );
}
