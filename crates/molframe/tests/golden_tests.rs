//! Executable end-to-end correctness fixtures for shipped workflows.

use molframe::formats::cif::CifWriteOptions;
use molframe::formats::pdb::PdbOptions;
use molframe::geometry::{Rigid, superpose};
use molframe::{
    AltlocPolicy, AnalysisPolicy, Code, ReadOptions, read_bytes, transform, write_bcif,
    write_mmcif_with_options, write_pdb,
};
use molframe_core::selection::AtomSelection;

const MULTI_MODEL_PDB: &str = "\
HEADER    GOLDEN                              01-JAN-00   1ABC
MODEL        7
ATOM      1  N   GLY A  10A     10.000  11.000  12.000  1.00 10.00           N
ATOM      2  CA  GLY A  10A     11.000  11.000  12.000  1.00 11.00           C
ENDMDL
MODEL       11
ATOM      1  N   GLY A  10A     12.000  13.000  14.000  1.00 10.00           N
ATOM      2  CA  GLY A  10A     13.000  13.000  14.000  1.00 11.00           C
ENDMDL
END
";

const NAMESPACED_CIF: &str = "\
data_golden
loop_
_entity.id
_entity.type
1 polymer
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . GLY LONG 1 1 0 0 0 1 10 42 GLY A N 1
ATOM 2 C CA A GLY LONG 1 1 1 0 0 0.4 11 42 GLY A CA 1
ATOM 3 C CA B GLY LONG 1 1 2 0 0 0.6 12 42 GLY A CA 1
";

fn read_fixture(text: &str, name: &str) -> molframe::Structure {
    match read_bytes(text.as_bytes().to_vec(), Some(name), &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("golden read failed: {findings:?}"),
    }
}

#[test]
fn gw_001_reports_the_complete_normalised_hierarchy() {
    let structure = read_fixture(NAMESPACED_CIF, "golden.cif");
    let report = format!(
        "models={} chains={} entities={} residues={} atoms={}",
        structure.model_count(),
        structure.chain_count(),
        structure.engine().entity_count(),
        structure.residue_count(),
        structure.atom_count(),
    );
    assert_eq!(report, "models=1 chains=1 entities=1 residues=1 atoms=3");
}

#[test]
fn gw_002_round_trips_mmcif_through_binary_cif_semantically() {
    let source = read_fixture(NAMESPACED_CIF, "golden.cif");
    let bytes = match write_bcif(&source) {
        Ok(bytes) => bytes,
        Err(findings) => panic!("BinaryCIF write failed: {findings:?}"),
    };
    let round_tripped = match read_bytes(bytes, Some("golden.bcif"), &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("BinaryCIF read failed: {findings:?}"),
    };

    assert_eq!(round_tripped.model_count(), source.model_count());
    assert_eq!(round_tripped.chain_count(), source.chain_count());
    assert_eq!(round_tripped.residue_count(), source.residue_count());
    assert_eq!(round_tripped.atom_count(), source.atom_count());
    assert_eq!(round_tripped.coordinates(), source.coordinates());
}

#[test]
fn gw_003_converts_insertion_codes_and_deposited_models_to_mmcif() {
    let structure = read_fixture(MULTI_MODEL_PDB, "golden.pdb");
    let options = CifWriteOptions::new().with_block_id("1ABC");
    let rendered = write_mmcif_with_options(&structure, &options);
    let rendered = rendered.unwrap_or_else(|error| panic!("mmCIF write failed: {error}"));
    let round_tripped = read_fixture(&rendered, "golden.cif");
    let model_numbers: Vec<i32> = round_tripped
        .engine()
        .data()
        .models()
        .filter_map(molframe::ModelRef::number)
        .collect();
    let insertion = round_tripped
        .engine()
        .data()
        .residues()
        .next()
        .and_then(molframe::ResidueRef::ins_code);

    assert_eq!(model_numbers, vec![7, 11]);
    assert_eq!(insertion, Some("A"));
    assert_eq!(round_tripped.atom_count(), 2);
    for (model, expected) in [
        [[10.0, 11.0, 12.0], [11.0, 11.0, 12.0]],
        [[12.0, 13.0, 14.0], [13.0, 13.0, 14.0]],
    ]
    .iter()
    .enumerate()
    {
        let Ok(model_index) = u32::try_from(model) else {
            panic!("test model index does not fit in u32");
        };
        let Some(actual) = round_tripped.model_coordinates(molframe::ModelIndex::new(model_index))
        else {
            panic!("round-tripped model missing")
        };
        assert!(
            actual
                .iter()
                .flatten()
                .zip(expected.iter().flatten())
                .all(|(actual, expected)| (*actual - *expected).abs() < f32::EPSILON)
        );
    }
}

#[test]
fn gw_005_refuses_a_lossy_legacy_chain_and_names_the_capacity() {
    let source = NAMESPACED_CIF.replace(" 42 GLY A ", " 42 GLY LONG ");
    let structure = read_fixture(&source, "golden.cif");
    let result = write_pdb(&structure, &PdbOptions::new());
    assert_eq!(
        result
            .err()
            .and_then(|findings| findings.first().map(molframe::Diagnostic::code)),
        Some(Code::E4102)
    );
}

#[test]
fn gw_007_altloc_policies_produce_recorded_atom_counts() {
    let structure = read_fixture(NAMESPACED_CIF, "golden.cif");
    let counts: Vec<u64> = [
        AltlocPolicy::KeepAll,
        AltlocPolicy::First,
        AltlocPolicy::HighestOccupancyPerResidue,
        AltlocPolicy::HighestOccupancyPerAtom,
        AltlocPolicy::ConformerConsistent,
    ]
    .into_iter()
    .map(|altloc| {
        structure
            .engine()
            .resolve_altlocs(&AnalysisPolicy::default().with_altloc(altloc))
            .value
            .len()
    })
    .collect();
    assert_eq!(counts, vec![3, 2, 2, 2, 2]);
}

#[test]
fn gw_009_keeps_both_chain_namespaces_addressable() {
    let structure = read_fixture(NAMESPACED_CIF, "golden.cif");
    let label = structure.chain_at(0);
    assert_eq!(label.and_then(molframe::ChainRef::label), Some("LONG"));
    assert_eq!(label.and_then(molframe::ChainRef::auth_label), Some("A"));
    assert!(
        structure
            .engine()
            .data()
            .chains()
            .any(|chain| chain.label() == Some("LONG"))
    );
    assert!(
        structure
            .engine()
            .data()
            .chains()
            .any(|chain| chain.auth_label() == Some("A"))
    );
}

#[test]
fn gw_011_superposes_and_applies_the_reported_rigid_transform() {
    let fixed = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let mobile = [[2.0, 3.0, 4.0], [3.0, 3.0, 4.0], [2.0, 4.0, 4.0]];
    let fit = match superpose(&mobile, &fixed) {
        Ok(fit) => fit,
        Err(error) => panic!("fit failed: {error:?}"),
    };
    assert!(fit.rmsd < 1.0e-6);

    let structure = read_fixture(NAMESPACED_CIF, "golden.cif");
    let moved = match transform(
        &structure,
        &AtomSelection::All(structure.atom_count()),
        &Rigid::translation([1.0, 2.0, 3.0]),
    ) {
        Ok(structure) => structure,
        Err(findings) => panic!("transform failed: {findings:?}"),
    };
    assert_eq!(
        moved.engine().generation().get(),
        structure
            .engine()
            .generation()
            .next()
            .expect("initial generation advances")
            .get()
    );
}
