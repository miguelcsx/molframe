use super::*;
use molframe_bench::{Sample, structure as sample};
use molframe_core::io::{InputBuffer, ReadOptions};

const HEADER: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

fn build(rows: &str) -> Result<Structure, String> {
    let text = format!("{HEADER}{rows}");
    let input = InputBuffer::from_bytes(text.into_bytes());
    molframe_cif::read(&input, &ReadOptions::new())
        .map(|(structure, _)| structure)
        .map_err(|findings| format!("{findings:?}"))
}

const N_H: &str = "ATOM 1 N N GLY A 1 0.0 0.0 0.0\nATOM 2 H H GLY A 1 1.0 0.0 0.0\n";

#[test]
fn polar_hydrogen_against_oxygen_is_a_hydrogen_bond() -> Result<(), Box<dyn std::error::Error>> {
    let structure = build(&format!("{N_H}HETATM 3 O O HOH B 2 2.8 0.0 0.0\n"))?;
    let result = clash_score(
        &structure,
        ClashScoreOptions::default(),
        &ExecutionContext::default(),
    )?;
    assert!(result.hydrogen_aware);
    assert_eq!(result.clashes, 0);
    assert_eq!(result.atoms, 3);
    Ok(())
}

#[test]
fn hydrogen_against_carbon_is_a_clash_only_when_hydrogens_count()
-> Result<(), Box<dyn std::error::Error>> {
    let structure = build(&format!("{N_H}ATOM 3 C CA ALA B 2 2.9 0.0 0.0\n"))?;
    let context = ExecutionContext::default();
    let aware = clash_score(&structure, ClashScoreOptions::default(), &context)?;
    assert_eq!(aware.clashes, 1);
    assert_eq!(aware.residue_clashes, vec![1, 1]);
    assert!((aware.score - 1000.0 / 3.0).abs() < 1e-9);
    let heavy = clash_score(
        &structure,
        ClashScoreOptions {
            hydrogens: HydrogenMode::Ignore,
            ..ClashScoreOptions::default()
        },
        &context,
    )?;
    assert!(!heavy.hydrogen_aware);
    assert_eq!((heavy.clashes, heavy.atoms), (0, 2));
    Ok(())
}

#[test]
fn heavy_atom_score_separates_clean_and_crowded_models() -> Result<(), BuiltinError> {
    let context = ExecutionContext::default();
    let ubiquitin = clash_score(
        &sample(Sample::Small),
        ClashScoreOptions::default(),
        &context,
    )?;
    let haemoglobin = clash_score(
        &sample(Sample::Medium),
        ClashScoreOptions::default(),
        &context,
    )?;
    assert!(!ubiquitin.hydrogen_aware);
    assert!(ubiquitin.score < 20.0, "{}", ubiquitin.score);
    assert!(
        haemoglobin.score > ubiquitin.score,
        "{} vs {}",
        haemoglobin.score,
        ubiquitin.score
    );
    Ok(())
}
