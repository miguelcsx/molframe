use super::*;
use super::{parse, written};
use crate::read;
use molframe_core::diagnostic::Code;
use molframe_core::index::AtomIndex;
use molframe_core::io::{InputBuffer, ReadOptions};

const TWO_MODELS: &str = "\
MODEL        1
ATOM      1  N   GLY A   1      27.340  24.430   2.614  1.00 10.00           N
ATOM      2  CA  GLY A   1      26.266  25.413   2.842  1.00 11.00           C
ENDMDL
MODEL        2
ATOM      1  N   GLY A   1      28.340  24.430   2.614  1.00 10.00           N
ATOM      2  CA  GLY A   1      27.266  25.413   2.842  1.00 11.00           C
ENDMDL
END
";

#[test]
fn every_model_is_written_between_model_and_endmdl() {
    let text = written(TWO_MODELS);
    assert_eq!(text.lines().filter(|l| l.starts_with("MODEL")).count(), 2);
    assert_eq!(text.lines().filter(|l| l.starts_with("ENDMDL")).count(), 2);
    let again = parse(&text);
    assert_eq!(again.model_count(), 2);
    let second = again.model_positions(molframe_core::index::ModelIndex::new(1));
    assert_eq!(second.map(|p| p[0]), Some([28.340, 24.430, 2.614]));
}

#[test]
fn a_later_model_coordinate_that_overflows_its_field_is_refused() {
    let context = molframe_core::ExecutionContext::builder()
        .scratch_policy(molframe_core::ScratchPolicy::new(0))
        .build()
        .expect("valid context");
    let mut editor = parse(TWO_MODELS)
        .edit_coordinates(&context)
        .expect("edit fits");
    editor
        .positions_mut(molframe_core::index::ModelIndex::new(1))
        .expect("dense frame")[0] = [123_456.0, 0.0, 0.0];
    let structure = editor.commit().expect("commit");
    let result = write(&structure, &PdbOptions::new());
    let Err(refusals) = result else {
        panic!("a model that cannot be written must be refused, not truncated")
    };
    assert!(refusals.iter().any(|d| d.code() == Code::E4104));
}

#[test]
fn models_with_different_atoms_are_each_written() {
    let text = "\
MODEL        1
ATOM      1  N   GLY A   1      1.000   1.000   1.000  1.00  0.00           N
ATOM      2  CA  GLY A   1      2.000   2.000   2.000  1.00  0.00           C
ENDMDL
MODEL        2
ATOM      1  N   GLY A   1      1.500   1.500   1.500  1.00  0.00           N
ATOM      2  O   GLY A   1      2.500   2.500   2.500  1.00  0.00           O
ENDMDL
END
";
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let options = ReadOptions::new().mode(molframe_core::io::ParseMode::Recover);
    let Ok((structure, _)) = read(&input, &options) else {
        panic!("ragged read")
    };
    assert!(structure.is_ragged_ensemble());
    let out = write(&structure, &PdbOptions::new()).expect("ragged write");
    assert_eq!(out.lines().filter(|l| l.starts_with("ENDMDL")).count(), 2);
    assert!(
        out.lines()
            .any(|l| l.starts_with("ATOM") && l.contains(" O   GLY"))
    );
}

const RAGGED: &str = "\
MODEL        1
ATOM      1  N   GLY A   1      1.000   1.000   1.000  1.00  0.00           N
ATOM      2  CA  GLY A   1      2.000   2.000   2.000  1.00  0.00           C
ENDMDL
MODEL        2
ATOM      1  N   GLY A   1      1.500   1.500   1.500  1.00  0.00           N
ATOM      2  O   GLY A   1      2.500   2.500   2.500  1.00  0.00           O
ENDMDL
END
";

fn ragged() -> Structure {
    let input = InputBuffer::from_bytes(RAGGED.as_bytes().to_vec());
    let options = ReadOptions::new().mode(molframe_core::io::ParseMode::Recover);
    let Ok((structure, _)) = read(&input, &options) else {
        panic!("ragged read")
    };
    assert!(structure.is_ragged_ensemble());
    structure
}

#[test]
fn an_atom_filter_is_refused_on_a_ragged_ensemble_instead_of_misapplied() {
    struct FirstAtomOnly;
    impl Select for FirstAtomOnly {
        fn accept_atom(&self, atom: AtomIndex) -> bool {
            atom.get() == 0
        }
    }
    let result = write_selected(&ragged(), &PdbOptions::new(), &FirstAtomOnly);
    let Err(refusals) = result else {
        panic!("a model-0 index filter must not be applied to other models")
    };
    assert!(refusals.iter().any(|d| d.code() == Code::E6003));
}

#[test]
fn a_model_filter_selects_whole_ragged_models() {
    struct SecondModel;
    impl Select for SecondModel {
        fn accept_model(&self, model: molframe_core::index::ModelIndex) -> bool {
            model.get() == 1
        }
    }
    let out = write_selected(&ragged(), &PdbOptions::new(), &SecondModel).expect("write");
    assert_eq!(out.lines().filter(|l| l.starts_with("ENDMDL")).count(), 1);
    assert!(
        out.lines()
            .any(|l| l.starts_with("ATOM") && l.contains(" O   GLY"))
    );
    assert!(
        !out.lines()
            .any(|l| l.starts_with("ATOM") && l.contains(" CA  GLY"))
    );
}

#[test]
fn conect_records_are_written_for_bare_bonds_in_both_directions() {
    let source = "\
HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C
HETATM    2  C2  LIG A   1       1.500   0.000   0.000  1.00  0.00           C
HETATM    3  O1  LIG A   1       2.500   0.000   0.000  1.00  0.00           O
CONECT    1    2
CONECT    2    1    3
CONECT    3    2
END
";
    let out = written(source);
    let conect: Vec<&str> = out.lines().filter(|l| l.starts_with("CONECT")).collect();
    assert_eq!(
        conect,
        [
            "CONECT    1    2",
            "CONECT    2    1    3",
            "CONECT    3    2"
        ]
    );
}

#[test]
fn disulfide_and_metal_bonds_become_ssbond_and_link_records() {
    use molframe_core::bond::{BondOrder, BondProvenance, BondRecord, BondTableBuilder};
    let source = "\
ATOM      1  SG  CYS A   3       0.000   0.000   0.000  1.00  0.00           S
ATOM      2  SG  CYS A  10       2.050   0.000   0.000  1.00  0.00           S
HETATM    3 ZN    ZN A 101       0.000   2.300   0.000  1.00  0.00          ZN
END
";
    let mut data = parse(source).data().clone();
    let disulf = data.dictionary.intern("disulf").expect("intern");
    let metalc = data.dictionary.intern("metalc").expect("intern");
    let mut bonds = BondTableBuilder::new();
    for (a, b, kind) in [(0, 1, disulf), (0, 2, metalc)] {
        bonds.push(BondRecord {
            atom_a: AtomIndex::new(a),
            atom_b: AtomIndex::new(b),
            order: BondOrder::Single,
            provenance: BondProvenance::File,
            kind: molframe_core::optional::OptionalSymbol::some(kind),
        });
    }
    data.bonds = bonds.finish();
    let out = write(&Structure::new(data), &PdbOptions::new()).expect("write");
    let ss = out
        .lines()
        .find(|l| l.starts_with("SSBOND"))
        .expect("SSBOND");
    assert_eq!(&ss[11..14], "CYS");
    assert_eq!(ss[17..21].trim(), "3");
    assert_eq!(ss[31..35].trim(), "10");
    assert_eq!(ss[73..78].trim(), "2.05");
    let link = out.lines().find(|l| l.starts_with("LINK")).expect("LINK");
    assert_eq!(link[12..16].trim(), "SG");
    assert_eq!(link[42..46].trim(), "ZN");
    assert_eq!(link[73..78].trim(), "2.30");
    let again = parse(&out);
    let kinds: Vec<_> = again
        .data()
        .bonds
        .iter()
        .filter_map(|bond| bond.kind.get().and_then(|kind| again.resolve(kind)))
        .collect();
    assert_eq!(
        kinds,
        ["disulf", "metalc"],
        "records read back as typed bonds"
    );
    assert!(!out.contains("CONECT"), "typed bonds are not CONECT");
}
