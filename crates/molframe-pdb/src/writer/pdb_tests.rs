use super::*;
use crate::read;
use molframe_core::diagnostic::Code;
use molframe_core::index::AtomIndex;
use molframe_core::io::{InputBuffer, ReadOptions};
use proptest::prelude::*;
use std::fmt::Write as _;
use std::io::{self, Write};

const DIPEPTIDE: &str = "\
ATOM      1  N   GLY A   1      27.340  24.430   2.614  1.00 10.00           N
ATOM      2  CA  GLY A   1      26.266  25.413   2.842  1.00 11.00           C
ATOM      3  N   ASN A   2      26.335  27.770   3.258  1.00 14.00           N
END
";

fn parse(text: &str) -> Structure {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    match read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("read failed: {findings:?}"),
    }
}

fn written(text: &str) -> String {
    match write(&parse(text), &PdbOptions::new()) {
        Ok(out) => out,
        Err(refusals) => panic!("write refused: {refusals:?}"),
    }
}

#[test]
fn an_anisou_record_parses_into_the_atom_it_precedes() {
    let text = concat!(
        "ATOM      1  N   GLY A   1      27.340  24.430   2.614  1.00 10.00           N\n",
        "ANISOU    1  N   GLY A   1   100000 200000 300000  40000  50000  60000       N\n",
        "END\n",
    );
    let structure = parse(text);
    let anisotropy = &structure.data().anisotropy;
    assert!(anisotropy.is_available());
    assert_eq!(anisotropy.len(), 1);
    let Some(record) = anisotropy.get(molframe_core::AnisotropyIndex::new(0)) else {
        panic!("the parsed ellipsoid was lost");
    };
    assert_eq!(record.atom, AtomIndex::new(0));
    for (stored, expected) in record.u.iter().zip([10.0, 20.0, 30.0, 4.0, 5.0, 6.0]) {
        assert!((stored - expected).abs() <= f32::EPSILON);
    }
}

#[test]
fn an_anisou_serial_its_atom_never_carried_is_ignored() {
    let text = concat!(
        "ATOM      1  N   GLY A   1      27.340  24.430   2.614  1.00 10.00           N\n",
        "ANISOU    7  N   GLY A   1   100000 200000 300000  40000  50000  60000       N\n",
        "END\n",
    );
    let structure = parse(text);
    let anisotropy = &structure.data().anisotropy;
    assert_eq!(anisotropy.len(), 0);
}

#[test]
fn an_anisou_record_round_trips_through_the_writer_at_stated_precision() {
    let text = concat!(
        "ATOM      1  N   GLY A   1      27.340  24.430   2.614  1.00 10.00           N\n",
        "ANISOU    1  N   GLY A   1   123456 234567 345678  45678  56789  67890       N\n",
        "END\n",
    );
    let out = written(text);
    let round = parse(&out);
    let anisotropy = &round.data().anisotropy;
    let record = anisotropy
        .get(molframe_core::AnisotropyIndex::new(0))
        .expect("round-tripped");
    // Each field was an integer in units of 0.0001, so dividing and multiplying
    // by 10 000 lands on the same six values within float rounding.
    let expected = [12.3456, 23.4567, 34.5678, 4.5678, 5.6789, 6.789];
    for (position, value) in expected.iter().enumerate() {
        assert!((record.u[position] - *value).abs() < 1e-5);
    }
    // And the written record itself carries 7-wide integer columns.
    let anisou = out
        .lines()
        .find(|line| line.starts_with("ANISOU"))
        .expect("the writer emitted an ANISOU record");
    assert_eq!(anisou[27..35].trim(), "123456");
    assert!(anisou.len() >= 70);
}

#[test]
fn a_written_file_reads_back_to_the_same_structure() {
    let original = parse(DIPEPTIDE);
    let round_tripped = parse(&written(DIPEPTIDE));

    assert_eq!(round_tripped.atom_count(), original.atom_count());
    assert_eq!(round_tripped.residue_count(), original.residue_count());
    assert_eq!(round_tripped.chain_count(), original.chain_count());
}

#[test]
fn positions_survive_a_round_trip_to_the_precision_the_format_records() {
    let original = parse(DIPEPTIDE);
    let round_tripped = parse(&written(DIPEPTIDE));

    for (before, after) in original.positions().iter().zip(round_tripped.positions()) {
        for axis in 0..3 {
            assert!(
                (before[axis] - after[axis]).abs() < 1e-3,
                "{before:?} became {after:?}",
            );
        }
    }
}

#[test]
fn names_and_elements_land_in_the_columns_a_reader_expects() {
    let out = written(DIPEPTIDE);
    let Some(first) = out.lines().find(|line| line.starts_with("ATOM")) else {
        panic!("expected an atom record")
    };
    assert_eq!(crate::fixed::text(first, 13, 16), "N");
    assert_eq!(crate::fixed::text(first, 18, 20), "GLY");
    assert_eq!(crate::fixed::text(first, 22, 22), "A");
    assert_eq!(crate::fixed::text(first, 23, 26), "1");
    assert_eq!(crate::fixed::text(first, 77, 78), "N");
}

#[test]
fn namespace_selects_one_complete_identifier_family() {
    const CIF: &str = "data_ids\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.auth_atom_id\n_atom_site.label_comp_id\n\
_atom_site.auth_comp_id\n_atom_site.label_asym_id\n_atom_site.auth_asym_id\n\
_atom_site.label_seq_id\n_atom_site.auth_seq_id\n_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n\
ATOM 1 C LC AC LIG DRG L A 7 42 0 0 0 0.75 12.5\n";
    let input = InputBuffer::from_bytes(CIF.as_bytes().to_vec());
    let Ok((structure, _)) = molframe_cif::read(&input, &ReadOptions::new()) else {
        panic!("namespace fixture should parse");
    };
    let Ok(label) = write(
        &structure,
        &PdbOptions::new().namespace(PdbIdentifierNamespace::Label),
    ) else {
        panic!("label namespace is complete");
    };
    let Ok(auth) = write(
        &structure,
        &PdbOptions::new().namespace(PdbIdentifierNamespace::Auth),
    ) else {
        panic!("auth namespace is complete");
    };
    let Some(label_atom) = label.lines().find(|line| line.starts_with("ATOM")) else {
        panic!("label atom line absent");
    };
    let Some(auth_atom) = auth.lines().find(|line| line.starts_with("ATOM")) else {
        panic!("auth atom line absent");
    };
    assert_eq!(crate::fixed::text(label_atom, 13, 16), "LC");
    assert_eq!(crate::fixed::text(label_atom, 18, 20), "LIG");
    assert_eq!(crate::fixed::text(label_atom, 22, 22), "L");
    assert_eq!(crate::fixed::text(label_atom, 23, 26), "7");
    assert_eq!(crate::fixed::text(auth_atom, 13, 16), "AC");
    assert_eq!(crate::fixed::text(auth_atom, 18, 20), "DRG");
    assert_eq!(crate::fixed::text(auth_atom, 22, 22), "A");
    assert_eq!(crate::fixed::text(auth_atom, 23, 26), "42");
}

#[test]
fn missing_occupancy_and_b_factor_are_refused_before_output() {
    const MISSING: &str =
        "ATOM      1  N   GLY A   1      27.340  24.430   2.614                      N\nEND\n";
    let structure = parse(MISSING);
    let Err(findings) = write(&structure, &PdbOptions::new()) else {
        panic!("missing mandatory numeric fields must be refused");
    };
    let contexts = findings
        .iter()
        .flat_map(|finding| finding.context().iter())
        .map(molframe_core::ContextItem::value)
        .collect::<Vec<_>>();
    assert!(contexts.contains(&"occupancy"));
    assert!(contexts.contains(&"B factor"));
}

#[test]
fn a_chain_label_the_format_cannot_hold_is_refused_rather_than_truncated() {
    let mut data = parse(DIPEPTIDE).data().clone();
    let Ok(long) = data.dictionary.intern("AA") else {
        panic!("interning failed")
    };
    data.topology.chains = rename_first_chain(&data.topology.chains, long);
    let structure = Structure::new(data);

    let refused = write(&structure, &PdbOptions::new());
    let codes: Vec<_> = refused
        .err()
        .unwrap_or_default()
        .iter()
        .map(molframe_core::diagnostic::Diagnostic::code)
        .collect();
    assert!(
        codes.contains(&Code::E4102),
        "expected a refusal, got {codes:?}"
    );
}

#[test]
fn an_explicit_chain_mapping_is_the_way_past_a_refusal() {
    let mut data = parse(DIPEPTIDE).data().clone();
    let Ok(long) = data.dictionary.intern("AA") else {
        panic!("interning failed")
    };
    data.topology.chains = rename_first_chain(&data.topology.chains, long);
    let structure = Structure::new(data);

    let options = PdbOptions::new().chain_map("AA", "B");
    let out = write(&structure, &options);
    assert!(out.is_ok(), "a named replacement should be accepted");
    assert!(out.unwrap_or_default().contains(" B "));
}

#[test]
fn a_coordinate_too_large_for_its_field_is_refused() {
    // The x field is columns 31 to 38; this one holds a value the format
    // cannot represent to three decimal places.
    let text = concat!(
        "ATOM      1  N   GLY A   1    99999.00   1.000   1.000  1.00  0.00           N\n",
        "END\n",
    );
    let structure = parse(text);
    let refused = write(&structure, &PdbOptions::new());
    let codes: Vec<_> = refused
        .err()
        .unwrap_or_default()
        .iter()
        .map(molframe_core::diagnostic::Diagnostic::code)
        .collect();
    assert!(
        codes.contains(&Code::E4104),
        "expected a refusal, got {codes:?}"
    );
}

#[test]
fn a_filter_writes_only_what_it_accepts() {
    struct FirstAtomOnly;
    impl Select for FirstAtomOnly {
        fn accept_atom(&self, atom: AtomIndex) -> bool {
            atom.get() == 0
        }
    }

    let structure = parse(DIPEPTIDE);
    let out = write_selected(&structure, &PdbOptions::new(), &FirstAtomOnly);
    let atoms = out
        .unwrap_or_default()
        .lines()
        .filter(|line| line.starts_with("ATOM"))
        .count();
    assert_eq!(atoms, 1);
}

#[test]
fn a_written_file_ends_with_the_record_that_says_so() {
    assert!(written(DIPEPTIDE).ends_with("END\n"));
}

#[derive(Default)]
struct CountingWriter {
    bytes: usize,
    writes: usize,
    largest_write: usize,
}

impl Write for CountingWriter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.bytes += bytes.len();
        self.writes += 1;
        self.largest_write = self.largest_write.max(bytes.len());
        Ok(bytes.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[test]
fn pdb_writer_emits_records_without_a_full_file_buffer() {
    let structure = parse(DIPEPTIDE);
    let expected = written(DIPEPTIDE);
    let mut writer = CountingWriter::default();
    if let Err(findings) = crate::write_to(&structure, &PdbOptions::new(), &mut writer) {
        panic!("stream write failed: {findings:?}");
    }
    assert_eq!(writer.bytes, expected.len());
    assert!(writer.writes > 1);
    assert!(writer.largest_write < expected.len());
}

proptest! {
    #[test]
    fn generated_representable_structures_survive_a_pdb_round_trip(
        positions in prop::collection::vec(
            (-999.0_f32..9_999.0, -999.0_f32..9_999.0, -999.0_f32..9_999.0),
            1..64,
        )
    ) {
        let source = generated_pdb(&positions);
        let original = parse(&source);
        let round_tripped = parse(&written(&source));
        prop_assert_eq!(round_tripped.atom_count(), original.atom_count());
        prop_assert_eq!(round_tripped.residue_count(), original.residue_count());
        for (before, after) in original.positions().iter().zip(round_tripped.positions()) {
            for axis in 0..3 {
                // PDB coordinates carry three decimals. Half one unit in the
                // last place plus f32 conversion error bounds the round trip.
                prop_assert!((before[axis] - after[axis]).abs() <= 5.1e-4);
            }
        }
    }
}

#[test]
fn a_negative_coordinate_that_loses_a_digit_to_the_sign_is_refused() {
    let original = parse(DIPEPTIDE);
    let context = molframe_core::ExecutionContext::default();
    let mut editor = original.edit_coordinates(&context).expect("edit fits");
    let Some(positions) = editor.positions_mut(molframe_core::index::ModelIndex::new(0)) else {
        panic!("first model missing")
    };
    positions[0][0] = -1_226.086;
    let edited = match editor.commit() {
        Ok(structure) => structure,
        Err(findings) => panic!("coordinate edit failed: {findings:?}"),
    };
    let codes: Vec<_> = write(&edited, &PdbOptions::new())
        .err()
        .unwrap_or_default()
        .iter()
        .map(molframe_core::diagnostic::Diagnostic::code)
        .collect();
    assert!(
        codes.contains(&Code::E4104),
        "expected E4104, got {codes:?}"
    );
}

fn generated_pdb(positions: &[(f32, f32, f32)]) -> String {
    let mut source = String::new();
    for (index, (x, y, z)) in positions.iter().enumerate() {
        let serial = index + 1;
        writeln!(
            source,
            "ATOM  {serial:>5} C{serial:<3} GLY A   1    \
             {x:>8.3}{y:>8.3}{z:>8.3}{occ:>6.2}{b:>6.2}           C  ",
            occ = 1.0,
            b = 10.0,
        )
        .expect("writing to a string is infallible");
    }
    source.push_str("END\n");
    source
}

/// Rebuilds a chain table with the first chain relabelled.
fn rename_first_chain(
    chains: &molframe_core::topology::ChainTable,
    label: molframe_core::symbol::SymbolId,
) -> molframe_core::topology::ChainTable {
    use molframe_core::optional::OptionalSymbol;
    use molframe_core::topology::{ChainRecord, ChainTable, PolymerKind};

    let mut rebuilt = ChainTable::default();
    for chain in chains.iter() {
        let Some(residues) = chains.residues(chain) else {
            continue;
        };
        let Some(entity) = chains.entity(chain) else {
            continue;
        };
        let relabelled = chain.get() == 0;
        rebuilt
            .push(
                ChainRecord {
                    label_asym_id: if relabelled {
                        label
                    } else {
                        match chains.label_asym_id(chain) {
                            Some(existing) => existing,
                            None => label,
                        }
                    },
                    auth_asym_id: OptionalSymbol::some(if relabelled {
                        label
                    } else {
                        match chains.auth_asym_id(chain) {
                            Some(existing) => existing,
                            None => label,
                        }
                    }),
                    entity,
                    polymer_kind: match chains.polymer_kind(chain) {
                        Some(kind) => kind,
                        None => PolymerKind::None,
                    },
                },
                residues,
            )
            .expect("small chain table");
    }
    rebuilt
}

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
    assert_eq!(
        again.data().bonds.len(),
        0,
        "LINK and SSBOND are headers, not CONECT"
    );
}
