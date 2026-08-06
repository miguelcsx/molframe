use super::*;
use crate::read;
use pdbiox_core::index::AtomIndex;
use pdbiox_core::io::{InputBuffer, ReadOptions};

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
        .map(pdbiox_core::diagnostic::Diagnostic::code)
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
        .map(pdbiox_core::diagnostic::Diagnostic::code)
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

/// Rebuilds a chain table with the first chain relabelled.
fn rename_first_chain(
    chains: &pdbiox_core::topology::ChainTable,
    label: pdbiox_core::symbol::SymbolId,
) -> pdbiox_core::topology::ChainTable {
    use pdbiox_core::optional::OptionalSymbol;
    use pdbiox_core::topology::{ChainRecord, ChainTable, PolymerKind};

    let mut rebuilt = ChainTable::default();
    for chain in chains.iter() {
        let Some(residues) = chains.residues(chain) else {
            continue;
        };
        let Some(entity) = chains.entity(chain) else {
            continue;
        };
        let relabelled = chain.get() == 0;
        rebuilt.push(
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
        );
    }
    rebuilt
}
