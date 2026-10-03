use crate::{PdbHeaders, PdbHeadersExt};
use molframe_bench::{Sample, input};
use molframe_core::io::ReadOptions;

fn headers_of(sample: Sample) -> PdbHeaders {
    let bytes = sample.pdb().expect("the sample ships a PDB");
    let (structure, _) = crate::read(&input(bytes), &ReadOptions::new()).expect("sample reads");
    structure
        .pdb_headers()
        .cloned()
        .expect("a PDB read keeps its headers")
}

#[test]
fn seqres_gives_each_chain_its_declared_sequence() {
    let chains = headers_of(Sample::Medium).seqres();
    assert_eq!(chains.len(), 4);
    for chain in &chains {
        assert_eq!(
            chain.declared,
            u32::try_from(chain.residues.len()).ok(),
            "chain {}",
            chain.chain
        );
    }
    assert_eq!(&*chains[0].chain, "A");
    assert_eq!(chains[0].residues.len(), 141);
    assert_eq!(&*chains[0].residues[0], "VAL");
}

#[test]
fn ssbonds_name_both_cysteines() {
    let bonds = headers_of(Sample::Tiny).ssbonds();
    let pairs: Vec<_> = bonds
        .iter()
        .map(|bond| (bond.first.sequence, bond.second.sequence))
        .collect();
    assert_eq!(pairs, [(3, 40), (4, 32), (16, 26)]);
    assert_eq!(&*bonds[0].first.name, "CYS");
    assert_eq!(&*bonds[0].second.chain, "A");
}

#[test]
fn links_carry_atoms_residues_and_lengths() {
    let links = headers_of(Sample::Medium).links();
    assert_eq!(links.len(), 4);
    assert_eq!(&*links[0].first_atom, "NE2");
    assert_eq!(&*links[0].first.name, "HIS");
    assert_eq!(links[0].first.sequence, 87);
    assert_eq!(&*links[0].second_atom, "FE");
    assert_eq!(&*links[0].second.name, "HEM");
    assert!(
        links[0]
            .distance
            .is_some_and(|value| (value - 2.0).abs() < 0.5)
    );
}

#[test]
fn the_typed_views_leave_the_preserved_lines_untouched() {
    let headers = headers_of(Sample::Medium);
    let before: Vec<_> = headers
        .records()
        .iter()
        .map(|record| record.line().to_owned())
        .collect();
    let _ = (headers.seqres(), headers.ssbonds(), headers.links());
    let after: Vec<_> = headers
        .records()
        .iter()
        .map(|record| record.line().to_owned())
        .collect();
    assert_eq!(before, after);
}
