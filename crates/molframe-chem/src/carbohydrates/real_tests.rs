use super::*;
use crate::{CarbohydrateLink, CarbohydrateReport, read_ccd};
use molframe_cif::Category;
use std::collections::{BTreeMap, BTreeSet};

type LinkKey = (String, i32, String, i32, String);

fn text<'a>(category: &'a Category, item: &str, row: usize) -> &'a str {
    category
        .text(item, row)
        .unwrap_or_else(|| panic!("missing {item} row {row}"))
}

fn identifier(category: &Category, item: &str, row: usize) -> String {
    category.identifier(item, row).unwrap().into_owned()
}

fn chain_label(structure: &Structure, residue: molframe_core::ResidueIndex) -> &str {
    structure
        .data()
        .chains()
        .find(|chain| chain.residues().any(|item| item.index() == residue))
        .unwrap()
        .label()
        .unwrap()
}

fn observed_links(structure: &Structure, links: &[CarbohydrateLink]) -> BTreeSet<LinkKey> {
    links
        .iter()
        .map(|link| {
            let donor = structure.atom(link.donor_atom).unwrap();
            let acceptor = structure.atom(link.acceptor_atom).unwrap();
            (
                chain_label(structure, donor.residue().unwrap().index()).to_owned(),
                donor.residue().unwrap().auth_seq_id().unwrap(),
                donor.name().unwrap().to_owned(),
                acceptor.residue().unwrap().auth_seq_id().unwrap(),
                acceptor.name().unwrap().to_owned(),
            )
        })
        .collect()
}

fn raw_links(scheme: &Category, links: &Category) -> BTreeSet<LinkKey> {
    let mut expected = BTreeSet::new();
    for row in 0..links.row_count() {
        let entity = identifier(links, "entity_id", row);
        let scheme_row = |item| {
            (0..scheme.row_count())
                .find(|r| {
                    identifier(scheme, "entity_id", *r) == entity
                        && identifier(scheme, "num", *r) == identifier(links, item, row)
                })
                .expect("branch residue in raw scheme")
        };
        let donor = scheme_row("entity_branch_list_num_1");
        let acceptor = scheme_row("entity_branch_list_num_2");
        let chain = text(scheme, "asym_id", donor);
        assert_eq!(chain, text(scheme, "asym_id", acceptor));
        expected.insert((
            chain.to_owned(),
            identifier(scheme, "pdb_seq_num", donor).parse().unwrap(),
            text(links, "atom_id_1", row).to_owned(),
            identifier(scheme, "pdb_seq_num", acceptor).parse().unwrap(),
            text(links, "atom_id_2", row).to_owned(),
        ));
    }
    expected
}

fn assert_composition(structure: &Structure, report: &CarbohydrateReport, scheme: &Category) {
    let expected: BTreeSet<_> = (0..scheme.row_count())
        .map(|row| {
            (
                text(scheme, "asym_id", row).to_owned(),
                identifier(scheme, "pdb_seq_num", row)
                    .parse::<i32>()
                    .unwrap(),
                text(scheme, "mon_id", row).to_owned(),
            )
        })
        .collect();
    let observed: BTreeSet<_> = report
        .monosaccharides
        .iter()
        .map(|sugar| {
            let residue = structure
                .data()
                .residues()
                .find(|r| r.index() == sugar.residue)
                .unwrap();
            assert!(sugar.geometry.is_some());
            assert_eq!(sugar.ring_atoms.len(), 6);
            (
                chain_label(structure, sugar.residue).to_owned(),
                residue.auth_seq_id().unwrap(),
                residue.name().unwrap().to_owned(),
            )
        })
        .collect();
    assert_eq!(observed, expected);
    assert_eq!(observed.len(), report.monosaccharides.len());
    let mut composition = BTreeMap::new();
    for sugar in &report.monosaccharides {
        *composition.entry(sugar.symbol.abbreviation).or_insert(0) += 1;
    }
    assert_eq!(
        composition,
        BTreeMap::from([("GlcNAc", 8), ("Man", 6), ("Gal", 3), ("Fuc", 1)])
    );
}

fn assert_no_deposited_or_close_protein_attachment(conn: &Category, atoms: &Category) {
    let sugars = ["NAG", "MAN", "BMA", "GAL", "FUC"];
    let mut sugar_links = 0;
    for row in 0..conn.row_count() {
        let a = sugars.contains(&text(conn, "ptnr1_label_comp_id", row));
        let b = sugars.contains(&text(conn, "ptnr2_label_comp_id", row));
        assert_eq!(a, b, "no sugar/protein struct_conn attachment");
        sugar_links += usize::from(a && b);
    }
    assert_eq!(sugar_links, 16);
    let position = |row| {
        ["Cartn_x", "Cartn_y", "Cartn_z"]
            .map(|item| atoms.value(item, row).unwrap().as_float().unwrap())
    };
    for (chain, expected) in [("E", 2.643_074_913_8), ("F", 2.453_818_656_7)] {
        let root = (0..atoms.row_count())
            .find(|row| {
                text(atoms, "label_asym_id", *row) == chain
                    && identifier(atoms, "auth_seq_id", *row) == "1"
                    && text(atoms, "label_atom_id", *row) == "C1"
            })
            .expect("root NAG C1");
        let root = position(root);
        let closest = (0..atoms.row_count())
            .filter(|row| {
                text(atoms, "label_comp_id", *row) == "ASN"
                    && text(atoms, "label_atom_id", *row) == "ND2"
            })
            .map(|row| {
                root.iter()
                    .zip(position(row))
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    .sqrt()
            })
            .min_by(f64::total_cmp)
            .expect("ASN acceptor");
        assert!(closest > 2.0);
        assert!((closest - expected).abs() < 1e-8, "{chain}: {closest}");
    }
}

fn without_sugar_links(structure: &Structure, report: &CarbohydrateReport) -> Structure {
    let mut data = structure.data().clone();
    let sugars: BTreeSet<_> = report
        .monosaccharides
        .iter()
        .map(|sugar| sugar.residue)
        .collect();
    let mut bonds = BondTableBuilder::new();
    for bond in data.bonds.iter() {
        let a = structure
            .atom(bond.atom_a)
            .unwrap()
            .residue()
            .unwrap()
            .index();
        let b = structure
            .atom(bond.atom_b)
            .unwrap()
            .residue()
            .unwrap()
            .index();
        if a == b || !sugars.contains(&a) || !sugars.contains(&b) {
            bonds.push(bond);
        }
    }
    data.bonds = bonds.finish();
    Structure::new(data)
}

#[test]
fn deposited_1hzh_composition_and_every_branch_link_match_the_raw_cif() {
    let input = molframe_bench::input_gzip(include_bytes!("../../data/1HZH.cif.gz"));
    let (document, structure, _) =
        molframe_cif::read_with_document(&input, &ReadOptions::new()).expect("actual 1HZH CIF");
    let (provider, _) = read_ccd(
        &InputBuffer::from_bytes(include_bytes!("../../data/CCD-saccharides.cif").to_vec()),
        DictionaryVersion::new("wwPDB-2026-10-02"),
    )
    .expect("actual CCD");
    let report = carbohydrates(&structure, Some(&provider), CarbohydrateOptions::default())
        .expect("1HZH inventory");
    assert_eq!(report.monosaccharides.len(), 18);
    assert_eq!(report.links.len(), 16);
    assert!(report.terminal_links.is_empty());
    assert!(report.incomplete_residues.is_empty());
    assert!(
        report
            .links
            .iter()
            .all(|link| link.provenance == BondProvenance::File)
    );
    let block = document.blocks().next().expect("block");
    let scheme = block
        .category("pdbx_branch_scheme")
        .expect("raw branch scheme");
    assert_composition(&structure, &report, scheme);
    let expected = raw_links(
        scheme,
        block
            .category("pdbx_entity_branch_link")
            .expect("raw branch links"),
    );
    assert_eq!(observed_links(&structure, &report.links), expected);
    assert_no_deposited_or_close_protein_attachment(
        block.category("struct_conn").unwrap(),
        block.category("atom_site").unwrap(),
    );
    let stripped = without_sugar_links(&structure, &report);
    let disabled = carbohydrates(
        &stripped,
        Some(&provider),
        CarbohydrateOptions {
            spatial_fallback: false,
        },
    )
    .unwrap();
    assert_eq!(disabled.monosaccharides.len(), 18);
    assert!(disabled.links.is_empty());
    let fallback = carbohydrates(&stripped, Some(&provider), CarbohydrateOptions::default())
        .expect("missing linkage fallback");
    assert_eq!(observed_links(&stripped, &fallback.links), expected);
    assert_eq!(fallback.links.len(), 16);
    assert!(fallback.terminal_links.is_empty());
    assert!(
        fallback
            .links
            .iter()
            .all(|link| link.provenance == BondProvenance::InferredDistance)
    );
    eprintln!(
        "1HZH: 18 rings (8 GlcNAc, 6 Man, 3 Gal, 1 Fuc), all 16 raw branch links recovered; no attachment declared or within 2 Å"
    );
}
