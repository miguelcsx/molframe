use super::structure_side_chain_torsions;
use crate::{
    AnalysisPolicy, BondOrder, Component, ComponentAtom, ComponentBond, ComponentKind,
    DictionaryVersion, Element, MemoryProvider, ReadOptions,
};
use std::sync::Arc;

const ENTRY: &str = "data_s\n\
loop_\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n\
_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n\
_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
1 N N LYS A 1 0 0 0\n2 C CA LYS A 1 1 0 0\n3 C CB LYS A 1 1 1 0\n\
4 C CG LYS A 1 1 1 1\n5 C CD LYS A 1 2 1 1\n6 C CE LYS A 1 2 2 1\n\
7 N NZ LYS A 1 2 2 2\n";

#[test]
fn structure_projection_delegates_all_chi_geometry_to_the_path_kernel() {
    let structure = match crate::read_bytes(
        ENTRY.as_bytes().to_vec(),
        Some("x.cif"),
        &ReadOptions::new(),
    ) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    let report =
        match structure_side_chain_torsions(&structure, &provider(), &AnalysisPolicy::default()) {
            Ok(report) => report,
            Err(finding) => panic!("torsions failed: {finding}"),
        };
    assert_eq!(report.records.len(), 1);
    assert_eq!(report.records[0].torsions.len(), 4);
    assert!(report.records[0].torsions.iter().all(Option::is_some));
}

fn provider() -> MemoryProvider {
    let names = ["N", "CA", "CB", "CG", "CD", "CE", "NZ"];
    let bonds = names.windows(2).map(|pair| ComponentBond {
        atom_a: pair[0].into(),
        atom_b: pair[1].into(),
        order: BondOrder::Single,
        aromatic: false,
        stereo: None,
    });
    MemoryProvider::new(
        DictionaryVersion::new("test"),
        [Component {
            id: "LYS".into(),
            name: "lysine".into(),
            kind: ComponentKind::AminoAcid,
            parent: None,
            formula: None,
            atoms: names
                .iter()
                .map(|name| ComponentAtom {
                    name: (*name).into(),
                    alternate_name: None,
                    element: if name.starts_with('N') {
                        Element::NITROGEN
                    } else {
                        Element::CARBON
                    },
                    charge: 0,
                    aromatic: false,
                    leaving: false,
                    stereo: None,
                })
                .collect(),
            bonds: bonds.collect(),
            ideal_coordinates: Some(Arc::from([])),
            model_coordinates: None,
        }],
    )
}
