//! The `assembly` selector over a materialised biological assembly.

use molframe::crystal::AssemblyExt;
use molframe::{AnalysisPolicy, ReadOptions};

#[test]
fn assembly_selects_each_materialised_chain_copy_of_1aon() {
    let input = molframe_bench::input_gzip(molframe_bench::large_cif_gz());
    let (structure, _) = molframe::read_buffer(&input, Some("1aon.cif"), &ReadOptions::new())
        .expect("the embedded 1AON reads");
    let materialized: molframe::Structure = structure
        .assembly("1")
        .expect("1AON declares assembly 1")
        .materialize()
        .expect("assembly 1 materialises")
        .into();
    let policy = AnalysisPolicy::default();

    let mut covered = 0_u64;
    for copy in 0..materialized.chain_count() {
        let chain = materialized
            .chain_at(copy)
            .expect("every copy index names a chain");
        let selected = materialized
            .select(&format!("assembly {copy}"), &policy)
            .expect("the instance annotation is present");
        let expected: Vec<u32> = chain
            .residues()
            .flat_map(molframe::ResidueRef::atoms)
            .map(|atom| atom.index().get())
            .collect();
        assert_eq!(
            selected
                .atoms()
                .map(|atom| atom.index().get())
                .collect::<Vec<_>>(),
            expected,
            "copy {copy} selects exactly its chain"
        );
        covered += selected.len();
    }
    assert_eq!(covered, u64::from(materialized.atom_count()));

    let refused = structure.select("assembly 0", &policy);
    assert!(
        refused.is_err(),
        "the asymmetric unit has no instance annotation"
    );
}
