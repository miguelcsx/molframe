//! A store with no profile and no cache: what the compiled-in layer answers.

use crate::{
    BondReach, CovalentRadii, HeteroPolicy, IonicRadii, MetalFallback, MetalSites, RadiusPolicy,
    RadiusSelection, Source, TableError, TableKind, TableRequest, TableStore, VdwRadii,
};
use molframe_core::contract::RadiiSet;

/// A store over an empty directory of its own, with no profile file.
fn clean(name: &str) -> TableStore {
    let dir = std::env::temp_dir().join(format!("molframe-bundled-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let _ = std::fs::create_dir_all(&dir);
    TableStore::new(dir)
}

fn chosen<T: TableKind + Send + Sync + 'static>(store: &TableStore) -> (String, Source) {
    match store.select::<T>(None, "--table") {
        Ok(selected) => (selected.provided.identity.name.to_string(), selected.source),
        Err(error) => panic!("{} resolves on a clean store: {error}", T::KIND),
    }
}

#[test]
fn a_clean_store_selects_every_default_kind_from_memory() {
    let store = clean("kinds");
    assert_eq!(
        chosen::<BondReach>(&store),
        ("macromolecular-reach".to_owned(), Source::Memory)
    );
    assert_eq!(
        chosen::<CovalentRadii>(&store),
        ("pyykko-covalent".to_owned(), Source::Memory)
    );
    assert_eq!(
        chosen::<VdwRadii>(&store),
        ("bondi".to_owned(), Source::Memory)
    );
    assert_eq!(
        chosen::<MetalSites>(&store),
        ("macromolecular-metals".to_owned(), Source::Memory)
    );
}

#[test]
fn the_default_policy_radius_set_resolves_on_a_clean_store() {
    let store = clean("policy");
    let selection = RadiusSelection::choose(None, true);
    match selection.resolve(RadiiSet::default(), &store) {
        Ok(source) => assert_eq!(&*source.identity().name, "bondi"),
        Err(error) => panic!("the default set resolves: {error}"),
    }
}

#[test]
fn the_default_metal_fallback_resolves_on_a_clean_store() {
    let store = clean("metals");
    let Ok(policy) = RadiusPolicy::resolve(
        &TableRequest::profile(),
        MetalFallback::default(),
        HeteroPolicy::default(),
        &store,
    ) else {
        panic!("the default radius policy resolves");
    };
    assert!(policy.radius(molframe_core::Element::IRON).is_some());
}

#[test]
fn a_table_the_flows_do_not_default_to_is_not_compiled_in() {
    let store = clean("ionic");
    assert!(matches!(
        store.select::<IonicRadii>(None, "--ionic-radii"),
        Err(TableError::Missing { .. })
    ));
    assert!(matches!(
        RadiusPolicy::resolve(
            &TableRequest::profile(),
            MetalFallback::Ionic,
            HeteroPolicy::default(),
            &store,
        ),
        Err(TableError::Missing { .. })
    ));
}

#[test]
fn a_cache_profile_and_file_win_over_the_compiled_in_choice() {
    let store = clean("cache-wins");
    let bytes = crate::test_data::shipped_bytes("alvarez.radii");
    let _ = std::fs::write(store.directory().join("profile.txt"), "vdw-radii: lab\n");
    let _ = std::fs::write(store.directory().join("lab.radii"), bytes);
    assert_eq!(
        chosen::<VdwRadii>(&store),
        ("alvarez".to_owned(), Source::Cache)
    );
    assert_eq!(
        chosen::<CovalentRadii>(&store),
        ("pyykko-covalent".to_owned(), Source::Memory)
    );
}

#[test]
fn a_memory_registration_wins_over_the_compiled_in_choice() {
    let store = clean("memory-wins");
    store.insert_profile("vdw-radii: held\n");
    store.insert_table::<VdwRadii>("held", crate::test_data::shipped_bytes("alvarez.radii"));
    assert_eq!(
        chosen::<VdwRadii>(&store),
        ("alvarez".to_owned(), Source::Memory)
    );
}
