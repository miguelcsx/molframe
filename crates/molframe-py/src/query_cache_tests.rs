use super::*;

#[test]
fn a_repeated_query_is_compiled_once_and_shared() {
    let mut cache = QueryCache::default();
    let (Ok(first), Ok(second)) = (
        cache.get_or_compile("protein"),
        cache.get_or_compile("protein"),
    ) else {
        panic!("the query compiles");
    };
    assert!(Arc::ptr_eq(&first, &second));
    assert_eq!(cache.len(), 1);
}

#[test]
fn the_cache_is_bounded_and_evicts_the_least_recently_used_query() {
    let mut cache = QueryCache::default();
    let Ok(kept) = cache.get_or_compile("resid 0") else {
        panic!("compiles");
    };
    for index in 1..QueryCache::CAPACITY {
        assert!(cache.get_or_compile(&format!("resid {index}")).is_ok());
    }
    // Touch the oldest entry, then overflow: the second-oldest is evicted.
    let Ok(touched) = cache.get_or_compile("resid 0") else {
        panic!("compiles");
    };
    assert!(Arc::ptr_eq(&kept, &touched));
    assert!(cache.get_or_compile("resid 999").is_ok());
    assert_eq!(cache.len(), QueryCache::CAPACITY);
    let Ok(again) = cache.get_or_compile("resid 0") else {
        panic!("compiles");
    };
    assert!(Arc::ptr_eq(&kept, &again));
}

#[test]
fn a_malformed_query_is_not_cached() {
    let mut cache = QueryCache::default();
    assert!(cache.get_or_compile("resname (").is_err());
    assert_eq!(cache.len(), 0);
}

#[test]
fn pending_rows_answer_only_the_query_that_produced_them() {
    let evaluation = || molframe::query::Evaluation {
        selection: std::iter::empty::<u32>().collect(),
        warnings: Vec::new(),
    };
    let mut pending = PendingRows::default();
    pending.keep("water", evaluation());
    assert!(pending.take("protein").is_none());
    pending.keep("water", evaluation());
    assert!(pending.take("water").is_some());
    assert!(pending.take("water").is_none());
}
