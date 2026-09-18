use super::*;

#[test]
fn the_registry_is_sorted_so_lookup_can_binary_search_it() {
    for pair in ENTRIES.windows(2) {
        let [earlier, later] = pair else { continue };
        assert!(
            earlier.code < later.code,
            "{} is not before {}",
            earlier.code,
            later.code
        );
    }
}

#[test]
fn every_registered_code_carries_a_cause_and_a_remedy() {
    for entry in ENTRIES {
        assert!(!entry.cause.is_empty(), "{} has no cause", entry.code);
        assert!(!entry.remedy.is_empty(), "{} has no remedy", entry.code);
    }
}

#[test]
fn every_registered_code_resolves_to_its_own_entry() {
    for entry in ENTRIES {
        assert!(entry.code.is_registered());
        assert_eq!(entry.code.severity(), entry.severity);
        assert_eq!(entry.code.cause(), entry.cause);
    }
}

#[test]
fn refusing_a_lossy_write_is_an_error_and_never_merely_a_warning() {
    for code in [Code::E4101, Code::E4102, Code::E4103, Code::E4104] {
        assert!(code.severity() >= Severity::Invalidating);
    }
}
