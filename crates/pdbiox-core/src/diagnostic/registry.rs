//! The diagnostic code registry.
//!
//! Declaring codes as data rather than as scattered string literals buys three
//! things: a remedy is structurally required rather than remembered, the whole
//! set can be enumerated for documentation, and a code's severity is decided in
//! one place instead of at each site that raises it.
//!
//! Entries are declared in sorted order — errors before warnings, then by
//! number — because lookup binary-searches this table. A test enforces it.

use super::code::{Code, Entry, Kind, Severity};

macro_rules! registry {
    ($($name:ident = $kind:ident $number:literal, $severity:ident, $cause:literal, $remedy:literal;)*) => {
        impl Code {
            $(
                #[doc = $cause]
                pub const $name: Self = Self::new(Kind::$kind, $number);
            )*
        }

        pub(super) static ENTRIES: &[Entry] = &[$(
            Entry {
                code: Code::new(Kind::$kind, $number),
                severity: Severity::$severity,
                cause: $cause,
                remedy: $remedy,
            },
        )*];
    };
}

registry! {
    E1001 = Error 1001, Breaking,
        "format could not be determined",
        "pass the format explicitly, or give the file an extension the reader recognises";
    E1101 = Error 1101, Breaking,
        "unterminated quoted value",
        "close the quote; a value containing a quote character must use the other quote style or a multi-line text field";
    E1102 = Error 1102, Breaking,
        "unterminated multi-line text field",
        "close the field with a semicolon at the start of a line";
    E1103 = Error 1103, Breaking,
        "loop_ row count is not a multiple of the column count",
        "count the values in the final row; an unquoted value containing a space parses as two values";
    E1104 = Error 1104, Invalidating,
        "duplicate tag within a category",
        "remove the repeated tag; the two values cannot both be the value of one item";
    E1105 = Error 1105, Breaking,
        "value outside any category",
        "precede the value with the tag it belongs to, or with a loop_ header";
    E1106 = Error 1106, Breaking,
        "missing data_ block header",
        "add a data_ line before the first category; a file with no block cannot be addressed";
    E1201 = Error 1201, Invalidating,
        "malformed fixed-column record",
        "check the record against its column layout; fields are positional, so one misplaced column shifts every field after it";
    E1202 = Error 1202, Invalidating,
        "unparseable numeric field",
        "correct the field, or read in recover mode to keep the rest of the file";
    E1203 = Error 1203, Invalidating,
        "atom serial exceeds format capacity and is not hybrid-36",
        "re-encode the serial in hybrid-36, or use a format without a five-column serial field";
    E1401 = Error 1401, Breaking,
        "declared length disagrees with the encoded payload",
        "the file is truncated or corrupt; fetch it again";
    E1402 = Error 1402, Breaking,
        "unknown column encoding",
        "the file uses an encoding this version does not implement; upgrade pdbiox";
    E1403 = Error 1403, Breaking,
        "codec chain type mismatch",
        "the encoding chain is inconsistent with itself; report it to whatever wrote the file";
    E1901 = Error 1901, Breaking,
        "resource limit exceeded",
        "raise the limit deliberately if the input is genuinely this large; the default guards against decompression bombs and absurd declared row counts";
    E2001 = Error 2001, Invalidating,
        "required category absent",
        "add the category, or read with a strictness that tolerates its absence";
    E2002 = Error 2002, Invalidating,
        "required item absent",
        "add the item; without it the affected rows cannot be interpreted";
    E2003 = Error 2003, Strict,
        "value outside the declared enumeration",
        "use one of the enumerated values, or relax strictness to accept the file as deposited";
    E2004 = Error 2004, Invalidating,
        "value has the wrong type for its item",
        "correct the value; a numeric item cannot hold text";
    E2005 = Error 2005, Strict,
        "parent-child relationship violated",
        "add the missing parent row, or correct the child's reference";
    E3001 = Error 3001, Invalidating,
        "column length disagrees with the chunk length",
        "report this as a pdbiox bug; a structure cannot be built this way from valid input";
    E3002 = Error 3002, Invalidating,
        "child range exceeds the child table",
        "report this as a pdbiox bug";
    E3003 = Error 3003, Invalidating,
        "overlapping child ranges",
        "report this as a pdbiox bug";
    E3004 = Error 3004, Invalidating,
        "atom belongs to no residue",
        "report this as a pdbiox bug";
    E3005 = Error 3005, Invalidating,
        "chain references a nonexistent entity",
        "check the entity identifiers the chains refer to";
    E3006 = Error 3006, Invalidating,
        "bond endpoint out of range",
        "check the connectivity records; they name an atom the file does not contain";
    E3007 = Error 3007, Invalidating,
        "self-bond",
        "remove the connectivity record bonding an atom to itself";
    E3008 = Error 3008, Invalidating,
        "non-finite coordinate without a validity mark",
        "correct the coordinate, or mark the atom's position as unrecorded";
    E3009 = Error 3009, Strict,
        "occupancy outside the zero-to-one range",
        "correct the occupancy; values outside the range make occupancy-weighted results meaningless";
    E3010 = Error 3010, Invalidating,
        "models share a coordinate set but not an atom count",
        "the models differ in composition and must be held as a ragged ensemble rather than as frames";
    E3011 = Error 3011, Invalidating,
        "custom annotation length disagrees with the atom count",
        "supply exactly one annotation row for every atom row";
    E3012 = Error 3012, Invalidating,
        "structures cannot share one merged model axis",
        "merge structures with equal dense model counts, or keep unlike models as a ragged ensemble";
    E3013 = Error 3013, Invalidating,
        "custom annotation types disagree across structures",
        "convert same-named annotation columns to one physical type before merging";
    E3014 = Error 3014, Invalidating,
        "domain extensions would be invalidated by a topology change",
        "explicitly clear extensions before editing topology or merging structures";
    E3015 = Error 3015, Invalidating,
        "residue boundary is not defined by the deposited identifiers",
        "repair the residue identifiers or explicitly select file-order boundary inference";
    E3201 = Error 3201, Invalidating,
        "bundled chemistry reference data is invalid",
        "reinstall pdbiox from a verified distribution";
    E3401 = Error 3401, Invalidating,
        "topology and coordinate atom counts disagree",
        "check that the coordinate source and the topology source describe the same system";
    E4001 = Error 4001, Invalidating,
        "selection syntax is invalid",
        "correct the selection at the reported token";
    E4002 = Error 4002, Invalidating,
        "selection value is incompatible with its column",
        "use a value of the column's declared type";
    E4003 = Error 4003, Invalidating,
        "selection requires data this structure does not carry",
        "provide the required annotations or choose a selector available for this structure";
    E4004 = Error 4004, Invalidating,
        "unknown selection keyword",
        "check the spelling against the keyword list; an annotation column can also be selected by its own name";
    E4101 = Error 4101, Invalidating,
        "atom count exceeds legacy format capacity",
        "write mmCIF or BinaryCIF, or enable hybrid-36 serials explicitly";
    E4102 = Error 4102, Invalidating,
        "chain identifier cannot be represented in the target format",
        "write mmCIF or BinaryCIF, or supply an explicit chain mapping";
    E4103 = Error 4103, Invalidating,
        "residue number exceeds legacy format capacity",
        "write mmCIF or BinaryCIF, or enable hybrid-36 residue numbers explicitly";
    E4104 = Error 4104, Invalidating,
        "coordinate magnitude exceeds the target field width",
        "translate the structure towards the origin, or write a format with wider coordinate fields";
    E4105 = Error 4105, Invalidating,
        "the target format cannot express a required category",
        "write a format that models the data, or drop the category deliberately";
    E5001 = Error 5001, Invalidating,
        "superposition is undefined for fewer than three non-collinear atoms",
        "widen the selection until it spans three dimensions";
    E5002 = Error 5002, Invalidating,
        "degenerate frame from collinear reference atoms",
        "choose reference atoms that are not in a straight line";
    E5003 = Error 5003, Invalidating,
        "empty selection where atoms are required",
        "check the selection against the structure; a keyword that matches nothing here may match elsewhere";
    E5004 = Error 5004, Invalidating,
        "periodic query without a unit cell",
        "supply a cell, or ask for a non-periodic query";
    E5005 = Error 5005, Invalidating,
        "cutoff exceeds half the minimum box dimension under the minimum image convention",
        "reduce the cutoff, or use a periodic mode that considers more than the nearest image";
    E6001 = Error 6001, Invalidating,
        "unqualified identifier while the policy requires explicit namespaces",
        "qualify the selector, for instance auth_chain rather than chain";
    E6002 = Error 6002, Invalidating,
        "requested assembly does not exist",
        "list the assemblies the entry defines and choose one of them";
    E6003 = Error 6003, Invalidating,
        "requested model does not exist",
        "check the model count; models are addressed by position, not by deposited number";
    E6004 = Error 6004, Invalidating,
        "policy field combination is contradictory",
        "reconcile the two fields; the analysis cannot honour both";
    E6005 = Error 6005, Invalidating,
        "verdict profile not found",
        "list the registered profiles; profile names carry a version";
    E6006 = Error 6006, Invalidating,
        "requested chain does not exist",
        "check the chain count or resolve the chain label before editing";
    E6007 = Error 6007, Invalidating,
        "identifier is empty",
        "supply a non-empty identifier";
    E6008 = Error 6008, Invalidating,
        "a whole-ensemble coordinate edit requires shared topology",
        "edit each ragged model snapshot independently";
    E6009 = Error 6009, Invalidating,
        "an atom selected for editing does not exist",
        "restrict the selection to the structure atom range";
    E6010 = Error 6010, Invalidating,
        "assembly operator expression is malformed",
        "use comma-separated identifiers and integer ranges, with each product factor in parentheses";
    E6011 = Error 6011, Invalidating,
        "assembly operator expression exceeds the instance limit",
        "reduce the operator product or raise the explicit assembly limit";
    E6012 = Error 6012, Invalidating,
        "assembly category contains a malformed required value",
        "correct the named assembly category item and row";
    E6013 = Error 6013, Invalidating,
        "assembly generator references an unknown definition",
        "correct the referenced assembly, operator, or label-asym identifier";
    E6014 = Error 6014, Invalidating,
        "non-crystallographic symmetry operator is malformed",
        "correct the NCS identifier, given/generate code, matrix, or translation vector";
    E6015 = Error 6015, Invalidating,
        "crystallographic symmetry metadata or algebraic operation is malformed",
        "correct the space-group number, operation identifier, or fractional xyz expression";
    E6016 = Error 6016, Invalidating,
        "crystal neighbour search configuration is incomplete or invalid",
        "supply explicit symmetry operators, a valid model, and a positive finite cutoff";
    E6017 = Error 6017, Invalidating,
        "crystal neighbour search exceeds its candidate-image limit",
        "reduce the cutoff or raise the explicit search limit after reviewing the expected cost";
    E6018 = Error 6018, Invalidating,
        "space-group setting is unknown or the bundled catalogue is invalid",
        "supply a valid Hall symbol or International Tables number, or reinstall pdbiox";
    E7901 = Error 7901, Invalidating,
        "output could not be written",
        "check the destination, permissions and available storage, then retry";
    E9001 = Error 9001, Breaking,
        "internal invariant violated",
        "report this as a pdbiox bug with the input that produced it; no input should be able to cause this";
    W1001 = Warning 1001, Loose,
        "trailing content after the final block",
        "remove the trailing content, or ignore this if the file is deliberately concatenated";
    W1002 = Warning 1002, Info,
        "non-standard whitespace or line ending",
        "normalise the line endings if the file is meant to be portable";
    W2001 = Warning 2001, Info,
        "unknown category retained but not interpreted",
        "no action needed; the category survives a preserving write untouched";
    W2002 = Warning 2002, Info,
        "unknown item retained but not interpreted",
        "no action needed; the item survives a preserving write untouched";
    W2003 = Warning 2003, Loose,
        "deprecated item used",
        "migrate to the replacement item named in the dictionary";
    W3011 = Warning 3011, Info,
        "residue boundary was inferred from file order by explicit policy",
        "check the residue identifiers in this region; consecutive residues carry annotations that do not distinguish them";
    W3012 = Warning 3012, Loose,
        "an alternate location carries a different component identity",
        "no action needed; this is one residue modelled with two chemical identities, and both are retained";
    W3201 = Warning 3201, Loose,
        "component not found in any chemical component provider",
        "supply a provider that knows this component if its expected atoms or bonds matter";
    W3202 = Warning 3202, Strict,
        "observed atoms are inconsistent with the component definition; file connectivity was used",
        "check the component's atom names against its definition";
    W3203 = Warning 3203, Loose,
        "element inferred from the atom name",
        "declare the element explicitly; the name-based rule is a convention that files are free to violate";
    W3204 = Warning 3204, Loose,
        "unusual valence",
        "check the connectivity and formal charge of the atom named";
    W3205 = Warning 3205, Strict,
        "chirality disagrees with the component definition",
        "check the geometry of the centre named; this is often a genuine modelling error";
    W3301 = Warning 3301, Info,
        "atoms expected by the component definition are absent from the model",
        "no action needed; the shortfall is reported as coverage rather than hidden";
    W3302 = Warning 3302, Loose,
        "chain break detected within a polymer",
        "no action needed; the gap is unmodelled density and is reported rather than bridged";
    W4001 = Warning 4001, Loose,
        "unparenthesised mix of and and or",
        "parenthesise the expression; and binds more tightly than or, which may not be what was meant";
    W4002 = Warning 4002, Loose,
        "float equality comparison",
        "compare against a range instead; exact equality on a measured quantity rarely matches";
    W4003 = Warning 4003, Loose,
        "selection symbol not present in this structure",
        "check the spelling; the selection is valid but matches nothing here";
    W4101 = Warning 4101, Strict,
        "metadata dropped by the target format",
        "write a format that models the metadata if it matters";
    W5001 = Warning 5001, Info,
        "neighbour list rebuilt more often than expected",
        "increase the skin distance to trade memory for rebuild frequency";
    W5002 = Warning 5002, Info,
        "distance matrix requested for a large selection",
        "use a neighbour query instead; a full matrix grows with the square of the selection";
    W6001 = Warning 6001, Strict,
        "per-atom highest occupancy may produce a conformation that was never modelled",
        "prefer a conformer-consistent altloc policy unless per-atom selection is genuinely intended";
    W6002 = Warning 6002, Info,
        "assembly generated multiple instances of a chain",
        "no action needed; instances carry their own identifiers and are distinguishable";
    W6003 = Warning 6003, Strict,
        "analysis returned an indeterminate result",
        "read the assumptions and coverage; no defensible single answer exists under this policy";
    W6004 = Warning 6004, Info,
        "a defaulted policy field materially affects this result",
        "set the field explicitly so the choice is recorded as yours rather than inherited";
}

#[cfg(test)]
#[path = "registry_tests.rs"]
mod tests;
