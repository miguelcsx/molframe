"""Hydrogen display classes use topology through every Python query surface."""

import molframe


def hydrogen_structure(edges):
    elements = ["N", "O", "S", "C", "P", "F"] + ["H"] * 10
    names = [
        "N",
        "O",
        "S",
        "C",
        "P",
        "F",
        "HN",
        "HO",
        "HS",
        "HC",
        "FREE",
        "MIX",
        "HP",
        "HF",
        "HH1",
        "HH2",
    ]
    atoms = [
        f"HETATM{index + 1:5d} {name:>4s} UNK A   1    "
        f"{index * 10:8.3f}{0:8.3f}{0:8.3f}{1:6.2f}{0:6.2f}          {element:>2s}"
        for index, (name, element) in enumerate(zip(names, elements, strict=True))
    ]
    bonds = [f"CONECT{first + 1:5d}{second + 1:5d}" for first, second in edges]
    return molframe.read("\n".join([*atoms, *bonds, "END", ""]).encode(), name="hydrogen.pdb")


def test_hydrogen_text_compiled_and_typed_queries_agree_and_roundtrip():
    structure = hydrogen_structure(
        [(0, 6), (1, 7), (2, 8), (3, 9), (3, 11), (0, 11), (4, 12), (5, 13), (14, 15)]
    )
    for source, selector, expected in [
        ("polar_hydrogen", molframe.sel.polar_hydrogen, [6, 7, 8, 11]),
        ("nonpolar_hydrogen", molframe.sel.nonpolar_hydrogen, [9, 10, 12, 13, 14, 15]),
    ]:
        query = selector()
        reparsed = molframe.Query(query.source)
        assert query.source == source
        assert query.fingerprint == reparsed.fingerprint
        assert structure.select(source).indices.tolist() == expected
        assert molframe.Query(source).select(structure).indices.tolist() == expected
        assert query.select(structure).indices.tolist() == expected
        assert structure.select(reparsed).indices.tolist() == expected
        assert structure.select(molframe.sel.hydrogen() & query).indices.tolist() == expected
    polar = structure.select(molframe.sel.polar_hydrogen())
    nonpolar = structure.select(molframe.sel.nonpolar_hydrogen())
    assert len(polar & nonpolar) == 0
    assert (polar | nonpolar).indices.tolist() == structure.select("hydrogen").indices.tolist()
    assert structure.select(~molframe.sel.nonpolar_hydrogen()).indices.tolist() == [
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        11,
    ]


def test_hydrogen_parent_survives_hydrogen_only_conjunctions():
    structure = hydrogen_structure([(0, 6), (3, 9)])
    for query in [
        "index 6 9 and polar_hydrogen",
        "hydrogen and polar_hydrogen",
        "polar_hydrogen and hydrogen",
    ]:
        assert structure.select(query).indices.tolist() == [6]
    assert structure.select("index 6 9 and nonpolar_hydrogen").indices.tolist() == [9]


def test_hydrogen_classification_changes_with_bonds_not_names_or_coordinates():
    nitrogen_parent = hydrogen_structure([(0, 9)])
    carbon_parent = hydrogen_structure([(3, 9)])
    query = molframe.sel.polar_hydrogen()
    assert nitrogen_parent.select(query).indices.tolist() == [9]
    assert carbon_parent.select(query).indices.tolist() == []


def test_known_empty_topology_selects_all_hydrogen_as_nonpolar():
    structure = hydrogen_structure([])
    assert structure.select("polar_hydrogen").indices.tolist() == []
    assert structure.select("nonpolar_hydrogen").indices.tolist() == list(range(6, 16))


def test_hydrogen_macros_are_native_completion_candidates():
    for prefix, label in [("polar_", "polar_hydrogen"), ("nonpolar_", "nonpolar_hydrogen")]:
        start, end, candidates = molframe.query.complete(prefix, len(prefix))
        assert (start, end) == (0, len(prefix))
        assert (label, "macro") in candidates
