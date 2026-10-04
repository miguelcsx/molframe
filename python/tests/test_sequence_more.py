"""Matrix alignment, multiple alignment and trees, each checked independently."""

import numpy as np
import pytest

import molframe
from molframe import sequence

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
rng = np.random.default_rng(11)


def random_sequence(length):
    return "".join(rng.choice(list(ALPHABET), size=length))


def reference_score(  # noqa: PLR0913, PLR0917
    left, right, matrix, open_cost, extend_cost, local
):
    """Gotoh affine-gap dynamic programme; a gap of length k costs open + (k - 1) * extend."""
    n, m = len(left), len(right)
    floor = -(10**9)
    best = [[0 if local else floor] * (m + 1) for _ in range(n + 1)]
    gap_left = [[floor] * (m + 1) for _ in range(n + 1)]
    gap_right = [[floor] * (m + 1) for _ in range(n + 1)]
    best[0][0] = 0
    result = 0 if local else floor
    for i in range(n + 1):
        for j in range(m + 1):
            if i == 0 and j == 0:
                continue
            if i > 0:
                gap_left[i][j] = max(best[i - 1][j] + open_cost, gap_left[i - 1][j] + extend_cost)
            if j > 0:
                gap_right[i][j] = max(best[i][j - 1] + open_cost, gap_right[i][j - 1] + extend_cost)
            diagonal = (
                best[i - 1][j - 1] + matrix.score(left[i - 1], right[j - 1])
                if i > 0 and j > 0
                else floor
            )
            value = max(diagonal, gap_left[i][j], gap_right[i][j])
            best[i][j] = max(value, 0) if local else value
            result = max(result, best[i][j]) if local else best[n][m]
    return result if local else best[n][m]


def test_blosum62_holds_the_published_scores():
    matrix = sequence.SubstitutionMatrix.blosum(62)
    assert (matrix.score("A", "A"), matrix.score("W", "W"), matrix.score("C", "C")) == (4, 11, 9)
    assert (matrix.score("A", "R"), matrix.score("W", "C"), matrix.score("I", "V")) == (-1, -2, 3)
    assert matrix.name == "BLOSUM62"
    assert matrix.score("a", "A") == 4
    with pytest.raises(ValueError, match="single-letter"):
        matrix.score("AA", "A")


def test_the_bundled_matrices_load_and_unknown_ones_are_refused():
    assert sequence.SubstitutionMatrix.blosum(50).name == "BLOSUM50"
    assert sequence.SubstitutionMatrix.pam(250).name == "PAM250"
    assert sequence.SubstitutionMatrix.nuc44().score("A", "A") > 0
    assert sequence.SubstitutionMatrix.identity().score("A", "A") == 1
    with pytest.raises(molframe.MolframeError) as raised:
        sequence.SubstitutionMatrix.blosum(63)
    assert raised.value.code is not None


@pytest.mark.parametrize("local", [False, True])
def test_matrix_alignment_scores_equal_an_independent_gotoh_programme(local):
    matrix = sequence.SubstitutionMatrix.blosum(62)
    for _ in range(25):
        left, right = random_sequence(rng.integers(6, 24)), random_sequence(rng.integers(6, 24))
        found = sequence.align(
            left,
            right,
            mode="local" if local else "global",
            matrix=matrix,
            gap_open=-10,
            gap_extend=-1,
        )
        assert found.score == reference_score(left, right, matrix, -10, -1, local), (left, right)


def test_a_semi_global_alignment_never_scores_below_the_global_one():
    matrix = sequence.SubstitutionMatrix.blosum(62)
    left, right = "WWWAGHKLM", "AGHK"
    glob = sequence.align(left, right, matrix=matrix, gap_open=-10, gap_extend=-1)
    semi = sequence.align(
        left, right, mode="semi_global", matrix=matrix, gap_open=-10, gap_extend=-1
    )
    assert semi.score >= glob.score
    assert semi.score == sum(matrix.score(c, c) for c in "AGHK")


def test_a_matrix_alignment_states_its_gaps_and_refuses_a_second_scoring():
    matrix = sequence.SubstitutionMatrix.blosum(62)
    with pytest.raises(molframe.MolframeValueError, match="gap_open and gap_extend"):
        sequence.align("AC", "AC", matrix=matrix)
    with pytest.raises(molframe.MolframeValueError, match="not both"):
        sequence.align(
            "AC", "AC", matrix=matrix, scoring=sequence.Scoring(), gap_open=-1, gap_extend=-1
        )
    with pytest.raises(molframe.MolframeValueError, match="zero or negative"):
        sequence.align("AC", "AC", matrix=matrix, gap_open=3, gap_extend=-1)
    with pytest.raises(molframe.MolframeValueError, match="go with a matrix"):
        sequence.align("AC", "AC", gap_open=-1)


def test_a_multiple_alignment_keeps_every_sequence_and_aligns_the_conserved_core():
    base = "MKTAYIAKQRQISFVKSHFSRQ"
    variants = [base, base[:5] + base[7:], base[:10] + "GG" + base[10:], base.replace("K", "R")]
    rows = sequence.msa(variants, scoring=sequence.Scoring(match_score=2, mismatch_score=-1))
    assert len({len(row) for row in rows}) == 1
    assert [row.replace("-", "") for row in rows] == variants
    columns = list(zip(*rows, strict=True))
    conserved = sum(1 for column in columns if len(set(column)) == 1 and column[0] != "-")
    assert conserved >= 12
    refined = sequence.msa(
        variants, scoring=sequence.Scoring(match_score=2, mismatch_score=-1), refinement_passes=2
    )
    assert [row.replace("-", "") for row in refined] == variants


def tip_to_tip(tree, first, second):  # noqa: C901
    """Return the path length between two tips of a Newick tree, from a small recursive parse."""
    text = tree.newick().rstrip(";")
    state = {"position": 0, "next": 0}

    def branch():
        position = state["position"]
        if position < len(text) and text[position] == ":":
            position += 1
            start = position
            while position < len(text) and text[position] not in ",)":
                position += 1
            state["position"] = position
            return float(text[start:position])
        return 0.0

    def parse():
        identity = state["next"]
        state["next"] += 1
        if text[state["position"]] == "(":
            state["position"] += 1
            children = [parse()]
            while text[state["position"]] == ",":
                state["position"] += 1
                children.append(parse())
            state["position"] += 1
            return (identity, None, children, branch())
        start = state["position"]
        while text[state["position"]] not in ":,)":
            state["position"] += 1
        return (identity, text[start : state["position"]], [], branch())

    ancestry = {}

    def walk(node, trail):
        identity, name, children, length = node
        trail = [*trail, (identity, length)]
        if name is not None:
            ancestry[name] = trail
        for child in children:
            walk(child, trail)

    walk(parse(), [])
    a, b = ancestry[first], ancestry[second]
    shared = 0
    while shared < min(len(a), len(b)) and a[shared][0] == b[shared][0]:
        shared += 1
    return sum(length for _, length in a[shared:]) + sum(length for _, length in b[shared:])


def test_neighbor_joining_recovers_an_additive_tree_exactly():
    # An unrooted tree with leaf branches A:2, B:3, C:4, D:5 and an internal branch of 1.5
    # (A, B on one side; C, D on the other).
    labels = ["A", "B", "C", "D"]
    distances = np.array(
        [
            [0.0, 5.0, 7.5, 8.5],
            [5.0, 0.0, 8.5, 9.5],
            [7.5, 8.5, 0.0, 9.0],
            [8.5, 9.5, 9.0, 0.0],
        ]
    )
    tree = sequence.neighbor_joining(labels, distances)
    assert sorted(tree.leaves()) == labels
    assert tree.leaf_count == 4
    for i in range(4):
        for j in range(i + 1, 4):
            assert tip_to_tip(tree, labels[i], labels[j]) == pytest.approx(
                distances[i, j], abs=1e-9
            )


def test_upgma_places_the_root_at_half_the_largest_distance_for_ultrametric_input():
    labels = ["A", "B", "C"]
    distances = np.array([[0.0, 4.0, 10.0], [4.0, 0.0, 10.0], [10.0, 10.0, 0.0]])
    tree = sequence.upgma(labels, distances)
    assert tip_to_tip(tree, "A", "B") == pytest.approx(4.0)
    assert tip_to_tip(tree, "A", "C") == pytest.approx(10.0)
    assert tip_to_tip(tree, "B", "C") == pytest.approx(10.0)


def test_a_tree_round_trips_through_newick_and_can_be_rerooted_and_ladderized():
    tree = sequence.Tree.from_newick("((A:1,B:2):1,(C:3,D:4):2);")
    assert tree.newick() == "((A:1,B:2):1,(C:3,D:4):2);"
    assert tree.leaves() == ["A", "B", "C", "D"]
    assert tree.ladderize("descending").leaf_count == 4
    rerooted = tree.reroot_at_leaf("C", 0.5)
    assert sorted(rerooted.leaves()) == ["A", "B", "C", "D"]
    for a, b in (("A", "B"), ("A", "C"), ("C", "D")):
        assert tip_to_tip(rerooted, a, b) == pytest.approx(tip_to_tip(tree, a, b))
    with pytest.raises(molframe.MolframeError):
        tree.reroot_at_leaf("Z", 0.5)
    with pytest.raises(molframe.MolframeError):
        sequence.Tree.from_newick("((A,B")
    with pytest.raises(molframe.MolframeValueError):
        sequence.neighbor_joining(["A", "B"], np.zeros((3, 3)))
