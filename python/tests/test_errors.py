import pickle
import warnings

import numpy as np
import pytest

import molframe
from molframe import errors

PDB = """\
ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N
ATOM      2  CA  ALA A   1      12.560   6.195  -6.504  1.00  0.00           C
END
"""

KINDS = [
    errors.ParseError,
    errors.SchemaError,
    errors.ConsistencyError,
    errors.ConversionError,
    errors.GeometryError,
    errors.PolicyError,
    errors.QueryError,
    errors.MolframeValueError,
    errors.MolframeKeyError,
    errors.MolframeIndexError,
    errors.MolframeTypeError,
    errors.ResourceError,
    errors.MolframeIOError,
    errors.MemoryBudgetError,
    errors.Cancelled,
    errors.InternalError,
]


@pytest.mark.parametrize("kind", KINDS)
def test_every_error_is_a_molframe_error_and_names_its_builtin(kind):
    assert issubclass(kind, errors.MolframeError)
    assert issubclass(errors.MolframeError, Exception)


def test_builtin_bases_keep_existing_handlers_working():
    assert issubclass(errors.MolframeValueError, ValueError)
    assert issubclass(errors.ParseError, ValueError)
    assert issubclass(errors.GeometryError, ValueError)
    assert issubclass(errors.MolframeKeyError, KeyError)
    assert issubclass(errors.MolframeIndexError, IndexError)
    assert issubclass(errors.MolframeTypeError, TypeError)
    assert issubclass(errors.MemoryBudgetError, MemoryError)
    assert issubclass(errors.MolframeIOError, OSError)
    assert issubclass(errors.ResourceError, RuntimeError)
    assert issubclass(errors.QueryError, ValueError)


def test_the_package_re_exports_every_class():
    for name in ("MolframeError", "QueryError", "QueryWarning", "MemoryBudgetError", "Diagnostic"):
        assert getattr(molframe, name) is getattr(errors, name)
    assert molframe.errors is errors


def test_errors_survive_pickling_with_their_structured_fields():
    original = errors.MolframeKeyError(
        "no such chain",
        code="MOLFRAME-E6006",
        remedy="check the chain count",
        span=(3, 9),
        findings=[errors.Diagnostic("MOLFRAME-E6006", "no such chain")],
    )
    copy = pickle.loads(pickle.dumps(original))  # noqa: S301
    assert type(copy) is type(original)
    assert (copy.code, copy.remedy, copy.span) == (
        "MOLFRAME-E6006",
        "check the chain count",
        (3, 9),
    )
    assert copy.findings[0].message == "no such chain"
    assert str(copy) == "no such chain [MOLFRAME-E6006]"


def test_a_missing_file_is_an_os_error_with_its_code():
    with pytest.raises(OSError, match="E7101") as raised:
        molframe.read("no/such/file.cif")
    assert isinstance(raised.value, errors.MolframeIOError)
    assert raised.value.code == "MOLFRAME-E7101"
    assert raised.value.remedy
    assert raised.value.findings[0].code == "MOLFRAME-E7101"


def test_a_bad_selection_is_a_query_error_that_carries_every_finding():
    structure = molframe.read(PDB.encode(), name="t.pdb")
    with pytest.raises(molframe.QueryError) as raised:
        structure.select("chain A and nonsense_keyword")
    error = raised.value
    assert isinstance(error, ValueError)
    assert error.code is not None
    assert error.findings
    assert error.span is not None


def test_a_bad_chain_label_is_a_key_error_and_a_bad_index_an_index_error():
    structure = molframe.read(PDB.encode(), name="t.pdb")
    with pytest.raises(KeyError):
        structure.chains["Z"]
    with pytest.raises(IndexError):
        structure.chains[99]


def test_numeric_failures_carry_their_diagnostic_code_and_remedy():
    with pytest.raises(errors.GeometryError) as raised:
        molframe.geometry.rmsd(
            molframe_array([[0.0, 0.0, 0.0]]), molframe_array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        )
    assert raised.value.code == "MOLFRAME-E5102"
    assert raised.value.remedy


def molframe_array(rows):
    return np.array(rows, dtype=np.float32)


def test_warnings_are_warnings_that_filters_can_see():
    with pytest.warns(molframe.QueryWarning) as caught:
        warnings.warn(errors.QueryWarning("odd", code="MOLFRAME-W4001"), stacklevel=1)
    assert caught[0].message.code == "MOLFRAME-W4001"


def test_an_indeterminate_error_is_a_molframe_error_that_carries_its_reason():
    error = molframe.IndeterminateError("the policy rejects 3 missing and 0 ambiguous inputs")
    assert isinstance(error, molframe.MolframeError)
    assert "missing" in str(error)
    assert not isinstance(error, ValueError)
