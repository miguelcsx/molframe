#!/usr/bin/env python3
"""Mechanical gates for the logic-free Python binding layer."""

from __future__ import annotations

import ast
import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "crates" / "pdbiox-py" / "Cargo.toml"
PACKAGE = ROOT / "python" / "pdbiox"
NATIVE_DIR = ROOT / "crates" / "pdbiox-py" / "src"
ALLOWED_DEPENDENCIES = {"arrow", "numpy", "pdbiox", "pyo3"}
REQUIRED_FILES = {
    "__init__.py",
    "__init__.pyi",
    "py.typed",
    "selections.pyi",
    "typing.pyi",
}


def check_dependencies(errors: list[str]) -> None:
    with MANIFEST.open("rb") as source:
        manifest = tomllib.load(source)
    actual = set(manifest.get("dependencies", {}))
    if actual != ALLOWED_DEPENDENCIES:
        errors.append(
            "pdbiox-py dependencies differ: "
            f"expected {sorted(ALLOWED_DEPENDENCIES)}, found {sorted(actual)}"
        )


def check_python_source(errors: list[str]) -> None:
    missing = REQUIRED_FILES.difference(path.name for path in PACKAGE.iterdir())
    if missing:
        errors.append(f"Python package files missing: {sorted(missing)}")
    for path in PACKAGE.glob("*.py"):
        tree = ast.parse(path.read_text(), filename=str(path))
        for index, node in enumerate(tree.body):
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                continue
            if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant):
                if index == 0 and isinstance(node.value.value, str):
                    continue
            if isinstance(node, ast.Assign) and all(
                isinstance(target, ast.Name) and target.id == "__all__"
                for target in node.targets
            ):
                continue
            errors.append(f"executable Python statement in {path}:{node.lineno}")


def check_stubs(errors: list[str]) -> None:
    native = "\n".join(
        path.read_text()
        for path in NATIVE_DIR.glob("*.rs")
        if not path.name.endswith("_tests.rs")
    )
    stubs = (PACKAGE / "__init__.pyi").read_text()
    functions = re.findall(r"#\[pyfunction\]\s*fn\s+(\w+)", native)
    for function in functions:
        if re.search(rf"^def\s+{re.escape(function)}\s*\(", stubs, re.MULTILINE) is None:
            errors.append(f"Python stub missing for pyfunction {function}")


def check_library_root(errors: list[str]) -> None:
    library = (NATIVE_DIR / "lib.rs").read_text().splitlines()
    for line_number, line in enumerate(library, start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("//!"):
            continue
        if re.fullmatch(r"(?:pub\s+)?mod\s+\w+;", stripped):
            continue
        if stripped.startswith("pub use "):
            continue
        errors.append(f"implementation in pdbiox-py/src/lib.rs:{line_number}")


def main() -> int:
    errors: list[str] = []
    check_dependencies(errors)
    check_python_source(errors)
    check_stubs(errors)
    check_library_root(errors)
    if errors:
        for error in errors:
            print(f"FAIL  {error}")
        return 1
    print("Python binding gates passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
