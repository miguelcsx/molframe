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
    for path in PACKAGE.rglob("*.py"):
        tree = ast.parse(path.read_text(), filename=str(path))
        imported: list[str] = []
        exported: list[str] = []
        for index, node in enumerate(tree.body):
            if isinstance(node, ast.Import):
                imported.extend(alias.asname or alias.name for alias in node.names)
                continue
            if isinstance(node, ast.ImportFrom):
                imported.extend(alias.asname or alias.name for alias in node.names)
                continue
            if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant):
                if index == 0 and isinstance(node.value.value, str):
                    continue
            if isinstance(node, ast.Assign) and all(
                isinstance(target, ast.Name) and target.id == "__all__"
                for target in node.targets
            ):
                if not isinstance(node.value, (ast.List, ast.Tuple)):
                    errors.append(f"non-literal __all__ in {path}:{node.lineno}")
                    continue
                for element in node.value.elts:
                    if not isinstance(element, ast.Constant) or not isinstance(
                        element.value, str
                    ):
                        errors.append(f"non-string __all__ entry in {path}:{element.lineno}")
                        continue
                    exported.append(element.value)
                continue
            errors.append(f"executable Python statement in {path}:{node.lineno}")
        check_exports(path, imported, exported, errors)
        check_stub_module_exports(path, exported, errors)


def check_stub_module_exports(
    path: Path, exported: list[str], errors: list[str]
) -> None:
    """Require every public symbol declared by a module stub at runtime."""
    if path.name == "__init__.py":
        return
    stub = path.with_suffix(".pyi")
    if not stub.exists():
        return
    tree = ast.parse(stub.read_text(), filename=str(stub))
    declared = {
        node.name
        for node in tree.body
        if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef))
        and not node.name.startswith("_")
    }
    missing = sorted(declared.difference(exported))
    if missing:
        errors.append(f"stub symbols absent from runtime module {path}: {missing}")


def check_exports(
    path: Path, imported: list[str], exported: list[str], errors: list[str]
) -> None:
    duplicate_exports = sorted(
        name for name in set(exported) if exported.count(name) > 1
    )
    missing_exports = sorted(set(imported).difference(exported))
    unknown_exports = sorted(set(exported).difference(imported))
    if duplicate_exports:
        errors.append(f"duplicate exports in {path}: {duplicate_exports}")
    if missing_exports:
        errors.append(f"imports absent from __all__ in {path}: {missing_exports}")
    if unknown_exports:
        errors.append(f"__all__ entries not imported in {path}: {unknown_exports}")


def check_stubs(errors: list[str]) -> None:
    native = "\n".join(
        path.read_text()
        for path in NATIVE_DIR.rglob("*.rs")
        if not path.name.endswith("_tests.rs")
    )
    stubs = "\n".join(path.read_text() for path in PACKAGE.rglob("*.pyi"))
    stub_classes = collect_stub_classes(errors)
    functions = re.findall(
        r'#\[pyfunction(?:\(name\s*=\s*"(\w+)"\))?\]\s*'
        r'(?:#\[pyo3\([^\]]*\)\]\s*)?(?:pub\(crate\)\s+)?fn\s+(\w+)',
        native,
    )
    for exposed_name, rust_name in functions:
        function = exposed_name or rust_name
        if re.search(rf"^def\s+{re.escape(function)}\s*\(", stubs, re.MULTILINE) is None:
            errors.append(f"Python stub missing for pyfunction {function}")
    registered = set(re.findall(r"add_class::<(\w+)>", native))
    rust_classes = {
        rust: python
        for rust, python in collect_rust_classes(native).items()
        if rust in registered
    }
    for rust_name, python_name in rust_classes.items():
        if python_name not in stub_classes:
            errors.append(f"Python stub missing for pyclass {python_name}")
    for rust_name, body in pymethod_blocks(native):
        python_name = rust_classes.get(rust_name)
        if python_name is None:
            continue
        members = stub_classes.get(python_name)
        if members is None:
            continue
        for method in exposed_methods(body):
            if method not in members:
                errors.append(f"Python stub missing {python_name}.{method}")


def collect_stub_classes(errors: list[str]) -> dict[str, set[str]]:
    classes: dict[str, set[str]] = {}
    for path in PACKAGE.rglob("*.pyi"):
        try:
            tree = ast.parse(path.read_text(), filename=str(path))
        except SyntaxError as error:
            errors.append(f"invalid Python stub {path}:{error.lineno}: {error.msg}")
            continue
        for node in tree.body:
            if not isinstance(node, ast.ClassDef):
                continue
            members = classes.setdefault(node.name, set())
            for item in node.body:
                if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    members.add(item.name)
                elif isinstance(item, ast.AnnAssign) and isinstance(item.target, ast.Name):
                    members.add(item.target.id)
                elif isinstance(item, ast.Assign):
                    members.update(
                        target.id for target in item.targets if isinstance(target, ast.Name)
                    )
    return classes


def collect_rust_classes(native: str) -> dict[str, str]:
    classes: dict[str, str] = {}
    pattern = re.compile(
        r"#\[pyclass\((?P<options>.*?)\)\]\s*"
        r"(?:#\[[^\]]*\]\s*)*"
        r"(?:pub(?:\(crate\))?\s+)?(?:struct|enum)\s+(?P<rust>\w+)",
        re.DOTALL,
    )
    for match in pattern.finditer(native):
        options = match.group("options")
        named = re.search(r'name\s*=\s*"(\w+)"', options)
        classes[match.group("rust")] = named.group(1) if named else match.group("rust")
    return classes


def pymethod_blocks(native: str) -> list[tuple[str, str]]:
    blocks: list[tuple[str, str]] = []
    header = re.compile(r"#\[pymethods\]\s*impl\s+(\w+)\s*\{")
    for match in header.finditer(native):
        opening = native.find("{", match.start())
        closing = matching_brace(native, opening)
        if closing is not None:
            blocks.append((match.group(1), native[opening + 1 : closing]))
    return blocks


def matching_brace(source: str, opening: int) -> int | None:
    depth = 0
    for index in range(opening, len(source)):
        character = source[index]
        if character == "{":
            depth += 1
        elif character == "}":
            depth -= 1
            if depth == 0:
                return index
    return None


def exposed_methods(body: str) -> set[str]:
    methods: set[str] = set()
    pattern = re.compile(
        r"(?P<attributes>(?:\s*#\[[^\]]*\]\s*)*)"
        r"(?:pub(?:\(crate\))?\s+)?(?:const\s+)?fn\s+(?P<name>\w+)",
    )
    for match in pattern.finditer(body):
        attributes = match.group("attributes")
        name = match.group("name")
        if "#[new]" in attributes:
            methods.add("__init__")
            continue
        renamed = re.search(r'#\[pyo3\([^\]]*name\s*=\s*"(\w+)"', attributes)
        methods.add(renamed.group(1) if renamed else name)
    return methods


def check_library_root(errors: list[str]) -> None:
    library = (NATIVE_DIR / "lib.rs").read_text().splitlines()
    for line_number, line in enumerate(library, start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("//!"):
            continue
        if stripped == "#![deny(unsafe_op_in_unsafe_fn)]":
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
