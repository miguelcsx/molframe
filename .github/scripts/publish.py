#!/usr/bin/env python3
"""Publish every member of the workspace to crates.io, in dependency order.

`cargo publish --workspace` cannot be resumed. It publishes in dependency order
and waits for each crate to reach the index, but a version that already exists
aborts the whole run, and there is no `--idempotent` to ask otherwise
(rust-lang/cargo#13397). With twenty-five irreversible uploads, one flaky
upload mid-run would leave the workspace half-published and the retry would
refuse to start.

This script does the same thing with the property that matters: it is safe to
re-run. Each crate is published on its own, and a crate already on the registry
at this version, with the checksum the locally built archive has, counts as
done. Anything else -- a different checksum at that version, a manifest that
disagrees with the workspace version, an upload that fails -- stops the run
before the next crate, so the failure is narrow and the resume is exact.
"""

from __future__ import annotations

import graphlib
import hashlib
import json
import subprocess
import sys
import time
import tomllib
import urllib.error
import urllib.request
from pathlib import Path

INDEX = "https://index.crates.io"
# crates.io asks that automated clients identify themselves.
HEADERS = {"User-Agent": "molframe-release-workflow (github.com/miguelcsx/molframe)"}

# How long to wait for a just-published crate to appear in the index. A upload
# is not instant, and the crates that depend on it cannot be packaged until it
# is there.
INDEX_TIMEOUT = 600
INDEX_POLL = 5


def run(*args: str) -> str:
    """Runs a command, returning stdout, and fails loudly on a bad exit."""
    result = subprocess.run(args, capture_output=True, text=True)
    if result.returncode != 0:
        sys.stderr.write(result.stdout)
        sys.stderr.write(result.stderr)
        raise SystemExit(f"command failed ({result.returncode}): {' '.join(args)}")
    return result.stdout


def index_path(name: str) -> str:
    """The sparse-index path for a crate name.

    One- and two-character names are sharded by length, three by their first
    character, and anything longer by the first two and the next two.
    """
    lowered = name.lower()
    if len(lowered) == 1:
        return f"1/{lowered}"
    if len(lowered) == 2:
        return f"2/{lowered}"
    if len(lowered) == 3:
        return f"3/{lowered[0]}/{lowered}"
    return f"{lowered[:2]}/{lowered[2:4]}/{lowered}"


def fetch(url: str) -> bytes | None:
    """Fetches a URL, returning `None` for a 404 rather than raising."""
    request = urllib.request.Request(url, headers=HEADERS)
    try:
        with urllib.request.urlopen(request) as response:
            return response.read()
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return None
        raise


def published_checksum(name: str, version: str) -> str | None:
    """The checksum crates.io holds for exactly this version, or `None`.

    Read from the sparse index rather than the web API, because the index is
    what `cargo` itself resolves against: agreement with it is agreement with
    what a consumer's build will see. Its `cksum` is the SHA-256 of the `.crate`
    archive, which is the same digest this script computes locally.
    """
    body = fetch(f"{INDEX}/{index_path(name)}")
    if body is None:
        return None
    for line in body.decode().splitlines():
        if not line.strip():
            continue
        entry = json.loads(line)
        if entry.get("vers") == version:
            checksum = entry.get("cksum")
            return checksum if isinstance(checksum, str) else None
    return None


def await_published(name: str, version: str) -> str:
    """Waits for a version to reach the index, returning its checksum."""
    deadline = time.monotonic() + INDEX_TIMEOUT
    while time.monotonic() < deadline:
        checksum = published_checksum(name, version)
        if checksum is not None:
            return checksum
        time.sleep(INDEX_POLL)
    raise SystemExit(f"{name} {version} did not reach the index within {INDEX_TIMEOUT}s")


def local_checksums(version: str) -> dict[str, str]:
    """The SHA-256 of each freshly built archive, keyed by crate name."""
    checksums: dict[str, str] = {}
    suffix = f"-{version}.crate"
    for archive in sorted(Path("target/package").glob("*.crate")):
        digest = hashlib.sha256(archive.read_bytes()).hexdigest()
        checksums[archive.name[: -len(suffix)]] = digest
    return checksums


def publish_order() -> list[str]:
    """Package names to publish, dependencies first.

    `publish = false` members are dropped, and so are dev-dependency edges: a
    dev-dependency is stripped from the packaged manifest, so it cannot
    constrain the order of two published crates. That is what keeps
    `molframe-query` and `molframe-spatial` from being a cycle here.
    """
    metadata = json.loads(run("cargo", "metadata", "--format-version", "1", "--no-deps"))
    packages = {package["name"]: package for package in metadata["packages"]}
    publishable = {name for name, package in packages.items() if package["publish"] != []}
    graph: dict[str, set[str]] = {}
    for name in publishable:
        dependencies = set()
        for dependency in packages[name]["dependencies"]:
            if dependency.get("kind") == "dev":
                continue
            if dependency["name"] in publishable and dependency.get("path"):
                dependencies.add(dependency["name"])
        dependencies.discard(name)
        graph[name] = dependencies
    return [name for name in graphlib.TopologicalSorter(graph).static_order() if name in publishable]


def main() -> int:
    with Path("Cargo.toml").open("rb") as handle:
        version = tomllib.load(handle)["workspace"]["package"]["version"]
    print(f"publishing workspace version {version}")

    checksums = local_checksums(version)
    if not checksums:
        raise SystemExit("no archives in target/package; the preflight must package first")

    order = publish_order()
    missing = sorted(set(order) - set(checksums))
    if missing:
        raise SystemExit(f"no archive was built for {', '.join(missing)}")

    published = 0
    skipped = 0
    for position, name in enumerate(order, start=1):
        expected = checksums[name]
        prefix = f"[{position}/{len(order)}]"
        existing = published_checksum(name, version)
        if existing is not None:
            if existing != expected:
                raise SystemExit(
                    f"{name} {version} is already on crates.io with checksum {existing}, "
                    f"but the local archive hashes to {expected}. The version is taken and "
                    f"cannot be replaced; bump the workspace version instead."
                )
            print(f"{prefix} {name} {version} already published, skipping")
            skipped += 1
            continue

        print(f"{prefix} publishing {name} {version}")
        subprocess.run(["cargo", "publish", "-p", name, "--locked"], check=True)
        actual = await_published(name, version)
        if actual != expected:
            raise SystemExit(
                f"{name} {version} reached the index with checksum {actual}, but the local "
                f"archive hashes to {expected}"
            )
        published += 1

    print(f"done: {published} published, {skipped} already present")
    return 0


if __name__ == "__main__":
    sys.exit(main())
