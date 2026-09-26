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
import re
import subprocess
import sys
import tarfile
import time
import tomllib
import urllib.error
import urllib.request
from email.utils import parsedate_to_datetime
from datetime import datetime, timezone
from pathlib import Path

INDEX = "https://index.crates.io"
# crates.io asks that automated clients identify themselves.
HEADERS = {"User-Agent": "molframe-release-workflow (github.com/miguelcsx/molframe)"}

# How long to wait for a just-published crate to appear in the index. A upload
# is not instant, and the crates that depend on it cannot be packaged until it
# is there.
INDEX_TIMEOUT = 600
INDEX_POLL = 5

# crates.io throttles new-crate uploads with a 429 that names the time to retry
# after. Waiting that out is worth it when short; a longer wait is better left
# to a re-run, which skips everything already published.
RATE_LIMIT_MAX_WAIT = 1800
RATE_LIMIT_SLACK = 5
RATE_LIMIT_ATTEMPTS = 5
DEPENDENCY_TABLES = ("dependencies", "build-dependencies", "dev-dependencies")


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


def packaged_dependencies(archive: Path, version: str) -> set[str]:
    """The crates a `.crate` archive depends on, per its normalised manifest.

    Cargo rewrites `Cargo.toml` when it packages: a path-only dev-dependency is
    dropped, and one that carries a version is kept and must resolve against the
    registry at publish time. Reading the result, rather than re-deriving
    Cargo's rules from `cargo metadata`, makes Cargo's own decision the source
    of truth for the order.
    """
    root = archive.name[: -len(".crate")]
    with tarfile.open(archive, "r:gz") as tar:
        member = tar.extractfile(f"{root}/Cargo.toml")
        if member is None:
            raise SystemExit(f"{archive.name} has no Cargo.toml")
        manifest = tomllib.loads(member.read().decode())

    tables = [manifest.get(table, {}) for table in DEPENDENCY_TABLES]
    for target in manifest.get("target", {}).values():
        tables.extend(target.get(table, {}) for table in DEPENDENCY_TABLES)

    names: set[str] = set()
    for table in tables:
        for key, spec in table.items():
            # A renamed dependency is keyed by its alias; `package` is the crate.
            names.add(spec.get("package", key) if isinstance(spec, dict) else key)
    return names


def publish_order(version: str) -> list[str]:
    """Package names to publish, dependencies first, deterministically.

    The graph comes from the manifests inside the packaged archives, so an edge
    exists exactly when Cargo kept the dependency. That is what keeps a
    version-carrying dev-dependency ahead of its dependent, and the path-only
    `molframe-query` -> `molframe-spatial` one out of a cycle.
    """
    suffix = f"-{version}.crate"
    archives = {
        archive.name[: -len(suffix)]: archive
        for archive in sorted(Path("target/package").glob("*.crate"))
    }
    graph = {
        name: sorted((packaged_dependencies(archive, version) & archives.keys()) - {name})
        for name, archive in sorted(archives.items())
    }
    return list(graphlib.TopologicalSorter(graph).static_order())


def rate_limit_wait(output: str) -> float | None:
    """Seconds to wait if `output` is a crates.io 429, else `None`."""
    if "429" not in output and "Too Many Requests" not in output:
        return None
    match = re.search(r"try again (?:after|at) ([A-Za-z]{3}, [^\n]*?GMT)", output)
    if match is None:
        return float(RATE_LIMIT_MAX_WAIT)
    retry_at = parsedate_to_datetime(match.group(1))
    return max(0.0, (retry_at - datetime.now(timezone.utc)).total_seconds())


def publish(name: str) -> None:
    """Publishes one crate, waiting out a short crates.io rate limit."""
    for attempt in range(1, RATE_LIMIT_ATTEMPTS + 1):
        result = subprocess.run(
            ["cargo", "publish", "-p", name, "--locked"], capture_output=True, text=True
        )
        sys.stdout.write(result.stdout)
        sys.stderr.write(result.stderr)
        if result.returncode == 0:
            return
        wait = rate_limit_wait(result.stdout + result.stderr)
        if wait is None:
            raise SystemExit(f"cargo publish failed for {name} ({result.returncode})")
        if wait > RATE_LIMIT_MAX_WAIT or attempt == RATE_LIMIT_ATTEMPTS:
            raise SystemExit(
                f"rate limited by crates.io while publishing {name}; retry in {wait:.0f}s. "
                f"Re-run the workflow later: crates already published are skipped."
            )
        print(f"rate limited; waiting {wait + RATE_LIMIT_SLACK:.0f}s before retrying {name}")
        time.sleep(wait + RATE_LIMIT_SLACK)


def main() -> int:
    with Path("Cargo.toml").open("rb") as handle:
        version = tomllib.load(handle)["workspace"]["package"]["version"]
    print(f"publishing workspace version {version}")

    checksums = local_checksums(version)
    if not checksums:
        raise SystemExit("no archives in target/package; the preflight must package first")

    order = publish_order(version)
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
        publish(name)
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
