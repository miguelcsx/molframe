# Releasing

Two registries, two workflows, both dispatched by hand and both authenticated by
OIDC rather than a stored token. Nothing in this repository holds a crates.io or
PyPI credential.

| Workflow | Registry | Trusted publisher |
|---|---|---|
| `.github/workflows/release.yml` | crates.io | owner `miguelcsx`, repo `molframe`, workflow `release.yml`, environment `release` |
| `.github/workflows/release-python.yml` | PyPI | owner `miguelcsx`, repo `molframe`, workflow `release-python.yml`, environment `release` |

Both require a GitHub environment named `release`. Add it under **Settings →
Environments**, and require a reviewer on it if a human should approve every
upload; the workflows already declare `environment: release`, and crates.io
checks that name as part of the publisher identity.

The environment also carries two secrets, which are the only credentials this
repository holds and the only ones a release ever uses:

| Secret | For |
|---|---|
| `RELEASE_GPG_PRIVATE_KEY` | the exported private key, ASCII-armoured |
| `RELEASE_GPG_KEY_ID` | that key's id, so the tag is signed by the intended key |

Export with `gpg --armor --export-secret-keys <KEY_ID>`. The key must have no
passphrase: `git tag -s` would otherwise block at the agent prompt, and there is
no terminal in a workflow to answer it. That is the same key already used for
commits, and it is protected by the environment rather than by a passphrase —
which is why requiring a reviewer on the `release` environment matters more here
than it would otherwise.

The tag step fails rather than falling back to an unsigned tag if the key is
missing. An unverifiable release tag is worse than a release that stops.

## The version is the input

Dispatch **Actions → Release → Run workflow** and give the version without the
leading `v`. The first step refuses to continue unless that string is exactly
`[workspace.package].version` in `Cargo.toml`, so a typo cannot publish under a
version nobody chose. `release-python.yml` additionally requires
`pyproject.toml`'s `version` to agree.

Bump first, in one commit that touches `Cargo.toml`, `Cargo.lock` and
`pyproject.toml`, and merge it to `main`. The release must come from a commit
that is an ancestor of `origin/main`; the workflow checks that itself.

## What the crates.io release does, in order

Every gate runs before the first upload, so a red run publishes nothing and
leaves the version number free.

1. Version input, packaging metadata, commit-on-main.
2. `cargo fmt --check`, `clippy -D warnings`, `cargo test --workspace`,
   doc tests, `RUSTDOCFLAGS="-D warnings" cargo doc --all-features`.
3. The full facade feature matrix.
4. `cargo publish --workspace --dry-run --locked`.
5. `cargo package`, then every archive is checked for size (under the 10 MB
   crates.io limit), for unexpected files, and for a packed `LICENSE`.
6. The Python test suite.
7. OIDC exchange for a scoped, short-lived crates.io token.
8. `.github/scripts/publish.py` publishes the 25 crates in dependency order.
9. Signed tag `vX.Y.Z`, pushed.
10. GitHub Release with generated notes.

The tag is created **last**. A tag asserts that a commit is the release, and a
half-published version is not one.

## Why publishing is a script and not `cargo publish --workspace`

`cargo publish --workspace` publishes in dependency order and waits for each
crate to reach the index, but it aborts the entire run on a version that already
exists — there is no skip and no `--idempotent`, which is still an open request
([rust-lang/cargo#13397](https://github.com/rust-lang/cargo/issues/13397)).

With twenty-five irreversible uploads that is the wrong failure mode. One flaky
upload halfway through would leave the workspace half-published, and the retry
would refuse to start because the crates that did upload now exist. A published
version cannot be replaced, only yanked.

`.github/scripts/publish.py` does the same work with the property that matters:
it is safe to re-run. It publishes one crate at a time, in an order derived from
`cargo metadata`, and treats a crate that is already on the registry at this
version **with the checksum the local archive has** as done. Anything else — a
different checksum at that version, an archive that was not built, an upload
that fails — stops the run before the next crate. So a failed release is resumed
by dispatching the same version again.

The order comes from the dependency edges that survive packaging. A
`dev-dependency` is stripped from a published manifest, which is why
`molframe-query` and `molframe-spatial` are not a cycle here; see the comment on
that edge in `crates/molframe-query/Cargo.toml`.

## First release only: the bootstrap

Trusted publishing is configured on crates that already exist, so `0.1.0` has to
be published once with a real token before any of the above applies.

```
cargo login                      # paste the token, stdin only; never as an argument
cargo publish --workspace --dry-run --locked

cargo package --workspace \
  --exclude molframe-bench --exclude molframe-py --exclude molframe-resource-bench \
  --no-verify --locked
python3 .github/scripts/publish.py
```

Use the script here too, not `cargo publish --workspace`. The bootstrap is the
run most likely to be interrupted, and it is the one where a half-published
workspace is most expensive to reason about. The script publishes the same 25
crates in the same order and can simply be run again.

Two prerequisites before the first upload, both website actions:

1. **A verified email address on the crates.io account.** crates.io refuses
   every publish with `A verified email address is required to publish crates to
   crates.io` until this is set, under https://crates.io/settings/profile.
2. **The bootstrap token itself**, from https://crates.io/settings/tokens. It
   needs `publish-new` for the first release, because none of the 25 names exist
   yet; `publish-update` alone cannot create them.

Then, for each of the 25 crates, add the trusted publisher on crates.io. Enable
**Trusted Publishing Only** for each once the setting exists for it, so a leaked
token cannot publish them.

When every crate has its publisher configured, retire the bootstrap: revoke the
token under https://crates.io/settings/tokens and run `cargo logout`. Revoking
is what invalidates it; `cargo logout` only forgets it locally.

## After publishing

Check the release from outside the workspace, where no `path = "…"` dependency
can hide a packaging mistake:

```bash
tmp=$(mktemp -d) && cd "$tmp"
cargo init --lib --name registry_smoke
cargo add molframe@0.1.0 && cargo check
cargo install molframe-cli --version 0.1.0 --locked && molframe --help
```

## Not yet wired: release notes and the benchmark register

`docs/`, `idea/` and `inspo/` are excluded from git, so no workflow reads them.
`molframe-bench`, `molframe-py` and `molframe-resource-bench` are
`publish = false` and are skipped by every publish path above.
