# Releasing PanIsoGuard

How to cut a release and ship it to bioconda. The version string is the single source of
truth in **three** files, kept in lockstep by CI (`scripts/check_version_sync.py`):

| file | field |
|------|-------|
| `CMakeLists.txt` | `project(panisoguard VERSION X.Y.Z ...)` |
| `python/pyproject.toml` | `version = "X.Y.Z"` |
| `recipes/bioconda/meta.yaml` | `{% set version = "X.Y.Z" %}` + `sha256:` |

## First-time status (v0.0.3)

The **first** bioconda submission is [bioconda-recipes PR #65953](https://github.com/bioconda/bioconda-recipes/pull/65953)
(`Add panisoguard 0.0.3`). It is OPEN / `REVIEW_REQUIRED` / `BLOCKED` — CI is green and the
`please review & merge` label is applied; it is waiting for a bioconda maintainer to
approve+merge (first submissions are queue-gated; days–weeks is normal). **Nothing is
required from us** until it merges. Etiquette: only ping after ~a week of no activity, with
a short friendly note — do not bump the recipe version or open a competing PR while it is
in review.

## Cutting the next release (e.g. v0.0.4) — do this AFTER v0.0.3 merges

> The next release will carry the **consensus axis** and everything else committed since
> the `v0.0.3` git tag (see [CHANGELOG.md](../CHANGELOG.md)). Hold it until v0.0.3 is on the
> bioconda channel so the first review isn't thrashed.

```bash
V=0.0.4

# 0) confirm v0.0.3 is actually published on the channel
conda search -c bioconda panisoguard

# 1) finalize the changelog: rename "[Unreleased] — targets 0.0.4" to "[0.0.4] - <date>"
#    (edit CHANGELOG.md)

# 2) bump the version in all THREE sources (sha256 is fixed in step 6)
sed -i -E "s/(project\(panisoguard VERSION )[0-9.]+/\1$V/" CMakeLists.txt
sed -i -E "s/^version = \"[0-9.]+\"/version = \"$V\"/" python/pyproject.toml
sed -i -E "s/(set version = \")[0-9.]+/\1$V/" recipes/bioconda/meta.yaml
python3 scripts/check_version_sync.py        # must print OK (checks the version string, not the sha)

# 3) commit + push the bump
git add CMakeLists.txt python/pyproject.toml recipes/bioconda/meta.yaml CHANGELOG.md
git commit -m "release: v$V"
git push origin main

# 4) tag + GitHub Release (this is what produces the source tarball the recipe points at)
git tag -a "v$V" -m "PanIsoGuard v$V"
git push origin "v$V"
gh release create "v$V" --title "v$V" --notes-from-tag   # or paste the CHANGELOG section

# 5) get the release tarball's sha256
curl -sL "https://github.com/jibeomko/PanIsoGuard/archive/refs/tags/v$V.tar.gz" | sha256sum

# 6) put that sha256 into recipes/bioconda/meta.yaml (replace the sha256: line), commit, push
#    git commit -am "recipe: pin v$V source sha256" && git push origin main
```

## Getting v0.0.4 onto bioconda

Once the GitHub Release exists, **either**:

- **Autobump (preferred, automatic after the first merge).** The bioconda autobump bot
  watches GitHub releases and opens a version-bump PR on `bioconda-recipes` for you. Review
  it, comment `@BiocondaBot please add label` once CI is green, and wait for a maintainer.
- **Manual PR.** In your `bioconda-recipes` fork, edit `recipes/panisoguard/meta.yaml`
  (version + sha256 to match step 5–6), open a PR, and add the `please review & merge`
  label via the bot.

## Optional: publish `panisoguard-report` to PyPI

The C++ binary on bioconda does **not** include the Python report tool. To distribute it:

```bash
python -m build python/                 # sdist + wheel -> python/dist/
python -m twine upload python/dist/*    # needs a PyPI account/token
```

A bioconda recipe for `panisoguard-report` (noarch python, `matplotlib` dep) can follow once
it is on PyPI; keep its version in lockstep with the C++ package.

## Checklist

- [ ] v0.0.3 merged & visible via `conda search -c bioconda panisoguard`
- [ ] CHANGELOG section finalized with a date
- [ ] version bumped in all three files; `check_version_sync.py` OK
- [ ] `ctest --test-dir build` green (incl. `integration_config_equivalence`)
- [ ] tag pushed + GitHub Release created
- [ ] recipe `sha256` updated to the release tarball
- [ ] bioconda PR (autobump or manual) opened & labelled
- [ ] (optional) report tool published to PyPI
