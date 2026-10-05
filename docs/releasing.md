# Releasing PanIsoGuard

How to cut a release and ship it to bioconda. The version string is the single source of
truth in **three** files, kept in lockstep by CI (`scripts/check_version_sync.py`):

| file | field |
|------|-------|
| `CMakeLists.txt` | `project(panisoguard VERSION X.Y.Z ...)` |
| `python/pyproject.toml` | `version = "X.Y.Z"` |
| `recipes/bioconda/meta.yaml` | `{% set version = "X.Y.Z" %}` + `sha256:` |

## GitHub releases do not wait for bioconda

PanIsoGuard has been on bioconda since 0.0.4: the first submission,
[bioconda-recipes PR #65953](https://github.com/bioconda/bioconda-recipes/pull/65953), was
merged on 2026-09-30 after four months in review. A GitHub release is cut whenever a change is
ready; the bioconda package follows through a version-bump PR (see below).

Two things went wrong with the first submission. Watch for them with any bioconda PR, and with
the first submission of any new package (for example `panisoguard-report`):

- **A push after approval takes the PR out of the merge queue.** Mergify drops a queued PR
  when its branch is "manually updated", and only a maintainer can requeue it. #65953 was
  queued two minutes before a push of an optional review suggestion and lost a day. Before
  pushing to an approved PR, check whether it is already queued.
- **The first upload of a new package can miss a platform.** After the merge, each platform
  uploads on its own, and they race to create the package on anaconda.org; a loser fails with
  `Conflict: ('Owner bioconda already have a package named ...', 409)`. For 0.0.4 the
  linux-64 upload failed this way while both macOS builds landed. The usual fix is a PR that
  only bumps `build: number`; we opened one
  ([#69756](https://github.com/bioconda/bioconda-recipes/pull/69756)), but the 0.0.5 autobump
  ([#69871](https://github.com/bioconda/bioconda-recipes/pull/69871)) rebuilt every platform
  before it was reviewed, so it was closed. A new version fixes it the same way, and once the
  package exists the race cannot recur. After a first merge, check the upload checks on the
  merge commit and <https://anaconda.org/bioconda/panisoguard/files>.

## Cutting a release

```bash
V=0.0.4

# 1) finalize the changelog: rename "[Unreleased]" to "[X.Y.Z] - <date>" (edit CHANGELOG.md);
#    step 3's workflow copies this section into the release notes

# 2) bump the version in all THREE sources (sha256 is fixed in step 4)
sed -i -E "s/(project\(panisoguard VERSION )[0-9.]+/\1$V/" CMakeLists.txt
sed -i -E "s/^version = \"[0-9.]+\"/version = \"$V\"/" python/pyproject.toml
sed -i -E "s/(set version = \")[0-9.]+/\1$V/" recipes/bioconda/meta.yaml
python3 scripts/check_version_sync.py        # must print OK (checks the version string, not the sha)
git add CMakeLists.txt python/pyproject.toml recipes/bioconda/meta.yaml CHANGELOG.md
git commit -m "release: v$V"
git push origin main

# 3) tag. Pushing it runs .github/workflows/release.yml, which creates the GitHub Release:
#    notes = bioconda source block + this version's CHANGELOG section, plus the
#    bioconda-source.yaml asset. Do NOT `gh release create` by hand (it will already exist).
git tag -a "v$V" -m "PanIsoGuard v$V"
git push origin "v$V"
sleep 10   # let the Release run register before watching it
gh run watch "$(gh run list --workflow release.yml --limit 1 --json databaseId --jq '.[0].databaseId')"

# 4) pin the sha256 the workflow computed into the recipe
SHA=$(gh release download "v$V" -p bioconda-source.yaml -O - | awk '/sha256/{print $2}')
sed -i -E "s/^(  sha256: )[0-9a-f]{64}/\1$SHA/" recipes/bioconda/meta.yaml
git commit -am "recipe: pin v$V source sha256" && git push origin main
```

## Getting a release onto bioconda

- **Autobump (preferred).** The bioconda autobump bot
  watches GitHub releases and opens a version-bump PR on `bioconda-recipes` for you. Review
  it, comment `@BiocondaBot please add label` once CI is green, and wait for a maintainer.
- **Manual PR.** In your `bioconda-recipes` fork, edit `recipes/panisoguard/meta.yaml`
  (version + sha256 to match step 4, and `build: number` back to 0), open a PR, and add the
  `please review & merge` label via the bot.

## Optional: publish `panisoguard-report` to PyPI

The C++ binary on bioconda does **not** include the Python report tool. To distribute it:

```bash
python -m build python/                 # sdist + wheel -> python/dist/
python -m twine upload python/dist/*    # needs a PyPI account/token
```

A bioconda recipe for `panisoguard-report` (noarch python, `matplotlib` dep) can follow once
it is on PyPI; keep its version in lockstep with the C++ package.

## Checklist

- [ ] CHANGELOG section finalized with a date
- [ ] version bumped in all three files; `check_version_sync.py` OK
- [ ] `ctest --test-dir build` green (incl. `integration_config_equivalence`)
- [ ] tag pushed + GitHub Release created
- [ ] recipe `sha256` updated to the release tarball
- [ ] bioconda: autobump/manual PR opened & labelled; after the merge, the new version is on
      <https://anaconda.org/bioconda/panisoguard/files> for every platform
- [ ] (optional) report tool published to PyPI
