#!/usr/bin/env python3
"""Regenerate and validate PanIsoGuard's committed benchmark result artifacts.

The numeric claims about PanIsoGuard's accuracy used to live only in README /
docs/validation.md prose. This script makes the *self-contained* ones tracked and
reproducible: it re-runs each protocol's `ci_check.py --emit-metrics`, wraps the
result in the schema envelope (benchmark/results/schema.json), and writes
benchmark/results/<protocol>/metrics.json. Heavy protocols that need external
tooling (sqanti_sim, end2end, hg002) keep a hand-committed, schema-valid envelope
with status="transcribed_pending_tracked_run" until a tracked run replaces it.

  collect.py                 regenerate self-contained metrics + validate all envelopes
  collect.py --check         CI mode: fail if a regenerated metric drifts from the
                             committed file, or if any envelope is schema-invalid
  collect.py --pig PATH      path to the panisoguard binary (default ./build/panisoguard)

Exit code is non-zero on drift or schema failure, so it can be wired into CI/CTest.
"""
import json
import math
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent          # benchmark/results
BENCH = HERE.parent                             # benchmark
REPO = BENCH.parent                             # repo root

# Self-contained protocols: regenerated here, byte-checkable. (scope, provenance.)
# `requires` lists external tools whose absence makes the *tracked* metric
# irreproducible locally -> the protocol is skipped (not drifted/clobbered) so a
# contributor without that tool neither sees false drift nor downgrades the file.
SELF_CONTAINED = {
    "controlled_truth": {
        "scope": "ci_fixture",
        "data_provenance": "synthetic disjoint-window reference GTF (benchmark/controlled_truth/ci_check.py)",
        "command": "python3 benchmark/results/collect.py",
        "requires": [],
    },
    "synthetic_axes": {
        "scope": "ci_fixture",
        "data_provenance": "fully synthetic contig (benchmark/synthetic_axes/gen.py)",
        "command": "python3 benchmark/results/collect.py",
        "requires": ["samtools"],  # committed metric includes the BAM/mapping axis
    },
}


def tool_version(pig: str) -> str:
    try:
        out = subprocess.run([pig, "version"], capture_output=True, text=True)
        first = out.stdout.splitlines()[0].strip()      # "panisoguard 0.0.2"
        return first.split()[-1] if first else None
    except Exception:
        return None


def metrics_equal(a, b) -> bool:
    """Deep-compare metric values, treating NaN == NaN as equal so a degenerate
    (e.g. precision=nan) metric does not report spurious drift on every run."""
    if isinstance(a, float) and isinstance(b, float):
        return a == b or (math.isnan(a) and math.isnan(b))
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(metrics_equal(a[k], b[k]) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(metrics_equal(x, y) for x, y in zip(a, b))
    return a == b


ROOT = HERE.parent.parent  # the repository


def ruleset_version(pig: str) -> str:
    """The ruleset the binary applies, read from the provenance log of a run on the tiny example."""
    tiny = ROOT / "examples" / "tiny" / "data"
    with tempfile.TemporaryDirectory() as t:
        r = subprocess.run([pig, "adjudicate", "--classification", str(tiny / "classification.tsv"),
                            "--isoforms-gtf", str(tiny / "caller.gtf"), "--out-prefix", f"{t}/p"],
                           capture_output=True, text=True)
        if r.returncode != 0:
            sys.stderr.write(f"cannot read the ruleset from {pig}:\n{r.stderr}\n")
            sys.exit(1)
        for line in open(f"{t}/p.provenance.log"):
            if line.startswith("ruleset_version\t"):
                return line.rstrip("\n").split("\t", 1)[1]
    sys.exit("provenance.log has no ruleset_version line")


def regen_metrics(proto: str, pig: str) -> dict:
    """Run the protocol's ci_check.py with --emit-metrics and return the bare dict."""
    ci = BENCH / proto / "ci_check.py"
    with tempfile.TemporaryDirectory() as t:
        out = Path(t) / "m.json"
        r = subprocess.run([sys.executable, str(ci), pig, "--emit-metrics", str(out)],
                           capture_output=True, text=True)
        if r.returncode != 0 or not out.exists():
            sys.stderr.write(f"[{proto}] regeneration failed:\n{r.stdout}\n{r.stderr}\n")
            sys.exit(1)
        return json.loads(out.read_text())


_JSON_TYPES = {"string": str, "object": dict, "array": list, "boolean": bool, "null": type(None),
               "number": (int, float), "integer": int}


def validate_envelope(path: Path, schema: dict) -> list:
    """Dependency-free check of schema.json: required keys, no extra keys, each key's type (or
    list of types), enum, and minProperties."""
    errs = []
    try:
        env = json.loads(path.read_text())
    except Exception as e:
        return [f"{path}: not valid JSON ({e})"]
    if not isinstance(env, dict):
        return [f"{path}: envelope must be a JSON object"]
    for key in schema["required"]:
        if key not in env:
            errs.append(f"{path}: missing required key '{key}'")
    props = schema["properties"]
    if schema.get("additionalProperties") is False:
        for key in env:
            if key not in props:
                errs.append(f"{path}: key '{key}' is not in the schema")
    for key, spec in props.items():
        if key not in env:
            continue
        val = env[key]
        types = spec["type"] if isinstance(spec["type"], list) else [spec["type"]]
        ok = any(isinstance(val, _JSON_TYPES[t]) and not (t in ("number", "integer") and isinstance(val, bool))
                 for t in types)
        if not ok:
            errs.append(f"{path}: '{key}' must be of type {types}, got {type(val).__name__}")
            continue
        if "enum" in spec and val not in spec["enum"]:
            errs.append(f"{path}: '{key}' = {val!r} not in {spec['enum']}")
        if "minProperties" in spec and isinstance(val, dict) and len(val) < spec["minProperties"]:
            errs.append(f"{path}: '{key}' must have at least {spec['minProperties']} entries")
    return errs


def main() -> int:
    args = sys.argv[1:]
    check = "--check" in args
    pig = "./build/panisoguard"
    if "--pig" in args:
        pig = args[args.index("--pig") + 1]
    if not Path(pig).exists():
        # allow a bare name resolvable on PATH
        from shutil import which
        if which(pig) is None:
            sys.stderr.write(f"panisoguard binary not found: {pig} (build it, or pass --pig)\n")
            return 2

    schema = json.loads((HERE / "schema.json").read_text())
    ver = tool_version(pig)
    ruleset = ruleset_version(pig)
    drift, changed, skipped = [], [], []
    from shutil import which

    for proto, meta in SELF_CONTAINED.items():
        missing = [t for t in meta.get("requires", []) if which(t) is None]
        if missing:
            skipped.append(f"{proto} (missing {','.join(missing)})")
            continue
        metrics = regen_metrics(proto, pig)
        out = HERE / proto / "metrics.json"
        prev = json.loads(out.read_text()) if out.exists() else None

        if check:
            if prev is None:
                drift.append(f"{proto}: no committed metrics.json (run collect.py)")
            elif not metrics_equal(prev.get("metrics"), metrics):
                drift.append(f"{proto}: metrics drifted\n    committed={prev.get('metrics')}\n    fresh    ={metrics}")
            continue

        # write mode: preserve generated_utc when the numbers are unchanged
        same = prev is not None and metrics_equal(prev.get("metrics"), metrics)
        stamp = prev.get("generated_utc") if same else datetime.now(timezone.utc).replace(microsecond=0).isoformat()
        envelope = {
            "protocol": proto,
            "scope": meta["scope"],
            "status": "tracked",
            "tool_version": ver,
            "ruleset_version": ruleset,
            "generated_utc": stamp,
            "data_provenance": meta["data_provenance"],
            "command": meta["command"],
            "source": None,
            "metrics": metrics,
            "notes": "Self-contained correctness check; regenerated and drift-checked by collect.py.",
        }
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(envelope, indent=2, sort_keys=True) + "\n")
        if not same:
            changed.append(proto)

    # validate every committed envelope (self-contained + transcribed/pending)
    errs = []
    for mj in sorted(HERE.glob("*/metrics.json")):
        errs += validate_envelope(mj, schema)

    print(f"tool_version={ver}  protocols(self-contained)={list(SELF_CONTAINED)}  "
          f"committed_envelopes={len(list(HERE.glob('*/metrics.json')))}"
          + (f"  skipped={skipped}" if skipped else ""))
    if check:
        if drift:
            for d in drift:
                sys.stderr.write(f"DRIFT: {d}\n")
        if errs:
            for e in errs:
                sys.stderr.write(f"SCHEMA: {e}\n")
        if drift or errs:
            return 1
        print("collect.py --check: OK (no drift, all envelopes schema-valid)")
        return 0

    if errs:
        for e in errs:
            sys.stderr.write(f"SCHEMA: {e}\n")
        return 1
    print(f"wrote/updated: {changed or 'none (no metric changes)'}; all envelopes schema-valid")
    return 0


if __name__ == "__main__":
    sys.exit(main())
