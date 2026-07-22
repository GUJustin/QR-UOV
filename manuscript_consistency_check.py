#!/usr/bin/env python3
"""Mechanical consistency checks for the corrected QR-UOV manuscript."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TEX = (ROOT / "main.tex").read_text()

spec = importlib.util.spec_from_file_location(
    "bounds", ROOT / "qruov_random_projective_slice_bounds_v13.py"
)
assert spec and spec.loader
bounds = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = bounds
spec.loader.exec_module(bounds)
rows = [bounds.row(p) for p in bounds.PARAMS]

# Principal numerical tables.
for r in rows:
    required = [
        f"{r['fast_direct']:.2f}",
        f"{r['fast_uniform']:.2f}",
        f"{r['fallback_direct']:.2f}",
        f"{r['fallback_uniform']:.2f}",
        f"{r['fallback_cubic_uniform']:.2f}",
        f"{r['pseudo_direct']:.2f}",
        f"{r['pseudo_uniform']:.2f}",
    ]
    for value in required:
        assert value in TEX, (r["level"], value)

# No stale unsupported sensitivity columns or old fixed-slice headline caveat.
for forbidden in (
    "Best $\\omega$",
    "Schoolbook",
    "Key outside a proper algebraic exceptional set",
    "unquantified density",
    "regularized residual",
):
    assert forbidden not in TEX, forbidden

# Main labels must be unique and all bibliography keys used.
import re
from collections import Counter
labels = re.findall(r"\\label\{([^}]+)\}", TEX)
assert not [x for x, n in Counter(labels).items() if n > 1]

cites: list[str] = []
for match in re.finditer(r"\\cite(?:\[[^]]*\])?\{([^}]+)\}", TEX):
    cites.extend(x.strip() for x in match.group(1).split(","))
bib = set(re.findall(r"\\bibitem\{([^}]+)\}", TEX))
assert set(cites) == bib, (set(cites) - bib, bib - set(cites))

print("PASS: manuscript numbers, labels, citations, and stale-claim checks agree.")
