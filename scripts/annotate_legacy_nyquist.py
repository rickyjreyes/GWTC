#!/usr/bin/env python3
"""Read-only audit of historical binned GWTC peaks against the sampling Nyquist limit.

Do NOT update historical science CSVs. The archived n_star is related to k_best
by n_star = k_best * active_log_width / (2*pi). Reconstructing that width
lets historical rows be audited without rerunning or silently editing scans.
This is sampling aliasing only; no correction to historical p-values is made.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

INPUTS = (
    "tables/gwtc_scan_primary.csv",
    "tables/gwtc_bin_stress.csv",
    "tables/gwtc_master_results.csv",
)
FIELDS = ("source_file", "source_row", "stage", "variable", "subset", "bin_count",
          "historical_k_best", "historical_n_star", "recovered_active_log_width",
          "nyquist_k", "sampling_classification", "historical_verdict", "annotation")


def classify(row: dict, source: str, row_number: int) -> dict:
    out = dict(source_file=source, source_row=row_number, stage=row.get("stage", ""),
               variable=row.get("variable", ""), subset=row.get("subset", ""),
               bin_count=row.get("bin_count", ""), historical_k_best=row.get("k_best", ""),
               historical_n_star=row.get("n_star", ""), recovered_active_log_width="",
               nyquist_k="", sampling_classification="INSUFFICIENT_METADATA",
               historical_verdict=row.get("verdict_label", ""),
               annotation="No reinterpretation of archived statistics")
    try:
        bins = int(row["bin_count"])
        k = float(row["k_best"])
        n = float(row["n_star"])
        if bins < 2 or not all(math.isfinite(v) and v > 0 for v in (k, n)):
            return out
        width = 2 * math.pi * n / k
        kn = math.pi * bins / width
        # Match strict code guard: k < k_N*(1-1e-10).
        classification = "RESOLVED" if k < kn * (1 - 1e-10) else "ALIASED_OR_NYQUIST"
        out.update(recovered_active_log_width=format(width, ".15g"),
                   nyquist_k=format(kn, ".15g"),
                   sampling_classification=classification,
                   annotation=("Historical maximum is above/at the binning Nyquist; NOT a resolved mode"
                               if classification != "RESOLVED"
                               else "Historical peak below Nyquist; not a new statistical validation"))
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        pass
    return out


def audit(root: Path, files=INPUTS) -> list[dict]:
    rows = []
    for rel in files:
        with (root / rel).open(newline="", encoding="utf-8") as fh:
            for i, row in enumerate(csv.DictReader(fh), 2):
                rows.append(classify(row, rel, i))
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--out", type=Path, default=None,
                   help="Separate derived annotation CSV; historical inputs remain unchanged")
    args = p.parse_args(argv)
    rows = audit(args.root)
    dest = args.out or args.root / "results" / "legacy_nyquist_annotations.csv"
    if dest.resolve() in {(args.root / s).resolve() for s in INPUTS}:
        p.error("Refusing to overwrite historical input")
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    counts = {name: sum(r["sampling_classification"] == name for r in rows)
              for name in ("RESOLVED", "ALIASED_OR_NYQUIST", "INSUFFICIENT_METADATA")}
    print(f"Read-only Nyquist annotation: {dest} / {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
