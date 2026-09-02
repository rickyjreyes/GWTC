#!/usr/bin/env python3
"""Generate publication figures for the strict GWTC frozen-mode evidence chain.

The figures are derived only from committed frozen/result artifacts. The script
intentionally does not use the historical exploratory scan as primary evidence.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

mpl.rcParams["svg.hashsalt"] = "gwtc-publication-figures-v1"
mpl.rcParams["svg.fonttype"] = "none"

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "tables"


def _save(fig: plt.Figure, stem: Path, formats: list[str]) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        kwargs = {"bbox_inches": "tight"}
        if fmt == "svg":
            kwargs["metadata"] = {"Date": None}
        fig.savefig(stem.with_suffix(f".{fmt}"), **kwargs)
    plt.close(fig)


def _load_json(name: str) -> dict:
    with (TABLES / name).open("r", encoding="utf-8") as fh:
        return json.load(fh)


def make_protocol(out: Path, formats: list[str]) -> None:
    fig, ax = plt.subplots(figsize=(10, 2.4))
    ax.axis("off")
    steps = [
        "Freeze split",
        "Select background\n(training only)",
        "Select residual mode\n(training only)",
        "Freeze model",
        "Evaluate GWTC-5\nholdout once",
    ]
    xs = [0.06, 0.27, 0.49, 0.70, 0.90]
    for x, label in zip(xs, steps):
        ax.text(
            x,
            0.5,
            label,
            ha="center",
            va="center",
            fontsize=11,
            bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="black"),
        )
    for a, b in zip(xs[:-1], xs[1:]):
        ax.annotate(
            "",
            xy=(b - 0.065, 0.5),
            xytext=(a + 0.065, 0.5),
            arrowprops=dict(arrowstyle="->"),
        )
    ax.text(
        0.5,
        0.92,
        "Strict frozen holdout protocol",
        ha="center",
        va="center",
        fontsize=14,
        fontweight="bold",
    )
    _save(fig, out / "publication_protocol", formats)


def make_v2_frozen_model(out: Path, formats: list[str]) -> None:
    v2 = _load_json("gwtc_frozen_mode.json")
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(v2["bin_centers"], v2["baseline_probabilities"], label="Degree-7 smooth baseline")
    ax.plot(
        v2["bin_centers"],
        v2["residual_probabilities"],
        label=f"Frozen residual model (k = {v2['k_star']:.4f})",
    )
    ax.set_xlabel(r"$\ell = \ln(M_{\rm chirp})$")
    ax.set_ylabel("frozen bin probability")
    ax.set_title("V2 training-frozen chirp-mass model")
    ax.legend(frameon=False)
    fig.tight_layout()
    _save(fig, out / "publication_v2_frozen_model", formats)


def make_v3_frequency_scan(out: Path, formats: list[str]) -> None:
    v3 = _load_json("gwtc_v3_frozen_unbinned_kde_mode.json")
    scan = pd.DataFrame(v3["scan"])
    zoom = scan[(scan["k"] >= 7.5) & (scan["k"] <= 12.5)].copy()
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    ax.plot(zoom["k"], zoom["train_delta_2logl"], marker="o")
    ax.axvspan(9.5, 10.0, alpha=0.15, label="Prospective band 9.5-10.0")
    ax.axvline(9.7, linestyle="--", label="Prospective k = 9.7")
    ax.scatter(
        [v3["k_star"]],
        [float(scan.loc[scan["k"].sub(v3["k_star"]).abs().idxmin(), "train_delta_2logl"])],
        s=70,
        zorder=5,
        label="V3 training maximum",
    )
    ax.set_xlabel("log-frequency k")
    ax.set_ylabel(r"training $\Delta 2\log L$")
    ax.set_title("V3 training-only frequency scan near the frozen mode")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    _save(fig, out / "publication_v3_frequency_scan", formats)


def make_cross_method_frequency(out: Path, formats: list[str]) -> None:
    v2 = _load_json("gwtc_frozen_mode.json")
    v3 = _load_json("gwtc_v3_frozen_unbinned_kde_mode.json")
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    ax.axvspan(9.5, 10.0, alpha=0.15, label="Prospective interval")
    ax.axvline(9.7, linestyle="--", label="Prospective center")
    ks = [v2["k_star"], v3["k_star"]]
    ax.scatter(ks, [1, 0], s=90)
    ax.set_yticks([1, 0], ["V2 binned Poisson", "V3 unbinned KDE"])
    ax.set_xlim(9.3, 10.1)
    ax.set_xlabel("selected log-frequency k")
    ax.set_title("Cross-method frequency agreement")
    ax.annotate(f"{ks[0]:.4f}", (ks[0], 1), xytext=(6, 8), textcoords="offset points")
    ax.annotate(f"{ks[1]:.4f}", (ks[1], 0), xytext=(6, 8), textcoords="offset points")
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    fig.tight_layout()
    _save(fig, out / "publication_cross_method_frequency", formats)


def make_v4_structured_nulls(out: Path, formats: list[str]) -> None:
    v4 = pd.read_csv(TABLES / "gwtc_v4_structured_null_result.csv")
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot(v4["selection_gamma"], v4["p_structured_null"], marker="o")
    ax.axhline(0.05, linestyle="--", label="Predeclared p = 0.05 threshold")
    ax.set_yscale("log")
    ax.set_xlabel(r"selection exponent $\gamma$")
    ax.set_ylabel("empirical structured-null p")
    ax.set_title("V4 frozen structured non-periodic null challenge")
    for _, row in v4.iterrows():
        ax.annotate(
            f"{row['p_structured_null']:.4g}",
            (row["selection_gamma"], row["p_structured_null"]),
            xytext=(5, 7),
            textcoords="offset points",
        )
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    _save(fig, out / "publication_v4_structured_nulls", formats)


def make_evidence_summary(out: Path, formats: list[str]) -> None:
    v2 = pd.read_csv(TABLES / "gwtc_frozen_holdout_result.csv").iloc[0]
    v3 = pd.read_csv(TABLES / "gwtc_v3_unbinned_kde_holdout_result.csv").iloc[0]
    v4 = pd.read_csv(TABLES / "gwtc_v4_structured_null_result.csv")
    names = ["V2 frozen holdout", "V3 KDE robustness", "V4 structured null\n(worst case)"]
    ps = [float(v2["holdout_p"]), float(v3["p_holdout"]), float(v4["p_structured_null"].max())]
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    ax.scatter(range(3), ps, s=90)
    ax.axhline(0.05, linestyle="--", label="p = 0.05")
    ax.set_yscale("log")
    ax.set_xticks(range(3), names)
    ax.set_ylabel("empirical p under declared null")
    ax.set_title("Frozen-mode evidence across successive tests")
    for i, p in enumerate(ps):
        ax.annotate(f"{p:.4g}", (i, p), xytext=(5, 7), textcoords="offset points")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    _save(fig, out / "publication_evidence_summary", formats)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "figures" / "publication",
        help="output directory (default: figures/publication)",
    )
    parser.add_argument(
        "--formats",
        default="svg,png",
        help="comma-separated output formats (default: svg,png)",
    )
    args = parser.parse_args()
    formats = [x.strip().lower() for x in args.formats.split(",") if x.strip()]
    unsupported = set(formats) - {"svg", "png", "pdf"}
    if unsupported:
        raise SystemExit(f"unsupported formats: {sorted(unsupported)}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    make_protocol(args.output_dir, formats)
    make_v2_frozen_model(args.output_dir, formats)
    make_v3_frequency_scan(args.output_dir, formats)
    make_cross_method_frequency(args.output_dir, formats)
    make_v4_structured_nulls(args.output_dir, formats)
    make_evidence_summary(args.output_dir, formats)
    print(f"wrote publication figures to {args.output_dir}")


if __name__ == "__main__":
    main()
