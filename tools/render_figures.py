#!/usr/bin/env python3
"""Render saved native Fortran fields; no scientific calculation or selection.

Requires matplotlib. Run from any directory. Existing original figures are never
read or written. --overwrite permits revising only this new derived figure set.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D

P = Path(__file__).resolve().parents[1]
E = P
R = E / "runs/correctness_01/work/results"
A = E / "runs/paper_minimum_aggregate_02"
OUT = P / "build/rendered_figures"
MAP = OUT / "figure_source_map.json"
BLUE, ORANGE, TEAL = "#0077BB", "#EE7733", "#009988"
DIRECTIONS = (("nominal", "-", "Nominal"),
              ("SAS_conservative", "--", "SAS high / 2D low"),
              ("SPACE_conservative", ":", "SAS low / 2D high"))
TRACE = {"operation": "Document rendering of saved scalar fields only",
         "limitations": ["No new science run or error/frontier recomputation.",
                         "Budget traces cover only the explicitly measured finite candidate lists; omitted intermediate nodes can change crossings and gains."],
         "inputs": {}, "figures": {}}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    TRACE["inputs"][str(path.relative_to(P))] = sha(path)
    with path.open(newline="") as f:
        return [dict(row, _csv_line=i) for i, row in enumerate(csv.DictReader(f), 2)]


def source(path):
    TRACE["inputs"][str(path.relative_to(P))] = sha(path)


def numbers(rows, key):
    return [float(row[key]) for row in rows]


def save(fig, name, info):
    fig.savefig(OUT / (name + ".png"), dpi=300)
    fig.savefig(OUT / (name + ".pdf"), metadata={"CreationDate": None, "ModDate": None})
    TRACE["figures"][name] = dict(info, outputs={
        str((OUT / (name + ext)).relative_to(P)): sha(OUT / (name + ext))
        for ext in (".png", ".pdf")})
    plt.close(fig)


def simple_axes(ax):
    ax.grid(alpha=.18, linewidth=.5)
    ax.tick_params(direction="out")


def figure1():
    rows = read(R / "case1_profiles.csv")
    assert len(rows) == 150
    fig, ax = plt.subplots(figsize=(6.9, 4.0), layout="constrained")
    for x, color, marker in ((0, BLUE, "o"), (2, ORANGE, "s"), (4, TEAL, "^")):
        rr = sorted([r for r in rows if float(r["x"]) == x], key=lambda r: float(r["z"]))
        assert len(rr) == 50 and all(float(r["y"]) == 0 for r in rr)
        ax.plot(numbers(rr, "z"), numbers(rr, "SAS"), color=color, lw=1.6,
                label=fr"SAS, $x={x}$")
        ax.plot(numbers(rr, "z"), numbers(rr, "partition120"), color=color,
                ls="--", lw=.9, marker=marker, markevery=5, ms=4, mfc="white", mew=.9,
                label=fr"Independent reference, $x={x}$")
    ax.set(xlabel=r"Depth $z$", ylabel=r"Pulse-integrated deficit $\Theta_{\mathrm{def}}$",
           xlim=(.1, 1.5), title=r"Uniform patch: $y=0$, $t_p=100$, $\alpha=1$")
    ax.legend(frameon=False, ncol=2, fontsize=8)
    simple_axes(ax)
    save(fig, "figure1", {"source": "R/case1_profiles.csv", "rows": 150,
         "fields": ["z", "SAS", "partition120"], "reference_markers": "Every fifth saved sample; all 50 saved samples per line plotted in both traces",
         "caption": "Pulse-integrated thermal deficit on the three original vertical sampling lines. SAS curves and partition-refined independent spatial references use the same 150 saved positions."})


def figure2():
    rows = read(R / "case2_profiles.csv")
    fig, ax = plt.subplots(figsize=(6.9, 4.0), layout="constrained")
    for m, color, style in ((2, ORANGE, "--"), (4, TEAL, "-."), (8, BLUE, "-")):
        rr = sorted([r for r in rows if int(r["M"]) == m], key=lambda r: float(r["r"]))
        assert len(rr) == 100
        ax.plot(numbers(rr, "r"), numbers(rr, "SAS"), color=color, ls=style, lw=1.6,
                label=fr"SAS, $M={m}$")
        if m == 8:
            ax.plot(numbers(rr, "r"), numbers(rr, "Gaussian_reference"), color="black",
                    lw=1.1, ls=":", label="Refined Gaussian reference")
    ax.set(xlabel=r"Radial coordinate $r$", ylabel=r"Steady temperature $T_s$",
           xlim=(0, 1.5), title=r"Gaussian boundary: $x=y=r/\sqrt{2}$, $z=1$")
    ax.legend(frameon=False)
    simple_axes(ax)
    save(fig, "figure2", {"source": "R/case2_profiles.csv", "M": [2, 4, 8],
         "rows_per_M": 100, "fields": ["r", "SAS", "Gaussian_reference"],
         "caption": "Steady Gaussian-boundary ray: saved Chebyshev-superposition approximations at M=2,4,8 and the independently refined true Gaussian reference."})


def figure3():
    rows = [r for i in range(1, 9) for r in read(R / f"case3_grid_{i}.csv")]
    assert len(rows) == 802
    assert len({(int(r["node"]), r["profile"]) for r in rows}) == 802
    fig, ax = plt.subplots(figsize=(6.9, 4.0), layout="constrained")
    for profile, color, style, label in (("U", BLUE, "-", "Uniform (U)"),
                                         ("Q", ORANGE, "--", "Biquadratic (Q)")):
        rr = sorted([r for r in rows if r["profile"] == profile], key=lambda r: int(r["node"]))
        assert len(rr) == 401 and [int(r["node"]) for r in rr] == list(range(401))
        ax.plot(numbers(rr, "ln_z"), numbers(rr, "SAS_t90_over_z2"), color=color,
                ls=style, label=label, lw=1.6)
    ax.set(xlabel=r"$\ln z$", ylabel=r"$t_{90}/z^2$",
           title=r"Centerline response: $L_x=2$, $L_y=1.5$, $\alpha=1$")
    ax.legend(frameon=False)
    simple_axes(ax)
    save(fig, "figure3", {"sources": [f"R/case3_grid_{i}.csv" for i in range(1, 9)],
         "rows": 802, "depths_per_profile": 401, "fields": ["ln_z", "SAS_t90_over_z2"],
         "caption": "Stored response-time ratios for uniform and biquadratic surface temperatures. Both coordinates are native saved scalar fields, without root or ratio recomputation."})


def figure4():
    rows = read(A / "cpu_main100_comparisons.csv")
    source(E / "vendor/performance/replay_case_bank.f90")
    # Literal case membership from build_cases('paper'): profile, geometry,
    # Fo, rho loops. Geometry 2 is (a,b)=(1,.75); retain Q,N1,B1 only.
    groups = (("Q", range(11, 21)), ("N1", range(31, 41)), ("B1", range(71, 81)))
    cells = []
    fig, axes = plt.subplots(1, 3, figsize=(6.9, 2.65), layout="constrained")
    norm = LogNorm(vmin=1/256, vmax=256)
    for ax, (profile, ids) in zip(axes, groups):
        selected = []
        for case_id in ids:
            matches = [r for r in rows if r["case"] == f"paper_{case_id:04d}"
                       and r["family"] == "DIRECT_2D" and r["cost"] == "cpu"
                       and r["mode"] == "ACTUAL_ERROR" and float(r["cap"]) == 1e-6]
            assert len(matches) == 1
            r = matches[0]
            assert r["profile"] == profile and r["status"] == "ACTUAL_ERROR_MATCH"
            assert r["separated"] == "T" and r["resolution"] == "BOTH_ABOVE_RELATIVE_FLOOR"
            assert float(r["SPACE_over_SAS_cost"]) > 0
            selected.append(r)
        grid = []
        for fo in (.75, 1.):
            line = []
            for rho in (.025, .1, .2, .4, .8):
                rr = [r for r in selected if float(r["Fo"]) == fo and float(r["rho"]) == rho]
                assert len(rr) == 1
                r = rr[0]
                line.append(float(r["SPACE_over_SAS_cost"]))
                cells.append(dict(r, a=1, b=.75, x_rho=rho, y_Fo=fo))
            grid.append(line)
        image = ax.imshow(grid, origin="lower", cmap="RdBu_r", norm=norm, aspect="auto")
        for y, line in enumerate(grid):
            for x, value in enumerate(line):
                ax.text(x, y, f"{value:.1f}", ha="center", va="center", fontsize=8,
                        color="white" if value > 32 else "black")
        ax.set(xticks=range(5), xticklabels=[".025", ".1", ".2", ".4", ".8"],
               yticks=[0, 1], yticklabels=[".75", "1"], xlabel=r"$\rho_d$", title=profile)
        ax.set_xticks([.5, 1.5, 2.5, 3.5], minor=True)
        ax.set_yticks([.5], minor=True)
        ax.grid(which="minor", color="white", lw=.55, alpha=.65)
        ax.tick_params(length=0)
        ax.tick_params(which="minor", length=0)
        ax.set_ylabel(r"$\mathrm{Fo}$")
    bar = fig.colorbar(image, ax=axes, shrink=.92, ticks=[.01, .1, 1, 10, 100], pad=.02)
    bar.set_label("2D / SAS CPU query cost")
    assert len(cells) == 30
    save(fig, "figure4", {"source": "A/cpu_main100_comparisons.csv", "selector":
         {"family": "DIRECT_2D", "cost": "cpu", "mode": "ACTUAL_ERROR", "cap": "1e-6", "a": 1, "b": .75},
         "geometry_mapping": "vendor/performance/replay_case_bank.f90 build_cases paper and geometry(2)",
         "colorbar": {"normalization": "log", "vmin": 1/256, "vmax": 256, "center": 1},
         "cells": cells,
         "caption": "CPU query cost ratios at matched observed error for the 30 Q,N1,B1 cells on (a,b)=(1,.75), within the complete 100-case collection. A ratio above one favors SAS."})


def budget_panel(ax, rows, title):
    assert rows
    floors = {float(r["relative_floor"]) for r in rows}
    assert len(floors) == 1
    floor = floors.pop()
    endpoints = [float(r[k]) for r in rows for k in ("start", "end")]
    for direction, style, _ in DIRECTIONS:
        rr = sorted([r for r in rows if r["direction"] == direction], key=lambda r: (float(r["start"]), float(r["end"])))
        assert rr and sum(r["finite_bank_right_endpoint"] == "T" for r in rr) == 1
        for route, color, endpoint_marker in (("SAS", BLUE, "s"), ("SPACE", ORANGE, "o")):
            previous = None
            for r in rr:
                start, end = float(r["start"]), float(r["end"])
                assert end >= start > 0
                if not r[route + "_key"] or not r[route + "_ER"]:
                    previous = None
                    continue
                # Clip only the displayed error at the reporting floor.
                y = max(float(r[route + "_ER"]), floor)
                ax.plot([start, end], [y, y], color=color, ls=style, lw=1.25)
                if previous is not None and previous[0] == start:
                    ax.plot([start, start], [previous[1], y], color=color, ls=style, lw=.8)
                if r["finite_bank_right_endpoint"] == "T":
                    ax.plot(end, y, marker=endpoint_marker, color=color, ms=3.5, clip_on=False)
                if r["zero_width"] == "T" or start == end:
                    ax.plot(start, y, marker="D", color=color, ms=4)
                previous = (end, y)
    ax.axhline(floor, color=".45", lw=.65, ls=(0, (6, 4)), zorder=0)
    ax.set(xscale="log", yscale="log", xlim=(min(endpoints), max(endpoints)),
           ylim=(floor * .6, 2), xlabel=r"CPU query budget $B$ (s)",
           ylabel=r"Relative error $E_R$", title=title)
    ax.grid(which="major", alpha=.18, linewidth=.5)
    ax.text(.98, .02, "Reporting floor", transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color=".35")
    return {"rows": rows, "display_floor": floor, "x_bounds_seconds": [min(endpoints), max(endpoints)],
            "zero_width_rows": [r["_csv_line"] for r in rows if r["zero_width"] == "T" or float(r["start"]) == float(r["end"])],
            "availability": "Each route drawn independently when its saved key and error are present"}


def budget_legend(fig):
    handles = [Line2D([], [], color=BLUE, label="SAS"), Line2D([], [], color=ORANGE, label="2D spatial kernel")]
    handles += [Line2D([], [], color=".25", ls=s, label=label) for _, s, label in DIRECTIONS]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .045), ncol=3, frameon=False, fontsize=8)
    fig.text(.5, .012, "Finite measured lists; blank traces denote unavailable candidates. Symbols mark final endpoints.",
             ha="center", fontsize=8)


def figure5():
    rows = read(A / "cpu_budget3_budgets.csv")
    fig, axes = plt.subplots(1, 2, figsize=(6.9, 4.2), sharey=True)
    fig.subplots_adjust(left=.095, right=.98, top=.88, bottom=.26, wspace=.23)
    panels = {}
    for ax, case, rho in zip(axes, ("paper_0076", "paper_0080"), (.025, .8)):
        rr = [r for r in rows if r["case"] == case and r["family"] == "DIRECT_2D" and r["cost"] == "cpu"]
        assert all(r["profile"] == "B1" and float(r["Fo"]) == 1 and float(r["rho"]) == rho for r in rr)
        panels[case] = budget_panel(ax, rr, fr"B1: $\rho_d={rho}$, $\mathrm{{Fo}}=1$")
    budget_legend(fig)
    save(fig, "figure5", {"source": "A/cpu_budget3_budgets.csv", "panels": panels,
         "fields": ["start", "end", "SAS_ER", "SPACE_ER", "relative_floor", "direction", "finite_bank_right_endpoint", "zero_width"],
         "caption": "Saved affordable representative errors over all finite CPU query-budget intervals for B1 at two depths. Nominal and both conservative timing directions, route-specific availability, reporting-floor segments and finite endpoints are retained. Intermediate unmeasured nodes can alter crossings and gains."})


def figure_s1():
    rows = read(A / "cpu_deep_Q_budgets.csv")
    rr = [r for r in rows if r["case"] == "temperature_scale_0096" and r["family"] == "DIRECT_2D" and r["cost"] == "cpu"]
    assert all(r["profile"] == "Q" and float(r["Fo"]) == 1000 for r in rr)
    fig, ax = plt.subplots(figsize=(6.9, 4.2))
    fig.subplots_adjust(left=.11, right=.98, top=.88, bottom=.26)
    panel = budget_panel(ax, rr, r"Deep Q: $\mathrm{Fo}=1000$, $\rho_d=3.3333$")
    budget_legend(fig)
    save(fig, "figure_s1", {"source": "A/cpu_deep_Q_budgets.csv", "panels": {"temperature_scale_0096": panel},
         "caption": "Deep, long-time Q control over the complete measured finite CPU query-budget range. Both accuracy-advantage directions, both conservative timing directions, route-specific unavailable intervals, reporting-floor values and final finite endpoints are retained."})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    expected = {name + ext for name in ("figure1", "figure2", "figure3", "figure4", "figure5", "figure_s1") for ext in (".png", ".pdf")}
    if OUT.exists():
        assert args.overwrite, "New derived output directory already exists; use --overwrite only for this generated set"
        assert {p.name for p in OUT.iterdir()} <= expected | {"figure_source_map.json"}, "Unexpected files in output directory"
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
        "axes.labelsize": 10, "axes.titlesize": 10, "xtick.labelsize": 8,
        "ytick.labelsize": 8, "legend.fontsize": 8, "axes.spines.top": False,
        "axes.spines.right": False, "pdf.fonttype": 42, "savefig.dpi": 300})
    for render in (figure1, figure2, figure3, figure4, figure5, figure_s1):
        render()
    assert all(sha(P / path) == digest for path, digest in TRACE["inputs"].items()), "An input changed during rendering"
    TRACE["renderer"] = {"path": str(Path(__file__).relative_to(P)), "sha256": sha(Path(__file__)), "matplotlib": matplotlib.__version__}
    MAP.write_text(json.dumps(TRACE, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"figures": len(TRACE["figures"]), "output_directory": str(OUT), "source_map": str(MAP)}))


if __name__ == "__main__":
    main()
