#!/usr/bin/env python3
"""Result figures for docs/paper/IHM-1.md, read only from committed result files in out/.

  stance_kernel_controls.png   recovery per policy x push set; peak COM per condition
  skin_stance_support.png      pelvis drift per contact arm at 1 s and 2 s (if present)

Writes PNG + JPG into docs/figures/results/.
"""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/figures/results"
INK, SOFT, GRID = "#222226", "#6b6b72", "#e4e2dc"
LEARNED, CONTROL = "#b8452f", "#6f8fb3"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": SOFT, "axes.labelcolor": INK, "xtick.color": SOFT,
                     "ytick.color": SOFT, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "white", "savefig.dpi": 200})


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.jpg", bbox_inches="tight", pil_kwargs={"quality": 88})
    plt.close(fig)
    print(OUT / f"{name}.png")


def stance():
    o = json.loads((ROOT / "out/stance_kernel_controls.json").read_text())
    R, learned, controls = o["rows"], o["learned"], o["controls"]
    order = learned + controls
    sets = "ABCD"
    sizes = {s: sum(1 for r in R[order[0]].values() if r["set"] == s) for s in sets}
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.2), gridspec_kw={"width_ratios": [1, 1.5]})
    bottom = np.zeros(len(order))
    shades = {"A": 1.0, "B": 0.78, "C": 0.56, "D": 0.34}
    for s in sets:
        n = np.array([sum(1 for r in R[p].values() if r["set"] == s and r["full"]["completed"]) for p in order])
        cols = [LEARNED if p in learned else CONTROL for p in order]
        a.bar(range(len(order)), n, bottom=bottom, color=cols, alpha=shades[s], edgecolor="white", linewidth=0.6,
              label=f"set {s} ({sizes[s]})")
        bottom += n
    for i, t in enumerate(bottom):
        a.text(i, t + 0.3, f"{int(t)}", ha="center", color=INK, fontsize=9)
    a.set_xticks(range(len(order)), order, rotation=40, ha="right")
    a.set_ylabel(f"conditions recovered (of {o['recovered_of']})")
    a.set_ylim(0, o["recovered_of"])
    a.set_title("Recovery: learned (red) vs control kernels (blue)", loc="left", color=INK, fontsize=11)
    a.legend(frameon=False, fontsize=8, title="push set (n)", title_fontsize=8, loc="upper right")
    conds = sorted(R[order[0]], key=lambda c: ("ABCD".index(R[order[0]][c]["set"]), c))
    x = np.arange(len(conds))
    for p in order:
        y = [R[p][c]["full"]["peak_m"] * 1000 if R[p][c]["full"]["completed"] else np.nan for c in conds]
        b.plot(x, y, "o", ms=4, color=LEARNED if p in learned else CONTROL, alpha=0.9 if p in learned else 0.55)
    b.set_yscale("log")
    nice = lambda c: c.split("_")[0] + " (" + ", ".join(f"{float(v):+.0f}" for v in c.split("_")[1:]) + ")"
    b.set_xticks(x, [nice(c) for c in conds], rotation=70, fontsize=7)
    b.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    b.yaxis.set_minor_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}" if v in (1.5, 2, 3, 5) else ""))
    b.set_ylabel("peak pelvis COM excursion, mm (completed runs only)")
    b.set_title(f"Peak excursion per push (x, z) N. Pre-registered Δ ln(peak), {o['complete_conditions']} all-complete conditions:\n"
                f"{o['mean_delta_ln_peak']:+.3f} [{o['ci95'][0]:+.3f}, {o['ci95'][1]:+.3f}] — {o['verdict']}",
                loc="left", color=INK, fontsize=9)
    b.grid(axis="y", color=GRID, lw=0.6)
    fig.text(0, -0.16, "Missing points are falls. Native OpenSim scaffold plant, 77.6 kg; 100 ms pelvis push at t = 1 s, "
             "10 s horizon. docs/STANCE_KERNEL_CONTROLS.md.", color=SOFT, fontsize=8)
    save(fig, "stance_kernel_controls")


def skin():
    files = [(s, ROOT / f"out/skin_stance_support_{s}.json") for s in (100, 200)]
    files = [(s, f) for s, f in files if f.exists()]
    if not files:
        return
    rows = {}
    for steps, f in files:
        for r in json.loads(f.read_text())["arms"]:
            if "failed" in r:
                rows.setdefault(r["arm"], {})[steps] = (np.nan, True)
            else:
                drift = (r["final"]["pelvis_ty_m"] - r["initial"]["pelvis_ty_m"]) * 1000
                rows.setdefault(r["arm"], {})[steps] = (drift, False)
    arms = list(rows)
    fig, ax = plt.subplots(figsize=(8, 3.6))
    for j, steps in enumerate(s for s, _ in files):
        vals = [rows[a].get(steps, (np.nan, False)) for a in arms]
        y = np.arange(len(arms)) + j * 0.38
        ax.barh(y, [max(v, -900) for v, _ in vals], height=0.36, color=["#c9c4b8", "#8a8478"][j],
                label=f"{steps / 100:.0f} s")
        for yy, (v, failed) in zip(y, vals):
            if failed:
                ax.text(-130, yy, "integrator failed at 1.08 s", va="center", ha="right", fontsize=8, color=LEARNED)
            else:
                ax.text(min(v, 0) - 8, yy, f"{v:.1f}", va="center", ha="right", fontsize=7, color=INK)
    ax.axvline(-20, color=SOFT, lw=0.8, ls="--"); ax.axvline(-100, color=LEARNED, lw=0.8, ls="--")
    ax.set_yticks(np.arange(len(arms)) + 0.19, arms)
    ax.set_xlabel("pelvis drift, mm (stand > −20; fall < −100); passive plant, joint stops, no controller")
    ax.set_xlim(-1000, 40)
    ax.legend(frameon=False)
    ax.set_title("Where the support is: skin contact with foot-sphere subsets", loc="left", color=INK, fontsize=11)
    save(fig, "skin_stance_support")


if __name__ == "__main__":
    stance()
    skin()
