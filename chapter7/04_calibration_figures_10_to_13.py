"""
PhD C7 - REDRAW THE CALIBRATION FIGURES (Figures 10-13) from the refit's saved test-set predictions (no refitting).
Deciles of predicted 5-year risk ordered from lowest (1) to highest (10), so risk increases left to right and the legend sits
clear of the bars. Same Hosmer-Lemeshow calculation as phd_c7_refit.py, so the p-values are unchanged.
Output: S:\\...\\PURE_v2\\ch7_refit\\figures\\Figure10..13_calibration_*.png (300 dpi). ~5 min (reading the outcome columns).
Run:  python "OneDrive - University College London\\phd_c7_calibration_plots.py"
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2

K = r"S:\FPHS_THIN_Projects\Karan"
DATA = os.path.join(K, "BZI_PhD_T2D", "PhD_descriptive_data", "mice_data", "impv2_phddata_230224.dta")
OUTD = os.path.join(K, "PURE_v2", "ch7_refit")
FIGD = os.path.join(OUTD, "figures")
os.makedirs(FIGD, exist_ok=True)
SEED, TEST_SHARE, T_STAR = 2026, 0.30, 5.0
GRID = np.linspace(0.05, 9.99, 100)
MODELS = [("CPH_original", 10, "Cox (original variables)"), ("CPH_nonlinear_interactions", 11, "Cox (non-linear and interaction terms)"),
          ("RSF", 12, "Random survival forest"), ("DeepSurv", 13, "DeepSurv")]

parts = [p.copy() for p in pd.read_stata(DATA, columns=["exposure", "mortality_status"], convert_categoricals=False, chunksize=20000)]
D = pd.concat(parts, ignore_index=True)
rng = np.random.default_rng(SEED)
te_idx = np.sort(rng.choice(len(D), int(round(TEST_SHARE * len(D))), replace=False))
T = D["exposure"].values[te_idx].astype(float)
E = D["mortality_status"].values[te_idx].astype(bool)


def km_at(t, e, t_star):
    o = np.argsort(t)
    t, e = t[o], e[o].astype(float)
    u, idx = np.unique(t, return_index=True)
    d = np.add.reduceat(e, idx)
    s = np.cumprod(1 - d / (len(t) - idx))
    k = np.searchsorted(u, t_star, side="right") - 1
    return float(s[k]) if k >= 0 else 1.0


j = int(np.argmin(np.abs(GRID - T_STAR)))
for name, fig_no, label in MODELS:
    p = os.path.join(OUTD, f"pred_{name}.npz")
    if not os.path.exists(p):
        print(f"missing {p}")
        continue
    S_at = np.load(p)["S"][:, j].astype(float)
    q = pd.qcut(pd.Series(-S_at).rank(method="first"), 10, labels=False).values      # decile 0 = lowest predicted risk
    obs, exp_, stat = [], [], 0.0
    for b in range(10):
        m = q == b
        o = km_at(T[m], E[m], GRID[j])
        pr = S_at[m].mean()
        obs.append(1 - o)
        exp_.append(1 - pr)
        stat += m.sum() * (o - pr) ** 2 / max(pr * (1 - pr), 1e-9)
    pval = 1 - chi2.cdf(stat, 8)
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    x = np.arange(1, 11)
    ax.bar(x - 0.18, obs, 0.36, label="Observed", color="#1f77b4")
    ax.bar(x + 0.18, exp_, 0.36, label="Predicted", color="#ff7f0e")
    ax.set_xticks(x)
    ax.set_xlabel("Decile of predicted risk (1 = lowest, 10 = highest)")
    ax.set_ylabel(f"Probability of death by {GRID[j]:.0f} years")
    ax.legend(loc="upper left", frameon=False)
    ax.text(0.02, 0.80, f"Hosmer–Lemeshow p = {pval:.3f}" if pval >= 0.001 else "Hosmer–Lemeshow p < 0.001", transform=ax.transAxes)
    ax.set_title(label, fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out = os.path.join(FIGD, f"Figure{fig_no}_calibration_{name}.png")
    fig.savefig(out)
    plt.close(fig)
    print(f"saved {out}  (p = {pval:.3f})")
