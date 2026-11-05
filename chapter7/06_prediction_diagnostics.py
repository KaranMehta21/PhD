"""
PhD C7 - DIAGNOSTIC for the refit: why does the CPH with splines + interactions have a worse IBS than the original CPH?
Uses the saved test-set predictions (ch7_refit\\pred_*.npz) and models; re-creates the identical split (same seed).
Reports per model: distribution of predicted 5- and 9.99-year survival, calibration-in-the-large (mean predicted vs Kaplan-Meier
observed), number of extreme predictions, Brier score over time, and IBS restricted to people without extreme predictions;
plus the largest coefficients in the fitted Cox models.  ~5-10 min (most of it reading the outcome columns).
Run AFTER phd_c7_refit.py has finished:  python "OneDrive - University College London\\phd_c7_diagnose.py"
Output: Desktop\\pure_inventory\\c7_diagnose.txt
"""
import os, gc, pickle, warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")
K = r"S:\FPHS_THIN_Projects\Karan"
DATA = os.path.join(K, "BZI_PhD_T2D", "PhD_descriptive_data", "mice_data", "impv2_phddata_230224.dta")
OUTD = os.path.join(K, "PURE_v2", "ch7_refit")
OUT = os.path.join(os.path.expanduser("~"), "Desktop", "pure_inventory", "c7_diagnose.txt")
SEED, TEST_SHARE = 2026, 0.30
GRID = np.linspace(0.05, 9.99, 100)
L = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    L.append(s)


parts = [p.copy() for p in pd.read_stata(DATA, columns=["exposure", "mortality_status"], convert_categoricals=False, chunksize=20000)]
D = pd.concat(parts, ignore_index=True)
T, E = D["exposure"].values.astype(float), D["mortality_status"].values.astype(bool)
rng = np.random.default_rng(SEED)
te_idx = np.sort(rng.choice(len(D), int(round(TEST_SHARE * len(D))), replace=False))
te = np.zeros(len(D), bool); te[te_idx] = True
Tt, Et, Ttr, Etr = T[te], E[te], T[~te], E[~te]


def km(t, e):
    o = np.argsort(t); t, e = t[o], e[o].astype(float)
    u, idx = np.unique(t, return_index=True)
    d = np.add.reduceat(e, idx)
    return u, np.cumprod(1 - d / (len(t) - idx))


Gu, Gs = km(Ttr, ~Etr)
G = lambda x: np.maximum(np.where(np.searchsorted(Gu, x, "right") - 1 >= 0, Gs[np.clip(np.searchsorted(Gu, x, "right") - 1, 0, None)], 1.0), 0.05)
ku, ks = km(Tt, Et)
obs = lambda t: float(ks[np.searchsorted(ku, t, "right") - 1])


def brier_curve(S, t=None, e=None):
    t = Tt if t is None else t
    e = Et if e is None else e
    Gt, Gg = G(t), G(GRID)
    died = (t[:, None] <= GRID[None, :]) & e[:, None]
    alive = t[:, None] > GRID[None, :]
    return (np.where(died, S ** 2 / Gt[:, None], 0) + np.where(alive, (1 - S) ** 2 / Gg[None, :], 0)).mean(axis=0)


trap = getattr(np, "trapezoid", None) or np.trapz
log(f"test set {te.sum():,} people, {int(Et.sum()):,} deaths; observed KM survival at 5 y {obs(5):.4f}, at 9.99 y {obs(9.99):.4f}")
for nm in ("CPH_original", "CPH_nonlinear_interactions", "DeepSurv", "RSF"):
    p = os.path.join(OUTD, f"pred_{nm}.npz")
    if not os.path.exists(p):
        continue
    z = np.load(p); S = z["S"].astype(float)
    j5, j9 = int(np.argmin(abs(GRID - 5))), len(GRID) - 1
    s5, s9 = S[:, j5], S[:, j9]
    bc = brier_curve(S)
    extreme = (s9 < 0.01) | ~np.isfinite(s9)
    ib_all = trap(bc, GRID) / (GRID[-1] - GRID[0])
    bc_ok = brier_curve(S[~extreme], Tt[~extreme], Et[~extreme]) if extreme.any() else bc
    log(f"\n{nm}")
    log(f"  predicted S(5y): mean {s5.mean():.4f} vs observed {obs(5):.4f} | percentiles 1/50/99: {np.percentile(s5, [1, 50, 99]).round(4).tolist()}")
    log(f"  predicted S(9.99y): mean {s9.mean():.4f} vs observed {obs(9.99):.4f} | min {np.nanmin(s9):.2e}")
    log(f"  people with predicted 10-y survival < 1%: {int(extreme.sum()):,} ({100 * extreme.mean():.2f}%) - of whom died: {int(Et[extreme].sum()):,}")
    if extreme.sum() > 20:
        eu, es = km(Tt[extreme], Et[extreme])
        log(f"    their OBSERVED survival: 1 y {float(es[np.searchsorted(eu, 1, 'right') - 1]) if (eu <= 1).any() else 1:.3f}, "
            f"5 y {float(es[np.searchsorted(eu, 5, 'right') - 1]) if (eu <= 5).any() else 1:.3f} (vs predicted < 0.01 at 10 y)")
    log(f"  Brier score at 1/3/5/7/9.99 y: " + ", ".join(f"{bc[int(np.argmin(abs(GRID - t)))]:.4f}" for t in (1, 3, 5, 7, 9.99)))
    log(f"  IBS all {ib_all:.4f} | IBS excluding the extreme predictions {trap(bc_ok, GRID) / (GRID[-1] - GRID[0]):.4f}")
for nm in ("cph_original.pkl", "cph_nonlinear_train.pkl"):
    p = os.path.join(OUTD, "models", nm)
    if os.path.exists(p):
        b = pickle.load(open(p, "rb"))
        co = b["model"].params_.abs().sort_values(ascending=False)
        log(f"\n{nm}: largest |coefficients| (per SD): " + "; ".join(f"{k} {b['model'].params_[k]:+.2f}" for k in co.index[:10]))
        log(f"  converged: {getattr(b['model'], '_convergence_success', 'n/a')}; log-likelihood {b['model'].log_likelihood_:.1f}")
os.makedirs(os.path.dirname(OUT), exist_ok=True)
open(OUT, "w", encoding="utf-8").write("\n".join(L))
