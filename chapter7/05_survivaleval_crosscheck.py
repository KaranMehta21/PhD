"""
PhD C7 - CROSS-CHECK THE REFIT'S PERFORMANCE MEASURES WITH SurvivalEVAL (the package used in the thesis), using only the saved
test-set predictions (ch7_refit\\pred_*.npz) - no models are loaded, so memory use is low. Same seeded split as the refit.
For each model, side by side:  refit (lifelines Harrell's C; IPCW IBS over 0-9.99 y)  vs  SurvivalEVAL (Harrell's C; IBS default
and over 0-9.99 y; 1-calibration at 5 years; D-calibration).
Output: Desktop\\pure_inventory\\c7_survivaleval_check.txt. ~10 min.
Run:  python "OneDrive - University College London\\phd_c7_survivaleval_check.py"
"""
import os, sys, gc, warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")

K = r"S:\FPHS_THIN_Projects\Karan"
DATA = os.path.join(K, "BZI_PhD_T2D", "PhD_descriptive_data", "mice_data", "impv2_phddata_230224.dta")
OUTD = os.path.join(K, "PURE_v2", "ch7_refit")
OD = os.path.join(os.path.expanduser("~"), "OneDrive - University College London")
OUT = os.path.join(os.path.expanduser("~"), "Desktop", "pure_inventory", "c7_survivaleval_check.txt")
SEED, TEST_SHARE = 2026, 0.30
GRID = np.linspace(0.05, 9.99, 100)
L = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    L.append(s)


# SurvivalEVAL: pip-installed package or the local copy used by the thesis notebooks
SE = None
for p in (os.path.join(OD, "SurvivalEVAL"), os.path.join(OD, "SurvivalEVAL-main"), OD):
    if os.path.isdir(p) and p not in sys.path:
        sys.path.append(p)
for imp in ("from SurvivalEVAL import SurvivalEvaluator as SE", "from SurvivalEVAL.Evaluator import SurvivalEvaluator as SE",
            "from Evaluator import SurvivalEvaluator as SE"):
    try:
        exec(imp)
        log(f"SurvivalEVAL imported via: {imp}")
        break
    except Exception:
        continue
if SE is None:
    sys.exit("SurvivalEVAL could not be imported (pip install SurvivalEVAL, or check the local copy in OneDrive)")

parts = [p.copy() for p in pd.read_stata(DATA, columns=["exposure", "mortality_status"], convert_categoricals=False, chunksize=20000)]
D = pd.concat(parts, ignore_index=True)
rng = np.random.default_rng(SEED)
te_idx = np.sort(rng.choice(len(D), int(round(TEST_SHARE * len(D))), replace=False))
te = np.zeros(len(D), bool)
te[te_idx] = True
T, E = D["exposure"].values.astype(float), D["mortality_status"].values.astype(int)
Tt, Et, Ttr, Etr = T[te], E[te], T[~te], E[~te]
del D, parts
log(f"test set {te.sum():,} ({Et.sum():,} deaths); train {(~te).sum():,}")

from lifelines.utils import concordance_index


def km(t, e):
    o = np.argsort(t)
    t, e = t[o], e[o].astype(float)
    u, idx = np.unique(t, return_index=True)
    d = np.add.reduceat(e, idx)
    return u, np.cumprod(1 - d / (len(t) - idx))


Gu, Gs = km(Ttr, 1 - Etr)
G = lambda x: np.maximum(np.where(np.searchsorted(Gu, x, "right") - 1 >= 0, Gs[np.clip(np.searchsorted(Gu, x, "right") - 1, 0, None)], 1.0), 0.05)


def own_ibs(S):
    Gt, Gg = G(Tt), G(GRID)
    died = (Tt[:, None] <= GRID[None, :]) & (Et[:, None] == 1)
    alive = Tt[:, None] > GRID[None, :]
    bs = (np.where(died, S ** 2 / Gt[:, None], 0) + np.where(alive, (1 - S) ** 2 / Gg[None, :], 0)).mean(axis=0)
    trap = getattr(np, "trapezoid", None) or np.trapz
    return float(trap(bs, GRID) / (GRID[-1] - GRID[0]))


def call(f, *variants):
    for kw in variants:
        try:
            return f(**kw)
        except TypeError:
            continue
    return f()


rows = []
coords = np.concatenate([[0.0], GRID])
for name in ("CPH_original", "CPH_nonlinear_interactions", "RSF", "DeepSurv"):
    p = os.path.join(OUTD, f"pred_{name}.npz")
    if not os.path.exists(p):
        log(f"{name}: predictions not found")
        continue
    z = np.load(p)
    S = z["S"].astype(float)
    risk = z["risk"].astype(float)
    r = {"model": name, "own_C": concordance_index(Tt, -risk, Et), "own_IBS": own_ibs(S)}
    ev = SE(np.hstack([np.ones((len(S), 1)), S]), coords, Tt, Et, Ttr, Etr)
    try:
        c = call(ev.concordance, {"ties": "None"}, {})
        r["SE_C"] = float(c[0] if isinstance(c, (tuple, list)) else c)
    except Exception as ex:
        r["SE_C"] = np.nan; log(f"  {name} concordance failed: {ex}")
    try:
        ib = call(ev.integrated_brier_score, {"num_points": None, "draw_figure": False}, {"draw_figure": False}, {})
        r["SE_IBS_default"] = float(ib[0] if isinstance(ib, tuple) else ib)
    except Exception as ex:
        r["SE_IBS_default"] = np.nan; log(f"  {name} IBS (default) failed: {ex}")
    try:
        ib = call(ev.integrated_brier_score, {"target_times": GRID, "draw_figure": False})
        r["SE_IBS_0_10y"] = float(ib[0] if isinstance(ib, tuple) else ib)
    except Exception:
        r["SE_IBS_0_10y"] = np.nan
    try:
        oc = call(ev.one_calibration, {"target_time": 5.0}, {"target_time": 5})
        r["SE_1cal_p_5y"] = float(oc[0] if isinstance(oc, (tuple, list)) else oc)
    except Exception as ex:
        r["SE_1cal_p_5y"] = np.nan; log(f"  {name} 1-calibration failed: {ex}")
    try:
        dc = ev.d_calibration()
        r["SE_Dcal_p"] = float(dc[0] if isinstance(dc, (tuple, list)) else dc)
    except Exception:
        r["SE_Dcal_p"] = np.nan
    rows.append(r)
    log(f"{name:<28} C: refit {r['own_C']:.4f} | SurvivalEVAL {r['SE_C']:.4f}    IBS: refit {r['own_IBS']:.4f} | SurvivalEVAL default "
        f"{r['SE_IBS_default']:.4f}, 0-10 y {r['SE_IBS_0_10y']:.4f}    1-cal p (5 y) {r['SE_1cal_p_5y']:.3g}    D-cal p {r['SE_Dcal_p']:.3g}")
    del ev, S, z
    gc.collect()
pd.DataFrame(rows).to_csv(os.path.join(OUTD, "c7_survivaleval_check.csv"), index=False)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
open(OUT, "w", encoding="utf-8").write("\n".join(L))
print(f"\nSaved {OUT}")
