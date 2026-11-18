"""
PhD Ch7 - EXACT SPLINE FORMULAS FOR THE PREDICTION TOOL. ~5-20 min. Run on the UCL PC (needs the thesis data).
Problem: the refit saved each mvrs spline as a lookup table (raw value -> spline terms, averaged over observed values) and
the tool interpolates between table points. For most variables that is near-exact; for diabetes duration (and SBP, albumin,
ALT) it is noticeably off.
Fix: mvrs splines are cubic regression splines, so every spline term is EXACTLY a linear combination of
    1, x, x^2, x^3, (x-k1)^3+, (x-k2)^3+, ...   for the knots k used by mvrs.
For each spline variable this script finds the knots (from the Stata labels, then mvrs's default centile knots, then a search),
fits that combination on people with an OBSERVED raw value, and keeps it only if it reproduces the stored terms
(max error < 1e-6 of the term's range). It also reports whether the term is a function of the raw value at all
(if not, the raw column is not the one mvrs used, and the script says so).
Output (aggregate only: knots + coefficients, no patient data):
    Desktop\\phd_model_export\\spline_maps_exact.pkl   (the old maps + exact formulas where found)
    Desktop\\phd_model_export\\spline_fix_report.txt
Copy both to Google Drive PhD_final_submission/outputs/ch7/models/phd_model_export/ . phd_c7_predict.py uses the exact
formulas automatically when spline_maps_exact.pkl is present.
Run:  python "OneDrive - University College London\\phd_c7_fix_splines.py"
"""
import os, re, gc, pickle, itertools, datetime
import numpy as np
import pandas as pd

K = r"S:\FPHS_THIN_Projects\Karan"
DATA = os.path.join(K, "BZI_PhD_T2D", "PhD_descriptive_data", "mice_data", "impv2_phddata_230224.dta")
MODD = os.path.join(K, "PURE_v2", "ch7_refit", "models")
OUT = os.path.join(os.path.expanduser("~"), "Desktop", "phd_model_export")
os.makedirs(OUT, exist_ok=True)
TOL = 1e-5          # relative to each term's range (Stata stores terms as float, ~1e-7 precision)
lines = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    lines.append(s)


def read_cols(cols, chunk=20000):
    parts = []
    for part in pd.read_stata(DATA, columns=list(cols), convert_categoricals=False, chunksize=chunk):
        parts.append(part.copy())
        gc.collect()
    return pd.concat(parts, ignore_index=True)


def basis(x, knots):
    x = np.asarray(x, float)
    cols = [np.ones_like(x), x, x ** 2, x ** 3] + [np.clip(x - k, 0, None) ** 3 for k in knots]
    return np.column_stack(cols)


def fit(x, Y, knots):
    """Least squares of every term on the cubic basis; returns coefficients and max error relative to each term's range."""
    B = basis(x, knots)
    scale = np.abs(B).max(axis=0)
    scale[scale == 0] = 1
    C, *_ = np.linalg.lstsq(B / scale, Y, rcond=None)
    C = C / scale[:, None]
    rng = np.ptp(Y, axis=0)
    rng[rng == 0] = 1
    err = float(np.max(np.abs(B @ C - Y) / rng))
    return C, err


def local_knots(x, Y, frac=0.08):
    """Windows of fixed x-span; a cubic fits exactly unless a knot lies inside. Returns knot brackets."""
    yr = np.ptp(Y, axis=0); yr[yr == 0] = 1
    lo, hi = x.min(), x.max(); span = frac * (hi - lo); stp = span / 4
    res = []
    a = lo
    while a < hi:
        m = (x >= a) & (x <= a + span)
        if m.sum() >= 8:
            xx = x[m]; YY = Y[m] / yr; c = (xx - xx.mean()) / (np.ptp(xx) or 1)
            B = np.column_stack([np.ones_like(c), c, c ** 2, c ** 3]); C, *_ = np.linalg.lstsq(B, YY, rcond=None)
            res.append((a, float(np.max(np.abs(B @ C - YY)))))
        a += stp
    r = np.array([v for _, v in res]); thr = max(50 * np.median(r), 1e-6)
    groups = []
    for a_, v in res:
        if v > thr:
            if groups and a_ <= groups[-1][1] + stp + 1e-12: groups[-1][1] = a_
            else: groups.append([a_, a_])
    # bracket = intersection-ish: from last start to first end of flagged windows
    return [(max(lo, g[1]), min(hi, g[0] + span)) for g in groups]


def refine(x, Y, kn, bounds, sweeps=40):
    """Move each knot in turn to minimise the RMS error, within its bracket."""
    from scipy.optimize import minimize_scalar
    yr = np.ptp(Y, axis=0)
    yr[yr == 0] = 1

    def rms(kk):
        C_, _ = fit(x, Y, sorted(kk))
        return float(np.sqrt(np.mean(((basis(x, sorted(kk)) @ C_ - Y) / yr) ** 2)))
    kn = list(kn)
    for _ in range(sweeps):
        before = rms(kn)
        for j, (a, z) in enumerate(bounds):
            r = minimize_scalar(lambda v: rms(kn[:j] + [v] + kn[j + 1:]), bounds=(a, z), method="bounded",
                                options={"xatol": 1e-11 * (np.ptp(x) + 1)})
            kn[j] = float(r.x)
        if before - rms(kn) < 1e-14:
            break
    return sorted(kn)


def label_knots(label):
    nums = [float(v) for v in re.findall(r"-?\d+(?:\.\d+)?", label or "")]
    return nums


log("=" * 90, f"\nPhD Ch7 spline fix   run {datetime.datetime.now():%Y-%m-%d %H:%M}\n" + "=" * 90)
SM = pickle.load(open(os.path.join(MODD, "spline_maps.pkl"), "rb"))
MAPS = SM["maps"]
with pd.io.stata.StataReader(DATA) as rdr:
    labels = rdr.variable_labels()
allcols = set(labels)

new_maps = {}
for b, m in sorted(MAPS.items()):
    terms = m["terms"]
    orig_terms = [t.replace("_1_", "", 1) if t.replace("_1_", "", 1) in allcols else t for t in terms]
    raw_candidates = [c for c in (b, f"_1_{b}") if c in allcols]
    log(f"\n{b}  ({len(terms)} terms: {', '.join(orig_terms)})")
    R = read_cols(sorted(set(raw_candidates + orig_terms)))
    best = None
    for raw in raw_candidates:
        sub = R[[raw] + orig_terms].dropna().astype(float)
        if len(sub) < 100:
            continue
        # is each term a function of the raw value at all?
        g = sub.groupby(sub[raw].round(9))[orig_terms]
        spread = float((g.max() - g.min()).max().max() / max(np.ptp(sub[orig_terms].values, axis=0).max(), 1e-12))
        u = g.mean().sort_index()
        x, Y = u.index.values.astype(float), u.values.astype(float)
        old_err = max(float(np.max(np.abs(np.interp(sub[raw].values, m["x"], m["y"][:, k]) - sub[orig_terms[k]].values)))
                      / max(np.ptp(sub[orig_terms[k]].values), 1e-12) for k in range(len(terms)))
        log(f"   raw column {raw}: {len(sub):,} observed rows, {len(x):,} distinct values; "
            f"same raw value -> different term values: {spread:.1e} (0 = a true function); old lookup-table error {old_err:.1e}")
        cands = []
        lk = sorted({v for t in orig_terms for v in label_knots(labels.get(t, "")) if x.min() < v < x.max()})
        if lk:
            cands.append(("Stata label", lk))
        xs_all = sub[raw].values
        for nk in (1, 2, 3, 4):
            q = np.linspace(0, 100, nk + 2)[1:-1]
            cands.append((f"{nk} knots at equal centiles (observed)", list(np.percentile(xs_all, q))))
            imp = f"_1_{b}"
            if imp in R.columns and imp != raw:
                cands.append((f"{nk} knots at equal centiles (imputed)", list(np.percentile(R[imp].dropna().values, q))))
        for name, kn in cands:
            C, err = fit(x, Y, kn)
            if best is None or err < best["err"]:
                best = {"err": err, "knots": kn, "coef": C, "how": name, "raw": raw, "old_err": old_err, "spread": spread}
        if best["err"] >= TOL:                                   # detect knots from where a single cubic stops fitting
            for frac in (0.08, 0.05, 0.12):
                br = local_knots(x, Y, frac)
                if 1 <= len(br) <= 4:
                    kn = refine(x, Y, [(a + z) / 2 for a, z in br], br)
                    C, err = fit(x, Y, kn)
                    if err < best["err"]:
                        best = {"err": err, "knots": kn, "coef": C, "how": f"detected ({len(kn)} knots)", "raw": raw,
                                "old_err": old_err, "spread": spread}
                    if err < TOL:
                        break
        if best["err"] >= TOL:
            # search: 1-3 knots over a centile grid (exhaustive), then refine each knot in turn (1-D searches)
            from scipy.optimize import minimize_scalar
            yr = np.ptp(Y, axis=0)
            yr[yr == 0] = 1

            def rms(kk):
                C_, _ = fit(x, Y, sorted(kk))
                return float(np.sqrt(np.mean(((basis(x, sorted(kk)) @ C_ - Y) / yr) ** 2)))
            lo, hi = float(x.min()), float(x.max())
            for nk, step in ((1, 1), (2, 1), (3, 2)):
                grid = np.unique(np.percentile(xs_all, np.arange(2, 99, step)))
                combos = sorted(itertools.combinations(grid, nk), key=rms)[:25]   # refine the 25 best starts
                found = []
                for start in combos:
                    kn = list(start)
                    for _sweep in range(30):
                        before = rms(kn)
                        for j in range(nk):
                            a = kn[j - 1] if j > 0 else lo
                            z = kn[j + 1] if j < nk - 1 else hi
                            r = minimize_scalar(lambda v: rms(kn[:j] + [v] + kn[j + 1:]), bounds=(a, z), method="bounded",
                                                options={"xatol": 1e-10 * (hi - lo + 1)})
                            kn[j] = float(r.x)
                        if before - rms(kn) < 1e-13:
                            break
                    found.append((fit(x, Y, kn)[1], kn))
                    if found[-1][0] < TOL:
                        break
                kn = min(found)[1]
                C, err = fit(x, Y, kn)
                if err < best["err"]:
                    best = {"err": err, "knots": kn, "coef": C, "how": f"searched ({nk} knots)", "raw": raw,
                            "old_err": old_err, "spread": spread}
                if err < TOL:
                    break
    del R
    gc.collect()
    nm = dict(m)
    if best and best["err"] < TOL:
        nm.update({"knots": [float(k) for k in best["knots"]], "coef": best["coef"], "exact": True})
        log(f"   EXACT: {best['how']}; knots {np.round(best['knots'], 4).tolist()}; max error {best['err']:.1e} "
            f"(old table {best['old_err']:.1e})")
    elif best:
        nm["exact"] = False
        why = ("the raw column is not what mvrs used (same raw value gives different terms)" if best["spread"] > 1e-5
               else "no cubic-spline knot set reproduces it")
        log(f"   NOT EXACT (best {best['err']:.1e} via {best['how']}): {why}. Tool keeps the lookup table for {b}.")
    else:
        nm["exact"] = False
        log("   raw column not found; lookup table kept")
    new_maps[b] = nm

out = dict(SM)
out["maps"] = new_maps
out["note"] = "exact cubic-spline formulas (knots + coefficients) where found; aggregate only"
pickle.dump(out, open(os.path.join(OUT, "spline_maps_exact.pkl"), "wb"))
n_ex = sum(v.get("exact", False) for v in new_maps.values())
log(f"\n{n_ex} of {len(new_maps)} spline variables now exact. Saved {os.path.join(OUT, 'spline_maps_exact.pkl')}")
open(os.path.join(OUT, "spline_fix_report.txt"), "w", encoding="utf-8").write("\n".join(lines))
print("Saved", os.path.join(OUT, "spline_fix_report.txt"))
