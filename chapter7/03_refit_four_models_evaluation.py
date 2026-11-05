"""
PhD CHAPTER 7 REFIT (C7.2, C7.4, C7.5) - all four Table 8 models on ONE seeded random 70/30 split.  Several hours (RSF longest).
Data      impv2_phddata_230224.dta, imputation 1 (as in the thesis notebooks).
Predictors read DIRECTLY from the thesis notebooks, so they are identical to the originals:
            CPH original (linear) = X_train columns of PhD_cph_scikitsurv_240324.ipynb (also used for RSF and DeepSurv)
            CPH non-linear + interactions = 'dataset' columns of PhD_CPHlifelines_nonandinter_230424.ipynb (mvrs splines + 13 interactions)
Split     seeded random 70% train / 30% test (SEED below); identical test set for every model. DeepSurv uses 20% of the training set
          (seeded) for early stopping. Continuous predictors standardised with TRAINING-set means/SDs.
Models    lifelines CoxPHFitter (both Cox models) | scikit-survival RandomSurvivalForest(n_estimators=32, min_samples_split=6,
          min_samples_leaf=3, max_depth=6, random_state=20) | pycox DeepSurv (MLPVanilla [32,32], batch norm, dropout 0.1, Adam lr 0.01,
          batch 256, up to 512 epochs, early stopping) - the thesis hyperparameters.
Outputs   (S:\\...\\PURE_v2\\ch7_refit and Desktop\\pure_inventory)
  C7.5    Table 8 corrected: Harrell's C-index and integrated Brier score (IBS, 0-9.99 y, IPCW) on the test set, with 95% bootstrap CIs
          (percentile, B resamples of the test set); SurvivalEVAL values printed as a cross-check when the package is importable.
  C7.4    final CPH (non-linear + interactions) refitted on ALL data: coefficients, HR (95% CI), p, standardisation mean/SD for every
          term, and baseline survival S0(t) at 1, 5 and 9 years -> c7_cox_final_coefficients.csv, c7_baseline_survival.csv
  C7.2    proportional hazards test (Schoenfeld residuals, rank time transform; on a seeded 30,000-person training subsample),
          and the spline variable labels from the Stata file (mvrs usually stores the knot positions there)
  Figures 1-calibration by decile of predicted survival at t* (observed KM vs mean predicted, Hosmer-Lemeshow p) for each model
Prediction bundle (ch7_refit\models): all four fitted models (+ the final full-data CPH), standardisation parameters, and the
          spline maps (raw value -> mvrs spline terms, learned from observed values and verified), for phd_c7_predict.py.
Checkpoints: each fitted model's test predictions are saved, so a rerun skips models already done (RSF is the slow one).
Run:  python "OneDrive - University College London\\phd_c7_refit.py"
"""
import os
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")      # Windows/Anaconda: PyTorch and numpy/scikit-learn each ship OpenMP
import re, sys, json, time, gc, datetime, warnings, pickle
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")

K = r"S:\FPHS_THIN_Projects\Karan"
DATA = os.path.join(K, "BZI_PhD_T2D", "PhD_descriptive_data", "mice_data", "impv2_phddata_230224.dta")
OD = os.path.join(os.path.expanduser("~"), "OneDrive - University College London")
NB_LINEAR = os.path.join(OD, "PhD_cph_scikitsurv_240324.ipynb")
NB_NONLIN = os.path.join(OD, "PhD_CPHlifelines_nonandinter_230424.ipynb")
OUTD = os.path.join(K, "PURE_v2", "ch7_refit")
DESK = os.path.join(os.path.expanduser("~"), "Desktop", "pure_inventory")
SEED, TEST_SHARE, B_BOOT, T_MAX, T_CAL = 2026, 0.30, 1000, 9.99, 5.0
GRID = np.linspace(0.05, T_MAX, 100)
PENALTIES = [0.0, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1]
MODD = os.path.join(OUTD, "models")
os.makedirs(OUTD, exist_ok=True)
os.makedirs(MODD, exist_ok=True)
os.makedirs(DESK, exist_ok=True)
INTERACTIONS = {"sex_neuro": ("sex", "hx_neuro"), "sex_age": ("sex", "age"), "age_albumin": ("age", "albumin"), "age_alp": ("age", "alp"),
                "age_sodium": ("age", "sodium"), "age_mcv": ("age", "mcv"), "age_smoking": ("age", "smoking_status2"), "age_bmi": ("age", "bmi"),
                "age_resp": ("age", "hx_resp"), "age_liver": ("age", "hx_liver"), "age_severekidney": ("age", "hx_severekidney"),
                "age_cancer": ("age", "hx_cancer"), "age_pvd": ("age", "hx_pvd")}     # phd_updatedmice_180224.do lines 22-34


def dump(obj, name):
    with open(os.path.join(MODD, name), "wb") as fh:
        pickle.dump(obj, fh)
lines = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    lines.append(s)
    open(os.path.join(OUTD, "c7_refit_summary.txt"), "w", encoding="utf-8").write("\n".join(lines))


def read_cols(cols, chunk=20000):
    """Read only `cols` from the very wide Stata file, chunk by chunk (keeps memory low)."""
    parts = []
    for part in pd.read_stata(DATA, columns=list(cols), convert_categoricals=False, chunksize=chunk):
        parts.append(part.copy())
        gc.collect()
    return pd.concat(parts, ignore_index=True)


def cell_lists(nb_path):
    """All column lists of the form df[[...]] or df[['...']] in a notebook, in order."""
    nb = json.load(open(nb_path, encoding="utf-8"))
    out = []
    for c in nb["cells"]:
        if c.get("cell_type") != "code":
            continue
        src = "".join(c.get("source", ""))
        for m in re.finditer(r"(\w+)\s*=\s*\w+\[\[(.*?)\]\]", src, re.S):
            out.append((m.group(1), re.findall(r"['\"]([^'\"]+)['\"]", m.group(2))))
        for m in re.finditer(r"(continuous_columns)\s*=\s*\[(.*?)\]", src, re.S):
            out.append((m.group(1), re.findall(r"['\"]([^'\"]+)['\"]", m.group(2))))
    return out


log("=" * 100, f"\nPhD CHAPTER 7 REFIT   run {datetime.datetime.now():%Y-%m-%d %H:%M}   seed {SEED}\n" + "=" * 100)
lin = cell_lists(NB_LINEAR)
non = cell_lists(NB_NONLIN)
X_LIN = next(v for k, v in lin if k == "X_train")
CONT_LIN = next(v for k, v in lin if k == "continuous_columns")
X_NON = [c for c in next(v for k, v in non if k == "dataset") if c not in ("exposure", "mortality_status")]
CONT_NON = next(v for k, v in non if k == "continuous_columns")
log(f"predictors read from notebooks: CPH original/RSF/DeepSurv {len(X_LIN)} (continuous {len(CONT_LIN)}); "
    f"CPH non-linear+interactions {len(X_NON)} (continuous {len(CONT_NON)})")

need = sorted(set(X_LIN) | set(X_NON) | {"exposure", "mortality_status"})
D = read_cols(need)
for c in need:
    D[c] = pd.to_numeric(D[c], errors="coerce").astype("float32" if c not in ("exposure",) else "float64")
n0 = len(D)
D = D.dropna(subset=need).reset_index(drop=True)
log(f"data: {n0:,} rows; complete on all predictors (imputation 1): {len(D):,}; deaths {int(D['mortality_status'].sum()):,}")
T_all, E_all = D["exposure"].values.astype(float), D["mortality_status"].values.astype(bool)

rng = np.random.default_rng(SEED)
test_idx = np.sort(rng.choice(len(D), int(round(TEST_SHARE * len(D))), replace=False))
is_test = np.zeros(len(D), bool)
is_test[test_idx] = True
tr, te = ~is_test, is_test
log(f"split: train {tr.sum():,} ({int(E_all[tr].sum()):,} deaths) | test {te.sum():,} ({int(E_all[te].sum()):,} deaths) - same test set for all models")
T_tr, E_tr, T_te, E_te = T_all[tr], E_all[tr], T_all[te], E_all[te]


def standardise(cols, cont):
    X = D[cols].copy()
    mu, sd = X.loc[tr, cont].mean(), X.loc[tr, cont].std().replace(0, 1)
    X[cont] = (X[cont] - mu) / sd
    return X, mu, sd


# ---------------------------------------------------------------- metrics (own implementation, used for point estimates AND bootstrap)
def km(times, events):
    o = np.argsort(times)
    t, e = times[o], events[o].astype(float)
    u, idx, cnt = np.unique(t, return_index=True, return_counts=True)
    d = np.add.reduceat(e, idx)
    at_risk = len(t) - idx
    return u, np.cumprod(1 - d / at_risk)


G_u, G_s = km(T_tr, ~E_tr)                                     # censoring distribution from the TRAINING set


def G(t):
    i = np.searchsorted(G_u, t, side="right") - 1
    return np.maximum(np.where(i >= 0, G_s[np.clip(i, 0, None)], 1.0), 0.05)


def ibs(S, t, e):
    """IPCW integrated Brier score over GRID (S: n x len(GRID) predicted survival)."""
    Gt, Ggrid = G(t), G(GRID)
    died = (t[:, None] <= GRID[None, :]) & e[:, None]
    alive = t[:, None] > GRID[None, :]
    bs = (np.where(died, S ** 2 / Gt[:, None], 0) + np.where(alive, (1 - S) ** 2 / Ggrid[None, :], 0)).mean(axis=0)
    trap = getattr(np, "trapezoid", None) or np.trapz
    return float(trap(bs, GRID) / (GRID[-1] - GRID[0]))


from lifelines.utils import concordance_index


def cidx(risk, t, e):
    return float(concordance_index(t, -risk, e))


def one_cal(S_at, t, e, t_star, bins=10):
    """1-calibration at t_star: deciles of predicted survival; observed KM vs mean predicted; Hosmer-Lemeshow chi2 (bins-2 df)."""
    from scipy.stats import chi2 as _chi2
    q = pd.qcut(pd.Series(S_at).rank(method="first"), bins, labels=False).values
    obs, exp_, stat = [], [], 0.0
    for b in range(bins):
        m = q == b
        u, s = km(t[m], e[m])
        o = s[np.searchsorted(u, t_star, side="right") - 1] if (u <= t_star).any() else 1.0
        p_ = S_at[m].mean()
        n = m.sum()
        obs.append(o); exp_.append(p_)
        stat += n * (o - p_) ** 2 / max(p_ * (1 - p_), 1e-9)
    return np.array(obs), np.array(exp_), float(1 - _chi2.cdf(stat, bins - 2))


# ---------------------------------------------------------------- models (checkpointed predictions)
def ckpt(name):
    return os.path.join(OUTD, f"pred_{name}.npz")


def save_pred(name, risk, S):
    np.savez_compressed(ckpt(name), risk=risk, S=S.astype(np.float32))


preds = {}

# 1. CPH original (linear)
name = "CPH_original"
if not os.path.exists(ckpt(name)):
    print(f"[{datetime.datetime.now():%H:%M}] starting {name} ...", flush=True)
    from lifelines import CoxPHFitter
    t0 = time.time()
    X, mu, sd = standardise(X_LIN, CONT_LIN)
    df = X.assign(T=T_all, E=E_all.astype(int))
    m = CoxPHFitter().fit(df[tr], "T", "E")
    risk = m.predict_partial_hazard(X[te]).values.ravel()
    S = m.predict_survival_function(X[te], times=GRID).T.values
    save_pred(name, risk, S)
    dump({"model": m, "columns": X_LIN, "continuous": CONT_LIN, "mean": mu, "sd": sd, "trained_on": "70% training set"}, "cph_original.pkl")
    log(f"fitted {name} in {(time.time() - t0) / 60:.1f} min (model saved); training log-likelihood {m.log_likelihood_:.1f}")
    del m, df, X
    gc.collect()

# 2. CPH non-linear + interactions (also the C7.4 model)
name = "CPH_nonlinear_interactions"
if not os.path.exists(ckpt(name)):
    print(f"[{datetime.datetime.now():%H:%M}] starting {name} (then the PH test) ...", flush=True)
    from lifelines import CoxPHFitter
    t0 = time.time()
    X, mu, sd = standardise(X_NON, CONT_NON)
    df = X.assign(T=T_all, E=E_all.astype(int))
    # The unpenalised fit fails to converge (near-collinear spline terms, e.g. triglycerides_0/_1). Choose a small ridge penalty
    # on a validation split WITHIN the training set (never the test set): best validation IBS among penalties that converge.
    tr_idx = np.where(tr)[0]
    rv = np.random.default_rng(SEED + 11)
    v_idx = np.sort(rv.choice(tr_idx, int(0.2 * len(tr_idx)), replace=False))
    f_idx = np.setdiff1d(tr_idx, v_idx)
    best = None
    for pen in PENALTIES:
        try:
            with warnings.catch_warnings(record=True) as wl:
                warnings.simplefilter("always")
                mm = CoxPHFitter(penalizer=pen).fit(df.iloc[f_idx], "T", "E")
            conv = not any("converge" in str(w_.message).lower() for w_ in wl)
            Sv = mm.predict_survival_function(X.iloc[v_idx], times=GRID).T.values
            vib = ibs(Sv, T_all[v_idx], E_all[v_idx])
            vc = cidx(mm.predict_partial_hazard(X.iloc[v_idx]).values.ravel(), T_all[v_idx], E_all[v_idx])
            log(f"  penalty {pen:<7g} converged {conv} | validation C {vc:.4f} IBS {vib:.4f} | max |coef| {mm.params_.abs().max():.2f}")
            if conv and (best is None or vib < best[1]):
                best = (pen, vib)
        except Exception as ex:
            log(f"  penalty {pen:<7g} failed ({type(ex).__name__})")
    PEN = best[0] if best else 0.01
    log(f"  chosen penalty: {PEN}")
    json.dump({"penalizer": PEN}, open(os.path.join(OUTD, "cph_nonlinear_penalty.json"), "w"))
    m = CoxPHFitter(penalizer=PEN).fit(df[tr], "T", "E")
    log(f"  training log-likelihood: {m.log_likelihood_:.1f}")
    risk = m.predict_partial_hazard(X[te]).values.ravel()
    S = m.predict_survival_function(X[te], times=GRID).T.values
    save_pred(name, risk, S)
    dump({"model": m, "columns": X_NON, "continuous": CONT_NON, "mean": mu, "sd": sd, "penalizer": PEN, "trained_on": "70% training set"}, "cph_nonlinear_train.pkl")
    log(f"fitted {name} in {(time.time() - t0) / 60:.1f} min (model saved)")
    # C7.2: proportional hazards test on a training subsample
    try:
        from lifelines.statistics import proportional_hazard_test
        sub = df[tr].sample(n=min(30000, int(tr.sum())), random_state=SEED)
        m_sub = CoxPHFitter(penalizer=PEN).fit(sub, "T", "E")
        ph = proportional_hazard_test(m_sub, sub, time_transform="rank").summary
        ph.to_csv(os.path.join(OUTD, "c7_ph_test.csv"))
        from scipy.stats import chi2 as _chi2
        glob_stat = ph["test_statistic"].sum()
        log(f"PH test (Schoenfeld, rank time; {len(sub):,}-person training subsample): {int((ph['p'] < 0.05).sum())} of {len(ph)} terms p<0.05, "
            f"{int((ph['p'] < 0.001).sum())} p<0.001; global chi2 {glob_stat:.1f} on {len(ph)} df, p={1 - _chi2.cdf(glob_stat, len(ph)):.3g}")
        log("  largest departures: " + "; ".join(f"{i} (p={r.p:.2g})" for i, r in ph.sort_values("test_statistic", ascending=False).head(8).iterrows()))
    except Exception as ex:
        log(f"PH test failed: {type(ex).__name__}: {str(ex)[:120]}")
    del m, df, X
    gc.collect()

if os.path.exists(os.path.join(OUTD, "c7_ph_test.csv")) and os.path.exists(ckpt("CPH_nonlinear_interactions")):
    ph = pd.read_csv(os.path.join(OUTD, "c7_ph_test.csv"), index_col=0)
    from scipy.stats import chi2 as _chi2
    gs = ph["test_statistic"].sum()
    log(f"PH test (from previous run): {int((ph['p'] < 0.05).sum())} of {len(ph)} terms p<0.05, {int((ph['p'] < 0.001).sum())} p<0.001; "
        f"global chi2 {gs:.1f} on {len(ph)} df, p={1 - _chi2.cdf(gs, len(ph)):.3g}")

# 3. DeepSurv (pycox)
name = "DeepSurv"
if not os.path.exists(ckpt(name)):
    print(f"[{datetime.datetime.now():%H:%M}] starting {name} (one line per training epoch) ...", flush=True)
    try:
        import torch, torchtuples as tt
        from pycox.models import CoxPH
        t0 = time.time()
        np.random.seed(1234); torch.manual_seed(123)
        X, mu, sd = standardise(X_LIN, CONT_LIN)
        Xa = X.values.astype("float32")
        tr_idx = np.where(tr)[0]
        r2 = np.random.default_rng(SEED + 1)
        val_idx = np.sort(r2.choice(tr_idx, int(0.2 * len(tr_idx)), replace=False))
        fit_idx = np.setdiff1d(tr_idx, val_idx)
        y = lambda ix: (T_all[ix].astype("float32"), E_all[ix].astype("float32"))
        net = tt.practical.MLPVanilla(Xa.shape[1], [32, 32], 1, batch_norm=True, dropout=0.1, output_bias=False)
        model = CoxPH(net, tt.optim.Adam)
        model.optimizer.set_lr(0.01)
        model.fit(Xa[fit_idx], y(fit_idx), batch_size=256, epochs=512, callbacks=[tt.callbacks.EarlyStopping()], verbose=True,
                  val_data=(Xa[val_idx], y(val_idx)), val_batch_size=256)
        model.compute_baseline_hazards()
        surv = model.predict_surv_df(Xa[te])                       # index = times, columns = test people
        S = np.vstack([np.interp(GRID, surv.index.values, surv.iloc[:, j].values) for j in range(surv.shape[1])])
        risk = model.predict(Xa[te]).ravel()
        save_pred(name, np.exp(risk), S)
        torch.save(net.state_dict(), os.path.join(MODD, "deepsurv_net.pt"))
        dump({"columns": X_LIN, "continuous": CONT_LIN, "mean": mu, "sd": sd, "layers": [32, 32], "batch_norm": True, "dropout": 0.1,
              "baseline_cumulative_hazards": model.baseline_cumulative_hazards_, "trained_on": "56% (80% of training set)"}, "deepsurv_meta.pkl")
        log(f"fitted {name} in {(time.time() - t0) / 60:.1f} min")
        del model, net, X, Xa
        gc.collect()
    except Exception as ex:
        log(f"DeepSurv failed: {type(ex).__name__}: {str(ex)[:160]}")

# 4. RSF (slowest)
name = "RSF"
if not os.path.exists(ckpt(name)):
    print(f"[{datetime.datetime.now():%H:%M}] starting {name} (32 trees; progress lines from scikit-learn) ...", flush=True)
    try:
        from sksurv.ensemble import RandomSurvivalForest
        from sksurv.util import Surv
        t0 = time.time()
        X = D[X_LIN]
        rsf = RandomSurvivalForest(n_estimators=32, min_samples_split=6, min_samples_leaf=3, max_depth=6, n_jobs=-1, random_state=20, verbose=2)
        rsf_path = os.path.join(MODD, "rsf.pkl")
        if os.path.exists(rsf_path):                               # trained in an earlier run whose predictions failed
            with open(rsf_path, "rb") as fh:
                rsf = pickle.load(fh)["model"]
            log(f"loaded the previously trained {name} from {rsf_path}")
        else:
            rsf.fit(X[tr].values, Surv.from_arrays(E_tr, T_tr))
            dump({"model": rsf, "columns": X_LIN, "trained_on": "70% training set"}, "rsf.pkl")   # saved BEFORE predicting
            log(f"fitted {name} in {(time.time() - t0) / 60:.1f} min (model saved)")
        rsf.verbose = 0
        Xt = X[te].values
        S_parts, r_parts = [], []
        ut = rsf.unique_times_
        for s0 in range(0, len(Xt), 500):                          # small batches keep memory low (~0.5 GB)
            arr = rsf.predict_survival_function(Xt[s0:s0 + 500], return_array=True)   # (batch, n unique times)
            S_parts.append(np.clip(np.vstack([np.interp(GRID, ut, a) for a in arr]), 0, 1))
            r_parts.append(rsf.predict(Xt[s0:s0 + 500]))
            del arr
            gc.collect()
            if (s0 // 500) % 20 == 0:
                print(f"    RSF predictions {min(s0 + 500, len(Xt)):,}/{len(Xt):,}", flush=True)
        S, risk = np.vstack(S_parts), np.concatenate(r_parts)
        save_pred(name, risk, S)
        del rsf, X, Xt
        gc.collect()
    except Exception as ex:
        log(f"RSF failed: {type(ex).__name__}: {str(ex)[:160]}")

# ---------------------------------------------------------------- Table 8 + bootstrap CIs (C7.5)
MODELS = ["CPH_original", "CPH_nonlinear_interactions", "RSF", "DeepSurv"]
for n_ in MODELS:
    if os.path.exists(ckpt(n_)):
        z = np.load(ckpt(n_))
        preds[n_] = (z["risk"], z["S"])
log(f"\nTABLE 8 (corrected) - test set n={te.sum():,}, deaths {int(E_te.sum()):,}; {B_BOOT} bootstrap resamples of the test set")
rows, boots = [], {}
rngb = np.random.default_rng(SEED + 7)
boot_idx = [rngb.integers(0, te.sum(), te.sum()) for _ in range(B_BOOT)]
for n_, (risk, S) in preds.items():
    t0 = time.time()
    c, b = cidx(risk, T_te, E_te), ibs(S, T_te, E_te)
    cs, bs = [], []
    print(f"[{datetime.datetime.now():%H:%M}] bootstrap {n_} ...", flush=True)
    for k, ix in enumerate(boot_idx, 1):
        cs.append(cidx(risk[ix], T_te[ix], E_te[ix]))
        bs.append(ibs(S[ix], T_te[ix], E_te[ix]))
        if k % 100 == 0:
            print(f"    {k}/{B_BOOT} resamples ({(time.time() - t0) / 60:.1f} min)", flush=True)
    boots[n_] = (np.array(cs), np.array(bs), c, b)
    cl, ch = np.percentile(cs, [2.5, 97.5])
    bl, bh = np.percentile(bs, [2.5, 97.5])
    rows.append(dict(model=n_, c_index=c, c_lo=cl, c_hi=ch, ibs=b, ibs_lo=bl, ibs_hi=bh))
    log(f"  {n_:<28} C-index {c:.3f} ({cl:.3f}-{ch:.3f})   IBS {b:.4f} ({bl:.4f}-{bh:.4f})   [{(time.time() - t0) / 60:.1f} min]")
pd.DataFrame(rows).to_csv(os.path.join(OUTD, "c7_table8.csv"), index=False)
# paired comparisons: every model was scored on the SAME bootstrap resamples, so differences have valid paired CIs
pr = []
names_ = list(boots.keys())
log("\nPAIRED DIFFERENCES, all pairs (same test resamples; 95% CI excluding 0 = a real difference)")
log("  C-index: positive = first model better | IBS: negative = first model better")
for i_ in range(len(names_)):
    for j_ in range(i_ + 1, len(names_)):
        a_, b_ = names_[i_], names_[j_]
        ca, ba, c0a, b0a = boots[a_]
        cb, bb, c0b, b0b = boots[b_]
        dcl, dch = np.percentile(ca - cb, [2.5, 97.5])
        dbl, dbh = np.percentile(ba - bb, [2.5, 97.5])
        pr.append(dict(comparison=f"{a_} minus {b_}", c_diff=c0a - c0b, c_lo=dcl, c_hi=dch, ibs_diff=b0a - b0b, ibs_lo=dbl, ibs_hi=dbh))
        log(f"  {a_} minus {b_}: C-index {c0a - c0b:+.4f} ({dcl:+.4f} to {dch:+.4f}) | IBS {b0a - b0b:+.4f} ({dbl:+.4f} to {dbh:+.4f})")
pd.DataFrame(pr).to_csv(os.path.join(OUTD, "c7_paired_differences.csv"), index=False)
del boot_idx
boots_small = {k: (v[2], v[3]) for k, v in boots.items()} if "boots" in globals() else {}
gc.collect()

# SurvivalEVAL cross-check (point estimates), if available
try:
    for p in (os.path.join(OD, "SurvivalEVAL"), os.path.join(OD, "SurvivalEVAL-main")):
        if os.path.isdir(p) and p not in sys.path:
            sys.path.insert(0, p)
    try:
        from SurvivalEVAL import SurvivalEvaluator
    except ImportError:
        from Evaluator import SurvivalEvaluator
    log("\nSurvivalEVAL cross-check (point estimates, same test set):")
    for n_, (risk, S) in preds.items():
        ev = SurvivalEvaluator(S, GRID, T_te, E_te, T_tr, E_tr)
        cc = ev.concordance()[0]
        ib = ev.integrated_brier_score(draw_figure=False) if "draw_figure" in ev.integrated_brier_score.__code__.co_varnames else ev.integrated_brier_score()
        log(f"  {n_:<28} SurvivalEVAL C-index {cc:.3f}   IBS {ib:.4f}")
except Exception as ex:
    log(f"\nSurvivalEVAL cross-check not run ({type(ex).__name__}: {str(ex)[:100]})")

# ---------------------------------------------------------------- calibration figures (1-calibration by decile at T_CAL years)
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    j = int(np.argmin(np.abs(GRID - T_CAL)))
    log(f"\n1-calibration at {GRID[j]:.2f} years (deciles of predicted survival; Hosmer-Lemeshow p):")
    for n_, (risk, S) in preds.items():
        obs, exp_, pval = one_cal(S[:, j], T_te, E_te, GRID[j])
        log(f"  {n_:<28} p = {pval:.3f}")
        fig, ax = plt.subplots(figsize=(8, 5), dpi=120)
        xb = np.arange(10)
        ax.bar(xb - 0.15, 1 - obs, 0.3, label="Observed", color="#1f77b4")
        ax.bar(xb + 0.15, 1 - exp_, 0.3, label="Expected", color="#ff7f0e")
        ax.set_xlabel("Decile of predicted risk"); ax.set_ylabel(f"Probability of death by {GRID[j]:.0f} years")
        ax.text(0.05, 0.85, f"p-value = {pval:.3f}", transform=ax.transAxes, fontsize=12)
        ax.legend(loc="upper left"); fig.tight_layout()
        fig.savefig(os.path.join(OUTD, f"c7_calibration_{n_}.png")); plt.close(fig)
except Exception as ex:
    log(f"calibration figures failed: {type(ex).__name__}: {str(ex)[:120]}")

preds.clear()
gc.collect()

# ---------------------------------------------------------------- C7.4: final CPH refitted on ALL data, coefficients + baseline survival
try:
    from lifelines import CoxPHFitter
    print(f"[{datetime.datetime.now():%H:%M}] C7.4: refitting the final CPH on all data ...", flush=True)
    t0 = time.time()
    X = D[X_NON].copy()
    mu, sd = X[CONT_NON].mean(), X[CONT_NON].std().replace(0, 1)
    X[CONT_NON] = (X[CONT_NON] - mu) / sd
    PEN = json.load(open(os.path.join(OUTD, "cph_nonlinear_penalty.json")))["penalizer"] if os.path.exists(os.path.join(OUTD, "cph_nonlinear_penalty.json")) else 0.0
    m = CoxPHFitter(penalizer=PEN).fit(X.assign(T=T_all, E=E_all.astype(int)), "T", "E")
    log(f"  final model penalty {PEN}; log-likelihood {m.log_likelihood_:.1f}")
    s = m.summary
    tab = pd.DataFrame({"coef": s["coef"], "HR": s["exp(coef)"], "HR_lo95": s["exp(coef) lower 95%"], "HR_hi95": s["exp(coef) upper 95%"],
                        "p": s["p"], "std_mean": [mu.get(i, np.nan) for i in s.index], "std_sd": [sd.get(i, np.nan) for i in s.index]})
    tab.to_csv(os.path.join(OUTD, "c7_cox_final_coefficients.csv"))
    dump({"model": m, "columns": X_NON, "continuous": CONT_NON, "mean": mu, "sd": sd, "trained_on": "all data (final model, C7.4)"}, "cph_nonlinear_final.pkl")
    bh = m.baseline_survival_                                    # at covariate means (continuous = 0 after standardisation)
    binary = [c for c in X_NON if c not in CONT_NON]
    shift = float((X[binary].mean() * m.params_[binary]).sum())
    out = []
    for t in (1, 5, 9):
        s_mean = float(np.interp(t, bh.index.values, bh.iloc[:, 0].values))
        out.append(dict(t_years=t, S0_at_means=s_mean, S0_reference=s_mean ** np.exp(-shift)))
    pd.DataFrame(out).to_csv(os.path.join(OUTD, "c7_baseline_survival.csv"), index=False)
    log(f"\nC7.4 final CPH (non-linear + interactions) on all {len(X):,} people, {len(tab)} terms [{(time.time() - t0) / 60:.1f} min]")
    for r in out:
        log(f"  baseline survival at {r['t_years']} y: {r['S0_reference']:.5f} (reference: continuous terms at their mean, binary terms 0); "
            f"{r['S0_at_means']:.5f} at all covariate means")
    log("  equation: S(t|x) = S0(t) ^ exp( sum_j b_j * z_j ), z_j = (x_j - mean_j)/sd_j for continuous terms, x_j for binary terms")
except Exception as ex:
    log(f"C7.4 refit failed: {type(ex).__name__}: {str(ex)[:160]}")

# ---------------------------------------------------------------- C7.2: spline knots from the Stata variable labels
try:
    from pandas.io.stata import StataReader
    with StataReader(DATA) as rdr:
        labels = rdr.variable_labels()
    spl = {k: v for k, v in labels.items() if re.search(r"_[0-3]$", k) and not k.startswith("_") and v}
    log(f"\nspline variable labels ({len(spl)}; mvrs usually records knot positions here):")
    for k, v in list(spl.items())[:60]:
        log(f"  {k}: {v}")
except Exception as ex:
    log(f"spline labels not read: {type(ex).__name__}")

# ---------------------------------------------------------------- spline maps: raw value -> mvrs spline terms (for new individuals)
try:
    from pandas.io.stata import StataReader
    with StataReader(DATA) as rdr:
        allcols = list(rdr.variable_labels().keys())
    bases = {}
    for c in X_NON:
        mm = re.match(r"^(_1_)?(.+)_(\d)$", c)
        if mm and not c.startswith(("sex_", "age_", "_1_age_")) or (mm and c in ("age_0", "age_1", "age_2", "age_3")):
            bases.setdefault(mm.group(2), []).append(c)
    want = set()
    for b, terms in bases.items():
        want |= {b} | {t.replace("_1_", "", 1) for t in terms} | set(terms) | {f"_1_{b}"}
    want = [c for c in want if c in allcols]
    R = read_cols(want)
    maps, rr = {}, np.random.default_rng(SEED + 3)
    log("\nSPLINE MAPS (raw value -> spline terms, learned from people with an OBSERVED raw value; max error on held-out rows):")
    for b, terms in sorted(bases.items()):
        raw = b if b in R.columns else f"_1_{b}"
        orig_terms = [t.replace("_1_", "", 1) if t.replace("_1_", "", 1) in R.columns else t for t in terms]
        sub = R[[raw] + orig_terms].dropna().astype(float)
        hold = rr.random(len(sub)) < 0.1
        def build(dd):
            g = dd.groupby(dd[raw].round(6))[orig_terms].mean().sort_index()
            return g.index.values.astype(float), g.values.astype(float)
        xs, ys = build(sub[~hold])
        err = max(float(np.nanmax(np.abs(np.interp(sub.loc[hold, raw].values, xs, ys[:, k]) - sub.loc[hold, orig_terms[k]].values))) for k in range(len(terms)))
        xs, ys = build(sub)
        maps[b] = {"raw": b, "terms": terms, "x": xs, "y": ys, "range": (float(xs.min()), float(xs.max()))}
        log(f"  {b:<26} {len(terms)} terms; from {len(sub):,} observed values; range {xs.min():.3g}-{xs.max():.3g}; max held-out error {err:.2e}")
    qcols = set(X_NON) | {"age", "sex"} | {v for pair in INTERACTIONS.values() for v in pair} | {f"_1_{v}" for pair in INTERACTIONS.values() for v in pair}
    Q = read_cols([c for c in qcols if c in allcols])
    log("  interaction terms (raw products, phd_updatedmice_180224.do), checked on people with OBSERVED values; scale learned from data:")
    scales = {}
    for name_, (a_, b_) in INTERACTIONS.items():
        col = next((c for c in (f"_1_{name_}", name_) if c in Q.columns), None)
        av = next((c for c in (f"_1_{a_}", a_) if c in Q.columns), None)
        bv = next((c for c in (f"_1_{b_}", b_) if c in Q.columns), None)
        if not (col and av and bv):
            log(f"    {name_:<18} columns not found - assumed scale 1")
            scales[name_] = 1.0
            continue
        obs = Q[col].notna() & Q[av].notna() & Q[bv].notna()
        for orig in (a_, b_):                                   # restrict to people whose raw values were observed (not imputed)
            if orig in Q.columns and orig != av and orig != bv:
                obs &= Q[orig].notna()
        x = (Q.loc[obs, av] * Q.loc[obs, bv]).astype(float).values
        y = Q.loc[obs, col].astype(float).values
        sc = float((x @ y) / (x @ x)) if (x @ x) > 0 else 1.0
        sc = 1.0 if abs(sc - 1) < 1e-6 else sc
        rel = float(np.max(np.abs(y - sc * x)) / max(np.max(np.abs(y)), 1e-9))
        scales[name_] = sc
        log(f"    {name_:<18} = {sc:g} x {a_} x {b_}   (n={int(obs.sum()):,}; max relative error {rel:.2e}{'  <-- CHECK' if rel > 1e-3 else ''})")
    dump({"maps": maps, "interactions": INTERACTIONS, "interaction_scales": scales, "sex_values": sorted(pd.unique(D["sex"].dropna()).tolist()),
          "binary_columns": [c for c in X_LIN if c not in CONT_LIN]}, "spline_maps.pkl")
    del R, Q
except Exception as ex:
    log(f"spline maps failed: {type(ex).__name__}: {str(ex)[:160]}")

import importlib
vers = []
for pkg in ("lifelines", "sksurv", "pycox", "torch", "torchtuples", "sklearn", "numpy", "pandas"):
    try:
        vers.append(f"{pkg} {importlib.import_module(pkg).__version__}")
    except Exception:
        pass
log(f"\nsoftware: Python {sys.version.split()[0]}; " + ", ".join(vers))
txt = "\n".join(lines)
open(os.path.join(DESK, "c7_refit_summary.txt"), "w", encoding="utf-8").write(txt)
print(f"\nDone. Outputs in {OUTD}; summary copied to Desktop\\pure_inventory\\c7_refit_summary.txt")
