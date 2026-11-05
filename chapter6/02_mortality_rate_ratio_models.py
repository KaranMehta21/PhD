"""
PhD CHAPTER 6 CORRECTIONS (C6.1, C6.2, C6.4) - mortality rate ratios, T2D vs non-T2D (Python replacement for the Stata do-file).
Poisson regression of deaths with log person-years as offset; practice clustering handled with CLUSTER-ROBUST standard errors
by practice (for log-link Poisson a practice random intercept only shifts the baseline rate, so rate ratios target the same
quantity as -mepoisson-). Data are aggregated into cells (identical likelihood to person-level data, much faster).
Input : S:\\FPHS_THIN_Projects\\Karan\\PURE_v2\\ch6_2009\\ch6_analysis.dta  (from phd_ch6_2009_build.py)
        original t2d_data.dta / nondiab_data.dta (to reproduce the original Table 5 for comparison)
Models: ORIG  original Table 5 as published (T2D vs non-T2D, no age adjustment overall/by sex; by age band x sex strata)
        A     age band x sex (overall); age band (by sex); stratum-specific (age band x sex)
        B     A + Townsend quintile (6 = missing level) + smoking (never/ex/current/missing)   [PRIMARY]
        Ccc   B + BMI category, complete cases          Cmi  B + BMI with a missing level        [BMI may partly mediate]
        D     B + CVD history + antihypertensive use    [mediator sensitivity]
        A, B without censoring controls at T2D diagnosis (original approach)                   [sensitivity]
        Interaction tests (Model A): T2D x sex, T2D x age band (cluster-robust Wald tests)
        E-values for Model B (VanderWeele & Ding 2017)
Also  : corrected Table 5 layout - mortality rate per 1,000 pyar (exact 95% CI) by age band and sex, each arm.
Output: ...\\ch6_2009\\ch6_results.csv, ch6_table5_rates.csv, ch6_models_summary.txt (copy on Desktop\\pure_inventory). ~5-15 min.
Run:  python "OneDrive - University College London\\phd_ch6_2009_models.py"
"""
import os, datetime, gc
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import chi2

K = r"S:\FPHS_THIN_Projects\Karan"
D = os.path.join(K, "PURE_v2", "ch6_2009")
P = os.path.join(K, "BZI_PhD_T2D", "PhD_descriptive_data")
DESK = os.path.join(os.path.expanduser("~"), "Desktop", "pure_inventory")
BANDS = {i: f"{30 + 5 * i}-{34 + 5 * i}" for i in range(1, 12)}          # 1 = 35-39 ... 11 = 85-89
lines, rows = [], []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    lines.append(s)


def fit(data, died, fu, rhs, model, stratum, analysis, cond=None):
    """Aggregate to cells, fit Poisson with offset and practice cluster-robust SEs; return the T2D rate ratio."""
    d = data if cond is None else data[cond(data)]
    n, deaths = len(d), int(d[died].sum())
    terms = sorted({t for t in ["diab", "prac"] + [v for v in ["ageband_cat", "sex", "townsend_cat", "smoking_cat", "bmi_cat", "cvd", "antihtn"]
                                                     if v in rhs]})
    g = d.groupby(terms, observed=True).agg(deaths=(died, "sum"), py=(fu, "sum")).reset_index()
    g = g[g["py"] > 0]
    try:
        m = smf.glm(f"deaths ~ diab{(' + ' + rhs) if rhs else ''}", data=g, family=sm.families.Poisson(),
                    offset=np.log(g["py"])).fit(cov_type="cluster", cov_kwds={"groups": g["prac"].values})
        b, se = m.params["diab"], m.bse["diab"]
        rows.append(dict(model=model, stratum=stratum, analysis=analysis, irr=np.exp(b), lo=np.exp(b - 1.959964 * se),
                         hi=np.exp(b + 1.959964 * se), n=n, deaths=deaths))
        del m, g
        gc.collect()
        return None
    except Exception as ex:
        rows.append(dict(model=model, stratum=stratum, analysis=analysis, irr=np.nan, lo=np.nan, hi=np.nan, n=n, deaths=deaths))
        log(f"  {model} {stratum} {analysis}: failed ({type(ex).__name__}: {str(ex)[:80]})")
        gc.collect()
        return None


def evalue(rr):
    if not np.isfinite(rr):
        return np.nan
    rr = rr if rr >= 1 else 1 / rr
    return rr + np.sqrt(rr * (rr - 1))


log("=" * 100, f"\nPhD CHAPTER 6 - MORTALITY RATE RATIOS (Python)   run {datetime.datetime.now():%Y-%m-%d %H:%M}\n" + "=" * 100)

# ---------------------------------------------------------------- ORIGINAL Table 5 reproduced
t2 = pd.read_stata(os.path.join(P, "t2d_data.dta"), convert_categoricals=False, columns=["prac", "sex", "age", "mortality_status", "exposure"])
nd = pd.read_stata(os.path.join(P, "nondiab_data.dta"), convert_categoricals=False)
nd = nd.rename(columns={"index_age": "age"}) if "age" not in nd.columns and "index_age" in nd.columns else nd
O = pd.concat([t2.assign(diab=1), nd[["prac", "sex", "age", "mortality_status", "exposure"]].assign(diab=0)], ignore_index=True)
O["prac"] = pd.factorize(O["prac"].astype(str))[0].astype(np.int32)
O["sex"] = pd.to_numeric(O["sex"], errors="coerce")
O["ageband_cat"] = np.clip(((pd.to_numeric(O["age"], errors="coerce") - 35) // 5 + 1), 1, 11).astype("Int64")
O = O[O["exposure"] > 0]
log(f"original data: T2D {int((O['diab'] == 1).sum()):,}, non-T2D {int((O['diab'] == 0).sum()):,}")
O = O[["diab", "prac", "sex", "ageband_cat", "mortality_status", "exposure"]].dropna()
O["ageband_cat"] = O["ageband_cat"].astype(np.int8)
fit(O, "mortality_status", "exposure", "", "ORIG", "overall", "orig")
for s, lab in ((1, "men"), (2, "women")):
    fit(O, "mortality_status", "exposure", "", "ORIG", lab, "orig", cond=lambda x, s=s: x["sex"] == s)
    for a in range(1, 12):
        fit(O, "mortality_status", "exposure", "", "ORIG", f"{lab} {BANDS[a]}", "orig", cond=lambda x, s=s, a=a: (x["sex"] == s) & (x["ageband_cat"] == a))
fit(O, "mortality_status", "exposure", "C(prac)", "ORIG-FE", "overall", "orig")
del O, t2, nd
gc.collect()

# ---------------------------------------------------------------- CORRECTED (aligned 2009-2018 cohort)
A = pd.read_stata(os.path.join(D, "ch6_analysis.dta"), convert_categoricals=False)
A["prac"] = pd.factorize(A["prac"].astype(str))[0].astype(np.int32)
for c in ("ageband_cat", "sex", "townsend_cat", "smoking_cat", "bmi_cat", "cvd", "antihtn", "diab", "died_main", "died_sens"):
    A[c] = A[c].astype(np.int8)
A = A[["diab", "prac", "sex", "ageband_cat", "townsend_cat", "smoking_cat", "bmi_cat", "cvd", "antihtn",
       "died_main", "fu_main", "died_sens", "fu_sens"]]
gc.collect()
log(f"corrected data: {len(A):,} people (T2D {int((A['diab'] == 1).sum()):,})")
CA, CS = "C(ageband_cat) * C(sex)", "C(ageband_cat)"
fit(A, "died_sens", "fu_sens", "", "NEWraw", "overall", "nocens")
for s, lab in ((1, "men"), (2, "women")):
    fit(A, "died_sens", "fu_sens", "", "NEWraw", lab, "nocens", cond=lambda x, s=s: x["sex"] == s)
fit(A, "died_main", "fu_main", f"{CA} + C(prac)", "A-FE", "overall", "main")
Bx = "C(townsend_cat) + C(smoking_cat)"
for v, an in (("main", "main"), ("sens", "nocens")):
    dd, ff = f"died_{v}", f"fu_{v}"
    mA = fit(A, dd, ff, CA, "A", "overall", an)
    fit(A, dd, ff, f"{CA} + {Bx}", "B", "overall", an)
    for s, lab in ((1, "men"), (2, "women")):
        fit(A, dd, ff, CS, "A", lab, an, cond=lambda x, s=s: x["sex"] == s)
        fit(A, dd, ff, f"{CS} + {Bx}", "B", lab, an, cond=lambda x, s=s: x["sex"] == s)
    if v == "main":
        fit(A, dd, ff, f"{CA} + {Bx} + C(bmi_cat)", "Ccc", "overall", an, cond=lambda x: x["bmi_cat"] < 7)
        fit(A, dd, ff, f"{CA} + {Bx} + C(bmi_cat)", "Cmi", "overall", an)
        fit(A, dd, ff, f"{CA} + {Bx} + cvd + antihtn", "D", "overall", an)
        for s, lab in ((1, "men"), (2, "women")):
            fit(A, dd, ff, f"{CS} + {Bx} + C(bmi_cat)", "Ccc", lab, an, cond=lambda x, s=s: (x["sex"] == s) & (x["bmi_cat"] < 7))
            fit(A, dd, ff, f"{CS} + {Bx} + C(bmi_cat)", "Cmi", lab, an, cond=lambda x, s=s: x["sex"] == s)
            fit(A, dd, ff, f"{CS} + {Bx} + cvd + antihtn", "D", lab, an, cond=lambda x, s=s: x["sex"] == s)
            for a in range(1, 12):
                c = lambda x, s=s, a=a: (x["sex"] == s) & (x["ageband_cat"] == a)
                fit(A, dd, ff, "", "A", f"{lab} {BANDS[a]}", an, cond=c)
                fit(A, dd, ff, Bx, "B", f"{lab} {BANDS[a]}", an, cond=c)

# interaction tests (Model A, main analysis): does the T2D rate ratio differ by sex / by age band?
g = A.groupby(["diab", "prac", "ageband_cat", "sex"]).agg(deaths=("died_main", "sum"), py=("fu_main", "sum")).reset_index()
g = g[g["py"] > 0]
grp = g["prac"].values
for lab, f in (("T2D x sex", f"deaths ~ diab * C(sex) + {CA}"), ("T2D x age band", f"deaths ~ diab * C(ageband_cat) + {CA}")):
    try:
        m = smf.glm(f, data=g, family=sm.families.Poisson(), offset=np.log(g["py"])).fit(cov_type="cluster", cov_kwds={"groups": grp})
        names = [n for n in m.params.index if n.startswith("diab:")]
        R = np.zeros((len(names), len(m.params)))
        for i, n in enumerate(names):
            R[i, list(m.params.index).index(n)] = 1
        w = m.wald_test(R, scalar=True)
        log(f"interaction {lab}: Wald chi2 = {float(w.statistic):.1f} on {len(names)} df, p = {float(w.pvalue):.2g}")
    except Exception as ex:
        log(f"interaction {lab}: failed ({type(ex).__name__})")

# corrected Table 5 layout: rates per 1,000 pyar by age band and sex
t5 = []
for s, lab in ((1, "Men"), (2, "Women")):
    for a in [None] + list(range(1, 12)):
        rec = {"sex": lab, "age": "Overall" if a is None else BANDS[a]}
        for arm, nm in ((1, "T2D"), (0, "non-T2D")):
            x = A[(A["sex"] == s) & (A["diab"] == arm) & ((A["ageband_cat"] == a) if a else True)]
            d, py = int(x["died_main"].sum()), x["fu_main"].sum()
            rec[f"{nm} deaths"], rec[f"{nm} pyar"] = d, round(py)
            rec[f"{nm} rate (95% CI)"] = (f"{1000 * d / py:.1f} ({chi2.ppf(0.025, 2 * d) / 2 / py * 1000 if d else 0:.1f} - "
                                          f"{chi2.ppf(0.975, 2 * (d + 1)) / 2 / py * 1000:.1f})") if py else ""
        t5.append(rec)
T5 = pd.DataFrame(t5)

R = pd.DataFrame(rows)
R["evalue"] = np.where(R["model"] == "B", R["irr"].map(evalue), np.nan)
R["evalue_ci"] = np.where(R["model"] == "B", np.where(R["lo"] > 1, R["lo"].map(evalue), np.where(R["hi"] < 1, R["hi"].map(evalue), 1.0)), np.nan)
R.to_csv(os.path.join(D, "ch6_results.csv"), index=False)
T5.to_csv(os.path.join(D, "ch6_table5_rates.csv"), index=False)

fmt = lambda r: f"{r.irr:.2f} ({r.lo:.2f}-{r.hi:.2f})" if np.isfinite(r.irr) else "n/a"
log("\nMRR (95% CI), T2D vs non-T2D - overall and by sex")
log(f"{'model':<7} {'analysis':<8} {'overall':<20} {'men':<20} {'women':<20} E-value (CI) overall")
for (mdl, an), grp_ in R[R["stratum"].isin(["overall", "men", "women"])].groupby(["model", "analysis"], sort=False):
    get = lambda s: grp_[grp_["stratum"] == s].iloc[0] if (grp_["stratum"] == s).any() else None
    o, m_, w_ = get("overall"), get("men"), get("women")
    ev = f"{o.evalue:.2f} ({o.evalue_ci:.2f})" if o is not None and np.isfinite(o.evalue) else ""
    log(f"{mdl:<7} {an:<8} {fmt(o) if o is not None else '':<20} {fmt(m_) if m_ is not None else '':<20} {fmt(w_) if w_ is not None else '':<20} {ev}")
def pick(mdl, st, an):
    r = R[(R["model"] == mdl) & (R["stratum"] == st) & (R["analysis"] == an)]
    return fmt(r.iloc[0]) if len(r) else "n/a"
log("\nSTEPWISE DECOMPOSITION - what changes the MRR (overall / men / women)")
log("  0 published Table 5 (Stata mepoisson, original data)        overall n/a | men 1.6 (1.5-1.6) | women 1.7 (1.6-1.7)")
for lab, mdl, an in (("1 original data, refitted in Python (method only)", "ORIG", "orig"),
                     ("2 NEW 2009 cohort, same unadjusted model, no censoring", "NEWraw", "nocens"),
                     ("3 + age adjustment (Model A), no censoring", "A", "nocens"),
                     ("4 + censoring controls at T2D diagnosis (Model A main)", "A", "main"),
                     ("5 + Townsend and smoking (Model B, primary)", "B", "main")):
    log(f"  {lab:<62} overall {pick(mdl, 'overall', an)} | men {pick(mdl, 'men', an)} | women {pick(mdl, 'women', an)}")
log(f"  CHECK practice fixed effects: original data {pick('ORIG-FE', 'overall', 'orig')} vs robust {pick('ORIG', 'overall', 'orig')}; "
    f"Model A {pick('A-FE', 'overall', 'main')} vs robust {pick('A', 'overall', 'main')}")
log("\nMRR by age band and sex - ORIGINAL vs corrected Model A and Model B (main)")
log(f"{'stratum':<16} {'ORIG':<18} {'A':<18} {'B':<18}")
for s in ("men", "women"):
    for a in range(1, 12):
        st = f"{s} {BANDS[a]}"
        cells = []
        for mdl, an in (("ORIG", "orig"), ("A", "main"), ("B", "main")):
            r = R[(R["model"] == mdl) & (R["stratum"] == st) & (R["analysis"] == an)]
            cells.append(fmt(r.iloc[0]) if len(r) else "")
        log(f"{st:<16} {cells[0]:<18} {cells[1]:<18} {cells[2]:<18}")
log("\nCorrected Table 5 rates per 1,000 pyar (main analysis): see ch6_table5_rates.csv")
log(T5[["sex", "age", "T2D rate (95% CI)", "non-T2D rate (95% CI)"]].to_string(index=False))
txt = "\n".join(lines)
open(os.path.join(D, "ch6_models_summary.txt"), "w", encoding="utf-8").write(txt)
os.makedirs(DESK, exist_ok=True)
open(os.path.join(DESK, "ch6_models_summary.txt"), "w", encoding="utf-8").write(txt)
print(f"\nSaved {os.path.join(D, 'ch6_results.csv')}, ch6_table5_rates.csv and Desktop\\pure_inventory\\ch6_models_summary.txt")
