"""
PhD CHAPTER 6 CORRECTIONS (C6.1, C6.2, C6.4) - 2009 DESIGN: thesis T2D cohort vs freshly matched controls.  ~4-5 h
(most time: reading the raw measurement and prescription files for BMI, smoking and antihypertensives).

T2D ARM   t2d_data.dta (155,399; index 1 Jan 2009); 'patprac' = patid+prac links to the registration files.
CONTROLS  drawn fresh at 1 Jan 2009 from patRecs_for_match_*.dta (one row per patient): practice with AMR and ACU dates,
          age 35-89, registered >= 9 months, alive and still registered at index, not diagnosed with T2D on/before index
          (complete_t2dm_cohort.dta). Frequency-matched 5:1 on 5-year age band x sex x practice (fixed seed); strata with
          fewer than 5 eligible controls per T2D person use all available.
BOTH      restricted to practices with control covariate files (AHD, antihypertensives, CVD records).
          Entry = latest(1 Jan 2009, registration + 9 months, AMR, ACU); exit = earliest(death, transfer out, practice
          last collection, 99th birthday, 31 Dec 2018). Main: controls censored at T2D diagnosis; sensitivity: not censored.
          Covariates at entry, identical code for both arms: BMI (recorded BMI, 15-60, latest within 24 months), smoking
          (latest status; never-after-ever = ex), Townsend (6/0/missing = missing level), CVD history (MI, CAD, angina,
          stroke, PVD, CABG/PCI, TIA, amputation; same code lists), antihypertensive use (5 classes, 12 months; same lists).
Outputs   PURE_v2\\ch6_2009\\ch6_analysis.dta (no patient ids) and ch6_build_summary.txt (copy on Desktop\\pure_inventory).
Then run phd_ch6_2009_models.do in Stata.
Run:  python "OneDrive - University College London\\phd_ch6_2009_build.py"
"""
import os, glob, datetime, gc
import numpy as np
import pandas as pd
from scipy.stats import chi2

K = r"S:\FPHS_THIN_Projects\Karan"
TP = os.path.join(K, "T2DM_python")
PV = os.path.join(K, "PURE_v2")
OUTD = os.path.join(PV, "ch6_2009")
DESK = os.path.join(os.path.expanduser("~"), "Desktop", "pure_inventory")
T2D_DATA = os.path.join(K, "BZI_PhD_T2D", "PhD_descriptive_data", "t2d_data.dta")
COMPLETE = os.path.join(K, "main Data", "complete_t2dm_cohort.dta")
POOL = os.path.join(TP, "nondiab_extract", "patRecs_for_match_*.dta")
PATREC_DIR = os.path.join(K, "BZI_PhD_T2D", "patrecs_marital")
TOWNSEND_DIR = os.path.join(K, "BZI_PhD_T2D", "pvirecs_townsend")
T2D_AHD = os.path.join(K, "BZI_PhD_T2D", "ahdfiles")
ND_AHD = os.path.join(TP, "nondiabahdfiles")
ND_RX = os.path.join(TP, "nondiab_antihypertensive_files")       # FULL pool (the 'filtered' folder holds old controls only)
ND_CVD = os.path.join(TP, "nondiab_cvdcomps_medrecs")             # FULL pool
CODES2A = os.path.join(PV, "step2a_parts")
RX4A = os.path.join(PV, "step4a_parts")
NEW = os.path.join(TP, "codelists", "new_codelists")
CM = os.path.join(NEW, "co-morbidities")
MED, MED_OLD = os.path.join(NEW, "medications"), os.path.join(NEW, "provisional_cover", "medications")
DICT = os.path.join(TP, "codelists", "readcodes1901.csv")
INDEX = pd.Timestamp("2009-01-01")
END = pd.Timestamp("2018-12-31")
REG_LIMIT = INDEX - pd.DateOffset(months=9)
RATIO, SEED = 5, 20090101
CHECKPOINT = os.path.join(OUTD, "checkpoint_matched_cohort.pkl")
HEART = {"mi": os.path.join(CM, "GenRe_CVD", "MI.csv"), "cad": os.path.join(CM, "GenRe_CVD", "CAD.csv"),
         "angina": os.path.join(CM, "angina_rclist.csv"), "cabg_pci": os.path.join(CM, "cabg_ptca_stent_rclist.csv"),
         "stroke": os.path.join(CM, "stroke_rclist.csv"), "tia": os.path.join(CM, "TIA_rclist.csv"),
         "pvd": os.path.join(CM, "new_pvd.csv"), "amputation": os.path.join(CM, "amputation_rclist.csv")}
CABG_WORDS = ("bypass", "cabg", "angioplast", "stent", "ptca", "percutaneous", "revasculari", "coronary artery op",
              "coronary endarterectomy", "coronary artery graft", "graft of coronary")
BP = {"acei_arb": ["acei_arb_rclist", "acei_arb_list"], "alphablockers": ["alphablockers_rclist"],
      "betablockers": ["betablockers_rclist", "betablockers_list"], "ca_blockers": ["cachannel_blockers_rclist", "cachannel_blockers_list"],
      "thiazides": ["thiazides_rclist"]}
lines = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    lines.append(s)


def as_date(s):
    s = pd.Series(s)
    if pd.api.types.is_datetime64_any_dtype(s):
        return s
    if pd.api.types.is_numeric_dtype(s):
        v = s.dropna()
        if len(v) and v.between(-40000, 40000).mean() > 0.9:
            return pd.to_datetime("1960-01-01") + pd.to_timedelta(s, unit="D")
    return pd.to_datetime(s.astype(str).str.strip().replace({"": None, "nan": None, "NaT": None}), errors="coerce")


def npat(s):
    """normalised patient id (string, no leading zeros, no '.0')"""
    s = pd.Series(s).astype(str).str.strip().str.replace(r"\.0$", "", regex=True).str.lstrip("0")
    return s.replace("", "0")


def nprac(s):
    return pd.Series(s).astype(str).str.strip().str.lower()


def mkkey(pat, prac):
    return npat(pat).values + "|" + nprac(prac).values


def prac_of(path):
    return os.path.splitext(os.path.basename(path).replace(".gz", ""))[0].split("_")[-1].lower()


def yob(s):
    return pd.to_numeric(pd.Series(s).astype(str).str[:4], errors="coerce")


def attrition(rows, arm):
    log(f"\n  {arm}")
    prev = None
    for lab, n in rows:
        log(f"    {lab:<52} {n:>11,}" + ("" if prev is None else f"   (excluded {prev - n:,})"))
        prev = n


log("=" * 100, f"\nPhD CHAPTER 6 - 2009 MATCHED DESIGN   run {datetime.datetime.now():%Y-%m-%d %H:%M}\n" + "=" * 100)
os.makedirs(OUTD, exist_ok=True)
cov_pracs = {prac_of(f) for f in glob.glob(os.path.join(ND_AHD, "*.pkl"))} & {prac_of(f) for f in glob.glob(os.path.join(ND_RX, "*.pkl"))} & \
            {prac_of(f) for f in glob.glob(os.path.join(ND_CVD, "*.pkl"))}
log(f"practices with all three control covariate sources (AHD, antihypertensives, CVD): {len(cov_pracs)}")

if os.path.exists(CHECKPOINT):
    A = pd.read_pickle(CHECKPOINT)
    prev = open(os.path.join(OUTD, "ch6_build_summary_part1.txt"), encoding="utf-8").read() if os.path.exists(os.path.join(OUTD, "ch6_build_summary_part1.txt")) else ""
    lines.append(prev)
    print(f"resumed from checkpoint ({len(A):,} people); delete {CHECKPOINT} to rebuild the cohort", flush=True)
else:
    # ====================================================================== T2D arm
    t2 = pd.read_stata(T2D_DATA, convert_categoricals=False)
    t2["prac_n"] = nprac(t2["prac"])
    pp = t2["patprac"].astype(str).str.strip()
    pr = t2["prac"].astype(str).str.strip()
    patid = [p[:-len(r)] if a.lower().endswith(r.lower()) else (p[len(r):] if a.lower().startswith(r.lower()) else p) for p, r, a in zip(pp, pr, pp)]
    t2["key"] = mkkey(patid, t2["prac"])
    regs = []
    for f in sorted(glob.glob(os.path.join(PATREC_DIR, "*.pkl"))):
        p = pd.read_pickle(f)
        regs.append(p[[c for c in ("patid", "prac", "sex", "d_regdate", "d_xferdate", "d_deathdate", "d_yob", "amrdate", "acu_date", "lastdate") if c in p.columns]])
    reg = pd.concat(regs, ignore_index=True)
    reg["key"] = mkkey(reg["patid"], reg["prac"])
    reg = reg.drop_duplicates("key")
    cc = pd.read_stata(COMPLETE, convert_categoricals=False, columns=["patid", "prac", "diabetes_start"])
    cc["key"] = mkkey(cc["patid"], cc["prac"])
    cc["ds"] = as_date(cc["diabetes_start"])
    cc = cc.drop_duplicates("key")
    T = t2[["key", "prac_n", "sex"]].merge(reg.drop(columns=["patid", "prac", "sex"], errors="ignore"), on="key", how="left") \
        .merge(cc[["key", "ds"]], on="key", how="left")
    log(f"\nT2D arm: t2d_data {len(t2):,} | linked to registration records {int(T['d_regdate'].notna().sum()):,} "
        f"({100 * T['d_regdate'].notna().mean():.1f}%) | to diagnosis dates {int(T['ds'].notna().sum()):,}")
    T = T.rename(columns={"prac_n": "prac"})
    T["sex"] = pd.to_numeric(T["sex"], errors="coerce")
    T["yob"] = yob(T["d_yob"])
    T["regdate"], T["amr"], T["acu"] = as_date(T["d_regdate"]), as_date(T["amrdate"]), as_date(T["acu_date"])
    T["death"] = as_date(T["d_deathdate"])
    T["end"] = pd.concat([as_date(T["d_xferdate"]), as_date(T["lastdate"])], axis=1).min(axis=1)
    T["diab"] = 1


    def eligible(df, arm_is_t2d):
        age = INDEX.year - df["yob"] - 0.5
        steps = [("linked to registration record", df["regdate"].notna() | df["yob"].notna()),
                 ("practice has AMR and ACU dates", df["amr"].notna() & df["acu"].notna()),
                 ("age 35-89 at index", age.between(35, 89.999)),
                 ("registered >= 9 months before index", df["regdate"] <= REG_LIMIT),
                 ("alive at index", df["death"].isna() | (df["death"] > INDEX)),
                 ("still registered at index", df["end"] > INDEX)]
        steps.append(("T2D diagnosed on or before index" if arm_is_t2d else "not diagnosed with T2D on or before index",
                      (df["ds"] <= INDEX) if arm_is_t2d else ~(df["ds"] <= INDEX)))
        keep, rows = pd.Series(True, index=df.index), [("source", len(df))]
        for lab, m in steps:
            keep &= m.fillna(False)
            rows.append((lab, int(keep.sum())))
        return keep, rows


    kT, rowsT = eligible(T, True)
    T = T[kT].copy()
    T["band"] = ((INDEX.year - T["yob"] - 0.5) // 5 * 5).astype(int)

    # ====================================================================== control pool and matching
    log("\nreading the control pool (patRecs_for_match) ...")
    diab_keys_before = set(cc.loc[cc["ds"] <= INDEX, "key"]) | set(t2["key"])
    pool, rowsC = [], None
    files = sorted(f for f in glob.glob(POOL) if "test" not in os.path.basename(f))
    cc_ds = cc.set_index("key")["ds"]
    for i, f in enumerate(files):
        try:
            r = pd.read_stata(f, convert_categoricals=False)
        except Exception as ex:
            log(f"  could not read {os.path.basename(f)}: {type(ex).__name__}")
            continue
        if "end" not in r.columns and "_end" in r.columns:
            r = r.rename(columns={"_end": "end"})
        r = r[[c for c in ("patid", "prac", "sex", "d_yob", "d_regdate", "d_deathdate", "amrdate", "acu_date", "end", "townsend") if c in r.columns]].copy()
        r["key"] = mkkey(r["patid"], r["prac"])
        r = r.drop_duplicates("key")
        r["prac"] = nprac(r["prac"])
        r["sex"] = pd.to_numeric(r["sex"], errors="coerce")
        r["yob"] = yob(r["d_yob"])
        r["regdate"], r["amr"], r["acu"] = as_date(r["d_regdate"]), as_date(r["amrdate"]), as_date(r["acu_date"])
        r["death"], r["end"] = as_date(r["d_deathdate"]), as_date(r["end"])
        r["ds"] = r["key"].map(cc_ds)
        r.loc[r["key"].isin(diab_keys_before), "ds"] = r.loc[r["key"].isin(diab_keys_before), "ds"].fillna(INDEX)
        k, rows = eligible(r, False)
        rowsC = rows if rowsC is None else [(a, x + y) for (a, x), (_, y) in zip(rowsC, rows)]
        pool.append(r.loc[k, ["key", "prac", "sex", "yob", "regdate", "amr", "acu", "death", "end", "ds", "townsend"]])
        if (i + 1) % 100 == 0:
            print(f"  {i + 1}/{len(files)} register files", flush=True)
    C = pd.concat(pool, ignore_index=True).drop_duplicates("key")
    del pool
    C["band"] = ((INDEX.year - C["yob"] - 0.5) // 5 * 5).astype(int)
    C["diab"] = 0

    # practice restriction (both arms)
    nT0, pT0 = len(T), T["prac"].nunique()
    T = T[T["prac"].isin(cov_pracs)]
    C = C[C["prac"].isin(cov_pracs)]
    rowsT.append(("in a practice with control covariate files", len(T)))
    rowsC.append(("in a practice with control covariate files", len(C)))
    log(f"\npractice restriction: T2D practices {pT0} -> {T['prac'].nunique()} | T2D people lost {nT0 - len(T):,}")

    # frequency matching 5:1 within band x sex x practice
    need = T.groupby(["band", "sex", "prac"]).size().rename("t2d")
    rng = np.random.default_rng(SEED)
    C = C.sort_values("key").reset_index(drop=True)
    picked = []
    for (b, s, p), g in C.groupby(["band", "sex", "prac"], sort=True):
        n = int(need.get((b, s, p), 0)) * RATIO
        if n == 0:
            continue
        idx = g.index.values
        picked.append(idx if len(idx) <= n else rng.choice(idx, n, replace=False))
    C = C.loc[np.concatenate(picked)].copy()
    rowsC.append(("frequency-matched sample (5:1, band x sex x practice)", len(C)))
    attrition(rowsT, "T2D ARM (thesis cohort)")
    attrition(rowsC, "CONTROLS (fresh 2009 pool)")
    got = C.groupby("band").size()
    tb = T.groupby("band").size()
    log(f"\n  achieved ratio overall {len(C) / len(T):.2f}:1")
    for b in sorted(tb.index):
        log(f"    age {b}-{b + 4}: T2D {tb[b]:>7,}  controls {int(got.get(b, 0)):>8,}  ratio {got.get(b, 0) / tb[b]:.2f}")
    short = (C.groupby(["band", "sex", "prac"]).size().reindex(need.index).fillna(0) < RATIO * need)
    log(f"  strata short of 5:1: {int(short.sum()):,} of {len(need):,} ({int(need[short].sum()):,} T2D people)")

    # ====================================================================== entry, exit, follow-up
    A = pd.concat([T, C], ignore_index=True)
    A["dob"] = pd.to_datetime(A["yob"].astype("Int64").astype(str) + "-07-01", errors="coerce")
    A["entry"] = pd.concat([pd.Series(INDEX, index=A.index), A["regdate"] + pd.DateOffset(months=9), A["amr"], A["acu"]], axis=1).max(axis=1)
    cap = pd.concat([A["death"], A["end"], A["dob"] + pd.DateOffset(years=99), pd.Series(END, index=A.index)], axis=1).min(axis=1)
    A["exit_sens"] = cap
    later_t2d = A["ds"].where((A["diab"] == 0) & (A["ds"] > A["entry"]))
    A["exit_main"] = pd.concat([cap, later_t2d], axis=1).min(axis=1)
    for v in ("main", "sens"):
        ok = A[f"exit_{v}"] > A["entry"]
        A[f"fu_{v}"] = ((A[f"exit_{v}"] - A["entry"]).dt.days / 365.25).where(ok, 0)
        A[f"died_{v}"] = (ok & A["death"].notna() & (A["death"] <= A[f"exit_{v}"]) & (A["death"] > A["entry"])).astype(int)
    A = A[A["fu_sens"] > 0]
    log(f"\ncontrols censored at a later T2D diagnosis (main analysis): {int(((A['diab'] == 0) & (A['exit_main'] < A['exit_sens'])).sum()):,}")
    reason = np.select([A["died_sens"] == 1, A["exit_sens"] == pd.Timestamp(END), A["exit_sens"] == A["dob"] + pd.DateOffset(years=99)],
                       ["died", "end of study", "99th birthday"], "transfer out / practice stopped")
    log("exit reasons (no censoring): " + "; ".join(f"{'T2D' if d == 1 else 'non-T2D'}: " +
        ", ".join(f"{k} {v:,}" for k, v in pd.Series(reason[A['diab'].values == d]).value_counts().items()) for d in (1, 0)))

    open(os.path.join(OUTD, "ch6_build_summary_part1.txt"), "w", encoding="utf-8").write("\n".join(lines))
    for _n in ("t2", "reg", "regs", "cc", "C", "T", "pool", "cc_ds", "diab_keys_before", "need"):
        globals().pop(_n, None)
    gc.collect()
if not os.path.exists(CHECKPOINT):
    A.to_pickle(CHECKPOINT)
    log(f"checkpoint saved: {CHECKPOINT}")
gc.collect()

# ====================================================================== covariates at entry (same code both arms)
log("\nCOVARIATES at entry")


STEP3A = os.path.join(PV, "step3a_parts")


def parts_fallback(prac, g):
    """T2D only: rebuild weight/BMI and smoking records for one practice from the PURE step3a parts (same AHD fields)."""
    mp = next((p for p in (os.path.join(STEP3A, f"{prac}_meas.pkl.gz"), os.path.join(STEP3A, f"{prac}_meas.pkl")) if os.path.exists(p)), None)
    sp = next((p for p in (os.path.join(STEP3A, f"{prac}_smoke.pkl.gz"), os.path.join(STEP3A, f"{prac}_smoke.pkl")) if os.path.exists(p)), None)
    if not mp and not sp:
        return None
    rows = []
    if mp:
        m = pd.read_pickle(mp)
        m = m[m["var"].astype(str) == "weight"]
        rows.append(pd.DataFrame({"patid": m["patid"].values, "ahdcode": "1005010200", "data1": "", "data3": m["v2"].values,
                                  "d_evntdate": pd.to_datetime(m["d"]).values}))
    if sp:
        sm = pd.read_pickle(sp)
        rows.append(pd.DataFrame({"patid": sm["patid"].values, "ahdcode": "1003040000", "data1": sm["status"].astype(str).values,
                                  "data3": np.nan, "d_evntdate": pd.to_datetime(sm["d"]).values}))
    return pd.concat(rows, ignore_index=True)


def ahd_covariates(folder, people):
    by_prac = {}
    for f in glob.glob(os.path.join(folder, "*.pkl")):
        by_prac.setdefault(prac_of(f), f)
    out, groups = [], people.groupby("prac")
    for i, (prac, g) in enumerate(groups):
        r = pd.DataFrame({"key": g["key"].values, "bmi": np.nan, "smoking": "missing"})
        f = by_prac.get(prac)
        a = None
        if f:
            try:
                a = pd.read_pickle(f)
            except MemoryError:
                gc.collect()
                a = None
                fb = parts_fallback(prac, g) if folder == T2D_AHD else None
                if fb is not None:
                    a = fb
                    log(f"    MEMORY: {os.path.basename(f)} too large - used the pre-parsed step3a parts for practice {prac} "
                        f"({a['patid'].nunique():,} of {len(g):,} people have records there)")
                else:
                    log(f"    MEMORY: could not load {os.path.basename(f)} - BMI/smoking left missing for {len(g):,} people in practice {prac}")
        if a is not None:
            a = a[["patid", "ahdcode", "data1", "data3", "d_evntdate"]]
            gc.collect()
            a["key"] = mkkey(a["patid"], pd.Series(prac, index=a.index))
            a["ahdcode"] = a["ahdcode"].astype(str).str.strip().str.replace(r"\.0$", "", regex=True)
            a = a[a["key"].isin(set(r["key"])) & a["ahdcode"].isin(["1005010200", "1003040000"])]
            a["d"] = pd.to_datetime(a["d_evntdate"], errors="coerce")
            e = g[["key", "entry"]]
            w = a[a["ahdcode"] == "1005010200"].assign(v=lambda z: pd.to_numeric(z["data3"], errors="coerce")).merge(e, on="key")
            w = w[w["v"].between(15, 60) & (w["d"] <= w["entry"]) & (w["d"] >= w["entry"] - pd.Timedelta(days=730))]
            r["bmi"] = r["key"].map(w.sort_values("d").groupby("key")["v"].last())
            s = a[a["ahdcode"] == "1003040000"].assign(st=lambda z: z["data1"].astype(str).str.strip().str.upper()).merge(e, on="key")
            s = s[s["d"] <= s["entry"]].sort_values("d")
            last = s.groupby("key")["st"].last().map({"Y": "current", "D": "ex", "N": "never"})
            ever = set(s.loc[s["st"].isin(["Y", "D"]), "key"])
            last = last.where(~((last == "never") & last.index.isin(ever)), "ex")
            r["smoking"] = r["key"].map(last).fillna("missing")
            del a, w, s, last
            gc.collect()
        out.append(r)
        if (i + 1) % 50 == 0:
            print(f"    {os.path.basename(folder)}: {i + 1}/{groups.ngroups} practices", flush=True)
    return pd.concat(out, ignore_index=True)


COVCHK = os.path.join(OUTD, "checkpoint_covariates.pkl")
if os.path.exists(COVCHK):
    del A
    gc.collect()
    A = pd.read_pickle(COVCHK)
    print(f"resumed covariates from {COVCHK} (BMI, smoking, Townsend, CVD done)", flush=True)
else:
    # ---- BMI/smoking: one short-lived subprocess per practice (memory freed after each); results cached per practice
    import subprocess, sys, tempfile
    HELPER = os.path.join(tempfile.gettempdir(), "phd_ch6_ahd_helper.py")
    import textwrap
    open(HELPER, "w", encoding="utf-8").write(textwrap.dedent(r"""
    import sys, re, gc
    import numpy as np, pandas as pd
    src, prac, people_p, out_p, mp, sp = sys.argv[1:7]
    g = pd.read_pickle(people_p)
    def npat(s):
        return pd.Series(s).astype(str).str.strip().str.replace(r"\.0$", "", regex=True).str.lstrip("0").replace("", "0")
    def key(pat):
        return npat(pat).values + "|" + prac
    r = pd.DataFrame({"key": g["key"].values, "bmi": np.nan, "smoking": "missing", "source": "none"})
    a = None
    if src != "-":
        try:
            a = pd.read_pickle(src)[["patid", "ahdcode", "data1", "data3", "d_evntdate"]]
            r["source"] = "raw"
        except MemoryError:
            a = None
            gc.collect()
    if a is None and (mp != "-" or sp != "-"):
        rows = []
        if mp != "-":
            m = pd.read_pickle(mp)
            m = m[m["var"].astype(str) == "weight"]
            rows.append(pd.DataFrame({"patid": m["patid"].values, "ahdcode": "1005010200", "data1": "", "data3": m["v2"].values, "d_evntdate": pd.to_datetime(m["d"]).values}))
            del m
        if sp != "-":
            sm = pd.read_pickle(sp)
            rows.append(pd.DataFrame({"patid": sm["patid"].values, "ahdcode": "1003040000", "data1": sm["status"].astype(str).values, "data3": np.nan, "d_evntdate": pd.to_datetime(sm["d"]).values}))
            del sm
        a = pd.concat(rows, ignore_index=True)
        r["source"] = "parts"
    if a is not None:
        a["key"] = key(a["patid"])
        a["ahdcode"] = a["ahdcode"].astype(str).str.strip().str.replace(r"\.0$", "", regex=True)
        a = a[a["key"].isin(set(r["key"])) & a["ahdcode"].isin(["1005010200", "1003040000"])].copy()
        a["d"] = pd.to_datetime(a["d_evntdate"], errors="coerce")
        e = g[["key", "entry"]]
        w = a[a["ahdcode"] == "1005010200"].assign(v=lambda z: pd.to_numeric(z["data3"], errors="coerce")).merge(e, on="key")
        w = w[w["v"].between(15, 60) & (w["d"] <= w["entry"]) & (w["d"] >= w["entry"] - pd.Timedelta(days=730))]
        r["bmi"] = r["key"].map(w.sort_values("d").groupby("key")["v"].last())
        s = a[a["ahdcode"] == "1003040000"].assign(st=lambda z: z["data1"].astype(str).str.strip().str.upper()).merge(e, on="key")
        s = s[s["d"] <= s["entry"]].sort_values("d")
        last = s.groupby("key")["st"].last().map({"Y": "current", "D": "ex", "N": "never"})
        ever = set(s.loc[s["st"].isin(["Y", "D"]), "key"])
        last = last.where(~((last == "never") & last.index.isin(ever)), "ex")
        r["smoking"] = r["key"].map(last).fillna("missing")
    r.to_pickle(out_p)
    """))
    COVDIR = os.path.join(OUTD, "cov_parts")
    os.makedirs(COVDIR, exist_ok=True)
    people_all = A[["key", "prac", "entry", "diab"]].copy()
    del A
    gc.collect()


    def run_arm(folder, arm, people):
        by_prac = {}
        for f in glob.glob(os.path.join(folder, "*.pkl")):
            by_prac.setdefault(prac_of(f), f)
        outs, groups = [], people.groupby("prac")
        for i, (prac, g) in enumerate(groups):
            out_p = os.path.join(COVDIR, f"{arm}_{prac}.pkl")
            if not os.path.exists(out_p):
                people_p = os.path.join(COVDIR, f"_people_{arm}_{prac}.pkl")
                g[["key", "entry"]].to_pickle(people_p)
                mp = sp = "-"
                if arm == "t2d":
                    mp = next((p for p in (os.path.join(STEP3A, f"{prac}_meas.pkl.gz"), os.path.join(STEP3A, f"{prac}_meas.pkl")) if os.path.exists(p)), "-")
                    sp = next((p for p in (os.path.join(STEP3A, f"{prac}_smoke.pkl.gz"), os.path.join(STEP3A, f"{prac}_smoke.pkl")) if os.path.exists(p)), "-")
                res = subprocess.run([sys.executable, HELPER, by_prac.get(prac, "-"), prac, people_p, out_p, mp, sp], capture_output=True, text=True)
                if res.returncode != 0 or not os.path.exists(out_p):
                    log(f"    {arm} practice {prac}: could not compute BMI/smoking ({len(g):,} people left missing) - {res.stderr.strip().splitlines()[-1][:120] if res.stderr.strip() else 'no message'}")
                    pd.DataFrame({"key": g["key"].values, "bmi": np.nan, "smoking": "missing", "source": "failed"}).to_pickle(out_p)
                try:
                    os.remove(people_p)
                except OSError:
                    pass
            outs.append(pd.read_pickle(out_p))
            if (i + 1) % 25 == 0:
                print(f"    {arm}: {i + 1}/{groups.ngroups} practices", flush=True)
        return pd.concat(outs, ignore_index=True)


    cv = pd.concat([run_arm(T2D_AHD, "t2d", people_all[people_all["diab"] == 1]), run_arm(ND_AHD, "nd", people_all[people_all["diab"] == 0])], ignore_index=True)
    log("  BMI/smoking source by arm: " + "; ".join(f"{('T2D' if d == 1 else 'non-T2D')}: " + ", ".join(f"{k} {v:,}" for k, v in cv[people_all.set_index('key').loc[cv['key'], 'diab'].values == d]['source'].value_counts().items()) for d in (1, 0)))
    cv = cv.drop(columns=["source"])
    del people_all
    gc.collect()
    A = pd.read_pickle(CHECKPOINT)
    A = A.merge(cv, on="key", how="left")

    # Townsend: controls from the register file; T2D from pvirecs_townsend (latest record on/before index, else earliest)
    tw = []
    for f in glob.glob(os.path.join(TOWNSEND_DIR, "*")):
        try:
            t = pd.read_pickle(f) if f.endswith((".pkl", ".gz")) else pd.read_stata(f)
        except Exception:
            continue
        t["key"] = mkkey(t["patid"], pd.Series(prac_of(f), index=t.index))
        t["d"] = pd.to_datetime(t.get("d_evntdate"), errors="coerce")
        tw.append(t[["key", "townsend", "d"]])
    if tw:
        tw = pd.concat(tw, ignore_index=True)
        tw["before"] = tw["d"] <= INDEX
        tw = tw.sort_values(["key", "before", "d"]).groupby("key").last()["townsend"]
        A["townsend"] = pd.to_numeric(A["townsend"], errors="coerce").astype(float)
        A.loc[A["diab"] == 1, "townsend"] = pd.to_numeric(A.loc[A["diab"] == 1, "key"].map(tw), errors="coerce")
    A["townsend"] = pd.to_numeric(A["townsend"], errors="coerce")

    # CVD history (same code lists)
    dd = pd.read_csv(DICT, dtype=str)
    dc = next((c for c in dd.columns if c.lower() in ("medcode", "readcode", "code")), dd.columns[0])
    desc = dict(zip(dd[dc].str.strip(), dd[next((c for c in dd.columns if "desc" in c.lower() or "term" in c.lower()), dd.columns[-1])].astype(str)))
    heart = set()
    for n, pth in HEART.items():
        if os.path.exists(pth):
            x = pd.read_csv(pth, dtype=str)
            col = next((c for c in x.columns if c.lower() in ("medcode", "readcode", "read_code", "code")), x.columns[0])
            codes = set(x[col].dropna().str.strip())
            heart |= {c for c in codes if any(w in desc.get(c, "").lower() for w in CABG_WORDS)} if n == "cabg_pci" else codes
    first = []
    for p in glob.glob(os.path.join(CODES2A, "*_codes.pkl*")):
        try:
            c = pd.read_pickle(p)
        except MemoryError:
            gc.collect()
            log(f"    MEMORY: skipped T2D code file {os.path.basename(p)}")
            continue
        c = c[c["medcode"].astype(str).str.strip().isin(heart)]
        if len(c):
            first.append(pd.Series(pd.to_datetime(c["first_d"]).values, index=mkkey(c["patid"], c["prac"])))
    for p in glob.glob(os.path.join(ND_CVD, "*.pkl")):
        try:
            c = pd.read_pickle(p)
        except MemoryError:
            gc.collect()
            log(f"    MEMORY: skipped CVD file {os.path.basename(p)}")
            continue
        c = c[c["medcode"].astype(str).str.strip().isin(heart)]
        if len(c):
            first.append(pd.Series(pd.to_datetime(c["d_evntdate"], errors="coerce").values, index=mkkey(c["patid"], pd.Series(prac_of(p), index=c.index))))
    first = pd.concat(first).groupby(level=0).min()
    A["cvd"] = (A["key"].map(first) <= A["entry"]).astype(int)

    A.to_pickle(COVCHK)
    log(f"covariate checkpoint saved: {COVCHK}")
for _n in ("first", "cv", "tw", "dd", "desc"):
    globals().pop(_n, None)
gc.collect()

# antihypertensives (same drug lists): T2D from the classified PURE scan, controls from the full raw files.
# Lean: each file is filtered at once to cohort members with a prescription in the 12 months before entry; only keys kept.
code2cls = {}
for cls, names in BP.items():
    for nm in names:
        for folder in (MED, MED_OLD):
            pth = os.path.join(folder, nm + ".csv")
            if os.path.exists(pth):
                x = pd.read_csv(pth, dtype=str)
                col = next((c for c in x.columns if c.lower().strip() in ("drugcode", "drug_code", "drugcodes", "code")), None)
                if col:
                    for c in x[col].dropna().astype(str).str.strip().str.replace(r"\.0$", "", regex=True):
                        code2cls.setdefault(c, set()).add(cls)
                break
entry_of = A.set_index("key")["entry"]
on_bp = set()


def flag(keys, dates):
    k = pd.Series(keys)
    d = pd.to_datetime(pd.Series(dates), errors="coerce")
    e = k.map(entry_of)
    ok = e.notna() & (d <= e) & (d > e - pd.Timedelta(days=365))
    on_bp.update(k[ok.values].tolist())


t2d_keys = set(A.loc[A["diab"] == 1, "key"])
for i, p in enumerate(glob.glob(os.path.join(RX4A, "*_rx.pkl.gz"))):
    try:
        r = pd.read_pickle(p)
        r = r[r["cls"].isin(list(BP))]
        prac = os.path.basename(p).split("_rx")[0]
        k = mkkey(r["patid"], pd.Series(prac, index=r.index))
        m = pd.Series(k).isin(t2d_keys).values
        flag(k[m], r["d"].values[m])
        del r
    except MemoryError:
        log(f"    MEMORY: skipped T2D prescription part {os.path.basename(p)}")
    gc.collect()
    if (i + 1) % 100 == 0:
        print(f"    T2D prescription parts {i + 1}", flush=True)
ctrl_keys = set(A.loc[A["diab"] == 0, "key"])
for i, p in enumerate(glob.glob(os.path.join(ND_RX, "*.pkl"))):
    try:
        t = pd.read_pickle(p)[["patid", "drugcode", "d_prscdate"]]
        t = t[t["drugcode"].astype(str).str.strip().str.replace(r"\.0$", "", regex=True).isin(code2cls)]
        k = mkkey(t["patid"], pd.Series(prac_of(p), index=t.index))
        m = pd.Series(k).isin(ctrl_keys).values
        flag(k[m], t["d_prscdate"].values[m])
        del t
    except MemoryError:
        log(f"    MEMORY: skipped antihypertensive file {os.path.basename(p)} (controls in practice {prac_of(p)} counted as not on antihypertensives)")
    gc.collect()
    if (i + 1) % 100 == 0:
        print(f"    control antihypertensive files {i + 1}", flush=True)
A["antihtn"] = A["key"].isin(on_bp).astype(int)
for d, lab in ((1, "T2D"), (0, "non-T2D")):
    g = A[A["diab"] == d]
    log(f"  {lab}: BMI recorded {100 * g['bmi'].notna().mean():.1f}% | smoking {g['smoking'].value_counts(normalize=True).round(3).to_dict()} | "
        f"Townsend 1-5 {100 * g['townsend'].between(1, 5).mean():.1f}% | CVD {100 * g['cvd'].mean():.1f}% | antihypertensive {100 * g['antihtn'].mean():.1f}%")

# ====================================================================== C6.4 table
log("\nC6.4 SUMMARY (crude rates per 1,000 person-years, exact Poisson 95% CI)")
for v, lab in (("main", "MAIN (controls censored at T2D diagnosis)"), ("sens", "SENSITIVITY (no censoring)")):
    log(f"\n  {lab}")
    log(f"  {'arm':<8} {'N':>9} {'deaths':>8} {'person-years':>13} {'rate (95% CI)':>24} {'median FU (IQR), years':>26}")
    for arm in (1, 0):
        g = A[A["diab"] == arm]
        d, py = int(g[f"died_{v}"].sum()), g[f"fu_{v}"].sum()
        lo = chi2.ppf(0.025, 2 * d) / 2 / py * 1000 if d > 0 else 0
        hi = chi2.ppf(0.975, 2 * (d + 1)) / 2 / py * 1000
        q = g[f"fu_{v}"].quantile([0.25, 0.5, 0.75])
        log(f"  {'T2D' if arm == 1 else 'non-T2D':<8} {len(g):>9,} {d:>8,} {py:>13,.0f} {1000 * d / py:>9.1f} ({lo:.1f}-{hi:.1f}) {q[0.5]:>12.2f} ({q[0.25]:.2f}-{q[0.75]:.2f})")

# ====================================================================== analysis file
A["ageband"] = ((INDEX.year - A["yob"] - 0.5) // 5 * 5).astype(int)
A["ageband_cat"] = ((A["ageband"] - 35) // 5 + 1).astype(int)
A["townsend_cat"] = A["townsend"].where(A["townsend"].between(1, 5), 6).fillna(6).astype(int)
A["smoking_cat"] = A["smoking"].map({"never": 1, "ex": 2, "current": 3}).fillna(4).astype(int)
A["bmi_cat"] = pd.cut(A["bmi"], [0, 18.5, 25, 30, 35, 40, 100], right=False, labels=[1, 2, 3, 4, 5, 6]).astype(float).fillna(7).astype(int)
out = A[["diab", "prac", "sex", "ageband", "ageband_cat", "died_main", "fu_main", "died_sens", "fu_sens",
         "townsend_cat", "smoking_cat", "bmi_cat", "cvd", "antihtn"]].copy()
out = out[out["sex"].isin([1, 2])]
out["prac"] = out["prac"].astype(str)
out.to_stata(os.path.join(OUTD, "ch6_analysis.dta"), write_index=False, version=118)
log(f"\nsaved {os.path.join(OUTD, 'ch6_analysis.dta')} ({len(out):,} rows; no patient ids)")
n_t, n_c = int((out["diab"] == 1).sum()), int((out["diab"] == 0).sum())
log("\nCOHORT DEFINITION (for the thesis)")
log(f"Mortality in people with T2D was compared with a comparison group without diabetes, both defined at a common index date of "
    f"1 January 2009 and followed until 31 December 2018. The T2D group comprised the thesis cohort ({n_t:,} people after "
    f"restriction to practices with complete comparison-group data). For each practice, sex and five-year age band, up to five "
    f"people without diabetes were sampled at random from those eligible at the index date ({n_c:,} people; achieved ratio "
    f"{n_c / n_t:.2f}:1); eligibility for both groups required age 35-89 years, registration for at least nine months at a practice "
    f"with acceptable mortality recording and computer usage dates, being alive and registered at the index date, and (for the "
    f"comparison group) no diagnosis of T2D on or before the index date. Follow-up started at the latest of the index date, nine months "
    f"after registration and the practice data-quality dates, and ended at the earliest of death, transfer out, the practice's last "
    f"data collection, the 99th birthday or 31 December 2018; in the main analysis, comparison-group members diagnosed with T2D during "
    f"follow-up were censored at diagnosis, and a sensitivity analysis retained them. Townsend quintile, smoking status, BMI, "
    f"cardiovascular history and antihypertensive use were measured at the start of follow-up in the same way for both groups "
    f"(most recent record in the preceding 24 months; 12 months for prescriptions).")
txt = "\n".join(lines)
open(os.path.join(OUTD, "ch6_build_summary.txt"), "w", encoding="utf-8").write(txt)
os.makedirs(DESK, exist_ok=True)
open(os.path.join(DESK, "ch6_2009_build_summary.txt"), "w", encoding="utf-8").write(txt)
print("\nDone. Next: run phd_ch6_2009_models.do in Stata.")
