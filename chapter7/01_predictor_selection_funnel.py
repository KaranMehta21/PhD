"""
PhD C7.2 - RECONSTRUCT THE PREDICTOR-SELECTION FUNNEL on the thesis data (t2d_data.dta, complete-case per variable, as in the
thesis notebooks), and compare with the final model's 58 predictors (read from PhD_cph_scikitsurv_240324.ipynb).
  1. candidate predictors (all columns except identifiers/outcome)
  2. Kaplan-Meier / log-rank screen: binary -> two groups; smoking -> its categories; continuous -> quartiles + missing as a
     separate group (as in 'Kaplan Meier Curves.ipynb'); multivariate log-rank test; threshold p < 0.005
  3. univariable Cox (lifelines; continuous linear per SD, categorical as dummies; complete cases for that variable); p < 0.005
  4. Pearson correlation among continuous candidates (pairwise complete): every pair with |r| > 0.8
  5. comparison with the final model
Output: Desktop\\pure_inventory\\c72_selection.txt (+ c72_selection_table.csv). ~15-25 min.
Run AFTER the refit has finished:  python "OneDrive - University College London\\phd_c72_selection.py"
"""
import os, re, json, warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")
from lifelines import CoxPHFitter
from lifelines.statistics import multivariate_logrank_test

DATA = r"S:\FPHS_THIN_Projects\Karan\BZI_PhD_T2D\PhD_descriptive_data\t2d_data.dta"
NB = os.path.join(os.path.expanduser("~"), "OneDrive - University College London", "PhD_cph_scikitsurv_240324.ipynb")
OUTD = os.path.join(os.path.expanduser("~"), "Desktop", "pure_inventory")
EXCLUDE = {"prac", "patprac", "patid", "mortality_status", "exposure"}
THR = 0.005
L = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    L.append(s)


D = pd.read_stata(DATA, convert_categoricals=False)
T, E = D["exposure"].astype(float), D["mortality_status"].astype(int)
cands = [c for c in D.columns if c not in EXCLUDE]
log(f"t2d_data: {len(D):,} people, {int(E.sum()):,} deaths | candidate predictors: {len(cands)}")
log("  " + ", ".join(cands))


def kind(c):
    v = pd.to_numeric(D[c], errors="coerce").dropna()
    if c == "smoking_status" or (v.nunique() <= 5 and set(v.unique()) - {0, 1}):
        return "categorical"
    if v.nunique() <= 2:
        return "binary"
    return "continuous"


rows = []
for i, c in enumerate(cands, 1):
    x = pd.to_numeric(D[c], errors="coerce")
    k = kind(c)
    if k == "continuous":
        grp = pd.qcut(x.rank(method="first"), 4, labels=["Q1", "Q2", "Q3", "Q4"]).astype(object)
        grp = grp.where(x.notna(), "missing")
    else:
        grp = x.astype(object).where(x.notna(), "missing")
    try:
        lr = multivariate_logrank_test(T, grp, E).p_value
    except Exception:
        lr = np.nan
    try:
        m = x.notna()
        if k == "continuous":
            df = pd.DataFrame({"z": (x[m] - x[m].mean()) / x[m].std(), "T": T[m], "E": E[m]})
        else:
            dm = pd.get_dummies(x[m].round().astype(int).astype(str), prefix="c", drop_first=True, dtype=float)
            df = pd.concat([dm, T[m].rename("T"), E[m].rename("E")], axis=1)
        cph = CoxPHFitter().fit(df, "T", "E")
        p_cox = float(cph.summary["p"].min())
        hr = float(np.exp(cph.params_.iloc[0]))
    except Exception:
        p_cox, hr = np.nan, np.nan
    rows.append(dict(variable=c, type=k, n_complete=int(x.notna().sum()), logrank_p=lr, cox_p=p_cox, cox_hr_first_term=hr,
                     pass_logrank=bool(np.isfinite(lr) and lr < THR), pass_cox=bool(np.isfinite(p_cox) and p_cox < THR)))
    if i % 10 == 0:
        print(f"  {i}/{len(cands)} variables screened", flush=True)
R = pd.DataFrame(rows)

cont = R.loc[R["type"] == "continuous", "variable"].tolist()
C = D[cont].apply(pd.to_numeric, errors="coerce").corr(method="pearson")
pairs = [(a, b, C.loc[a, b]) for i, a in enumerate(cont) for b in cont[i + 1:] if abs(C.loc[a, b]) > 0.8]

final = []
try:
    nb = json.load(open(NB, encoding="utf-8"))
    for cell in nb["cells"]:
        mm = re.search(r"X_train\s*=\s*\w+\[\[(.*?)\]\]", "".join(cell.get("source", "")), re.S)
        if mm:
            final = [re.sub(r"^_1_", "", v) for v in re.findall(r"['\"]([^'\"]+)['\"]", mm.group(1))]
            break
except Exception as ex:
    log(f"(final model list not read: {type(ex).__name__})")
final_set = {("smoking_status" if f == "smoking_status2" else f) for f in final}
R["in_final_model"] = R["variable"].isin(final_set)

s1 = R[R["pass_logrank"]]
s2 = R[R["pass_logrank"] & R["pass_cox"]]
log(f"\nFUNNEL (threshold p < {THR})")
log(f"  1. candidates                                     {len(R)}")
log(f"  2. pass KM/log-rank                               {len(s1)}   failed: {', '.join(R.loc[~R['pass_logrank'], 'variable']) or '-'}")
log(f"  3. also pass univariable Cox                      {len(s2)}   failed Cox: {', '.join(s1.loc[~s1['pass_cox'], 'variable']) or '-'}")
log(f"  4. correlated pairs |r| > 0.8 (continuous)        {len(pairs)}")
for a, b, r in sorted(pairs, key=lambda t: -abs(t[2])):
    log(f"       {a} - {b}: r = {r:.3f}   (in final model: {a in final_set} / {b in final_set})")
log(f"  5. final model predictors (notebook)              {len(final)}")
log(f"     passed both screens but NOT in the final model: {', '.join(s2.loc[~s2['in_final_model'], 'variable']) or '-'}")
log(f"     in the final model but failed a screen:         {', '.join(R.loc[R['in_final_model'] & ~(R['pass_logrank'] & R['pass_cox']), 'variable']) or '-'}")
log("\nPER-VARIABLE RESULTS")
log(R.assign(logrank_p=R["logrank_p"].map(lambda v: f"{v:.2g}"), cox_p=R["cox_p"].map(lambda v: f"{v:.2g}"),
             cox_hr_first_term=R["cox_hr_first_term"].round(3)).to_string(index=False))
os.makedirs(OUTD, exist_ok=True)
R.to_csv(os.path.join(OUTD, "c72_selection_table.csv"), index=False)
open(os.path.join(OUTD, "c72_selection.txt"), "w", encoding="utf-8").write("\n".join(L))
print(f"\nSaved {os.path.join(OUTD, 'c72_selection.txt')}")
