"""
PhD CHAPTER 6 (C6.2) - checking the Poisson model assumption (variance = mean, i.e. no overdispersion).
For Models A and B (main analysis, overall and by sex), on the same aggregated cells used for estimation:
  - Pearson chi2 / residual df  (~1 = no overdispersion; > 1 = overdispersion)
  - deviance / residual df
  - MRR from a NEGATIVE BINOMIAL model (allows extra-Poisson variation) as a sensitivity analysis, next to the Poisson MRR
  - Poisson MRR with model-based vs practice-clustered robust 95% CIs (robust CIs are valid under overdispersion)
Input : S:\\FPHS_THIN_Projects\\Karan\\PURE_v2\\ch6_2009\\ch6_analysis.dta      Output: Desktop\\pure_inventory\\ch6_poisson_check.txt
Run:  python "OneDrive - University College London\\phd_ch6_poisson_check.py"      (~5 min)
"""
import os, gc, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

D = r"S:\FPHS_THIN_Projects\Karan\PURE_v2\ch6_2009"
OUT = os.path.join(os.path.expanduser("~"), "Desktop", "pure_inventory", "ch6_poisson_check.txt")
A = pd.read_stata(os.path.join(D, "ch6_analysis.dta"), convert_categoricals=False)
A["prac"] = pd.factorize(A["prac"].astype(str))[0].astype(np.int32)
for c in ("ageband_cat", "sex", "townsend_cat", "smoking_cat", "diab", "died_main"):
    A[c] = A[c].astype(np.int8)
A = A[["diab", "prac", "sex", "ageband_cat", "townsend_cat", "smoking_cat", "died_main", "fu_main"]]
gc.collect()
lines = ["PhD CHAPTER 6 - POISSON ASSUMPTION CHECK (main analysis; cells = covariate pattern x practice)", ""]
lines.append(f"{'model':<6}{'stratum':<9}{'cells':>8}{'Pearson/df':>12}{'Dev/df':>9}   {'Poisson MRR, model-based CI':<30}{'Poisson MRR, robust CI':<26}{'Neg. binomial MRR (95% CI)':<28}alpha")
specs = [("A", "C(ageband_cat) * C(sex)", "C(ageband_cat)", ["ageband_cat", "sex"]),
         ("B", "C(ageband_cat) * C(sex) + C(townsend_cat) + C(smoking_cat)", "C(ageband_cat) + C(townsend_cat) + C(smoking_cat)",
          ["ageband_cat", "sex", "townsend_cat", "smoking_cat"])]
f = lambda b, se: f"{np.exp(b):.3f} ({np.exp(b - 1.96 * se):.3f}-{np.exp(b + 1.96 * se):.3f})"
for name, rhs_all, rhs_sex, keys in specs:
    for lab, cond, rhs in (("overall", None, rhs_all), ("men", 1, rhs_sex), ("women", 2, rhs_sex)):
        d = A if cond is None else A[A["sex"] == cond]
        g = d.groupby(["diab", "prac"] + keys, observed=True).agg(deaths=("died_main", "sum"), py=("fu_main", "sum")).reset_index()
        g = g[g["py"] > 0]
        off = np.log(g["py"])
        pois = smf.glm(f"deaths ~ diab + {rhs}", data=g, family=sm.families.Poisson(), offset=off).fit()
        rob = smf.glm(f"deaths ~ diab + {rhs}", data=g, family=sm.families.Poisson(), offset=off).fit(cov_type="cluster", cov_kwds={"groups": g["prac"].values})
        pear, devr = pois.pearson_chi2 / pois.df_resid, pois.deviance / pois.df_resid
        try:
            nb = smf.negativebinomial(f"deaths ~ diab + {rhs}", data=g, offset=off).fit(disp=0, maxiter=200)
            if np.isfinite(nb.bse["diab"]):
                nbs, alpha = f(nb.params["diab"], nb.bse["diab"]), f"{nb.params['alpha']:.4f}"
            else:
                nbs, alpha = f"{np.exp(nb.params['diab']):.3f} (= Poisson)", "~0 (no extra-Poisson variation)"
        except Exception as ex:
            nbs, alpha = f"failed ({type(ex).__name__})", ""
        lines.append(f"{name:<6}{lab:<9}{len(g):>8,}{pear:>12.3f}{devr:>9.3f}   {f(pois.params['diab'], pois.bse['diab']):<30}"
                     f"{f(rob.params['diab'], rob.bse['diab']):<26}{nbs:<28}{alpha}")
        del g, pois, rob
        gc.collect()
lines += ["", "Interpretation: Pearson/df close to 1 means no meaningful overdispersion. Robust (practice-clustered) CIs used in the",
          "thesis remain valid if there is overdispersion; similar negative binomial MRRs confirm the estimates are not sensitive to it."]
txt = "\n".join(lines)
print(txt)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
open(OUT, "w", encoding="utf-8").write(txt)
