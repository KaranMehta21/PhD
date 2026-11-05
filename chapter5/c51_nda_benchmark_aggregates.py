"""
PhD C5.1 - descriptive aggregates of the thesis T2D cohort (t2d_data.dta) for the comparison with the National Diabetes Audit.
No new analysis: counts, percentages, means, medians. Output: Desktop\\pure_inventory\\c51_nda_aggregates.txt (aggregates only).
Run:  python "OneDrive - University College London\\phd_c51_nda_aggregates.py"     (~1 min)
"""
import os
import numpy as np
import pandas as pd

DATA = r"S:\FPHS_THIN_Projects\Karan\BZI_PhD_T2D\PhD_descriptive_data\t2d_data.dta"
OUT = os.path.join(os.path.expanduser("~"), "Desktop", "pure_inventory", "c51_nda_aggregates.txt")
D = pd.read_stata(DATA, convert_categoricals=False)
L = []
n = len(D)
L.append(f"N = {n:,}")
sex = pd.to_numeric(D["sex"], errors="coerce")
L.append(f"sex: " + ", ".join(f"{int(k)} = {v:,} ({100 * v / n:.1f}%)" for k, v in sex.value_counts().sort_index().items()) + "  (1 = men expected)")
age = pd.to_numeric(D["age"], errors="coerce")
q = age.quantile([0.25, 0.5, 0.75])
L.append(f"age: mean {age.mean():.1f} (SD {age.std():.1f}); median {q[0.5]:.1f} (IQR {q[0.25]:.1f}-{q[0.75]:.1f}); range {age.min():.1f}-{age.max():.1f}")
bands = pd.cut(age, [0, 50, 60, 70, 80, 200], right=False, labels=["<50", "50-59", "60-69", "70-79", ">=80"])
L.append("age bands %: " + " / ".join(f"{k} {100 * v / n:.1f}" for k, v in bands.value_counts().reindex(bands.cat.categories).items()))
dur = pd.to_numeric(D["duration_diabetes"], errors="coerce")
qd = dur.quantile([0.25, 0.5, 0.75])
L.append(f"diabetes duration (years): mean {dur.mean():.1f} (SD {dur.std():.1f}); median {qd[0.5]:.1f} (IQR {qd[0.25]:.1f}-{qd[0.75]:.1f})")
h = pd.to_numeric(D["hba1c"], errors="coerce")
hr = h.dropna()
qh = hr.quantile([0.25, 0.5, 0.75])
L.append(f"HbA1c recorded: {len(hr):,} ({100 * len(hr) / n:.1f}%); mean {hr.mean():.1f} (SD {hr.std():.1f}); median {qh[0.5]:.1f} (IQR {qh[0.25]:.1f}-{qh[0.75]:.1f}) mmol/mol")
L.append("HbA1c among recorded, %: <=48 {:.1f} | <=53 {:.1f} | <=58 {:.1f} | 59-75 {:.1f} | >75 {:.1f} | >86 {:.1f}".format(
    100 * (hr <= 48).mean(), 100 * (hr <= 53).mean(), 100 * (hr <= 58).mean(), 100 * ((hr > 58) & (hr <= 75)).mean(),
    100 * (hr > 75).mean(), 100 * (hr > 86).mean()))
eth = [c for c in D.columns if "eth" in c.lower() or "race" in c.lower()]
L.append(f"ethnicity-like columns in the extract: {eth if eth else 'none'}")
txt = "\n".join(L)
print(txt)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
open(OUT, "w", encoding="utf-8").write(txt)
