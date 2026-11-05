"""
PhD correction C7.3 - exits from the Chapter 7 T2D cohort (t2d_data.dta, 155,399; index 1 Jan 2009) by reason:
died / reached 31 Dec 2018 / 99th birthday / transferred out / practice stopped contributing.
Dates come from the registration files (BZI_PhD_T2D\\patrecs_marital) linked via 'patprac'. Exit = earliest of death,
transfer out, practice last collection, 99th birthday, 31 Dec 2018 (ties: death > transfer > practice > birthday > end).
Check: recomputed follow-up and deaths are compared with t2d_data's own 'exposure' and 'mortality_status'.
~5-10 min. Output: Desktop\\pure_inventory\\c73_exits.txt
Run:  python "OneDrive - University College London\\phd_c73_exits.py"
"""
import os, glob
import numpy as np
import pandas as pd

K = r"S:\FPHS_THIN_Projects\Karan"
T2D_DATA = os.path.join(K, "BZI_PhD_T2D", "PhD_descriptive_data", "t2d_data.dta")
PATREC_DIR = os.path.join(K, "BZI_PhD_T2D", "patrecs_marital")
OUT = os.path.join(os.path.expanduser("~"), "Desktop", "pure_inventory", "c73_exits.txt")
INDEX, END = pd.Timestamp("2009-01-01"), pd.Timestamp("2018-12-31")


def as_date(s):
    s = pd.Series(s)
    if pd.api.types.is_numeric_dtype(s) and s.dropna().between(-40000, 40000).mean() > 0.9:
        return pd.to_datetime("1960-01-01") + pd.to_timedelta(s, unit="D")
    return pd.to_datetime(s.astype(str).str.strip().replace({"": None, "nan": None}), errors="coerce")


def npat(s):
    return pd.Series(s).astype(str).str.strip().str.replace(r"\.0$", "", regex=True).str.lstrip("0").replace("", "0")


t2 = pd.read_stata(T2D_DATA, convert_categoricals=False)
pp, pr = t2["patprac"].astype(str).str.strip(), t2["prac"].astype(str).str.strip()
pat = [p[:-len(r)] if p.lower().endswith(r.lower()) else (p[len(r):] if p.lower().startswith(r.lower()) else p) for p, r in zip(pp, pr)]
t2["key"] = npat(pat).values + "|" + pr.str.lower().values
reg = pd.concat([pd.read_pickle(f) for f in sorted(glob.glob(os.path.join(PATREC_DIR, "*.pkl")))], ignore_index=True)
reg["key"] = npat(reg["patid"]).values + "|" + reg["prac"].astype(str).str.strip().str.lower().values
d = t2[["key", "exposure", "mortality_status"]].merge(reg.drop_duplicates("key"), on="key", how="left")
lines = [f"C7.3 EXITS - Chapter 7 T2D cohort (t2d_data.dta, n={len(t2):,}); linked to registration records: "
         f"{int(d['d_regdate'].notna().sum()):,} ({100 * d['d_regdate'].notna().mean():.1f}%)"]
dob = pd.to_datetime(pd.to_numeric(d["d_yob"].astype(str).str[:4], errors="coerce").astype("Int64").astype(str) + "-07-01", errors="coerce")
cand = pd.DataFrame({"died": as_date(d["d_deathdate"]), "transferred out": as_date(d["d_xferdate"]),
                     "practice stopped contributing": as_date(d["lastdate"]), "99th birthday": dob + pd.DateOffset(years=99),
                     "reached 31 Dec 2018": pd.Series(END, index=d.index)})
exit_date = cand.min(axis=1)
reason = pd.Series("reached 31 Dec 2018", index=d.index)
for r in ["reached 31 Dec 2018", "99th birthday", "practice stopped contributing", "transferred out", "died"]:   # last wins = priority
    reason = reason.where(~(cand[r] == exit_date), r)
fu = ((exit_date - INDEX).dt.days / 365.25).clip(lower=0)
lines.append("\nExit reason            n         %")
for r in ["died", "reached 31 Dec 2018", "99th birthday", "transferred out", "practice stopped contributing"]:
    n = int((reason == r).sum())
    lines.append(f"  {r:<30} {n:>8,} {100 * n / len(d):>6.1f}%")
lines.append(f"\nCheck vs thesis file: deaths {int((reason == 'died').sum()):,} recomputed vs {int(d['mortality_status'].sum()):,} in t2d_data | "
             f"person-years {fu.sum():,.0f} vs {d['exposure'].sum():,.0f} | follow-up within 0.1 y of the file for "
             f"{100 * ((fu - d['exposure']).abs() < 0.1).mean():.1f}% of people")
aligned = reason.where(d["mortality_status"] != 1, "died")
aligned = aligned.where(~((d["mortality_status"] != 1) & (aligned == "died")), "transferred out")
lines.append("\nTABLE FOR THE THESIS (deaths = thesis mortality_status; other exits by the rule above)")
lines.append("Exit reason            n         %")
for r in ["died", "reached 31 Dec 2018", "99th birthday", "transferred out", "practice stopped contributing"]:
    n = int((aligned == r).sum())
    lines.append(f"  {r:<30} {n:>8,} {100 * n / len(d):>6.1f}%")
moved = (d["mortality_status"] == 1) & (reason != "died")
lines.append(f"  (thesis deaths whose recorded exit was not a dated death: {int(moved.sum()):,} - "
             + ", ".join(f"{k} {v:,}" for k, v in reason[moved].value_counts().items()) + ")")
lines.append("If the check differs materially, the thesis used a different exit rule (e.g. 1 Jan 2019 or a different "
             "transfer date); report the reasons using the thesis rule.")
txt = "\n".join(lines)
print(txt)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
open(OUT, "w", encoding="utf-8").write(txt)
