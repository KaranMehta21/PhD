"""
PhD CHAPTER 7 - PREDICT FOR ONE PERSON with the four saved models (run phd_c7_refit.py first; it creates ch7_refit\\models).
You give RAW values (age, sex, HbA1c, BMI, labs, drugs, conditions ...). The script builds everything the models need:
  - the mvrs spline terms (from the saved spline maps, learned from observed data and verified in the refit log)
  - the 13 interaction terms (raw products, as defined in phd_updatedmice_180224.do)
  - standardisation with each model's saved means/SDs
and returns, for every model: survival S(t), probability of death 1 - S(t) and cumulative hazard H(t) at all times 0-10 years,
plus a summary at 1, 5 and 9 years, a CSV of the curves and a plot.
Usage:
  python phd_c7_predict.py --template            writes case_template.json with every required input (edit the values)
  python phd_c7_predict.py case.json             predictions for that person -> case_predictions.csv, case_survival.png
Inputs use the variable names of the thesis data without the imputation prefix (e.g. hba1c, bmi, egfr; sex coded as in the data;
smoking_status2 = 1 ever-smoker, 0 never; drugs and hx_* conditions 0/1). All inputs are required.
Values outside the range seen in the data are clamped to that range for the spline terms (a warning is printed).
"""
import os, sys, json, pickle, warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")

MODD = os.path.join(r"S:\FPHS_THIN_Projects\Karan", "PURE_v2", "ch7_refit", "models")
TIMES = np.round(np.arange(0.0, 10.0001, 0.05), 2)


def load(name):
    p = os.path.join(MODD, name)
    if not os.path.exists(p):
        return None
    with open(p, "rb") as fh:
        return pickle.load(fh)


SM = load("spline_maps.pkl")
if SM is None:
    sys.exit(f"spline_maps.pkl not found in {MODD} - run phd_c7_refit.py first")
MAPS, INTER = SM["maps"], SM["interactions"]
SCALE = SM.get("interaction_scales", {})
TERM2BASE = {t: (b, k) for b, m in MAPS.items() for k, t in enumerate(m["terms"])}


def raw_name(col):
    return col[3:] if col.startswith("_1_") else col


def required_inputs():
    req = set()
    for nm in ("cph_original.pkl", "cph_nonlinear_final.pkl", "rsf.pkl", "deepsurv_meta.pkl"):
        b = load(nm)
        if b is None:
            continue
        for c in b["columns"]:
            if c in TERM2BASE:
                req.add(MAPS[TERM2BASE[c][0]]["raw"])
            elif raw_name(c) in INTER:
                req |= set(INTER[raw_name(c)])
            else:
                req.add(raw_name(c))
    return sorted(req)


def features(case, bundle):
    row, warn = {}, []
    for c in bundle["columns"]:
        if c in TERM2BASE:
            b, k = TERM2BASE[c]
            m = MAPS[b]
            x = float(case[m["raw"]])
            lo, hi = m["range"]
            if x < lo or x > hi:
                warn.append(f"{m['raw']}={x} outside data range {lo:.3g}-{hi:.3g}; clamped")
                x = min(max(x, lo), hi)
            row[c] = float(np.interp(x, m["x"], m["y"][:, k]))
        elif raw_name(c) in INTER:
            a, b = INTER[raw_name(c)]
            row[c] = SCALE.get(raw_name(c), 1.0) * float(case[a]) * float(case[b])
        else:
            row[c] = float(case[raw_name(c)])
    X = pd.DataFrame([row])[bundle["columns"]]
    if "continuous" in bundle:
        cont = bundle["continuous"]
        X[cont] = (X[cont] - bundle["mean"][cont]) / bundle["sd"][cont]
    return X, warn


def predict(case):
    out, warns = {}, set()
    for label, nm in (("CPH original", "cph_original.pkl"), ("CPH non-linear + interactions (final, all data)", "cph_nonlinear_final.pkl")):
        b = load(nm)
        if b is None:
            continue
        X, w = features(case, b)
        warns |= set(w)
        S = b["model"].predict_survival_function(X, times=TIMES).values.ravel()
        out[label] = S
    b = load("rsf.pkl")
    if b is not None:
        X, w = features(case, b)
        warns |= set(w)
        f = b["model"].predict_survival_function(X.values)[0]
        out["RSF"] = np.clip(f(np.clip(TIMES, f.x[0], f.x[-1])), 0, 1)
    b = load("deepsurv_meta.pkl")
    if b is not None and os.path.exists(os.path.join(MODD, "deepsurv_net.pt")):
        try:
            import torch, torchtuples as tt
            X, w = features(case, b)
            warns |= set(w)
            net = tt.practical.MLPVanilla(len(b["columns"]), b["layers"], 1, batch_norm=b["batch_norm"], dropout=b["dropout"], output_bias=False)
            net.load_state_dict(torch.load(os.path.join(MODD, "deepsurv_net.pt"), map_location="cpu"))
            net.eval()
            with torch.no_grad():
                lr = float(net(torch.tensor(X.values.astype("float32"))).numpy().ravel()[0])
            H0 = b["baseline_cumulative_hazards"]
            h0 = np.interp(TIMES, H0.index.values, H0.values, left=0.0)
            out["DeepSurv"] = np.exp(-h0 * np.exp(lr))
        except Exception as ex:
            print(f"DeepSurv not available: {type(ex).__name__}: {ex}")
    for w in sorted(warns):
        print("WARNING:", w)
    return out


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--template":
        tmpl = {k: 0 for k in required_inputs()}
        tmpl.update({k: v for k, v in {"age": 65, "sex": SM["sex_values"][0] if SM["sex_values"] else 1, "duration_diabetes": 8,
                                        "hba1c": 58, "bmi": 31}.items() if k in tmpl})
        json.dump(tmpl, open("case_template.json", "w"), indent=1)
        print(f"wrote case_template.json with {len(tmpl)} inputs (sex values in the data: {SM['sex_values']}); edit it, then run with it")
        sys.exit()
    case = json.load(open(sys.argv[1]))
    missing = [k for k in required_inputs() if k not in case]
    if missing:
        sys.exit(f"missing inputs: {missing}")
    res = predict(case)
    tab = pd.DataFrame({"time_years": TIMES})
    for k, S in res.items():
        tab[f"{k} | survival"] = S
        tab[f"{k} | death probability"] = 1 - S
        tab[f"{k} | cumulative hazard"] = -np.log(np.clip(S, 1e-12, 1))
    stem = os.path.splitext(sys.argv[1])[0]
    tab.to_csv(f"{stem}_predictions.csv", index=False)
    print("\nProbability of death by 1 / 5 / 9 years:")
    for k, S in res.items():
        p = [1 - float(np.interp(t, TIMES, S)) for t in (1, 5, 9)]
        print(f"  {k:<50} {p[0]:.3f} / {p[1]:.3f} / {p[2]:.3f}")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 5), dpi=120)
        for k, S in res.items():
            ax.plot(TIMES, S, label=k, linewidth=2)
        ax.set_xlabel("Years"); ax.set_ylabel("Survival probability"); ax.set_ylim(0, 1.02); ax.grid(alpha=0.3); ax.legend()
        fig.tight_layout(); fig.savefig(f"{stem}_survival.png")
        print(f"\nsaved {stem}_predictions.csv and {stem}_survival.png")
    except Exception:
        print(f"\nsaved {stem}_predictions.csv")
