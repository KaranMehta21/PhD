import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2
rates_csv, res_csv, out = sys.argv[1], sys.argv[2], sys.argv[3]
BLUE, RED = "#4285F4", "#EA4335"
AGES = ["35 - 39","40 - 44","45 - 49","50 - 54","55 - 59","60 - 64","65 - 69","70 - 74","75 - 79","80 - 84","85 - 89"]
T = pd.read_csv(rates_csv)
T = T[T["age"] != "Overall"].copy()
T["age"] = T["age"].str.replace("-", " - ", regex=False)
def rate_ci(d, py):
    d = np.asarray(d, float); py = np.asarray(py, float)
    r = 1000 * d / py
    lo = np.where(d > 0, chi2.ppf(0.025, 2 * d) / 2, 0) / py * 1000
    hi = chi2.ppf(0.975, 2 * (d + 1)) / 2 / py * 1000
    return r, lo, hi
def style(ax, ylabel):
    ax.set_xlabel("Age, years", fontsize=11); ax.set_ylabel(ylabel, fontsize=11)
    ax.grid(axis="y", color="#d9d9d9", linewidth=0.8); ax.set_axisbelow(True)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color("#666666"); ax.tick_params(axis="y", length=0)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=2, frameon=False, fontsize=11)
    ax.set_ylim(bottom=0)
def line(ax, x, r, lo, hi, color, label, band=True):
    ax.plot(x, r, color=color, linewidth=2, label=label)
    if band: ax.fill_between(range(len(x)), lo, hi, color=color, alpha=0.15, linewidth=0)
figs = []
for fid, title, series, ylab in [
    (5, "men", [("Men", "T2D", BLUE, "T2D Male"), ("Men", "non-T2D", RED, "Non-T2D Male")], "Mortality rate (per 1000 pyar)"),
    (6, "women", [("Women", "T2D", BLUE, "T2D Female"), ("Women", "non-T2D", RED, "Non-T2D Female")], "Mortality rate (per 1000 pyar)"),
    (7, "T2D men vs women", [("Men", "T2D", BLUE, "T2D Male"), ("Women", "T2D", RED, "T2D Female")], "Mortality rate (per 1000 pyar)")]:
    fig, ax = plt.subplots(figsize=(12, 7.42), dpi=100)
    for sex, arm, col, lab in series:
        x = T[T["sex"] == sex].set_index("age").reindex(AGES)
        r, lo, hi = rate_ci(x[f"{arm} deaths"], x[f"{arm} pyar"])
        line(ax, AGES, r, lo, hi, col, lab)
    style(ax, ylab); fig.tight_layout()
    p = f"{out}/Figure{fid}_ch6.png"; fig.savefig(p, dpi=100); plt.close(fig); figs.append(p)
R = pd.read_csv(res_csv)
fig, ax = plt.subplots(figsize=(12, 7.42), dpi=100)
for k, (sex, col, lab) in enumerate((("men", BLUE, "Men T2D vs Non T2D MRR"), ("women", RED, "Women T2D vs Non T2D MRR"))):
    rr = []
    for a in AGES:
        st = f"{sex} {a.replace(' - ', '-')}"
        r = R[(R["model"] == "A") & (R["analysis"] == "main") & (R["stratum"] == st)]
        rr.append(r[["irr", "lo", "hi"]].iloc[0].values if len(r) else (np.nan,) * 3)
    rr = np.array(rr, float); xs = np.arange(len(AGES)) + (k - 0.5) * 0.12
    ax.plot(xs, rr[:, 0], color=col, linewidth=2, label=lab, marker="o", markersize=4)
    ax.errorbar(xs, rr[:, 0], yerr=[rr[:, 0] - rr[:, 1], rr[:, 2] - rr[:, 0]], fmt="none", ecolor=col, alpha=0.6, capsize=3)
ax.axhline(1, color="#888888", linewidth=1, linestyle="--")
ax.set_xticks(range(len(AGES))); ax.set_xticklabels(AGES)
style(ax, "Mortality rate ratio"); fig.tight_layout()
p = f"{out}/Figure8_ch6.png"; fig.savefig(p, dpi=100); plt.close(fig); figs.append(p)
print("\n".join(figs))
