"""
Kaplan-Meier curves and log-rank tests for candidate predictors (Supplementary Methods 7A1, 7A2)
Original Jupyter notebook exported as a script (code cells only; outputs removed). Source last modified 2023-10-24.
"""

# %% cell 1
import pandas as pd
import numpy as np
import glob

# %% cell 2
import pysurvival
import numpy as np
from pysurvival.utils.metrics import concordance_index
from pysurvival.models.semi_parametric import NonLinearCoxPHModel
from sklearn.model_selection import train_test_split
from pysurvival.utils import save_model
from pysurvival.utils import load_model

# %% cell 3
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt

# %% cell 4
df = pd.read_csv('S:\FPHS_THIN_Projects\Karan\BZI_PhD_T2D\PhD_descriptive_data/desc_prac_data.csv')

# %% cell 5
df

# %% cell 6
df.columns

# %% cell 7
df[['patprac', 'exposure']]

# %% cell 8
#Convert years into days
df['exposure'] = df['exposure']*365.2422

# %% cell 9
df['exposure']

# %% cell 10
#Histogram of survival time
df['exposure'].hist(bins=50)

# %% cell 11
#Kaplan Meier curve for full cohort - lifelines
kmf = KaplanMeierFitter()

# %% cell 12
T = df["exposure"]
E = df["mortality_status"]

# %% cell 13
kmf.fit(durations = T, event_observed = E)

# %% cell 14
kmf.plot_survival_function()

# %% cell 15
kmf.plot_cumulative_density()

# %% cell 16
kmf.median_survival_time_

# %% cell 17
from lifelines.utils import median_survival_times
median_ = kmf.median_survival_time_
median_confidence_interval_ = median_survival_times(kmf.confidence_interval_)
print(median_)
print(median_confidence_interval_)

# %% cell 18
#Kaplan-Meier curve for full cohort - pysurvival
from pysurvival.models.non_parametric import KaplanMeierModel
from pysurvival.utils.display import display_non_parametric

# %% cell 19
km_model = KaplanMeierModel()
km1_model = KaplanMeierModel()

# %% cell 20
km_model.fit(T, E, alpha=0.95)

# %% cell 21
display_non_parametric(km_model)

# %% cell 22
#Survival of different gender groups
ax = plt.subplot(111)
m = (df["sex"] == 1)
kmf.fit(durations = T[m], event_observed = E[m], label = "Male")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~m], event_observed = E[~m], label = "Female")
kmf.plot_survival_function(ax = ax, at_risk_counts = True)
plt.title("Survival of different gender group")
results=logrank_test(T[m],T[~m], event_observed_A = E[m], event_observed_B = E[~m])
results.print_summary()

# %% cell 23
df.columns

# %% cell 24
#Survival of T2D on acei/arb versus not on it (hypertension + proteinuria)
ax = plt.subplot(111)
ace_arb = (df['acei/arb'] == 1)
kmf.fit(durations = T[ace_arb], event_observed = E[ace_arb], label = "on acei/arb")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~ace_arb], event_observed = E[~ace_arb], label = "not on acei/arb")
kmf.plot_survival_function(ax = ax, at_risk_counts = True)
plt.title("Survival of T2D on acei/arb versus not on it")
results=logrank_test(T[ace_arb],T[~ace_arb], event_observed_A = E[ace_arb], event_observed_B = E[~ace_arb])
results.print_summary()

# %% cell 25
#Survival of T2D on calcium channel blocker versus not on it (more severe hypertension: probably second line agent)
ax = plt.subplot(111)
ca_blocker = (df['ca_blocker'] == 1)
kmf.fit(durations = T[ca_blocker], event_observed = E[ca_blocker], label = "on calcium channel blocker")
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~ca_blocker], event_observed = E[~ca_blocker], label = "not on caclcium channel blocker")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D on calcium channel blocker versus not on it")
results=logrank_test(T[ca_blocker],T[~ca_blocker], event_observed_A = E[ca_blocker], event_observed_B = E[~ca_blocker])
results.print_summary()

# %% cell 26
#Survival of T2D on beta blocker versus not on it (severe hypertension, ?IHD)
ax = plt.subplot(111)
beta_blocker = (df['beta_blocker'] == 1)
kmf.fit(durations = T[beta_blocker], event_observed = E[beta_blocker], label = "on beta channel blocker")
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~beta_blocker], event_observed = E[~beta_blocker], label = "not on beta blocker")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D on beta blocker versus not on it")
results=logrank_test(T[beta_blocker],T[~beta_blocker], event_observed_A = E[beta_blocker], event_observed_B = E[~beta_blocker])
results.print_summary()

# %% cell 27
#Survival of T2D on thiazides versus not on it
ax = plt.subplot(111)
thiazides = (df['thiazides'] == 1)
kmf.fit(durations = T[thiazides], event_observed = E[thiazides], label = "on thiazides blocker")
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~thiazides], event_observed = E[~thiazides], label = "not on thiazides")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D on thiazides versus not on it")
results=logrank_test(T[thiazides],T[~thiazides], event_observed_A = E[thiazides], event_observed_B = E[~thiazides])
results.print_summary()

# %% cell 28
#Survival of T2D on alpha_blocker versus not on it
ax = plt.subplot(111)
alpha_blocker = (df['alpha_blocker'] == 1)
kmf.fit(durations = T[alpha_blocker], event_observed = E[alpha_blocker], label = "on alpha blocker")
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~alpha_blocker], event_observed = E[~alpha_blocker], label = "not on alpha blocker")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D on alpha blocker versus not on it")
results=logrank_test(T[alpha_blocker],T[~alpha_blocker], event_observed_A = E[alpha_blocker], event_observed_B = E[~alpha_blocker])
results.print_summary()

# %% cell 29
#Survival of T2D on metformin versus not on it
ax = plt.subplot(111)
metformin = (df['metformin'] == 1)
kmf.fit(durations = T[metformin], event_observed = E[metformin], label = "on metformin")
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~metformin], event_observed = E[~metformin], label = "not on metformin")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D on metformin versus not on it")
results=logrank_test(T[metformin],T[~metformin], event_observed_A = E[metformin], event_observed_B = E[~metformin])
results.print_summary()

# %% cell 30
#Survival of T2D on sulphonylureas versus not on it
ax = plt.subplot(111)
sulphonylureas = (df['sulphonylureas'] == 1)
kmf.fit(durations = T[sulphonylureas], event_observed = E[sulphonylureas], label = "on sulphonylureas")
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~sulphonylureas], event_observed = E[~sulphonylureas], label = "not on sulphonylureas")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D on sulphonylureas versus not on it")
results=logrank_test(T[sulphonylureas],T[~sulphonylureas], event_observed_A = E[sulphonylureas], event_observed_B = E[~sulphonylureas])
results.print_summary()

# %% cell 31
#Survival of T2D on meglitinides versus not on it
ax = plt.subplot(111)
meglitinides = (df['meglitinides'] == 1)
kmf.fit(durations = T[meglitinides], event_observed = E[meglitinides], label = "on meglitinides")
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~meglitinides], event_observed = E[~meglitinides], label = "not on meglitinides")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D on meglitinides versus not on it")
results=logrank_test(T[meglitinides],T[~meglitinides], event_observed_A = E[meglitinides], event_observed_B = E[~meglitinides])
results.print_summary()

# %% cell 32
#Survival of T2D on gliptins versus not on it
ax = plt.subplot(111)
gliptins = (df['gliptins'] == 1)
kmf.fit(durations = T[gliptins], event_observed = E[gliptins], label = "on gliptins")
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~gliptins], event_observed = E[~gliptins], label = "not on gliptins")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D on gliptins versus not on it")
results=logrank_test(T[gliptins],T[~gliptins], event_observed_A = E[gliptins], event_observed_B = E[~gliptins])
results.print_summary()

# %% cell 33
#Survival of T2D on thiazols versus not on it
ax = plt.subplot(111)
thiazols = (df['thiazols'] == 1)
kmf.fit(durations = T[thiazols], event_observed = E[thiazols], label = "on thiazols")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~thiazols], event_observed = E[~thiazols], label = "not on thiazols")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D on thiazols versus not on it")
results=logrank_test(T[thiazols],T[~thiazols], event_observed_A = E[thiazols], event_observed_B = E[~thiazols])
results.print_summary()

# %% cell 34
#Survival of T2D on glp1 versus not on it
ax = plt.subplot(111)
glp1 = (df['glp1'] == 1)
kmf.fit(durations = T[glp1], event_observed = E[glp1], label = "on glp1")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~glp1], event_observed = E[~glp1], label = "not on glp1")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D on glp1 versus not on it")
results=logrank_test(T[glp1],T[~glp1], event_observed_A = E[glp1], event_observed_B = E[~glp1])
results.print_summary()

# %% cell 35
#Survival of T2D on glp1 versus not on it
ax = plt.subplot(111)
acarbose = (df['acarbose'] == 1)
kmf.fit(durations = T[acarbose], event_observed = E[acarbose], label = "on acarbose")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~acarbose], event_observed = E[~acarbose], label = "not on acarbose")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D on acarbose versus not on it")
results=logrank_test(T[acarbose],T[~acarbose], event_observed_A = E[acarbose], event_observed_B = E[~acarbose])
results.print_summary()

# %% cell 36
#Survival of T2D on insulin versus not on it
ax = plt.subplot(111)
insulin = (df['insulin'] == 1)
kmf.fit(durations = T[insulin], event_observed = E[insulin], label = "on insulin")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~insulin], event_observed = E[~insulin], label = "not on insulin")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D on insulin versus not on it")
results=logrank_test(T[insulin],T[~insulin], event_observed_A = E[insulin], event_observed_B = E[~insulin])
results.print_summary()

# %% cell 37
#Survival of T2D on statins versus not on it
ax = plt.subplot(111)
statins = (df['statins'] == 1)
kmf.fit(durations = T[statins], event_observed = E[statins], label = "on statins")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~statins], event_observed = E[~statins], label = "not on statins")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D on statins versus not on it")
results=logrank_test(T[statins],T[~statins], event_observed_A = E[statins], event_observed_B = E[~statins])
results.print_summary()

# %% cell 38
#Survival of T2D with history of repiratory versus no history
ax = plt.subplot(111)
hx_resp = (df['hx_resp'] == 1)
kmf.fit(durations = T[hx_resp], event_observed = E[hx_resp], label = "history of respiratory disease")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~hx_resp], event_observed = E[~hx_resp], label = "no history of respiratory disease")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D with history of repiratory disease versus no history")
results=logrank_test(T[hx_resp],T[~hx_resp], event_observed_A = E[hx_resp], event_observed_B = E[~hx_resp])
results.print_summary()

# %% cell 39
#Survival of T2D with history of liver versus no history
ax = plt.subplot(111)
hx_liver = (df['hx_liver'] == 1)
kmf.fit(durations = T[hx_liver], event_observed = E[hx_liver], label = "history of liver disease")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~hx_liver], event_observed = E[~hx_liver], label = "no history of liver disease")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D with history of liver disease versus no history")
results=logrank_test(T[hx_liver],T[~hx_liver], event_observed_A = E[hx_liver], event_observed_B = E[~hx_liver])
results.print_summary()

# %% cell 40
#Survival of T2D with history of gi versus no history
ax = plt.subplot(111)
hx_gi = (df['hx_gi'] == 1)
kmf.fit(durations = T[hx_gi], event_observed = E[hx_gi], label = "history of gi disease")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~hx_gi], event_observed = E[~hx_gi], label = "no history of gi disease")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D with history of gi disease versus no history")
results=logrank_test(T[hx_gi],T[~hx_gi], event_observed_A = E[hx_gi], event_observed_B = E[~hx_gi])
results.print_summary()

# %% cell 41
#Survival of T2D with history of neuro versus no history
ax = plt.subplot(111)
hx_neuro = (df['hx_neuro'] == 1)
kmf.fit(durations = T[hx_neuro], event_observed = E[hx_neuro], label = "history of neuro disease")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~hx_neuro], event_observed = E[~hx_neuro], label = "no history of neuro disease")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D with history of neuro disease versus no history")
results=logrank_test(T[hx_neuro],T[~hx_neuro], event_observed_A = E[hx_neuro], event_observed_B = E[~hx_neuro])
results.print_summary()

# %% cell 42
#Survival of T2D with history of autoimmune versus no history
ax = plt.subplot(111)
hx_autoimmune = (df['hx_autoimmune'] == 1)
kmf.fit(durations = T[hx_autoimmune], event_observed = E[hx_autoimmune], label = "history of neuro disease")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~hx_autoimmune], event_observed = E[~hx_autoimmune], label = "no history of autoimmune disease")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D with history of autoimmune disease versus no history")
results=logrank_test(T[hx_autoimmune],T[~hx_autoimmune], event_observed_A = E[hx_autoimmune], event_observed_B = E[~hx_autoimmune])
results.print_summary()

# %% cell 43
#Survival of T2D with history of severe kidney disease versus no history
ax = plt.subplot(111)
hx_severekidney = (df['hx_severekidney'] == 1)
kmf.fit(durations = T[hx_severekidney], event_observed = E[hx_severekidney], label = "history of severe kidney disease")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~hx_severekidney], event_observed = E[~hx_severekidney], label = "no history of severe kidney disease")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D with history of severe kidney disease versus no history")
results=logrank_test(T[hx_severekidney],T[~hx_severekidney], event_observed_A = E[hx_severekidney], event_observed_B = E[~hx_severekidney])
results.print_summary()

# %% cell 44
#Survival of T2D with history of cancer versus no history
ax = plt.subplot(111)
hx_cancer = (df['hx_cancer'] == 1)
kmf.fit(durations = T[hx_cancer], event_observed = E[hx_cancer], label = "history of cancer disease")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~hx_cancer], event_observed = E[~hx_cancer], label = "no history of cancer disease")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D with history of cancer versus no history")
results=logrank_test(T[hx_cancer],T[~hx_cancer], event_observed_A = E[hx_cancer], event_observed_B = E[~hx_cancer])
results.print_summary()

# %% cell 45
#Survival of T2D with history of coronary artery disease versus no history
ax = plt.subplot(111)
hx_cad = (df['hx_cad'] == 1)
kmf.fit(durations = T[hx_cad], event_observed = E[hx_cad], label = "history of coronary artery disease")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~hx_cad], event_observed = E[~hx_cad], label = "no history of coronary artery disease")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D with history of coronary artery disease versus no history")
results=logrank_test(T[hx_cad],T[~hx_cad], event_observed_A = E[hx_cad], event_observed_B = E[~hx_cad])
results.print_summary()

# %% cell 46
#Survival of T2D with history of coronary artery disease versus no history
ax = plt.subplot(111)
hx_cerebro = (df['hx_cerebro'] == 1)
kmf.fit(durations = T[hx_cerebro], event_observed = E[hx_cerebro], label = "history of cerebrovascular disease")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~hx_cerebro], event_observed = E[~hx_cerebro], label = "no history of cerebrovascular disease")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D with history of cerebrovscular disease versus no history")
results=logrank_test(T[hx_cerebro],T[~hx_cerebro], event_observed_A = E[hx_cerebro], event_observed_B = E[~hx_cerebro])
results.print_summary()

# %% cell 47
#Survival of T2D with history of peripheral vascular disease versus no history
ax = plt.subplot(111)
hx_pvd = (df['hx_pvd'] == 1)
kmf.fit(durations = T[hx_pvd], event_observed = E[hx_pvd], label = "history of peripheral vascular disease")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~hx_pvd], event_observed = E[~hx_pvd], label = "no history of peripheral vascular disease")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D with history of peripheral vascular disease versus no history")
results=logrank_test(T[hx_pvd],T[~hx_pvd], event_observed_A = E[hx_pvd], event_observed_B = E[~hx_pvd])
results.print_summary()

# %% cell 48
#Survival of T2D with history of arrhythmia versus no history
ax = plt.subplot(111)
hx_arrhythmia = (df['hx_arrhythmia'] == 1)
kmf.fit(durations = T[hx_arrhythmia], event_observed = E[hx_arrhythmia], label = "history of arrhythmia")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~hx_arrhythmia], event_observed = E[~hx_arrhythmia], label = "no history of arrhythmia")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D with history of arrhythmia versus no history")
results=logrank_test(T[hx_arrhythmia],T[~hx_arrhythmia], event_observed_A = E[hx_arrhythmia], event_observed_B = E[~hx_arrhythmia])
results.print_summary()

# %% cell 49
#Survival of T2D with history of heart failure versus no history
ax = plt.subplot(111)
hx_heartfailure = (df['hx_heartfailure'] == 1)
kmf.fit(durations = T[hx_heartfailure], event_observed = E[hx_heartfailure], label = "history of heart failure")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~hx_heartfailure], event_observed = E[~hx_heartfailure], label = "no history of heart failure")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D with history of heart failure versus no history")
results=logrank_test(T[hx_heartfailure],T[~hx_heartfailure], event_observed_A = E[hx_heartfailure], event_observed_B = E[~hx_heartfailure])
results.print_summary()

# %% cell 50
#Survival of T2D with history of MI versus no history
ax = plt.subplot(111)
hx_mi = (df['hx_mi'] == 1)
kmf.fit(durations = T[hx_mi], event_observed = E[hx_mi], label = "history of MI")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~hx_mi], event_observed = E[~hx_mi], label = "no history of MI")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D with history of MI versus no history")
results=logrank_test(T[hx_mi],T[~hx_mi], event_observed_A = E[hx_mi], event_observed_B = E[~hx_mi])
results.print_summary()

# %% cell 51
#Survival of T2D with history of DVT versus no history
ax = plt.subplot(111)
hx_dvt = (df['hx_dvt'] == 1)
kmf.fit(durations = T[hx_dvt], event_observed = E[hx_dvt], label = "history of DVT")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~hx_dvt], event_observed = E[~hx_dvt], label = "no history of DVT")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D with history of DVT versus no history")
results=logrank_test(T[hx_dvt],T[~hx_dvt], event_observed_A = E[hx_dvt], event_observed_B = E[~hx_dvt])
results.print_summary()

# %% cell 52
#Survival of T2D with history of coma/shock versus no history
ax = plt.subplot(111)
hx_comashock = (df['hx_comashock'] == 1)
kmf.fit(durations = T[hx_comashock], event_observed = E[hx_comashock], label = "history of coma/shock")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~hx_comashock], event_observed = E[~hx_comashock], label = "no history of coma/shock")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D with history of coma/shock versus no history")
results=logrank_test(T[hx_comashock],T[~hx_comashock], event_observed_A = E[hx_comashock], event_observed_B = E[~hx_comashock])
results.print_summary()

# %% cell 53
#Survival of T2D with history of footleg ulcer versus no history
ax = plt.subplot(111)
hx_footlegulcer = (df['hx_footlegulcer'] == 1)
kmf.fit(durations = T[hx_footlegulcer], event_observed = E[hx_footlegulcer], label = "hx_footlegulcer")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~hx_footlegulcer], event_observed = E[~hx_footlegulcer], label = "no history of footlegulcer")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D with history of footleg ulcer versus no history")
results=logrank_test(T[hx_footlegulcer],T[~hx_footlegulcer], event_observed_A = E[hx_footlegulcer], event_observed_B = E[~hx_footlegulcer])
results.print_summary()

# %% cell 54
#Survival of T2D with history of protein in urine versus no history
ax = plt.subplot(111)
hx_proteinuria = (df['hx_proteinuria'] == 1)
kmf.fit(durations = T[hx_proteinuria], event_observed = E[hx_proteinuria], label = "histroy of protein in urine")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~hx_proteinuria], event_observed = E[~hx_proteinuria], label = "no history of protein in urine")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D with history of protein in urine versus no history")
results=logrank_test(T[hx_proteinuria],T[~hx_proteinuria], event_observed_A = E[hx_proteinuria], event_observed_B = E[~hx_proteinuria])
results.print_summary()

# %% cell 55
#Survival of T2D with history of retinopathy in urine versus no history
ax = plt.subplot(111)
hx_retinopathy = (df['hx_retinopathy'] == 1)
kmf.fit(durations = T[hx_retinopathy], event_observed = E[hx_retinopathy], label = "history of retinopathy")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~hx_retinopathy], event_observed = E[~hx_retinopathy], label = "no history of retinopathy")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D with history of retinopathy versus no history")
results=logrank_test(T[hx_retinopathy],T[~hx_retinopathy], event_observed_A = E[hx_retinopathy], event_observed_B = E[~hx_retinopathy])
results.print_summary()

# %% cell 56
#Survival of T2D with history of nephropathy in urine versus no history
ax = plt.subplot(111)
hx_nephropathy = (df['hx_nephropathy'] == 1)
kmf.fit(durations = T[hx_nephropathy], event_observed = E[hx_nephropathy], label = "history of nephropathy")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~hx_nephropathy], event_observed = E[~hx_nephropathy], label = "no history of nephropathy")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D with history of nephropathy versus no history")
results=logrank_test(T[hx_nephropathy],T[~hx_nephropathy], event_observed_A = E[hx_nephropathy], event_observed_B = E[~hx_nephropathy])
results.print_summary()

# %% cell 57
#Survival of T2D with history of nephropathy in urine versus no history
ax = plt.subplot(111)
hx_neuropathy = (df['hx_neuropathy'] == 1)
kmf.fit(durations = T[hx_neuropathy], event_observed = E[hx_neuropathy], label = "history of neuropathy")
plt.ylim(0, 1.0)
kmf.plot_survival_function(ax = ax)
kmf.fit(T[~hx_neuropathy], event_observed = E[~hx_neuropathy], label = "no history of neuropathy")
kmf.plot_survival_function(ax = ax)
plt.title("Survival of T2D with history of neuropathy versus no history")
results=logrank_test(T[hx_neuropathy],T[~hx_neuropathy], event_observed_A = E[hx_neuropathy], event_observed_B = E[~hx_neuropathy])
results.print_summary()

# %% cell 60
df.columns

# %% cell 61
df10_cont = df[['patprac', 'Duration_diabetes', 'BMI', 'Diastolic Blood Pressure',
       'Systolic Blood Pressure', 'Albumin','ALP', 'ALT', 'bilirubin',
       'total cholesterol', 'serum creatinine', 'Hb', 'HDL Cholesterol',
       'Sodium', 'LDL-C', 'Total protein', 'Triglycerides', 'MCV', 'TSH',
       'Urea', 'Neutrophil count', 'WCC', 'Lymphocyte count', 'Potassium',
       'Hba1c', 'eGFR', 'pulse_pressure',
       'mean arterial pressure', 'totalchol/hdl_ratio',
       'Neutrophil/lymphocyte', 'mortality_status', 'exposure']]

# %% cell 64
print(df10_cont['Hba1c'].describe())
df1 = df10_cont[['Hba1c', 'exposure', 'mortality_status']]
print(sum(df1['Hba1c'].isna()))
df1['Hba1c'] = df1['Hba1c'].fillna(0)
print(sum(df1['Hba1c'] == 0))

# %% cell 65
df1['hba1c_bin']=pd.cut(x = df1['Hba1c'],
                        bins = [-1,9,41.99,48,58,75, 192.896850], 
                        labels = [0,1,2,3,4,5])

# %% cell 66
df1[df1['hba1c_bin'].isna()]

# %% cell 67
df1['hba1c_bin'].unique()

# %% cell 68
result = multivariate_logrank_test(df1['exposure'], df1['hba1c_bin'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 69
# Define the categories
one_hba1c = (df1.hba1c_bin == 1) 
two_hba1c = (df1.hba1c_bin == 2)
three_hba1c = (df1.hba1c_bin == 3) 
four_hba1c = (df1.hba1c_bin == 4)
five_hba1c = (df1.hba1c_bin == 5)
missing_hba1c = (df1.hba1c_bin == 0)
print(missing_hba1c.value_counts())

# %% cell 70
# Plot the survival function
ax = plt.subplot()
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[one_hba1c].exposure,
event_observed=df1[one_hba1c].mortality_status, label='Low')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[two_hba1c].exposure, event_observed=df1[two_hba1c].mortality_status, label='Excellent control')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[three_hba1c].exposure, event_observed=df1[three_hba1c].mortality_status, label='Good control')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[four_hba1c].exposure, event_observed=df1[four_hba1c].mortality_status, label='High')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[five_hba1c].exposure, event_observed=df1[five_hba1c].mortality_status, label='Poor control')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing_hba1c].exposure, event_observed=df1[missing_hba1c].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.title('Survival Function based on HbA1c category')
plt.show()

# %% cell 71
df1 = df10_cont[['Hba1c', 'exposure', 'mortality_status']]

# %% cell 72
bin_labels = [1,2,3,4]
df1['Hba1c'] = pd.qcut(df1['Hba1c'], q=4, labels=bin_labels)

# %% cell 73
# Define the categories
Low = (df1.Hba1c == 1) 
Medium = (df1.Hba1c == 2)
High = (df1.Hba1c == 3) 
Very_high = (df1.Hba1c == 4)
missing = (df1.Hba1c.isna())

# %% cell 74
sum(df1['Hba1c'].isna())

# %% cell 75
missing.value_counts()

# %% cell 76
df1 = df10_cont[['Hba1c', 'exposure', 'mortality_status']]

# %% cell 77
sum(df1['Hba1c'].isna())

# %% cell 78
df1

# %% cell 79
from lifelines.statistics import multivariate_logrank_test

# %% cell 80
result = multivariate_logrank_test(df1['exposure'], df1['hba1c_bin'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 81
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on HbA1c category')
plt.show()

# %% cell 83
print(df10_cont['Diastolic Blood Pressure'].describe())
df1 = df10_cont[['Diastolic Blood Pressure', 'exposure', 'mortality_status']]
print(sum(df1['Diastolic Blood Pressure'].isna()))
df1['Diastolic Blood Pressure'] = df1['Diastolic Blood Pressure'].fillna(0)
print(sum(df1['Diastolic Blood Pressure'] == 0))

# %% cell 84
df1['dbp_bin']=pd.cut(x = df1['Diastolic Blood Pressure'],
                        bins = [-1,9,64.99,79.99,89.99,119.99,147], 
                        labels = [0,1,2,3,4,5])

# %% cell 85
df1['dbp_bin'].unique()

# %% cell 86
result = multivariate_logrank_test(df1['exposure'], df1['dbp_bin'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 87
# Define the categories
one_dbp = (df1.dbp_bin == 1) 
two_dbp = (df1.dbp_bin == 2)
three_dbp = (df1.dbp_bin == 3) 
four_dbp = (df1.dbp_bin == 4)
five_dbp = (df1.dbp_bin == 5)
missing_dbp = (df1.dbp_bin == 0)
print(missing_dbp.value_counts())

# %% cell 88
# Plot the survival function
ax = plt.subplot()
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[one_dbp].exposure,
event_observed=df1[one_dbp].mortality_status, label='Low')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[two_dbp].exposure, event_observed=df1[two_dbp].mortality_status, label='Good')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[three_dbp].exposure, event_observed=df1[three_dbp].mortality_status, label='Borderline')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[four_dbp].exposure, event_observed=df1[four_dbp].mortality_status, label='High')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[five_dbp].exposure, event_observed=df1[five_dbp].mortality_status, label='Malignant')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing_dbp].exposure, event_observed=df1[missing_dbp].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.title('Survival Function based on Diastolic Blood Pressure category')
plt.show()

# %% cell 89
df1 = df10_cont[['Diastolic Blood Pressure', 'exposure', 'mortality_status']]

# %% cell 90
bin_labels = [1,2,3,4]
df1['Diastolic Blood Pressure'] = pd.qcut(df1['Diastolic Blood Pressure'], q=4, labels=bin_labels)

# %% cell 91
result = multivariate_logrank_test(df1['exposure'], df1['Diastolic Blood Pressure'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 92
df1.columns

# %% cell 93
# Define the categories
Low = (df1['Diastolic Blood Pressure'] == 1) 
Medium = (df1['Diastolic Blood Pressure'] == 2)
High = (df1['Diastolic Blood Pressure'] == 3) 
Very_high = (df1['Diastolic Blood Pressure'] == 4)
missing = (df1['Diastolic Blood Pressure'].isna())
print(sum(df1['Diastolic Blood Pressure'].isna()))
print(missing.value_counts())
print(sum(df1['Diastolic Blood Pressure'].isna()))

# %% cell 94
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on Diastolic Blood Pressure category')
plt.show()

# %% cell 96
print(df10_cont['Systolic Blood Pressure'].describe())
df1 = df10_cont[['Systolic Blood Pressure', 'exposure', 'mortality_status']]
print(sum(df1['Systolic Blood Pressure'].isna()))
df1['Systolic Blood Pressure'] = df1['Systolic Blood Pressure'].fillna(0)
print(sum(df1['Systolic Blood Pressure'] == 0))

# %% cell 97
df1['sbp_bin']=pd.cut(x = df1['Systolic Blood Pressure'],
                        bins = [-1,79.99,89.99,129.99,139.99,159.99, 178.99, 217], 
                        labels = [0,1,2,3,4,5,6])

# %% cell 98
result = multivariate_logrank_test(df1['exposure'], df1['sbp_bin'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 99
# Define the categories
one_sbp = (df1.sbp_bin == 1) 
two_sbp = (df1.sbp_bin == 2)
three_sbp = (df1.sbp_bin == 3) 
four_sbp = (df1.sbp_bin == 4)
five_sbp = (df1.sbp_bin == 5)
six_sbp = (df1.sbp_bin == 6)
missing_sbp = (df1.sbp_bin == 0)
print(missing_sbp.value_counts())

# %% cell 101
# Plot the survival function
ax = plt.subplot()
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[one_sbp].exposure,
event_observed=df1[one_sbp].mortality_status, label='Low')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[two_sbp].exposure, event_observed=df1[two_sbp].mortality_status, label='Normal')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[three_sbp].exposure, event_observed=df1[three_sbp].mortality_status, label='Borderline')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[four_sbp].exposure, event_observed=df1[four_sbp].mortality_status, label='High')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[five_sbp].exposure, event_observed=df1[five_sbp].mortality_status, label='Very High')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[six_sbp].exposure, event_observed=df1[six_sbp].mortality_status, label='Malignant')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing_sbp].exposure, event_observed=df1[missing_sbp].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.title('Survival Function based on Systolic Blood Pressure category')
plt.show()

# %% cell 102
df1 = df10_cont[['Systolic Blood Pressure', 'exposure', 'mortality_status']]

# %% cell 103
bin_labels = [1,2,3,4]
df1['Systolic Blood Pressure'] = pd.qcut(df1['Systolic Blood Pressure'], q=4, labels=bin_labels)

# %% cell 104
result = multivariate_logrank_test(df1['exposure'], df1['Systolic Blood Pressure'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 105
# Define the categories
Low = (df1['Systolic Blood Pressure'] == 1) 
Medium = (df1['Systolic Blood Pressure'] == 2)
High = (df1['Systolic Blood Pressure'] == 3) 
Very_high = (df1['Systolic Blood Pressure'] == 4)
missing = (df1['Systolic Blood Pressure'].isna())
print(sum(df1['Systolic Blood Pressure'].isna()))
print(missing.value_counts())
print(sum(df1['Systolic Blood Pressure'].isna()))

# %% cell 106
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on Systolic Blood Pressure category')
plt.show()

# %% cell 108
- eGFR:
    - <15 : CKD 5
    - 15 - 29 : CKD 4
    - 29 - 44 : CKD 3b
    - 45 - 59 : CKD 3a
    - 60 - 89 : CKD 2
    - >89 : normal

# %% cell 109
df10_cont.columns

# %% cell 110
print(df10_cont['eGFR'].describe())
df1 = df10_cont[['eGFR', 'exposure', 'mortality_status']]
print(sum(df1['eGFR'].isna()))
df1['eGFR'] = df1['eGFR'].fillna(999)
print(sum(df1['eGFR'] == 999))

# %% cell 111
df1['gfr_bin']=pd.cut(x = df1['eGFR'],
                        bins = [-1,14.99,28.99,44.99,59.99,88.99,92, 999], 
                        labels = [1,2,3,4,5,6,7])

# %% cell 112
result = multivariate_logrank_test(df1['exposure'], df1['gfr_bin'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 113
# Define the categories
one_gfr = (df1.gfr_bin == 1) 
two_gfr = (df1.gfr_bin == 2)
three_gfr = (df1.gfr_bin == 3) 
four_gfr = (df1.gfr_bin == 4)
five_gfr = (df1.gfr_bin == 5)
six_gfr = (df1.gfr_bin == 6)
missing_gfr = (df1.gfr_bin == 7)
print(missing_gfr.value_counts())

# %% cell 114
# Plot the survival function
ax = plt.subplot()
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[one_gfr].exposure,
event_observed=df1[one_gfr].mortality_status, label='CKD 5')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[two_gfr].exposure, event_observed=df1[two_gfr].mortality_status, label='CKD 4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[three_gfr].exposure, event_observed=df1[three_gfr].mortality_status, label='CKD 3b')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[four_gfr].exposure, event_observed=df1[four_gfr].mortality_status, label='CKD 3a')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[five_gfr].exposure, event_observed=df1[five_gfr].mortality_status, label='CKD 2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[six_gfr].exposure, event_observed=df1[six_gfr].mortality_status, label='Normal')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing_gfr].exposure, event_observed=df1[missing_gfr].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.title('Survival Function based on eGFR category')
plt.show()

# %% cell 115
df1 = df10_cont[['eGFR', 'exposure', 'mortality_status']]

# %% cell 116
bin_labels = [1,2,3,4]
df1['eGFR'] = pd.qcut(df1['eGFR'], q=4, labels=bin_labels)

# %% cell 117
result = multivariate_logrank_test(df1['exposure'], df1['eGFR'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 118
# Define the categories
Low = (df1['eGFR'] == 1) 
Medium = (df1['eGFR'] == 2)
High = (df1['eGFR'] == 3) 
Very_high = (df1['eGFR'] == 4)
missing = (df1['eGFR'].isna())
print(sum(df1['eGFR'].isna()))
print(missing.value_counts())
print(sum(df1['eGFR'].isna()))

# %% cell 119
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on eGFR category')
plt.show()

# %% cell 121
df10_cont['Duration_diabetes'].describe()

# %% cell 122
df1 = df10_cont[['Duration_diabetes', 'exposure', 'mortality_status']]

# %% cell 123
bin_labels = [1,2,3,4]
df1['Duration_diabetes'] = pd.qcut(df1['Duration_diabetes'], q=4, labels=bin_labels)

# %% cell 124
# Define the categories
Low = (df1.Duration_diabetes == 1) 
Medium = (df1.Duration_diabetes == 2)
High = (df1.Duration_diabetes == 3) 
Very_high = (df1.Duration_diabetes == 4)
missing = (df1.Duration_diabetes.isna())

# %% cell 125
sum(df1['Duration_diabetes'].isna())

# %% cell 126
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
#kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
#kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on Duration of Diabetes category')
plt.show()

# %% cell 128
df10_cont.head()

# %% cell 129
df10_cont['Albumin'].describe()

# %% cell 130
df1 = df10_cont[['Albumin', 'exposure', 'mortality_status']]

# %% cell 131
bin_labels = [1,2,3,4]
df1['Albumin'] = pd.qcut(df1['Albumin'], q=4, labels=bin_labels)

# %% cell 132
df1

# %% cell 133
result = multivariate_logrank_test(df1['exposure'], df1['Albumin'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 134
df1

# %% cell 135
# Define the categories
Low = (df1.Albumin == 1) 
Medium = (df1.Albumin == 2)
High = (df1.Albumin == 3) 
Very_high = (df1.Albumin == 4)
missing = (df1.Albumin.isna())

# %% cell 136
sum(df['Albumin'].isna())

# %% cell 137
missing.value_counts()

# %% cell 138
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on Albumin category')
plt.show()

# %% cell 140
df10_cont['BMI'].describe()

# %% cell 141
df1 = df10_cont[['BMI', 'exposure', 'mortality_status']]

# %% cell 142
sum(df1['BMI'].isna())

# %% cell 143
df1['BMI'] = df1['BMI'].fillna(0)

# %% cell 144
sum(df1['BMI'] == 0)

# %% cell 145
df1['bmi_bin']=pd.cut(x = df1['BMI'],
                        bins = [-1,15,18.499,24.499, 29.999,60], 
                        labels = [0,1,2,3,4])

# %% cell 146
result = multivariate_logrank_test(df1['exposure'], df1['bmi_bin'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 147
# Define the categories
one_bmi = (df1.bmi_bin == 1) 
two_bmi = (df1.bmi_bin == 2)
three_bmi = (df1.bmi_bin == 3) 
four_bmi = (df1.bmi_bin == 4)
missing_bmi = (df1.bmi_bin == 0)

# %% cell 148
missing_bmi.value_counts()

# %% cell 149
# Plot the survival function
ax = plt.subplot()
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[one_bmi].exposure,
event_observed=df1[one_bmi].mortality_status, label='Underweight')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[two_bmi].exposure, event_observed=df1[two_bmi].mortality_status, label='Normal weight')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[three_bmi].exposure, event_observed=df1[three_bmi].mortality_status, label='Overweight')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[four_bmi].exposure, event_observed=df1[four_bmi].mortality_status, label='Obese')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing_bmi].exposure, event_observed=df1[missing_bmi].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on BMI category')
plt.show()

# %% cell 150
df1 = df10_cont[['BMI', 'exposure', 'mortality_status']]

# %% cell 151
bin_labels = [1,2,3,4]
df1['BMI'] = pd.qcut(df1['BMI'], q=4, labels=bin_labels)

# %% cell 152
result = multivariate_logrank_test(df1['exposure'], df1['BMI'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 153
# Define the categories
Low = (df1['BMI'] == 1) 
Medium = (df1['BMI'] == 2)
High = (df1['BMI'] == 3) 
Very_high = (df1['BMI'] == 4)
missing = (df1['BMI'].isna())
print(sum(df1['BMI'].isna()))
print(missing.value_counts())
print(sum(df1['BMI'].isna()))

# %% cell 154
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on BMI category')
plt.show()

# %% cell 156
df10_cont['ALP'].describe()

# %% cell 157
df1 = df10_cont[['ALP', 'exposure', 'mortality_status']]

# %% cell 158
bin_labels = [1,2,3,4]
df1['ALP'] = pd.qcut(df1['ALP'], q=4, labels=bin_labels)

# %% cell 159
df1

# %% cell 160
result = multivariate_logrank_test(df1['exposure'], df1['ALP'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 161
# Define the categories
Low = (df1.ALP == 1) 
Medium = (df1.ALP == 2)
High = (df1.ALP == 3) 
Very_high = (df1.ALP == 4)
missing = (df1.ALP.isna())

# %% cell 162
sum(df['ALP'].isna())

# %% cell 163
missing.value_counts()

# %% cell 164
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on ALP category')
plt.show()

# %% cell 166
df10_cont['ALT'].describe()

# %% cell 167
df1 = df10_cont[['ALT', 'exposure', 'mortality_status']]

# %% cell 168
bin_labels = [1,2,3,4]
df1['ALT'] = pd.qcut(df1['ALT'], q=4, labels=bin_labels)

# %% cell 169
df1

# %% cell 170
result = multivariate_logrank_test(df1['exposure'], df1['ALT'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 171
# Define the categories
Low = (df1['ALT'] == 1) 
Medium = (df1['ALT'] == 2)
High = (df1['ALT'] == 3) 
Very_high = (df1['ALT'] == 4)
missing = (df1['ALT'].isna())
print(sum(df1['ALT'].isna()))
print(missing.value_counts())
print(sum(df1['ALT'].isna()))

# %% cell 172
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on ALT category')
plt.show()

# %% cell 174
df10_cont['bilirubin'].describe()

# %% cell 175
df1 = df10_cont[['bilirubin', 'exposure', 'mortality_status']]

# %% cell 176
bin_labels = [1,2,3,4]
df1['bilirubin'] = pd.qcut(df1['bilirubin'], q=4, labels=bin_labels)

# %% cell 177
df1

# %% cell 178
result = multivariate_logrank_test(df1['exposure'], df1['bilirubin'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 179
# Define the categories
Low = (df1['bilirubin'] == 1) 
Medium = (df1['bilirubin'] == 2)
High = (df1['bilirubin'] == 3) 
Very_high = (df1['bilirubin'] == 4)
missing = (df1['bilirubin'].isna())
print(sum(df1['bilirubin'].isna()))
print(missing.value_counts())
print(sum(df1['bilirubin'].isna()))

# %% cell 180
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on Bilirubin category')
plt.show()

# %% cell 182
df10_cont['total cholesterol'].describe()

# %% cell 183
df1 = df10_cont[['total cholesterol', 'exposure', 'mortality_status']]

# %% cell 184
bin_labels = [1,2,3,4]
df1['total cholesterol'] = pd.qcut(df1['total cholesterol'], q=4, labels=bin_labels)

# %% cell 185
df1

# %% cell 186
result = multivariate_logrank_test(df1['exposure'], df1['total cholesterol'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 187
# Define the categories
Low = (df1['total cholesterol'] == 1) 
Medium = (df1['total cholesterol'] == 2)
High = (df1['total cholesterol'] == 3) 
Very_high = (df1['total cholesterol'] == 4)
missing = (df1['total cholesterol'].isna())
print(sum(df1['total cholesterol'].isna()))
print(missing.value_counts())
print(sum(df1['total cholesterol'].isna()))

# %% cell 188
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on Total cholesterol category')
plt.show()

# %% cell 190
df1 = df10_cont[['serum creatinine', 'exposure', 'mortality_status']]

# %% cell 191
bin_labels = [1,2,3,4]
df1['serum creatinine'] = pd.qcut(df1['serum creatinine'], q=4, labels=bin_labels)

# %% cell 192
df1

# %% cell 193
result = multivariate_logrank_test(df1['exposure'], df1['serum creatinine'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 194
# Define the categories
Low = (df1['serum creatinine'] == 1) 
Medium = (df1['serum creatinine'] == 2)
High = (df1['serum creatinine'] == 3) 
Very_high = (df1['serum creatinine'] == 4)
missing = (df1['serum creatinine'].isna())
print(sum(df1['serum creatinine'].isna()))
print(missing.value_counts())
print(sum(df1['serum creatinine'].isna()))

# %% cell 195
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on Serum Creatinine category')
plt.show()

# %% cell 197
df1 = df10_cont[['Hb', 'exposure', 'mortality_status']]

# %% cell 198
bin_labels = [1,2,3,4]
df1['Hb'] = pd.qcut(df1['Hb'], q=4, labels=bin_labels)

# %% cell 199
result = multivariate_logrank_test(df1['exposure'], df1['Hb'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 200
# Define the categories
Low = (df1['Hb'] == 1) 
Medium = (df1['Hb'] == 2)
High = (df1['Hb'] == 3) 
Very_high = (df1['Hb'] == 4)
missing = (df1['Hb'].isna())
print(sum(df1['Hb'].isna()))
print(missing.value_counts())
print(sum(df1['Hb'].isna()))

# %% cell 201
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on Hb category')
plt.show()
plt.savefig('S:\FPHS_THIN_Projects\Karan\BZI_PhD_T2D\PhD_descriptive_data\KM_curves\Hb.png', bbox_inches='tight')

# %% cell 202
plt.savefig('S:\FPHS_THIN_Projects\Karan\BZI_PhD_T2D\PhD_descriptive_data\KM_curves\Hb.png', bbox_inches='tight')

# %% cell 204
df10_cont.columns

# %% cell 205
df1 = df10_cont[['HDL Cholesterol', 'exposure', 'mortality_status']]

# %% cell 206
bin_labels = [1,2,3,4]
df1['HDL Cholesterol'] = pd.qcut(df1['HDL Cholesterol'], q=4, labels=bin_labels)

# %% cell 207
result = multivariate_logrank_test(df1['exposure'], df1['HDL Cholesterol'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 208
# Define the categories
Low = (df1['HDL Cholesterol'] == 1) 
Medium = (df1['HDL Cholesterol'] == 2)
High = (df1['HDL Cholesterol'] == 3) 
Very_high = (df1['HDL Cholesterol'] == 4)
missing = (df1['HDL Cholesterol'].isna())
print(sum(df1['HDL Cholesterol'].isna()))
print(missing.value_counts())
print(sum(df1['HDL Cholesterol'].isna()))

# %% cell 209
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on HDL-Cholesterol category')
plt.show()

# %% cell 211
df10_cont.columns

# %% cell 212
df1 = df10_cont[['LDL-C', 'exposure', 'mortality_status']]

# %% cell 213
bin_labels = [1,2,3,4]
df1['LDL-C'] = pd.qcut(df1['LDL-C'], q=4, labels=bin_labels)

# %% cell 214
result = multivariate_logrank_test(df1['exposure'], df1['LDL-C'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 215
# Define the categories
Low = (df1['LDL-C'] == 1) 
Medium = (df1['LDL-C'] == 2)
High = (df1['LDL-C'] == 3) 
Very_high = (df1['LDL-C'] == 4)
missing = (df1['LDL-C'].isna())
print(sum(df1['LDL-C'].isna()))
print(missing.value_counts())
print(sum(df1['LDL-C'].isna()))

# %% cell 216
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on LDL-C category')
plt.show()

# %% cell 218
df1 = df10_cont[['Sodium', 'exposure', 'mortality_status']]

# %% cell 219
bin_labels = [1,2,3,4]
df1['Sodium'] = pd.qcut(df1['Sodium'], q=4, labels=bin_labels)

# %% cell 220
df1

# %% cell 221
result = multivariate_logrank_test(df1['exposure'], df1['Sodium'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 222
# Define the categories
Low = (df1['Sodium'] == 1) 
Medium = (df1['Sodium'] == 2)
High = (df1['Sodium'] == 3) 
Very_high = (df1['Sodium'] == 4)
missing = (df1['Sodium'].isna())
print(sum(df1['Sodium'].isna()))
print(missing.value_counts())
print(sum(df1['Sodium'].isna()))

# %% cell 223
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on Sodium category')
plt.show()

# %% cell 225
df1 = df10_cont[['Total protein', 'exposure', 'mortality_status']]

# %% cell 226
bin_labels = [1,2,3,4]
df1['Total protein'] = pd.qcut(df1['Total protein'], q=4, labels=bin_labels)

# %% cell 227
df1

# %% cell 228
result = multivariate_logrank_test(df1['exposure'], df1['Total protein'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 229
# Define the categories
Low = (df1['Total protein'] == 1) 
Medium = (df1['Total protein'] == 2)
High = (df1['Total protein'] == 3) 
Very_high = (df1['Total protein'] == 4)
missing = (df1['Total protein'].isna())
print(sum(df1['Total protein'].isna()))
print(missing.value_counts())
print(sum(df1['Total protein'].isna()))

# %% cell 230
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on Total protein category')
plt.show()

# %% cell 232
df1 = df10_cont[['Triglycerides', 'exposure', 'mortality_status']]

# %% cell 233
bin_labels = [1,2,3,4]
df1['Triglycerides'] = pd.qcut(df1['Triglycerides'], q=4, labels=bin_labels)

# %% cell 234
result = multivariate_logrank_test(df1['exposure'], df1['Triglycerides'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 235
# Define the categories
Low = (df1['Triglycerides'] == 1) 
Medium = (df1['Triglycerides'] == 2)
High = (df1['Triglycerides'] == 3) 
Very_high = (df1['Triglycerides'] == 4)
missing = (df1['Triglycerides'].isna())
print(sum(df1['Triglycerides'].isna()))
print(missing.value_counts())
print(sum(df1['Triglycerides'].isna()))

# %% cell 236
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on Triglcyerides category')
plt.show()

# %% cell 238
df1 = df10_cont[['MCV', 'exposure', 'mortality_status']]

# %% cell 239
bin_labels = [1,2,3,4]
df1['MCV'] = pd.qcut(df1['MCV'], q=4, labels=bin_labels)

# %% cell 240
result = multivariate_logrank_test(df1['exposure'], df1['MCV'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 241
# Define the categories
Low = (df1['MCV'] == 1) 
Medium = (df1['MCV'] == 2)
High = (df1['MCV'] == 3) 
Very_high = (df1['MCV'] == 4)
missing = (df1['MCV'].isna())
print(sum(df1['MCV'].isna()))
print(missing.value_counts())
print(sum(df1['MCV'].isna()))

# %% cell 242
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on MCV category')
plt.show()

# %% cell 244
df10_cont.columns

# %% cell 245
df1 = df10_cont[['TSH', 'exposure', 'mortality_status']]

# %% cell 246
bin_labels = [1,2,3,4]
df1['TSH'] = pd.qcut(df1['TSH'], q=4, labels=bin_labels)

# %% cell 247
result = multivariate_logrank_test(df1['exposure'], df1['TSH'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 248
# Define the categories
Low = (df1['TSH'] == 1) 
Medium = (df1['TSH'] == 2)
High = (df1['TSH'] == 3) 
Very_high = (df1['TSH'] == 4)
missing = (df1['TSH'].isna())
print(sum(df1['TSH'].isna()))
print(missing.value_counts())
print(sum(df1['TSH'].isna()))

# %% cell 249
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on TSH category')
plt.show()

# %% cell 251
df10_cont.columns

# %% cell 252
df1 = df10_cont[['Urea', 'exposure', 'mortality_status']]

# %% cell 253
bin_labels = [1,2,3,4]
df1['Urea'] = pd.qcut(df1['Urea'], q=4, labels=bin_labels)

# %% cell 254
result = multivariate_logrank_test(df1['exposure'], df1['Urea'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 255
# Define the categories
Low = (df1['Urea'] == 1) 
Medium = (df1['Urea'] == 2)
High = (df1['Urea'] == 3) 
Very_high = (df1['Urea'] == 4)
missing = (df1['Urea'].isna())
print(sum(df1['Urea'].isna()))
print(missing.value_counts())
print(sum(df1['Urea'].isna()))

# %% cell 256
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on Urea category')
plt.show()

# %% cell 258
df1 = df10_cont[['Neutrophil count', 'exposure', 'mortality_status']]

# %% cell 259
bin_labels = [1,2,3,4]
df1['Neutrophil count'] = pd.qcut(df1['Neutrophil count'], q=4, labels=bin_labels)

# %% cell 260
result = multivariate_logrank_test(df1['exposure'], df1['Neutrophil count'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 261
# Define the categories
Low = (df1['Neutrophil count'] == 1) 
Medium = (df1['Neutrophil count'] == 2)
High = (df1['Neutrophil count'] == 3) 
Very_high = (df1['Neutrophil count'] == 4)
missing = (df1['Neutrophil count'].isna())
print(sum(df1['Neutrophil count'].isna()))
print(missing.value_counts())
print(sum(df1['Neutrophil count'].isna()))

# %% cell 262
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on Neutrophil count category')
plt.show()

# %% cell 264
df1 = df10_cont[['WCC', 'exposure', 'mortality_status']]

# %% cell 265
bin_labels = [1,2,3,4]
df1['WCC'] = pd.qcut(df1['WCC'], q=4, labels=bin_labels)

# %% cell 266
result = multivariate_logrank_test(df1['exposure'], df1['WCC'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 267
# Define the categories
Low = (df1['WCC'] == 1) 
Medium = (df1['WCC'] == 2)
High = (df1['WCC'] == 3) 
Very_high = (df1['WCC'] == 4)
missing = (df1['WCC'].isna())
print(sum(df1['WCC'].isna()))
print(missing.value_counts())
print(sum(df1['WCC'].isna()))

# %% cell 268
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on WCC category')
plt.show()

# %% cell 270
df1 = df10_cont[['Lymphocyte count', 'exposure', 'mortality_status']]

# %% cell 271
bin_labels = [1,2,3,4]
df1['Lymphocyte count'] = pd.qcut(df1['Lymphocyte count'], q=4, labels=bin_labels)

# %% cell 272
df1

# %% cell 273
result = multivariate_logrank_test(df1['exposure'], df1['Lymphocyte count'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 274
# Define the categories
Low = (df1['Lymphocyte count'] == 1) 
Medium = (df1['Lymphocyte count'] == 2)
High = (df1['Lymphocyte count'] == 3) 
Very_high = (df1['Lymphocyte count'] == 4)
missing = (df1['Lymphocyte count'].isna())
print(sum(df1['Lymphocyte count'].isna()))
print(missing.value_counts())
print(sum(df1['Lymphocyte count'].isna()))

# %% cell 275
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on lymphocyte count category')
plt.show()

# %% cell 277
df1 = df10_cont[['Potassium', 'exposure', 'mortality_status']]

# %% cell 278
bin_labels = [1,2,3,4]
df1['Potassium'] = pd.qcut(df1['Potassium'], q=4, labels=bin_labels)

# %% cell 279
df1

# %% cell 280
result = multivariate_logrank_test(df1['exposure'], df1['Potassium'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 281
# Define the categories
Low = (df1['Potassium'] == 1) 
Medium = (df1['Potassium'] == 2)
High = (df1['Potassium'] == 3) 
Very_high = (df1['Potassium'] == 4)
missing = (df1['Potassium'].isna())
print(sum(df1['Potassium'].isna()))
print(missing.value_counts())
print(sum(df1['Potassium'].isna()))

# %% cell 282
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on potassium category')
plt.show()

# %% cell 284
df10_cont.columns

# %% cell 285
df10_cont[['patprac', 'pulse_pressure']]

# %% cell 286
df1 = df10_cont[['pulse_pressure', 'exposure', 'mortality_status']]

# %% cell 287
bin_labels = [1,2,3,4]
df1['pulse_pressure'] = pd.qcut(df1['pulse_pressure'], q=4, labels=bin_labels)

# %% cell 288
# Define the categories
Low = (df1['pulse_pressure'] == 1) 
Medium = (df1['pulse_pressure'] == 2)
High = (df1['pulse_pressure'] == 3) 
Very_high = (df1['pulse_pressure'] == 4)
missing = (df1['pulse_pressure'].isna())
print(sum(df1['pulse_pressure'].isna()))
print(missing.value_counts())
print(sum(df1['pulse_pressure'].isna()))

# %% cell 289
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on Pulse Pressure category')
plt.show()

# %% cell 291
df1 = df10_cont[['mean arterial pressure', 'exposure', 'mortality_status']]

# %% cell 292
bin_labels = [1,2,3,4]
df1['mean arterial pressure'] = pd.qcut(df1['mean arterial pressure'], q=4, labels=bin_labels)

# %% cell 293
# Define the categories
Low = (df1['mean arterial pressure'] == 1) 
Medium = (df1['mean arterial pressure'] == 2)
High = (df1['mean arterial pressure'] == 3) 
Very_high = (df1['mean arterial pressure'] == 4)
missing = (df1['mean arterial pressure'].isna())
print(sum(df1['mean arterial pressure'].isna()))
print(missing.value_counts())
print(sum(df1['mean arterial pressure'].isna()))

# %% cell 294
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on Mean Arterial Pressure category')
plt.show()

# %% cell 296
df1 = df10_cont[['totalchol/hdl_ratio', 'exposure', 'mortality_status']]

# %% cell 297
bin_labels = [1,2,3,4]
df1['totalchol/hdl_ratio'] = pd.qcut(df1['totalchol/hdl_ratio'], q=4, labels=bin_labels)

# %% cell 298
# Define the categories
Low = (df1['totalchol/hdl_ratio'] == 1) 
Medium = (df1['totalchol/hdl_ratio'] == 2)
High = (df1['totalchol/hdl_ratio'] == 3) 
Very_high = (df1['totalchol/hdl_ratio'] == 4)
missing = (df1['totalchol/hdl_ratio'].isna())
print(sum(df1['totalchol/hdl_ratio'].isna()))
print(missing.value_counts())
print(sum(df1['totalchol/hdl_ratio'].isna()))

# %% cell 299
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on Total chol / HDL ratio category')
plt.show()

# %% cell 301
df1 = df10_cont[['Neutrophil/lymphocyte', 'exposure', 'mortality_status']]

# %% cell 302
bin_labels = [1,2,3,4]
df1['Neutrophil/lymphocyte'] = pd.qcut(df1['Neutrophil/lymphocyte'], q=4, labels=bin_labels)

# %% cell 303
# Define the categories
Low = (df1['Neutrophil/lymphocyte'] == 1) 
Medium = (df1['Neutrophil/lymphocyte'] == 2)
High = (df1['Neutrophil/lymphocyte'] == 3) 
Very_high = (df1['Neutrophil/lymphocyte'] == 4)
missing = (df1['Neutrophil/lymphocyte'].isna())
print(sum(df1['Neutrophil/lymphocyte'].isna()))
print(missing.value_counts())
print(sum(df1['Neutrophil/lymphocyte'].isna()))

# %% cell 304
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[Low].exposure, event_observed=df1[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Medium].exposure, event_observed=df1[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[High].exposure, event_observed=df1[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[Very_high].exposure, event_observed=df1[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on Neutrophil / lymphocyte ratio category')
plt.show()

# %% cell 306
data_acr = pd.read_csv('S:\FPHS_THIN_Projects\Karan\BZI_PhD_T2D\PhD_descriptive_data/desc_data_acr.csv') 

# %% cell 307
data_acr.columns

# %% cell 308
data_acr['acr'].describe()

# %% cell 309
data_acr.shape

# %% cell 310
(sum(data_acr['acr'].isna()) / 155399)*100

# %% cell 311
sum(data_acr['acr'].isna())

# %% cell 312
df_acr = data_acr[['acr', 'exposure', 'mortality_status']]

# %% cell 313
bin_labels = [1,2,3,4]
df_acr['acr_bin'] = pd.qcut(df_acr['acr'], q=4, labels=bin_labels)

# %% cell 314
df_acr

# %% cell 315
result = multivariate_logrank_test(df_acr['exposure'], df_acr['acr_bin'], df_acr['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 316
sum(df_acr['acr'].isna())

# %% cell 317
df_acr['acr_bin'].unique()

# %% cell 318
# Define the categories
Low = (df_acr['acr'] == 1) 
Medium = (df_acr['acr'] == 2)
High = (df_acr['acr'] == 3) 
Very_high = (df_acr['acr'] == 4)
missing = (df_acr['acr'].isna())
print(sum(df_acr['acr'].isna()))
print(missing.value_counts())
print(sum(df_acr['acr'].isna()))

# %% cell 319
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df_acr[Low].exposure, event_observed=df_acr[Low].mortality_status, label='Q1')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df_acr[Medium].exposure, event_observed=df_acr[Medium].mortality_status, label='Q2')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df_acr[High].exposure, event_observed=df_acr[High].mortality_status, label='Q3')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df_acr[Very_high].exposure, event_observed=df_acr[Very_high].mortality_status, label='Q4')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df_acr[missing].exposure, event_observed=df_acr[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on Albumin Creatinine ratio category')
plt.show()

# %% cell 320
result = multivariate_logrank_test(df_acr['exposure'], df_acr['acr_bin'], df_acr['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 322
data_acr.columns

# %% cell 323
df1 = data_acr[['Smoking_status', 'exposure', 'mortality_status']]

# %% cell 324
df1

# %% cell 325
df1['Smoking_status'].unique()

# %% cell 326
df1.fillna(3, inplace=True)

# %% cell 327
df1['Smoking_status'].unique()

# %% cell 328
df1

# %% cell 329
result = multivariate_logrank_test(df1['exposure'], df1['Smoking_status'], df1['mortality_status'])
result.test_statistic
result.p_value
result.print_summary()

# %% cell 330
    # Define the categories
non = (df1.Smoking_status == 0) 
ex = (df1.Smoking_status == 1)
current = (df1.Smoking_status == 2) 
missing = (df1.Smoking_status == 3)
# print(missing.value_counts())

# %% cell 331
# Plot the survival function
ax = plt.subplot(111)
kmf = KaplanMeierFitter()
kmf.fit(durations=df1[non].exposure, event_observed=df1[non].mortality_status, label='non-smoker')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[ex].exposure, event_observed=df1[ex].mortality_status, label='ex-smoker')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[current].exposure, event_observed=df1[current].mortality_status, label='current smoker')
kmf.survival_function_.plot(ax=ax)
kmf.fit(durations=df1[missing].exposure, event_observed=df1[missing].mortality_status, label='missing')
kmf.survival_function_.plot(ax=ax)
plt.ylim(0, 1.0)
plt.xlim(0, 10)
plt.xlabel('timeline (years)')
plt.title('Survival Function based on Smoking status')
plt.show()
