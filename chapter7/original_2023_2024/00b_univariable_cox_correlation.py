"""
Univariable Cox regression and Pearson correlation matrix (univariable table; collinearity screen)
Original Jupyter notebook exported as a script (code cells only; outputs removed). Source last modified 2024-01-30.
"""

# %% cell 1
import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
import seaborn as sb

# %% cell 2
from lifelines import CoxPHFitter
from lifelines.datasets import load_rossi

# %% cell 4
df = pd.read_stata('S:\FPHS_THIN_Projects\Karan\BZI_PhD_T2D\PhD_descriptive_data\mice_data\imp_phddata_170423.dta')

# %% cell 5
list(df.columns)

# %% cell 6
df1 = df[['patprac','sex','duration_diabetes', 'mortality_status',
 'age',
 'exposure',
 'aceiarb',
 'ca_blocker',
 'beta_blocker',
 'thiazides',
 'alpha_blocker',
 'metformin',
 'sulphonylureas',
 'meglitinides',
 'gliptins',
 'thiazols',
 'glp1',
 'acarbose',
 'insulin',
 'statins',
 'hx_resp',
 'hx_liver',
 'hx_gi',
 'hx_neuro',
 'hx_autoimmune',
 'hx_severekidney',
 'hx_cancer',
 'hx_proteinuria',
 'hx_neuropathy',
 'hx_retinopathy',
 'hx_nephropathy',
 'hx_footlegulcer',
 'hx_comashock',
 'hx_cad',
 'hx_cerebro',
 'hx_pvd',
 'hx_arrhythmia',
 'hx_heartfailure',
 'hx_mi',
 'hx_dvt', '_1_acr', '_1_albumin', '_1_alp', '_1_alt', '_1_bilirubin', '_1_totalcholesterol', '_1_serumcreatinine', '_1_hb', '_1_hdlcholesterol',
         '_1_sodium', '_1_ldlc', '_1_totalprotein', '_1_triglycerides', '_1_mcv', '_1_tsh', '_1_urea', '_1_neutrophilcount', '_1_wcc', '_1_lymphocytecount',
         '_1_potassium', '_1_hba1c', '_1_egfr', '_1_smoking_status', '_1_bmi', '_1_diastolicbloodpressure', '_1_systolicbloodpressure']]   

# %% cell 7
df1.shape

# %% cell 8
df1.columns

# %% cell 9
cph = CoxPHFitter()
cph.fit(df1[['_1_triglycerides', 'exposure', 'mortality_status']], 'exposure', 'mortality_status')
cph.print_summary()

# %% cell 11
df2 = df1[['sex', 'duration_diabetes','age',
       'aceiarb', 'ca_blocker', 'beta_blocker',
       'alpha_blocker', 'metformin', 'sulphonylureas', 
       'gliptins', 'thiazols', 'glp1', 'acarbose', 'insulin', 'statins',
       'hx_resp', 'hx_liver', 'hx_gi', 'hx_neuro', 'hx_autoimmune',
       'hx_severekidney', 'hx_cancer', 'hx_proteinuria', 'hx_neuropathy',
       'hx_retinopathy', 'hx_nephropathy', 'hx_footlegulcer', 'hx_comashock',
       'hx_cad', 'hx_cerebro', 'hx_pvd', 'hx_arrhythmia', 'hx_heartfailure',
       'hx_mi', 'hx_dvt', '_1_acr', '_1_albumin', '_1_alp', '_1_alt',
        '_1_totalcholesterol', '_1_serumcreatinine', '_1_hb',
       '_1_hdlcholesterol', '_1_sodium', '_1_ldlc', '_1_totalprotein',
       '_1_triglycerides', '_1_mcv', '_1_tsh', '_1_urea', '_1_neutrophilcount',
       '_1_wcc', '_1_lymphocytecount', '_1_potassium', '_1_hba1c', '_1_egfr',
       '_1_smoking_status', '_1_bmi', '_1_diastolicbloodpressure',
       '_1_systolicbloodpressure']]

# %% cell 12
df2.shape

# %% cell 13
pearsoncorr = df2.corr(method='pearson')

# %% cell 14
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

# %% cell 15
pearsoncorr.head()

# %% cell 16
import pysurvival

# %% cell 17
from pysurvival.utils.display import correlation_matrix

# %% cell 18
correlation_matrix(df2, figure_size=(50,30), text_fontsize=10)

# %% cell 19
dftest = pd.read_csv('S:\FPHS_THIN_Projects\Karan\BZI_PhD_T2D\PhD_descriptive_data\mice_data\phd_mice_data.csv')

# %% cell 20
dftest.shape

# %% cell 21
dftest.columns

# %% cell 22
phd_mice_data2 = dftest[['patprac', 'sex', 'Duration_diabetes', 'Albumin', 'ALP', 'ALT',
        'total cholesterol', 'serum creatinine', 'Hb',
       'HDL Cholesterol', 'Sodium','Total protein', 'Triglycerides',
       'MCV', 'TSH', 'Urea', 'Neutrophil count', 'Lymphocyte count',
       'Potassium', 'Hba1c', 'eGFR', 'Smoking_status', 'BMI',
       'Diastolic Blood Pressure', 'Systolic Blood Pressure',
       'mortality_status', 'age', 'exposure', 'acei/arb', 'ca_blocker',
       'beta_blocker', 'alpha_blocker', 'metformin',
       'sulphonylureas', 'gliptins', 'thiazols', 'glp1',
       'acarbose', 'insulin', 'statins', 'hx_resp', 'hx_liver', 'hx_gi',
       'hx_neuro', 'hx_autoimmune', 'hx_severekidney', 'hx_cancer',
       'hx_proteinuria', 'hx_neuropathy', 'hx_retinopathy', 'hx_nephropathy',
       'hx_footlegulcer', 'hx_comashock', 'hx_cad', 'hx_cerebro', 'hx_pvd',
       'hx_arrhythmia', 'hx_heartfailure', 'hx_mi', 'hx_dvt', 'acr']]

# %% cell 23
phd_mice_data2.to_csv('S:\FPHS_THIN_Projects\Karan\BZI_PhD_T2D\PhD_descriptive_data\mice_data\phd_mice_data_290124.csv', index=False)

# %% cell 24
nonlinear = pd.read_stata('S:/FPHS_THIN_Projects/Karan/BZI_PhD_T2D/PhD_descriptive_data/mice_data/nonlinearvariables_080923.dta')

# %% cell 25
nonlinear.shape

# %% cell 26
nonlinear.columns

# %% cell 27
nonlinear1 = nonlinear[['patprac', '_1_albumin_0', '_1_albumin_1', 'age_0', 'age_1',
       '_1_alp_0', '_1_alp_1', '_1_alp_2',
       '_1_alp_3', '_1_alt_0', '_1_alt_1', '_1_alt_2',
       '_1_totalcholesterol_0', '_1_totalcholesterol_1',
       '_1_serumcreatinine_0', '_1_serumcreatinine_1', '_1_serumcreatinine_2',
       '_1_serumcreatinine_3', '_1_hb_0', '_1_hb_1', '_1_hb_2', '_1_hb_3',
       '_1_hdlcholesterol_0', '_1_hdlcholesterol_1', '_1_hdlcholesterol_2',
       '_1_hdlcholesterol_3', '_1_sodium_0', '_1_sodium_1', '_1_sodium_2',
       '_1_sodium_3', '_1_totalprotein_0', '_1_totalprotein_1',
       '_1_totalprotein_2', '_1_totalprotein_3', '_1_triglycerides_0',
       '_1_triglycerides_1', '_1_triglycerides_2', '_1_mcv_0', '_1_mcv_1',
       '_1_mcv_2', '_1_mcv_3', '_1_tsh_0', '_1_tsh_1', '_1_tsh_2', '_1_urea_0',
       '_1_urea_1', '_1_urea_2', '_1_urea_3', '_1_neutrophilcount_0',
       '_1_neutrophilcount_1', '_1_neutrophilcount_2','_1_lymphocytecount_0', '_1_lymphocytecount_1',
       '_1_lymphocytecount_2', '_1_lymphocytecount_3', '_1_potassium_0',
       '_1_potassium_1', '_1_hba1c_0', '_1_hba1c_1', '_1_hba1c_2', '_1_egfr_0',
       '_1_egfr_1', '_1_egfr_2', '_1_egfr_3', '_1_acr_0', '_1_acr_1',
       '_1_acr_2', '_1_acr_3', '_1_bmi_0', '_1_bmi_1', '_1_bmi_2', '_1_bmi_3',
       '_1_diastolicbloodpressure_0', '_1_diastolicbloodpressure_1',
       '_1_systolicbloodpressure_0', '_1_systolicbloodpressure_1',
       '_1_systolicbloodpressure_2', '_1_systolicbloodpressure_3',
       'duration_diabetes_0', 'duration_diabetes_1']]

# %% cell 28
nonlinear1.to_csv('S:/FPHS_THIN_Projects/Karan/BZI_PhD_T2D/PhD_descriptive_data/mice_data/nonlinearvariables_290124.csv', index=False)

# %% cell 29
nonlinear1.shape

# %% cell 30
dfdf = pd.read_csv('S:/FPHS_THIN_Projects/Karan/BZI_PhD_T2D/PhD_descriptive_data/mice_data/nonlinearvariables_290124.csv')

# %% cell 31
dfdf.shape

# %% cell 32
dfdf.head()

# %% cell 33
df25 = pd.read_stata('S:/FPHS_THIN_Projects/Karan/BZI_PhD_T2D/PhD_descriptive_data/mice_data/nonlineartermsdata_310823.dta')

# %% cell 34
df25.shape

# %% cell 35
df25.columns
