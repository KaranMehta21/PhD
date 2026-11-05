*Set
	clear
	set more off
*cap log close*
*Import dataset
use "S:\FPHS_THIN_Projects\Karan\BZI_PhD_T2D\PhD_descriptive_data\mice_data\nonlineartermsdata_310823.dta"

stset exposure, failure(mortality_status)
	
*Model 1: simple cox model*

stcox _1_albumin _1_alp _1_alt _1_bilirubin _1_totalcholesterol _1_serumcreatinine _1_hb _1_hdlcholesterol _1_sodium _1_ldlc _1_totalprotein _1_triglycerides _1_mcv _1_tsh _1_urea _1_neutrophilcount _1_wcc _1_lymphocytecount _1_potassium _1_hba1c _1_egfr _1_acr _1_bmi _1_diastolicbloodpressure _1_systolicbloodpressure  duration_diabetes age i._1_smoking_status i.hx_neuro i.acarbose i.insulin i.statins i.hx_resp i.hx_liver i.hx_gi i.hx_severekidney i.hx_autoimmune i.thiazols i.hx_cancer i.hx_proteinuria i.hx_neuropathy i.hx_retinopathy i.hx_cad i.hx_cerebro i.hx_pvd i.hx_arrhythmia i.hx_heartfailure i.hx_mi i.hx_dvt i.hx_nephropathy i.hx_footlegulcer i.glp1 i.gliptins i.meglitinides i.sex i.aceiarb i.ca_blocker i.beta_blocker i.thiazides i.alpha_blocker i.metformin i.sulphonylureas i.hx_comashock

estimates store model1 

*Model 2: cox model with interaction terms*

stcox i.sex#i.hx_neuro i.sex#c.age i.sex#i.thiazides c.age#c._1_albumin c.age#c._1_alp c.age#c._1_bilirubin c.age#c._1_sodium c.age#c._1_mcv c.age#i._1_smoking_status c.age#c._1_bmi c.age#i.hx_resp c.age#i.hx_liver c.age#i.hx_severekidney c.age#i.hx_cancer c.age#i.hx_pvd _1_albumin _1_alp _1_alt _1_bilirubin _1_totalcholesterol _1_serumcreatinine _1_hb _1_hdlcholesterol _1_sodium _1_ldlc _1_totalprotein _1_triglycerides _1_mcv _1_tsh _1_urea _1_neutrophilcount _1_wcc _1_lymphocytecount _1_potassium _1_hba1c _1_egfr _1_acr _1_bmi _1_diastolicbloodpressure _1_systolicbloodpressure  duration_diabetes age i._1_smoking_status i.hx_neuro i.acarbose i.insulin i.statins i.hx_resp i.hx_liver i.hx_gi i.hx_severekidney i.hx_autoimmune i.thiazols i.hx_cancer i.hx_proteinuria i.hx_neuropathy i.hx_retinopathy i.hx_cad i.hx_cerebro i.hx_pvd i.hx_arrhythmia i.hx_heartfailure i.hx_mi i.hx_dvt i.hx_nephropathy i.hx_footlegulcer i.glp1 i.gliptins i.meglitinides i.sex i.aceiarb i.ca_blocker i.beta_blocker i.thiazides i.alpha_blocker i.metformin i.sulphonylureas i.hx_comashock

estimates store model2

*Model 3: cox model with nonlinearterms*

stcox _1_albumin_0 _1_albumin_1 _1_alp_0 _1_alp_1 _1_alp_2 _1_alp_3 _1_alt_0 _1_alt_1 _1_alt_2 _1_bilirubin_0 _1_bilirubin_1 _1_bilirubin_2 _1_bilirubin_3 _1_totalcholesterol_0 _1_totalcholesterol_1 _1_serumcreatinine_0 _1_serumcreatinine_1 _1_serumcreatinine_2 _1_serumcreatinine_3 _1_hb_0 _1_hb_1 _1_hb_2 _1_hb_3 _1_hdlcholesterol_0 _1_hdlcholesterol_1 _1_hdlcholesterol_2 _1_hdlcholesterol_3 _1_sodium_0 _1_sodium_1 _1_sodium_2 _1_sodium_3 _1_ldlc_0 _1_ldlc_1 _1_totalprotein_0 _1_totalprotein_1 _1_totalprotein_2 _1_totalprotein_3  _1_triglycerides_0 _1_triglycerides_1 _1_triglycerides_2 _1_mcv_0 _1_mcv_1 _1_mcv_2 _1_mcv_3 _1_tsh_0 _1_tsh_1 _1_tsh_2 _1_urea_0 _1_urea_1 _1_urea_2 _1_urea_3 _1_neutrophilcount_0 _1_neutrophilcount_1 _1_neutrophilcount_2 _1_wcc_0 _1_wcc_1 _1_wcc_2 _1_lymphocytecount_0 _1_lymphocytecount_1 _1_lymphocytecount_2 _1_lymphocytecount_3 _1_potassium_0 _1_potassium_1  _1_hba1c_0 _1_hba1c_1 _1_hba1c_2 _1_egfr_0 _1_egfr_1 _1_egfr_2 _1_egfr_3 _1_acr_1 _1_acr_2 _1_acr_3 _1_bmi_0 _1_bmi_1 _1_bmi_2 _1_bmi_3 _1_diastolicbloodpressure_0 _1_diastolicbloodpressure_1  _1_systolicbloodpressure_0 _1_systolicbloodpressure_1 _1_systolicbloodpressure_2 _1_systolicbloodpressure_3  duration_diabetes_0 duration_diabetes_1 age_0 age_1 i._1_smoking_status i.hx_neuro i.acarbose i.insulin i.statins i.hx_resp i.hx_liver i.hx_gi i.hx_severekidney i.hx_autoimmune i.thiazols i.hx_cancer i.hx_proteinuria i.hx_neuropathy i.hx_retinopathy i.hx_cad i.hx_cerebro i.hx_pvd i.hx_arrhythmia i.hx_heartfailure i.hx_mi i.hx_dvt i.hx_nephropathy i.hx_footlegulcer i.glp1 i.gliptins i.meglitinides i.sex i.aceiarb i.ca_blocker i.beta_blocker i.thiazides i.alpha_blocker i.metformin i.sulphonylureas i.hx_comashock

estimates store model3

*Model 4: cox model with nonlinearterms and interaction terms*

stcox i.sex#i.hx_neuro i.sex#c.age i.sex#i.thiazides c.age#c._1_albumin c.age#c._1_alp c.age#c._1_bilirubin c.age#c._1_sodium c.age#c._1_mcv c.age#i._1_smoking_status c.age#c._1_bmi c.age#i.hx_resp c.age#i.hx_liver c.age#i.hx_severekidney c.age#i.hx_cancer c.age#i.hx_pvd _1_albumin_0 _1_albumin_1 _1_alp_0 _1_alp_1 _1_alp_2 _1_alp_3 _1_alt_0 _1_alt_1 _1_alt_2 _1_bilirubin_0 _1_bilirubin_1 _1_bilirubin_2 _1_bilirubin_3 _1_totalcholesterol_0 _1_totalcholesterol_1 _1_serumcreatinine_0 _1_serumcreatinine_1 _1_serumcreatinine_2 _1_serumcreatinine_3 _1_hb_0 _1_hb_1 _1_hb_2 _1_hb_3 _1_hdlcholesterol_0 _1_hdlcholesterol_1 _1_hdlcholesterol_2 _1_hdlcholesterol_3 _1_sodium_0 _1_sodium_1 _1_sodium_2 _1_sodium_3 _1_ldlc_0 _1_ldlc_1 _1_totalprotein_0 _1_totalprotein_1 _1_totalprotein_2 _1_totalprotein_3  _1_triglycerides_0 _1_triglycerides_1 _1_triglycerides_2 _1_mcv_0 _1_mcv_1 _1_mcv_2 _1_mcv_3 _1_tsh_0 _1_tsh_1 _1_tsh_2 _1_urea_0 _1_urea_1 _1_urea_2 _1_urea_3 _1_neutrophilcount_0 _1_neutrophilcount_1 _1_neutrophilcount_2 _1_wcc_0 _1_wcc_1 _1_wcc_2 _1_lymphocytecount_0 _1_lymphocytecount_1 _1_lymphocytecount_2 _1_lymphocytecount_3 _1_potassium_0 _1_potassium_1  _1_hba1c_0 _1_hba1c_1 _1_hba1c_2 _1_egfr_0 _1_egfr_1 _1_egfr_2 _1_egfr_3 _1_acr_1 _1_acr_2 _1_acr_3 _1_bmi_0 _1_bmi_1 _1_bmi_2 _1_bmi_3 _1_diastolicbloodpressure_0 _1_diastolicbloodpressure_1  _1_systolicbloodpressure_0 _1_systolicbloodpressure_1 _1_systolicbloodpressure_2 _1_systolicbloodpressure_3  duration_diabetes_0 duration_diabetes_1 age_0 age_1 i._1_smoking_status i.hx_neuro i.acarbose i.insulin i.statins i.hx_resp i.hx_liver i.hx_gi i.hx_severekidney i.hx_autoimmune i.thiazols i.hx_cancer i.hx_proteinuria i.hx_neuropathy i.hx_retinopathy i.hx_cad i.hx_cerebro i.hx_pvd i.hx_arrhythmia i.hx_heartfailure i.hx_mi i.hx_dvt i.hx_nephropathy i.hx_footlegulcer i.glp1 i.gliptins i.meglitinides i.sex i.aceiarb i.ca_blocker i.beta_blocker i.thiazides i.alpha_blocker i.metformin i.sulphonylureas i.hx_comashock

estimates store model4

estimates stats _all

lrtest model4 model3

lrtest model4 model2

lrtest model4 model1


