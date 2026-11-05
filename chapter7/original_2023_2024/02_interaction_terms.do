*Set
	clear
	set more off
*cap log close*
*Import dataset
	import delimited "S:/FPHS_THIN_Projects/Karan/BZI_PhD_T2D/PhD_descriptive_data/mice_data/1_imputeddataset_250823.csv"
	
stset exposure, failure(mortality_status)
	/*CCA MODEL*/
	stcox _1_albumin _1_alp _1_alt _1_bilirubin _1_totalcholesterol _1_serumcreatinine _1_hb _1_hdlcholesterol _1_sodium _1_ldlc _1_totalprotein _1_triglycerides _1_mcv _1_tsh _1_urea _1_neutrophilcount _1_wcc _1_lymphocytecount _1_potassium _1_hba1c _1_egfr _1_acr i._1_smoking_status _1_bmi _1_diastolicbloodpressure _1_systolicbloodpressure i.hx_neuro i.acarbose i.insulin i.statins i.hx_resp i.hx_liver i.hx_gi i.hx_severekidney i.hx_autoimmune i.thiazols i.hx_cancer i.hx_proteinuria i.hx_neuropathy i.hx_retinopathy i.hx_cad i.hx_cerebro i.hx_pvd i.hx_arrhythmia i.hx_heartfailure i.hx_mi i.hx_dvt i.hx_nephropathy i.hx_footlegulcer i.glp1 i.gliptins i.meglitinides i.sex duration_diabetes age i.aceiarb i.ca_blocker i.beta_blocker i.thiazides i.alpha_blocker i.metformin i.sulphonylureas i.hx_comashock, nohr /*fully observed (n=)*/
	estimates store model1
	/*Interaction Models with sex*/
foreach i in c.duration_diabetes c.age i.aceiarb i.ca_blocker i.beta_blocker i.thiazides i.alpha_blocker i.metformin i.sulphonylureas i.hx_comashock {
	
	quietly: stcox i.sex#`i' _1_albumin _1_alp _1_alt _1_bilirubin _1_totalcholesterol _1_serumcreatinine _1_hb _1_hdlcholesterol _1_sodium _1_ldlc _1_totalprotein _1_triglycerides _1_mcv _1_tsh _1_urea _1_neutrophilcount _1_wcc _1_lymphocytecount _1_potassium _1_hba1c _1_egfr _1_acr i._1_smoking_status _1_bmi _1_diastolicbloodpressure _1_systolicbloodpressure i.hx_neuro i.acarbose i.insulin i.statins i.hx_resp i.hx_liver i.hx_gi i.hx_severekidney i.hx_autoimmune i.thiazols i.hx_cancer i.hx_proteinuria i.hx_neuropathy i.hx_retinopathy i.hx_cad i.hx_cerebro i.hx_pvd i.hx_arrhythmia i.hx_heartfailure i.hx_mi i.hx_dvt i.hx_nephropathy i.hx_footlegulcer i.glp1 i.gliptins i.meglitinides duration_diabetes age i.aceiarb i.ca_blocker i.beta_blocker i.thiazides i.alpha_blocker i.metformin i.sulphonylureas i.hx_comashock i.sex, nohr
	display "`i'"
	estimates store model_2
	lrtest model_2 model1
}	

stcox  i.sex#i.hx_neuro i.sex#c.age i.sex#i.thiazides _1_albumin _1_alp _1_alt _1_bilirubin _1_totalcholesterol _1_serumcreatinine _1_hb _1_hdlcholesterol _1_sodium _1_ldlc _1_totalprotein _1_triglycerides _1_mcv _1_tsh _1_urea _1_neutrophilcount _1_wcc _1_lymphocytecount _1_potassium _1_hba1c _1_egfr _1_acr i._1_smoking_status _1_bmi _1_diastolicbloodpressure _1_systolicbloodpressure i.hx_neuro i.acarbose i.insulin i.statins i.hx_resp i.hx_liver i.hx_gi i.hx_severekidney i.hx_autoimmune i.thiazols i.hx_cancer i.hx_proteinuria i.hx_neuropathy i.hx_retinopathy i.hx_cad i.hx_cerebro i.hx_pvd i.hx_arrhythmia i.hx_heartfailure i.hx_mi i.hx_dvt i.hx_nephropathy i.hx_footlegulcer i.glp1 i.gliptins i.meglitinides duration_diabetes age i.aceiarb i.ca_blocker i.beta_blocker i.thiazides i.alpha_blocker i.metformin i.sulphonylureas i.hx_comashock i.sex, nohr
	estimates store model3

lrtest model3 model1 		


/*interaction terms with age*/
		
foreach i in c._1_albumin c._1_alp c._1_alt c._1_bilirubin c._1_totalcholesterol c._1_serumcreatinine c._1_hb c._1_hdlcholesterol c._1_sodium c._1_ldlc c._1_totalprotein c._1_triglycerides c._1_mcv c._1_tsh c._1_urea c._1_neutrophilcount c._1_wcc c._1_lymphocytecount c._1_potassium c._1_hba1c c._1_egfr c._1_acr i._1_smoking_status c._1_bmi c._1_diastolicbloodpressure c._1_systolicbloodpressure i.hx_neuro i.acarbose i.insulin i.statins i.hx_resp i.hx_liver i.hx_gi i.hx_severekidney i.hx_autoimmune i.thiazols i.hx_cancer i.hx_proteinuria i.hx_neuropathy i.hx_retinopathy i.hx_cad i.hx_cerebro i.hx_pvd i.hx_arrhythmia i.hx_heartfailure i.hx_mi i.hx_dvt i.hx_nephropathy i.hx_footlegulcer i.glp1 i.gliptins i.meglitinides c.duration_diabetes i.aceiarb i.ca_blocker i.beta_blocker i.thiazides i.alpha_blocker i.metformin i.sulphonylureas i.hx_comashock {
	
	quietly: stcox c.age#`i' _1_albumin _1_alp _1_alt _1_bilirubin _1_totalcholesterol _1_serumcreatinine _1_hb _1_hdlcholesterol _1_sodium _1_ldlc _1_totalprotein _1_triglycerides _1_mcv _1_tsh _1_urea _1_neutrophilcount _1_wcc _1_lymphocytecount _1_potassium _1_hba1c _1_egfr _1_acr i._1_smoking_status _1_bmi _1_diastolicbloodpressure _1_systolicbloodpressure i.hx_neuro i.acarbose i.insulin i.statins i.hx_resp i.hx_liver i.hx_gi i.hx_severekidney i.hx_autoimmune i.thiazols i.hx_cancer i.hx_proteinuria i.hx_neuropathy i.hx_retinopathy i.hx_cad i.hx_cerebro i.hx_pvd i.hx_arrhythmia i.hx_heartfailure i.hx_mi i.hx_dvt i.hx_nephropathy i.hx_footlegulcer i.glp1 i.gliptins i.meglitinides c.duration_diabetes i.aceiarb i.ca_blocker i.beta_blocker i.thiazides i.alpha_blocker i.metformin i.sulphonylureas i.hx_comashock i.sex c.age, nohr
	display "`i'"
	estimates store model_2
	lrtest model_2 model1
}

stcox  c.age#c._1_albumin c.age#c._1_alp c.age#c._1_bilirubin c.age#c._1_sodium c.age#c._1_mcv c.age#c._1_smoking_status c.age#c._1_bmi c.age#i.hx_resp c.age#i.hx_liver c.age#i.hx_severekidney c.age#i.hx_cancer c.age#i.hx_pvd _1_albumin _1_alp _1_alt _1_bilirubin _1_totalcholesterol _1_serumcreatinine _1_hb _1_hdlcholesterol _1_sodium _1_ldlc _1_totalprotein _1_triglycerides _1_mcv _1_tsh _1_urea _1_neutrophilcount _1_wcc _1_lymphocytecount _1_potassium _1_hba1c _1_egfr _1_acr i._1_smoking_status _1_bmi _1_diastolicbloodpressure _1_systolicbloodpressure i.hx_neuro i.acarbose i.insulin i.statins i.hx_resp i.hx_liver i.hx_gi i.hx_severekidney i.hx_autoimmune i.thiazols i.hx_cancer i.hx_proteinuria i.hx_neuropathy i.hx_retinopathy i.hx_cad i.hx_cerebro i.hx_pvd i.hx_arrhythmia i.hx_heartfailure i.hx_mi i.hx_dvt i.hx_nephropathy i.hx_footlegulcer i.glp1 i.gliptins i.meglitinides duration_diabetes age i.aceiarb i.ca_blocker i.beta_blocker i.thiazides i.alpha_blocker i.metformin i.sulphonylureas i.hx_comashock i.sex, nohr
	estimates store model3
lrtest model3 model1

