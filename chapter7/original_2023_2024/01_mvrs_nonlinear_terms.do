*Set
	clear
	set more off
	findit mvrs
*cap log close*
*Import dataset
	import delimited "S:/FPHS_THIN_Projects/Karan/BZI_PhD_T2D/PhD_descriptive_data/mice_data/1_imputeddataset_250823.csv"
	
/*Cox model with splines*/	
	stset exposure, failure(mortality_status)
	mvrs stcox _1_albumin _1_alp _1_alt _1_bilirubin _1_totalcholesterol _1_serumcreatinine _1_hb _1_hdlcholesterol _1_sodium _1_ldlc _1_totalprotein _1_triglycerides _1_mcv _1_tsh _1_urea _1_neutrophilcount _1_wcc _1_lymphocytecount _1_potassium _1_hba1c _1_egfr _1_acr _1_bmi _1_diastolicbloodpressure _1_systolicbloodpressure  duration_diabetes age, select(0.05) //identifies non-linear relationship and the optimal nodes/knots to capture non-linear relationship//to consider adding interaction terms 
	estimates store cox_splines //storing the AIC and BIC estimates of the final model 
	
	i._1_smoking_status i.hx_neuro i.acarbose i.insulin i.statins i.hx_resp i.hx_liver i.hx_gi i.hx_severekidney i.hx_autoimmune i.thiazols i.hx_cancer i.hx_proteinuria i.hx_neuropathy i.hx_retinopathy i.hx_cad i.hx_cerebro i.hx_pvd i.hx_arrhythmia i.hx_heartfailure i.hx_mi i.hx_dvt i.hx_nephropathy i.hx_footlegulcer i.glp1 i.gliptins i.meglitinides i.sex i.aceiarb i.ca_blocker i.beta_blocker i.thiazides i.alpha_blocker i.metformin i.sulphonylureas i.hx_comashock
	
	
	
