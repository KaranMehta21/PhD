*Set
	clear
	set more off
*cap log close*

*cap log close*
import delimited "S:\FPHS_THIN_Projects\Karan\BZI_PhD_T2D\PhD_descriptive_data\mice_data\nonlinearvariables_290124.csv", clear

save "S:\FPHS_THIN_Projects\Karan\BZI_PhD_T2D\PhD_descriptive_data\mice_data\nonlinearvariables_290124.dta", replace

*Import dataset
	import delimited "S:\FPHS_THIN_Projects\Karan\BZI_PhD_T2D\PhD_descriptive_data\mice_data\phd_mice_data_290124.csv", clear
	merge 1:1 patprac using "S:\FPHS_THIN_Projects\Karan\BZI_PhD_T2D\PhD_descriptive_data\mice_data\nonlinearvariables_290124.dta"
	
gen smoking_status2 = .
replace smoking_status2 = 1 if smoking_status==1 | smoking_status==2
replace smoking_status2 = 0 if smoking_status==0
tab smoking_status smoking_status2, missing
	
	
	
generate sex_neuro = sex*hx_neuro
generate sex_age = sex*age 
generate age_albumin = age*albumin
generate age_alp = age*alp
generate age_sodium = age*sodium
generate age_mcv = age*mcv
generate age_smoking = age*smoking_status2
generate age_bmi = age*bmi
generate age_resp = age*hx_resp
generate age_liver = age*hx_liver
generate age_severekidney = age*hx_severekidney
generate age_cancer = age*hx_cancer
generate age_pvd = age*hx_pvd
	
gen miss_albumin = 1 if missing(albumin)
gen miss_alp = 1 if missing(alp)
gen miss_alt = 1 if missing(alt)
gen miss_totalcholesterol = 1 if missing(totalcholesterol)
gen miss_serumcreatinine = 1 if missing(serumcreatinine)
gen miss_hb = 1 if missing(hb)
gen miss_hdlcholesterol = 1 if missing(hdlcholesterol)
gen miss_sodium = 1 if missing(sodium)
gen miss_triglycerides = 1 if missing(triglycerides)
gen miss_totalprotein = 1 if missing(totalprotein)
gen miss_mcv = 1 if missing(mcv)
gen miss_tsh = 1 if missing(tsh)	
gen miss_urea = 1 if missing(urea)
gen miss_neutrophilcount = 1 if missing(neutrophilcount)
gen miss_lymphocytecount = 1 if missing(lymphocytecount)
gen miss_potassium = 1 if missing(potassium)
gen miss_hba1c = 1 if missing(hba1c)
gen miss_egfr = 1 if missing(egfr)
gen miss_acr = 1 if missing(acr)
gen miss_bmi = 1 if missing(bmi)
gen miss_diastolicbloodpressure = 1 if missing(diastolicbloodpressure)
gen miss_systolicbloodpressure = 1 if missing(systolicbloodpressure)

replace _1_albumin_0 = . if (miss_albumin==1)
replace _1_albumin_1 = . if (miss_albumin==1)

replace _1_alp_0 = . if (miss_alp==1)
replace _1_alp_1 = . if (miss_alp==1)
replace _1_alp_2 = . if (miss_alp==1)
replace _1_alp_3 = . if (miss_alp==1)

replace _1_alt_0 = . if (miss_alt==1)
replace _1_alt_1 = . if (miss_alt==1)
replace _1_alt_2 = . if (miss_alt==1)


replace _1_totalcholesterol_0 = . if (miss_totalcholesterol==1)
replace _1_totalcholesterol_1 = . if (miss_totalcholesterol==1)

replace _1_serumcreatinine_0 = . if (miss_serumcreatinine==1)
replace _1_serumcreatinine_1 = . if (miss_serumcreatinine==1)
replace _1_serumcreatinine_2 = . if (miss_serumcreatinine==1)
replace _1_serumcreatinine_3 = . if (miss_serumcreatinine==1)

replace _1_hb_0 = . if (miss_hb==1)
replace _1_hb_1 = . if (miss_hb==1)
replace _1_hb_2 = . if (miss_hb==1)
replace _1_hb_3 = . if (miss_hb==1)

replace _1_hdlcholesterol_0 = . if (miss_hdlcholesterol==1)
replace _1_hdlcholesterol_1 = . if (miss_hdlcholesterol==1)
replace _1_hdlcholesterol_2 = . if (miss_hdlcholesterol==1)
replace _1_hdlcholesterol_3 = . if (miss_hdlcholesterol==1)

replace _1_sodium_0 = . if (miss_sodium==1)
replace _1_sodium_1 = . if (miss_sodium==1)
replace _1_sodium_2 = . if (miss_sodium==1)
replace _1_sodium_3 = . if (miss_sodium==1)


replace _1_totalprotein_0 = . if (miss_totalprotein==1)
replace _1_totalprotein_1 = . if (miss_totalprotein==1)
replace _1_totalprotein_2 = . if (miss_totalprotein==1)
replace _1_totalprotein_3 = . if (miss_totalprotein==1)

replace _1_triglycerides_0 = . if (miss_triglycerides==1)
replace _1_triglycerides_1 = . if (miss_triglycerides==1)
replace _1_triglycerides_2 = . if (miss_triglycerides==1)

replace _1_mcv_0 = . if (miss_mcv==1)
replace _1_mcv_1 = . if (miss_mcv==1)
replace _1_mcv_2 = . if (miss_mcv==1)
replace _1_mcv_3 = . if (miss_mcv==1)

replace _1_tsh_0 = . if (miss_tsh==1)
replace _1_tsh_1 = . if (miss_tsh==1)
replace _1_tsh_2 = . if (miss_tsh==1)

replace _1_urea_0 = . if (miss_urea==1)
replace _1_urea_1 = . if (miss_urea==1)
replace _1_urea_2 = . if (miss_urea==1)
replace _1_urea_3 = . if (miss_urea==1)

replace _1_neutrophilcount_0 = . if (miss_neutrophilcount==1)
replace _1_neutrophilcount_1 = . if (miss_neutrophilcount==1)
replace _1_neutrophilcount_2 = . if (miss_neutrophilcount==1)

replace _1_lymphocytecount_0 = . if (miss_lymphocytecount==1)
replace _1_lymphocytecount_1 = . if (miss_lymphocytecount==1)
replace _1_lymphocytecount_2 = . if (miss_lymphocytecount==1)
replace _1_lymphocytecount_3 = . if (miss_lymphocytecount==1)

replace _1_potassium_0 = . if (miss_potassium==1)
replace _1_potassium_1 = . if (miss_potassium==1)

replace _1_hba1c_0 = . if (miss_hba1c==1)
replace _1_hba1c_1 = . if (miss_hba1c==1)

replace _1_egfr_0 = . if (miss_egfr==1)
replace _1_egfr_1 = . if (miss_egfr==1)
replace _1_egfr_2 = . if (miss_egfr==1)
replace _1_egfr_3 = . if (miss_egfr==1)

replace _1_acr_0 = . if (miss_acr==1)
replace _1_acr_1 = . if (miss_acr==1)
replace _1_acr_2 = . if (miss_acr==1)
replace _1_acr_3 = . if (miss_acr==1)

replace _1_bmi_0 = . if (miss_bmi==1)
replace _1_bmi_1 = . if (miss_bmi==1)
replace _1_bmi_2 = . if (miss_bmi==1)
replace _1_bmi_3 = . if (miss_bmi==1)

replace _1_diastolicbloodpressure_0 = . if (miss_diastolicbloodpressure==1)
replace _1_diastolicbloodpressure_1 = . if (miss_diastolicbloodpressure==1)
replace _1_systolicbloodpressure_0 = . if (miss_systolicbloodpressure==1)
replace _1_systolicbloodpressure_1 = . if (miss_systolicbloodpressure==1)
replace _1_systolicbloodpressure_2 = . if (miss_systolicbloodpressure==1)
replace _1_systolicbloodpressure_3 = . if (miss_systolicbloodpressure==1)



*Estimation with CCA 
	stset exposure, failure(mortality_status)
	sts gen na=na /*GENERATE NELSON-AALEN ESTIMATOR*/


	
forvalues a = 0(1)3{
	foreach x in albumin alp alt totalcholesterol serumcreatinine hb hdlcholesterol sodium totalprotein triglycerides mcv tsh urea neutrophilcount lymphocytecount potassium hba1c egfr acr bmi diastolicbloodpressure systolicbloodpressure{
		cap rename _1_`x'_`a' `x'_`a'
	}
}	

		
*Imputation
	mi set wide //can use flong instead of wide//
	mi register regular sex_neuro sex_age age_resp age_liver age_severekidney age_cancer age_pvd hx_neuro acarbose insulin statins hx_resp hx_liver hx_gi hx_severekidney hx_autoimmune thiazols hx_cancer hx_proteinuria hx_neuropathy hx_retinopathy hx_cad hx_cerebro hx_pvd hx_arrhythmia hx_heartfailure hx_mi hx_dvt hx_nephropathy hx_footlegulcer glp1 gliptins sex aceiarb ca_blocker beta_blocker alpha_blocker metformin sulphonylureas hx_comashock duration_diabetes_0 duration_diabetes_1 age_0 age_1 na
	mi register imputed age_albumin age_alp age_sodium age_mcv age_smoking age_bmi smoking_status2 albumin_0 albumin_1 alp_0 alp_1 alp_2 alp_3 alt_0 alt_1 alt_2 totalcholesterol_0 totalcholesterol_1 serumcreatinine_0 serumcreatinine_1 serumcreatinine_2 serumcreatinine_3 hb_0 hb_1 hb_2 hb_3 hdlcholesterol_0 hdlcholesterol_1 hdlcholesterol_2 hdlcholesterol_3 sodium_0 sodium_1 sodium_2 sodium_3 totalprotein_0 totalprotein_1 totalprotein_2 totalprotein_3 triglycerides_0 triglycerides_1 triglycerides_2 mcv_0 mcv_1 mcv_2 mcv_3 tsh_0 tsh_1 tsh_2 urea_0 urea_1 urea_2 urea_3 neutrophilcount_0 neutrophilcount_1 neutrophilcount_2 lymphocytecount_0 lymphocytecount_1 lymphocytecount_2 lymphocytecount_3 potassium_0 potassium_1 hba1c_0 hba1c_1 egfr_0 egfr_1 egfr_2 egfr_3 acr_0 acr_1 acr_2 acr_3 bmi_0 bmi_1 bmi_2 bmi_3 diastolicbloodpressure_0 diastolicbloodpressure_1 systolicbloodpressure_0 systolicbloodpressure_1 systolicbloodpressure_2 systolicbloodpressure_3 alp albumin alt totalcholesterol serumcreatinine hb hdlcholesterol sodium totalprotein triglycerides mcv tsh urea neutrophilcount lymphocytecount potassium hba1c egfr acr bmi diastolicbloodpressure systolicbloodpressure
	mi imp chain (logit, augment omit(age_smoking)) smoking_status2 (regress, omit(alp)) alp_0 alp_1 alp_2 alp_3 (regress, omit(alp_0 alp_1 alp_2 alp_3)) alp (regress, omit(albumin)) albumin_0 albumin_1 (regress, omit(albumin_0 albumin_1)) albumin (regress, omit(alt)) alt_0 alt_1 alt_2 (regress, omit(alt_0 alt_1 alt_2)) alt (regress, omit(totalcholesterol)) totalcholesterol_0 totalcholesterol_1 (regress, omit(totalcholesterol_0 totalcholesterol_1)) totalcholesterol (regress, omit(serumcreatinine)) serumcreatinine_0 serumcreatinine_1 serumcreatinine_2 serumcreatinine_3 (regress, omit(serumcreatinine_0 serumcreatinine_1 serumcreatinine_2 serumcreatinine_3)) serumcreatinine (regress, omit(hb)) hb_0 hb_1 hb_2 hb_3 (regress, omit(hb_0 hb_1 hb_2 hb_3)) hb (regress, omit(hdlcholesterol)) hdlcholesterol_0 hdlcholesterol_1 hdlcholesterol_2 hdlcholesterol_3 (regress, omit(hdlcholesterol_0 hdlcholesterol_1 hdlcholesterol_2 hdlcholesterol_3)) hdlcholesterol (regress, omit(sodium)) sodium_0 sodium_1 sodium_2 sodium_3 (regress, omit(sodium_0 sodium_1 sodium_2 sodium_3)) sodium (regress, omit(triglycerides)) triglycerides_0 triglycerides_1 triglycerides_2 (regress, omit(triglycerides_0 triglycerides_1 triglycerides_2)) triglycerides (regress, omit(totalprotein)) totalprotein_0 totalprotein_1 totalprotein_2 totalprotein_3 (regress, omit(totalprotein_0 totalprotein_1 totalprotein_2 totalprotein_3)) totalprotein (regress, omit(mcv)) mcv_0 mcv_1 mcv_2 mcv_3 (regress, omit(mcv_0 mcv_1 mcv_2 mcv_3)) mcv (regress, omit(tsh)) tsh_0 tsh_1 tsh_2 (regress, omit(tsh_0 tsh_1 tsh_2)) tsh (regress, omit(urea)) urea_0 urea_1 urea_2 urea_3 (regress, omit(urea_0 urea_1 urea_2 urea_3)) urea (regress, omit(neutrophilcount)) neutrophilcount_0 neutrophilcount_1 neutrophilcount_2 (regress, omit(neutrophilcount_0 neutrophilcount_1 neutrophilcount_2)) neutrophilcount (regress, omit(lymphocytecount)) lymphocytecount_0 lymphocytecount_1 lymphocytecount_2 lymphocytecount_3 (regress, omit(lymphocytecount_0 lymphocytecount_1 lymphocytecount_2 lymphocytecount_3)) lymphocytecount (regress, omit(potassium)) potassium_0 potassium_1 (regress, omit(potassium_0 potassium_1)) potassium (regress, omit(hba1c)) hba1c_0 hba1c_1 (regress, omit(hba1c_0 hba1c_1)) hba1c (regress, omit(egfr)) egfr_0 egfr_1 egfr_2 egfr_3 (regress, omit(egfr_0 egfr_1 egfr_2 egfr_3)) egfr (regress, omit(acr)) acr_0 acr_1 acr_2 acr_3 (regress, omit(acr_0 acr_1 acr_2 acr_3)) acr (regress, omit(bmi)) bmi_0 bmi_1 bmi_2 bmi_3 (regress, omit(bmi_0 bmi_1 bmi_2 bmi_3)) bmi (regress, omit(diastolicbloodpressure)) diastolicbloodpressure_0 diastolicbloodpressure_1 (regress, omit(diastolicbloodpressure_0 diastolicbloodpressure_1)) diastolicbloodpressure (regress, omit(systolicbloodpressure)) systolicbloodpressure_0 systolicbloodpressure_1 systolicbloodpressure_2 systolicbloodpressure_3 (regress, omit(systolicbloodpressure_0 systolicbloodpressure_1 systolicbloodpressure_2 systolicbloodpressure_3)) systolicbloodpressure (regress) age_albumin age_alp age_sodium age_mcv age_smoking age_bmi = sex_neuro sex_age age_resp age_liver age_severekidney age_cancer age_pvd i.hx_neuro i.acarbose i.insulin i.statins i.hx_resp i.hx_liver i.hx_gi i.hx_severekidney i.hx_autoimmune i.thiazols i.hx_cancer i.hx_proteinuria i.hx_neuropathy  i.hx_retinopathy i.hx_cad i.hx_cerebro i.hx_pvd i.hx_arrhythmia i.hx_heartfailure i.hx_mi i.hx_dvt i.hx_nephropathy i.hx_footlegulcer i.glp1 i.gliptins i.sex i.aceiarb i.ca_blocker i.beta_blocker i.alpha_blocker i.metformin i.sulphonylureas i.hx_comashock duration_diabetes_0 duration_diabetes_1 age_0 age_1 na, replace burn(10) add(20) dots rseed(2984792) noisily //dryrun //noisily*/*for now*//*replace burn (10) instead of 100*/
				
	
*Saving file
	save "S:\FPHS_THIN_Projects\Karan\BZI_PhD_T2D\PhD_descriptive_data\mice_data\impv2_phddata_180224.dta", replace
	
	
log close
exit
