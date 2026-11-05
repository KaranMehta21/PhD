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

//*Duplicate 50% of cohort*
//gen random = runiform()
//sort random
//gen dummy = 0
//count
//local a = round(r(N)/2)
//replace dummy = 1 in 1/`a'
//expand 2 if dummy==1, gen(expand) 



*Estimation with CCA 
	stset exposure, failure(mortality_status)
	sts gen na=na /*GENERATE NELSON-AALEN ESTIMATOR*/


*Transformation, log/reciprocal/ be prudent with transformation. don't overdo it. 

//qladder acr /*best is */
	//ladder acr
	//qladder acr
	//replace acr=sqrt(acr)



	//qladder bmi /*best is */
	//ladder bmi
	//replace bmi=log(bmi)
	
	//qladder albumin
	//hist albumin
	//replace albumin=(albumin)^2
	
	*Transformation
	//qladder alp /*best is */
	//hist alp
	//replace alp=1/sqrt(alp)
	
	*Transformation
	//qladder alt /*best is */
	//ladder systolic_bp
	//hist alt
	//replace alt=log(alt)
	
	//qladder bilirubin
	//hist bilirubin
	//replace bilirubin=(bilirubin)^2
	
	
	
*Transformation
	//qladder totalcholesterol /*best is */
	//ladder serum_cholesterol
	//hist totalcholesterol
	//replace totalcholesterol=1/sqrt(totalcholesterol)
	
	
	
	*Transformation
	//qladder serumcreatinine /*best is */
	//ladder systolic_bp
	//hist serumcreatinine
	//replace serumcreatinine=log(serumcreatinine)
	
	*Transformation
	//qladder hb /*best is */
	//ladder systolic_bp
	//hist hb
	//replace hb=(hb)^2
	
	*Transformation
	//qladder hdlcholesterol /*best is */
	//ladder hdl_chol
	//hist hdlcholesterol
	//replace hdlcholesterol=sqrt(hdlcholesterol)
	
	
*Transformation
	//qladder sodium /*best is */
	//ladder systolic_bp
	//hist sodium
	//replace sodium=1/sqrt(sodium)
	
		
*Transformation
	//qladder ldlc /*best is */
	//ladder non_hdlchol
	//hist ldlc
	//replace ldlc=sqrt(ldlc)
	
	*Transformation
	//qladder totalprotein /*best is */
	//hist totalprotein
	//replace totalprotein=(totalprotein)^2
	
	*Transformation
	//qladder triglycerides /*best is */
	//ladder triglycerides
	//hist triglycerides
	//replace triglycerides=log(triglycerides)
	//call log version another name 
	
	
*Transformation
	//qladder mcv /*best is */
	//ladder systolic_bp
	//hist mcv 
	//replace mcv=(mcv)^2
	
	*Transformation
	//qladder tsh /*best is */
	//ladder systolic_bp
	//hist tsh 
	//replace tsh=log(tsh)
	

*Transformation
	//qladder urea /*best is */
	//ladder systolic_bp
	//hist urea
	//replace urea=log(urea)
	
	
	*Transformation
	//qladder neutrophilcount /*best is */
	//ladder systolic_bp
	//hist neutrophilcount
	//replace neutrophilcount=1/sqrt(neutrophilcount)
	
	//qladder wcc /*best is */
	//ladder systolic_bp
	//hist wcc
	//replace wcc=log(wcc)

*Transformation
	//qladder lymphocytecount /*best is */
	//ladder systolic_bp
	//hist lymphocytecount
	//replace lymphocytecount=log(lymphocytecount)
	
	*Transformation
	//qladder potassium /*best is */
	//ladder systolic_bp
	//hist potassium
	//replace potassium=1/sqrt(potassium)
	
	*Transformation
	//qladder hba1c /*best is */
	//ladder systolic_bp
	//hist hba1c
	//replace hba1c=log(hba1c)
	
	*Transformation
	//qladder egfr /*best is */
	//ladder systolic_bp
	//hist egfr
	//replace egfr=(egfr)^2
	
	*Transformation
	//qladder diastolicbloodpressure /*best is */
	//ladder diastolic_bp
	//hist diastolicbloodpressure
	//replace diastolicbloodpressure=sqrt(diastolicbloodpressure)

*Transformation
	//qladder systolicbloodpressure /*best is */
	//ladder systolic_bp
	//hist systolicbloodpressure
	//replace systolicbloodpressure=sqrt(systolicbloodpressure)
	
forvalues a = 0(1)3{
	foreach x in albumin alp alt totalcholesterol serumcreatinine hb hdlcholesterol sodium totalprotein triglycerides mcv tsh urea neutrophilcount lymphocytecount potassium hba1c egfr acr bmi diastolicbloodpressure systolicbloodpressure{
		cap rename _1_`x'_`a' `x'_`a'
	}
}	



		
*Imputation
	mi set wide //can use flong instead of wide//
	mi register regular sex_neuro sex_age age_resp age_liver age_severekidney age_cancer age_pvd hx_neuro acarbose insulin statins hx_resp hx_liver hx_gi hx_severekidney hx_autoimmune thiazols hx_cancer hx_proteinuria hx_neuropathy hx_retinopathy hx_cad hx_cerebro hx_pvd hx_arrhythmia hx_heartfailure hx_mi hx_dvt hx_nephropathy hx_footlegulcer glp1 gliptins sex aceiarb ca_blocker beta_blocker alpha_blocker metformin sulphonylureas hx_comashock duration_diabetes_0 duration_diabetes_1 age_0 age_1 na
	mi register imputed age_albumin age_alp age_sodium age_mcv age_smoking age_bmi smoking_status2 albumin_0 albumin_1 alp_0 alp_1 alp_2 alp_3 alt_0 alt_1 alt_2 totalcholesterol_0 totalcholesterol_1 serumcreatinine_0 serumcreatinine_1 serumcreatinine_2 serumcreatinine_3 hb_0 hb_1 hb_2 hb_3 hdlcholesterol_0 hdlcholesterol_1 hdlcholesterol_2 hdlcholesterol_3 sodium_0 sodium_1 sodium_2 sodium_3 totalprotein_0 totalprotein_1 totalprotein_2 totalprotein_3 triglycerides_0 triglycerides_1 triglycerides_2 mcv_0 mcv_1 mcv_2 mcv_3 tsh_0 tsh_1 tsh_2 urea_0 urea_1 urea_2 urea_3 neutrophilcount_0 neutrophilcount_1 neutrophilcount_2 lymphocytecount_0 lymphocytecount_1 lymphocytecount_2 lymphocytecount_3 potassium_0 potassium_1 hba1c_0 hba1c_1 egfr_0 egfr_1 egfr_2 egfr_3 acr_0 acr_1 acr_2 acr_3 bmi_0 bmi_1 bmi_2 bmi_3 diastolicbloodpressure_0 diastolicbloodpressure_1 systolicbloodpressure_0 systolicbloodpressure_1 systolicbloodpressure_2 systolicbloodpressure_3 albumin alp alt totalcholesterol serumcreatinine hb hdlcholesterol sodium triglycerides totalprotein mcv tsh urea neutrophilcount lymphocytecount potassium hba1c egfr acr bmi diastolicbloodpressure systolicbloodpressure
	mi imp chain (logit, augment omit(age_smoking)) smoking_status2 (regress) albumin alp alt totalcholesterol serumcreatinine hb hdlcholesterol sodium triglycerides totalprotein mcv tsh urea neutrophilcount lymphocytecount potassium hba1c egfr acr bmi diastolicbloodpressure systolicbloodpressure age_albumin age_alp age_sodium age_mcv age_smoking age_bmi albumin_0 albumin_1 alp_0 alp_1 alp_2 alp_3 alt_0 alt_1 alt_2 totalcholesterol_0 totalcholesterol_1 serumcreatinine_0 serumcreatinine_1 serumcreatinine_2 serumcreatinine_3 hb_0 hb_1 hb_2 hb_3 hdlcholesterol_0 hdlcholesterol_1 hdlcholesterol_2 hdlcholesterol_3 sodium_0 sodium_1 sodium_2 sodium_3 totalprotein_0 totalprotein_1 totalprotein_2 totalprotein_3 triglycerides_0 triglycerides_1 triglycerides_2 mcv_0 mcv_1 mcv_2 mcv_3 tsh_0 tsh_1 tsh_2 urea_0 urea_1 urea_2 urea_3 neutrophilcount_0 neutrophilcount_1 neutrophilcount_2 lymphocytecount_0 lymphocytecount_1 lymphocytecount_2 lymphocytecount_3 potassium_0 potassium_1 hba1c_0 hba1c_1 egfr_0 egfr_1 egfr_2 egfr_3 acr_0 acr_1 acr_2 acr_3 bmi_0 bmi_1 bmi_2 bmi_3 diastolicbloodpressure_0 diastolicbloodpressure_1 systolicbloodpressure_0 systolicbloodpressure_1 systolicbloodpressure_2 systolicbloodpressure_3  = sex_neuro sex_age age_resp age_liver age_severekidney age_cancer age_pvd i.hx_neuro i.acarbose i.insulin i.statins i.hx_resp i.hx_liver i.hx_gi i.hx_severekidney i.hx_autoimmune i.thiazols i.hx_cancer i.hx_proteinuria i.hx_neuropathy  i.hx_retinopathy i.hx_cad i.hx_cerebro i.hx_pvd i.hx_arrhythmia i.hx_heartfailure i.hx_mi i.hx_dvt i.hx_nephropathy i.hx_footlegulcer i.glp1 i.gliptins i.sex i.aceiarb i.ca_blocker i.beta_blocker i.alpha_blocker i.metformin i.sulphonylureas i.hx_comashock duration_diabetes_0 duration_diabetes_1 age_0 age_1 na, replace burn(10) add(10) dots rseed(2984792) noisily //dryrun //noisily*/*for now*//*replace burn (10) instead of 100*/
		
		//Use MC error//
logit smoking_status2 age_albumin age_alp age_bilirubin age_sodium age_mcv age_bmi _1_albumin_0 _1_albumin_1 _1_alp_0 _1_alp_1 _1_alp_2 _1_alp_3 _1_alt_0 _1_alt_1 _1_alt_2 _1_bilirubin_0 _1_bilirubin_1 _1_bilirubin_2 _1_bilirubin_3 _1_totalcholesterol_0 _1_totalcholesterol_1 _1_serumcreatinine_0 _1_serumcreatinine_1 _1_serumcreatinine_2 _1_serumcreatinine_3 _1_hb_0 _1_hb_1 _1_hb_2 _1_hb_3 _1_hdlcholesterol_0 _1_hdlcholesterol_1 _1_hdlcholesterol_2 _1_hdlcholesterol_3 _1_sodium_0 _1_sodium_1 _1_sodium_2 _1_sodium_3 _1_ldlc_0 _1_ldlc_1 _1_totalprotein_0 _1_totalprotein_1 _1_totalprotein_2 _1_totalprotein_3 _1_triglycerides_0 _1_triglycerides_1 _1_triglycerides_2 _1_mcv_0 _1_mcv_1 _1_mcv_2 _1_mcv_3 _1_tsh_0 _1_tsh_1 _1_tsh_2 _1_urea_0 _1_urea_1 _1_urea_2 _1_urea_3 _1_neutrophilcount_0 _1_neutrophilcount_1 _1_neutrophilcount_2 _1_wcc_0 _1_wcc_1 _1_wcc_2 _1_lymphocytecount_0 _1_lymphocytecount_1 _1_lymphocytecount_2 _1_lymphocytecount_3 _1_potassium_0 _1_potassium_1 _1_hba1c_0 _1_hba1c_1 _1_egfr_0 _1_egfr_1 _1_egfr_2 _1_egfr_3 _1_acr_0 _1_acr_1 _1_acr_2 _1_acr_3 _1_bmi_0 _1_bmi_1 _1_bmi_2 _1_bmi_3 _1_diastolicbloodpressure_0 _1_diastolicbloodpressure_1 _1_systolicbloodpressure_0 _1_systolicbloodpressure_1 _1_systolicbloodpressure_2 _1_systolicbloodpressure_3 sex_neuro sex_thiazides sex_age age_resp age_liver age_severekidney age_cancer age_pvd i.hx_neuro i.acarbose i.insulin i.statins i.hx_resp i.hx_liver i.hx_gi i.hx_severekidney i.hx_autoimmune i.thiazols i.hx_cancer i.hx_proteinuria i.hx_neuropathy  i.hx_retinopathy i.hx_cad i.hx_cerebro i.hx_pvd i.hx_arrhythmia i.hx_heartfailure i.hx_mi i.hx_dvt i.hx_nephropathy i.hx_footlegulcer i.glp1 i.gliptins i.meglitinides i.sex duration_diabetes age i.aceiarb i.ca_blocker i.beta_blocker i.thiazides i.alpha_blocker i.metformin i.sulphonylureas i.hx_comashock duration_diabetes_0 duration_diabetes_1 age_0 age_1 na 

		
		
		
		
		
*Back transformation bmi
	forvalues i=1(1)20{
		replace _`i'_bmi=exp(_`i'_bmi) /*back transformation*/
	}
	replace bmi=exp(bmi) /*back transformation*/
	

*Back transformation alt
	forvalues i=1(1)20{
		replace _`i'_alt=exp(_`i'_alt) /*back transformation*/
	}
	replace alt=exp(alt) /*back transformation*/
	
*Back transformation serum creatinine
	forvalues i=1(1)20{
		replace _`i'_serumcreatinine=exp(_`i'_serumcreatinine) /*back transformation*/
	}
	replace serumcreatinine=exp(serumcreatinine) /*back transformation*/


*Back transformation triglcyerides
	forvalues i=1(1)20{
		replace _`i'_triglycerides=exp(_`i'_triglycerides) /*back transformation*/
	}
	replace triglycerides=exp(triglycerides) /*back transformation*/
	//replace trig = exp(logtrig) if missing(trig)//

	
*Back transformation tsh
	forvalues i=1(1)20{
		replace _`i'_tsh=exp(_`i'_tsh) /*back transformation*/
	}
	replace tsh=exp(tsh) /*back transformation*/


*Back transformation hba1c
	forvalues i=1(1)20{
		replace _`i'_hba1c=exp(_`i'_hba1c) /*back transformation*/
	}
	replace hba1c=exp(hba1c) /*back transformation*/

/*Check MC Error*/

		mi estimate, mcerror: stcox acr albumin alp alt bilirubin totalcholesterol serumcreatinine hb hdlcholesterol sodium ldlc totalprotein triglycerides mcv tsh urea neutrophilcount wcc lymphocytecount potassium hba1c egfr i.smoking_status bmi diastolicbloodpressure systolicbloodpressure /* 
		*/i.hx_neuro i.acarbose i.insulin i.statins i.hx_cad i.hx_cerebro i.hx_pvd i.hx_arrhythmia i.hx_heartfailure i.hx_mi i.hx_resp i.hx_liver i.hx_gi i.hx_severekidney i.hx_autoimmune i.thiazols i.hx_cancer i.hx_proteinuria i.hx_neuropathy i.hx_retinopathy i.hx_nephropathy i.hx_footlegulcer i.glp1 i.gliptins i.meglitinides i.sex duration_diabetes age i.aceiarb i.ca_blocker i.beta_blocker i.thiazides i.alpha_blocker i.metformin i.sulphonylureas i.hx_comashock, nohr /*fully observed (n=)*/
	
//local vars x1 x2 x3 x4
foreach var of local vars {
gen miss`var' = missing(`var')
hist `var'_1, by(miss`var', col(1)) name(fig`var', replace)
}	//
	
	
//hist acr
//hist acr_1 if missing(acr)//
//gen macr = missing(acr// this line and below
//hist acr_1, by(macr, col(1))//	
	
	
	
*Saving file
	save "S:\FPHS_THIN_Projects\Karan\BZI_PhD_T2D\PhD_descriptive_data\mice_data\impv2_phddata_310124.dta", replace
	
	
log close
exit
