***GREAT MIGRATION AND MOBILITY MASTER DO-FILE ***

*0. Set preconditions and directory definitions.
*--------------------------------------------
* Clear all and set max matsize. *
	capture log close
	clear all
	set more off
	set maxvar 30000
	
global XXX 		"C:\Users\sebsg\Documents\Para_la_U\Octavo\Urbana\Proyecto\Migration\replication_AER"

* Action required: Install packages 

						ssc install estout, replace
						ssc install maptile, replace
						ssc install spmap, replace
						ssc install shp2dta, replace
						ssc install parmest, replace
						ssc install ivreg2, replace
						ssc install ranktest, replace
						ssc install statastates, replace
						ssc install mdesc, replace
						ssc install coefplot, replace
						ssc install rsource, replace
						ssc install binscatter, replace
						ssc install keeporder, replace
						ssc install lincomest, replace
						ssc install egenmore, replace
						ssc install distinct, replace
						ssc install unique, replace 
						
						maptile_install using "http://files.michaelstepner.com/geo_cz1990.zip", replace
						
* Code dir and sub-dirs. *	
    global code "${XXX}/code"
	global lasso "$code/lasso"
	global bartik "$code/bartik"
	
* Data dir and sub-dirs. *
	global data "$XXX/data"
	global xwalks "$data/crosswalks"
	global urbrural "$xwalks/documentation/urban_rural_county_classification"
	global msanecma "$xwalks/documentation/msanecma_1999_codes_names"
	global city_sample "$data/city_sample"
	global mobdata "$data/mobility"
	global instrument "$data/instrument"
	global migshares "$instrument/shares"
	global migdata "$instrument/migration"
	global mechanisms "$data/mechanisms"
	global jobs "$mechanisms/jobs"
	global pf "$mechanisms/public_finance"
	global political "$mechanisms/political"
	global nbhds "$mechanisms/neighborhoods"
	global incarceration "$mechanisms/incarceration"
	global schools "$mechanisms/schools"
	global population "$mechanisms/population"
	global ri "$data/randomization_inference"	
	
*1. Build intermediate datasets*

***1.1 Build data/crosswalks data
	
	cd "$xwalks/documentation/urban_rural_county_classification"
	do clean_census_tract_2010_rural_urban_classification.do 

	cd "$msanecma"
	do msa_1990_codes_names.do 

	cd "$msanecma"
	do necma_1990_codes_names.do 

	cd "$urbrural"
	do urbrural_county.do 	

	cd "$xwalks/documentation"
	do county1940_crosswalks.do 

	cd "$xwalks/documentation"	
	do US_place_point_2010_crosswalks.do
	
***1.2. Build data/instrument data.

*	cd "$migdata/documentation"
*	do clean_south_county_white_nonwhite_mig_1940_1970.do /base de datos no disponible
	cd "$migshares/documentation"
	do clean_IPUMS_1935_1940_extract_to_construct_migration_weights.do
	cd "$migshares/documentation"
	do clean_IPUMS_1940_extract_to_construct_migration_weights.do
	cd "$instrument/rndmig"
	unzipfile instrument/rndmig/rndmig.zip
	cd "$migshares/rndmig"
	unzipfile migshares/rndmig/shares_rndmig.zip
	
***1.3. Build data/mechanisms/population data.

*	cd "$population/raw/documentation"
*	do interpolatedpopulations.do /data no disponible
	cd "$population/documentation"
	do clean_city_population_by_race_2000.do
*	cd "$population/documentation"
*	do clean_city_population_ccdb_1944_1977.do /data no disponible
	cd "$population/documentation"
	do clean_city_population_census_1940.do
*	cd "$population/documentation"
*	do clean_county_pop_1870_2000_bpopshare_1870_2000_cz.do /data no disponible
*	cd "$population/documentation"
*	do clean_cz_population_1940_1970.do /data no disponible
*	cd "$population/documentation"
*	do clean_cz_population_density_1940.do /data no disponible
*	cd "$population/documentation"
*	do clean_cz_snq_european_immigration_instrument.do /data no disponible
*	cd "$population/documentation"
*	do City_Book_1944_1977_clean.do /data no disponible
	
***1.4. Build data/mobility datasets.

	cd "$mobdata/documentation"
	do clean_cz_mobility_1900_2015.do
	cd "$mobdata/documentation"
	do clean_cz_p25_p75_race_share.do
	cd "$mobdata/documentation"
	do mobdata/documentation/clean_teen_school_attendance_mother_south_north.do	
	
***1.5. Build mechanisms/incarceration datasets.

	cd "$incarceration/documentation"
	do clean_cz_city_murder_rates_1931_1969.do
	cd "$incarceration/documentation"
	cap do clean_cz_riots_1964_1971.do
	cd "$incarceration/documentation"
	do clean_cz_jail_rates_1920_1960.do
	cd "$incarceration/documentation"
	do clean_cz_iob_crime_incarceration.do
	

***1.6. Build mechanisms/jobs datasets.

	cd "$jobs/documentation"
	do clean_cz_industry_employment_1940_1970.do

	
***1.7. Build mechanisms/neighborhoods datasets.

	cd "$nbhds/documentation"
	do clean_cz_marriage_income_occscore_1940.do
	do nbhds/documentation/clean_cz_neighborhoods.do
	
	
***1.8. Build mechanisms/political datasets.

	cd "$political/documentation"
	do clean_cz_wallace_share_1968.do		
	do political/documentation/clean_cz_weighted_racial_animus.do

	
***1.9. Build mechanisms/public_finance datasets.

	cd "$pf/raw/documentation"
	do CountyAreaFin_Eight_Var_InSampleCZ_1957_2002.do
	cd "$pf/documentation"
	do clean_cz_city_police_per_capita_1920_2007.do
	do pf/documentation/clean_cz_public_finance_1932_2012.do
	
	
***1.10. Build mechanisms/schools datasets.

	cd "$schools/documentation"
	do clean_cz_med_educd_25plus_1940.do
	do schools/documentation/clean_cz_prvschl_share_1920_2010.do
		
*2. Predict netmigration from the South between 1940 and 1970.
	
***2.1. Create and clean the migration datasets for prediction for each decade.


	/* Data on black netmigration for southern counties come from Boustan (2016):
	south_county.dta. These data were downloaded from the following link: 
	https://economics.princeton.edu/dl/Boustan/Chapter4.zip. */
	cd "$migdata/raw"
	use south_county.dta, clear
	drop if netbmig==.
	
	/* Instructions for cleaning the data from Boustan (2016) replication files
	are prefaced with "Boustan (2016)".
	
	Boustan (2016): This data set includes all of the southern data by county, from 
	CCDB and ICPSR Great Plains project.
	
	Boustan (2016): There are 350 or so counties with missing mining or 
	manufacturing information in 1950 & 1970. In this case, replace with the 
	1960 info. */
	
	sort state countyicp year
	replace permin=permin[_n+1] if year==1950 & permin==. & countyicp==countyicp[_n+1]
	replace permin=permin[_n-1] if year==1970 & permin==. & countyicp==countyicp[_n-1]
	replace perman=perman[_n+1] if year==1950 & perman==. & countyicp==countyicp[_n+1]
	replace perman=perman[_n-1] if year==1970 & perman==. & countyicp==countyicp[_n-1]
	
	/* Note that migration data is missing for several counties in Virginia. */
	
	/* Boustan (2016): Interact variables with % cotton and % agriculture. */
	
	replace perten=perten/100
	replace perag=perag/100 if year==1950
	replace perag=perag/10 if year==1960
	
	/* Boustan (2016): Create dummy for SA (GA, FL, VA, WV, SC, NC) and interact. */
	
	gen satl=(state==12 | state==13 | state==37 | state==45 | state==51 | state==54)
	gen pertensa=perten*satl
	gen permansa=perman*satl
	
	/* Boustan (2016): Create dummy for tobacco growing states and interact with 
	agriculture (NC, KY, TN). */
	
	gen tob=(state==37 | state==21 | state==47)
	gen peragtob=perag*tob
	
	/* Boustan (2016): Create dummy for mineral region (OK, TX). */
	
	gen ot=(state==40 | state==48)
	gen perminot=permin*ot
	
	/* Save the dataset that will be used for post LASSO */
	
	cd "$migdata"
	save clean_south_county.dta, replace
	
	/* Additional cleaning to prepare the data for R and running LASSO. */
	*local final_varlist netbmig perten perag permin perman aaa_pc warfac_pc warcon_pc avesz3 tmpav30 pcpav30 mxsw30s mxsd30s dustbowa summit swamp valley elevmax riv1120 riv2150 riv51up riv0510 elevrang awc clay kffact om perm thick minem35a Astate_5 Astate_12 Astate_13 Astate_21 Astate_22 Astate_28 Astate_37 Astate_40 Astate_45 Astate_47 Astate_48 Astate_51 Astate_54 satl pertensa permansa tob peragtob ot perminot
	local final_varlist netbmig percot perten perag peragtob tob warfac_pc permin perminot ot
	*perten perag permin perman warfac_pc satl pertensa permansa tob peragtob ot perminot
	
	/* Replace any vars that are still missing with the state-year mean value 
	for that var. */
	foreach var in `final_varlist'{
	mdesc `var'
	tab countyfips if `var'==.
	drop if `var'==.
	egen mean_`var'=mean(`var'), by(state year)
	replace `var'=mean_`var' if `var'==.
	}
	
	keep year `final_varlist'
	
	/* Create a separate dataset for each decade. */
	preserve
	cd "$lasso"
	keep if year==1950
	drop year
	save south_county_migration_dataset_for_prediction_1950.dta, replace
	restore
	
	preserve
	cd "$lasso"
	keep if year==1960
	drop year
	save south_county_migration_dataset_for_prediction_1960.dta, replace
	restore
	
	preserve
	cd "$lasso"
	keep if year==1970
	drop year
	save south_county_migration_dataset_for_prediction_1970.dta, replace
	restore
	
	clear
	
***2.2. Run lasso on each decade's migration dataset to obtain predictors.


	/* Initiate R and run LASSO using cv glmnet. */
	
/*
	global Rterm_path `"/usr/local/bin/r"'
	
	rsource, terminator(END_OF_R) roptions(--vanilla)
	
		// 'haven' is an R package for importing Stata '.dta' file
		library(haven)
		library(ggplot2)
		library(dplyr)
		library(tidyr)
		library(data.table)
		library(stargazer)
		library(randomForest)
		library(glmnet)
		library(rpart)
		library(parallel)
		library(stringr)
		
		// 0. Clear all
		rm(list = ls()) 
	
		// 1. Change to your directory
		setwd("/Users/elloraderenoncourt/Great_Migration_Mobility/code/lasso")
		
		// 2. Set code to run or not
		runLasso = FALSE
		
		if (runLasso) {
		// 3. Load the data, run lasso, get list of selected variables. Repeat for each year (1950, 1960, 1970)
		train = read_dta('south_county_migration_dataset_for_prediction_1950.dta')
	
		x = model.matrix(~., data=train %>% select(-netbmig))
		y = train$netbmig
		lassoPred=cv.glmnet(x=x, y=y,alpha=1,nfolds=5,standardize=TRUE)
		tmp_coeffs<-coef(lassoPred, s="lambda.min")
		lasso_list_of_vars_1950<-data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
	
	write.csv(lasso_list_of_vars_1950[2:dim(lasso_list_of_vars_1950)[1],1:dim(lasso_list_of_vars_1950)[2]], file='lasso_list_of_vars_1950.csv')
	
		train = read_dta('south_county_migration_dataset_for_prediction_1960.dta')
	
		x = model.matrix(~., data=train %>% select(-netbmig))
		y = train$netbmig
		lassoPred=cv.glmnet(x=x, y=y,alpha=1,nfolds=5,standardize=TRUE)
		tmp_coeffs<-coef(lassoPred, s="lambda.min")
		lasso_list_of_vars_1960<-data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
	
	write.csv(lasso_list_of_vars_1960[2:dim(lasso_list_of_vars_1960)[1],1:dim(lasso_list_of_vars_1960)[2]], file='lasso_list_of_vars_1960.csv')
	
		train = read_dta('south_county_migration_dataset_for_prediction_1970.dta')
	
		x = model.matrix(~., data=train %>% select(-netbmig))
		y = train$netbmig
		lassoPred=cv.glmnet(x=x, y=y,alpha=1,nfolds=5,standardize=TRUE)
		tmp_coeffs<-coef(lassoPred, s="lambda.min")
		lasso_list_of_vars_1970<-data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
	
	write.csv(lasso_list_of_vars_1970[2:dim(lasso_list_of_vars_1970)[1],1:dim(lasso_list_of_vars_1970)[2]], file='lasso_list_of_vars_1970.csv')
	}
	END_OF_R
*/
	
***2.3. Predict using original Boustan (2016) variables.

	/* Load full clean data. */
	
	cd "$migdata"
	use clean_south_county.dta, clear
	
	/* Predict county-level net migration rate, decade by decade with southern 
	variables chosen by LASSO. Predict net migration rate ("netbmig_pred") based on 
	these vars alone. */
	reg netbmig percot perten perag peragtob tob warfac_pc permin perminot ot if year==1950
	predict netbmig_pred if year==1950
	reg netbmig percot perten perag peragtob tob warfac_pc permin perminot ot if year==1960
	predict netbmig_pred01 if year==1960
	reg netbmig percot perten perag peragtob tob warfac_pc permin perminot ot if year==1970
	predict netbmig_pred02 if year==1970	
	
	replace netbmig_pred=netbmig_pred01 if year==1960
	replace netbmig_pred=netbmig_pred02 if year==1970
	drop netbmig_pred01 netbmig_pred02
	
	/* Boustan (2016): Total number leaving/coming to county: actual and predicted. Note 
	that netbmig is a migration rate (per 100 residents). So, the range is -100 to 
	+whatever. -100 because it is impossible for more than all of the residents to 
	leave. But, on the positive side, the rate is unrestricted, because the growth 
	could be quite high (for a county with 100 blacks in 1940, could have 100,000 
	blacks in 1950 which would be a rate of 1000. */ 
	
	gen totbmig=((bpop_l/100)*netbmig)
	gen totbmig_pred=((bpop_l/100)*netbmig_pred)
	gen weight=netbmig_pred*bpop_l
	
	/* One observation per county, year. */
	
	drop if year==year[_n-1]
	sort countyfips year
	drop if countyfips==.
	rename totbmig actoutmig
	rename totbmig_pred proutmig
	label var proutmig "predicted out migration, by county-year, south"
	drop _merge
	
	/* Merge with 1940 crosswalks data file. */
	
	cd "$xwalks"
	merge m:1 stateicp countyicp using county1940_crosswalks.dta, keepusing(fips state_name county_name)
	drop if _merge==2 
	g origin_fips=fips
	rename state_name origin_state_name
	rename county_name origin_county_name 
	
	/* Hand correct counties that didn't match using crosswalk file and internet search. */
	
	replace origin_fips = 51067 if countyfips==51620 & _merge==1
	replace origin_fips = 48203 if countyfips==48203 & _merge==1
	replace origin_fips = 51037 if countyfips==54039 & _merge==1
	replace origin_fips = 54041 if countyfips==54041 & _merge==1
	replace origin_fips = 51189 if countyfips==189 & _merge==1
	drop _merge
	
	tostring origin_fips, replace
	keep origin_fips year proutmig actoutmig netbmig_pred
	
	drop if netbmig_pred==. | proutmig==.
	
	bysort origin_fips year: gen dup= cond(_N==1,0,_n)
	tab dup
	
	drop dup
	
	cd "$instrument"
	save 1_boustan_predict_mig.dta, replace
	
	collapse (sum) netbmig_pred actoutmig proutmig, by(origin_fips year)
	save "$instrument/1_boustan_predict_mig_collapsed.dta", replace
	
	
***2.4. Run Post-LASSO to generate predicted migration figures for each county by decade.

	/* Load full clean data. */
	
	cd "$migdata"
	use clean_south_county.dta, clear
	
	/* Predict county-level net migration rate, decade by decade with southern 
	variables chosen by LASSO. Predict net migration rate ("netbmig_pred") based on 
	these vars alone. */

	reg netbmig perten perag warfac_pc percot peragtob ot perminot if year==1950
	predict netbmig_pred if year==1950
	reg netbmig percot perten perag peragtob tob warfac_pc permin perminot ot if year==1960
	predict netbmig_pred01 if year==1960
	reg netbmig percot perten perag peragtob tob warfac_pc permin perminot ot if year==1970
	predict netbmig_pred02 if year==1970	
	
	replace netbmig_pred=netbmig_pred01 if year==1960
	replace netbmig_pred=netbmig_pred02 if year==1970
	drop netbmig_pred01 netbmig_pred02
	
	/* Boustan (2016): Total number leaving/coming to county: actual and predicted. Note 
	that netbmig is a migration rate (per 100 residents). So, the range is -100 to 
	+whatever. -100 because it is impossible for more than all of the residents to 
	leave. But, on the positive side, the rate is unrestricted, because the growth 
	could be quite high (for a county with 100 blacks in 1940, could have 100,000 
	blacks in 1950 which would be a rate of 1000. */ 
	
	gen totbmig=((bpop_l/100)*netbmig)
	gen totbmig_pred=((bpop_l/100)*netbmig_pred)
	gen weight=netbmig_pred*bpop_l
	
	/* One observation per county, year. */
	
	drop if year==year[_n-1]
	sort countyfips year
	drop if countyfips==.
	rename totbmig actoutmig
	rename totbmig_pred proutmig
	label var proutmig "predicted out migration, by county-year, south"
	drop _merge
	
	/* Merge with 1940 crosswalks data file. */

	/*  Two methods for achieving consistent fips codes between the migration data, historical census data, 
	and the crosswalks file created for this project. Using county icp and state icp to match to the 
	crosswalk file yields the best results. Then one can merge the data with the migration weights from 
	census extract located here: data/shares/clean_IPUMS_1935_1940_extract_to_construct_migration_weights.dta.
	
	Alternatively, one can use state fip and county icp as used to produce the census 
	extract referenced above. The approximate code for this alternative method is below, but may need to be tweaked. 
	With either approach, a few counties (4-5) don't match and must be hand checked.
	
	tostring state, gen(southstatefip_str) 
	replace southstatefip_str=southstatefip_str+"0"
	gen southcounty=countyicp 
	replace southcounty=southcounty+20 if countyicp==24 & southcounty!=5100 & southcounty>50 // county ICP codes in the NHGIS file are shifted forward by 2 digits
	tostring southcounty, gen(southcountyicp_str)  
	replace southcountyicp_str="00"+southcountyicp_str if length(southcountyicp_str)==2 
	replace southcountyicp_str="0"+southcountyicp_str if length(southcountyicp_str)==3
	replace southcountyicp_str=substr(southcountyicp_str,1,length(southcountyicp_str)-2)+ "10" if countyicp==41 & southcountyicp_str=="0605" // Union county in Oregon is 605 in IPUMS census extract but 610 in NHGIS file
	replace southcountyicp_str =substr(southcountyicp_str,1,length(southcountyicp_str)-1)+ "0" if(regexm(southcountyicp_str, "[0-9][0-9][0-9][5]")) // IPUMS Census extract notes county code changes with 0 or 5 but all county codes end in 0 in NHGIS file
	replace southcountyicp_str="1860" if southcountyicp_str=="1930" & countyicp==29 // Discrepancy between Missouri county St Genevieve county code in IPUMS Census extract vs. NHGIS file
	replace southcountyicp_str="7805" if southcounty==7850 & southstatefip_str=="510" // Possible typo with Greenbrier county coded as 785 instead of 775 in IPUMS Census extract. Reassigned to South Norfolk's code from NHGIS file because both are part of Chesapeake (independent city) today.
	replace southcountyicp_str="0050" if southcountyicp_str=="0053" & countyicp==22 // Possible typo with Jefferson Davis county coded as 53 instead of 50 in IPUMS Census extract. Recoded as 50.
	gen gisjoin2_str = southstatefip_str + southcountyicp_str
	cd "$xwalks"
	merge m:1 gisjoin2_str using county1940_crosswalks.dta, keepusing(fips state_name county_name)
	*/
	
	/*
	
	Virginia counties for which migration data are missing:
		  51520 |          1        5.88        5.88
		  51540 |          1        5.88       11.76
		  51560 |          1        5.88       17.65
		  51590 |          1        5.88       23.53
		  51670 |          1        5.88       29.41
		  51680 |          1        5.88       35.29
		  51690 |          1        5.88       41.18
		  51740 |          1        5.88       47.06
		  51750 |          1        5.88       52.94
		  51760 |          1        5.88       58.82
		  51770 |          1        5.88       64.71
		  51790 |          1        5.88       70.59
		  51800 |          1        5.88       76.47
		  51830 |          1        5.88       82.35
		  51840 |          1        5.88       88.24
	
	*/

	cd "$xwalks"
	merge m:1 stateicp countyicp using county1940_crosswalks.dta, keepusing(fips state_name county_name)
	drop if _merge==2 
	g origin_fips=fips
	rename state_name origin_state_name
	rename county_name origin_county_name 
	
	/* Hand correct counties that didn't match using crosswalk file and internet search. */
	
	replace origin_fips = 51067 if countyfips==51620 & _merge==1
	replace origin_fips = 48203 if countyfips==48203 & _merge==1
	replace origin_fips = 51037 if countyfips==54039 & _merge==1
	replace origin_fips = 54041 if countyfips==54041 & _merge==1
	replace origin_fips = 51189 if countyfips==189 & _merge==1
	drop _merge
	
	tostring origin_fips, replace
	keep origin_fips origin_state_name year proutmig actoutmig netbmig_pred
	
	drop if netbmig_pred==. | proutmig==.
	
	bysort origin_fips year: gen dup= cond(_N==1,0,_n)
	tab dup
	drop dup
	
	preserve
	keep origin_fips year proutmig actoutmig netbmig_pred
	
	cd "$instrument"
	save 2_lasso_boustan_predict_mig.dta, replace
	restore
	
	preserve
	collapse (sum) netbmig_pred actoutmig proutmig, by(origin_fips year)
	save "$instrument/2_lasso_boustan_predict_mig_collapsed.dta", replace
	restore
	
	collapse (sum) netbmig_pred actoutmig proutmig, by(origin_state_name year)
	drop if origin_state_name==""
	cd "$xwalks"
	statastates, abbrev(origin_state_name)
	drop _merge 
	rename state_fips origin_state_fips
	tostring origin_state_fips, replace
	
	save "$instrument/3_lasso_boustan_predict_mig_state.dta", replace

	
***2.5. Within state variation in migration.


	/* Load full clean data. */
	
	cd "$migdata"
	use clean_south_county.dta, clear

	/* Predict county-level net migration rate within state, decade by decade. */
	reg netbmig  if year==1950, absorb(stateicp)
	predict netbmig_resid if year==1950, resid
	reg netbmig  if year==1960, absorb(stateicp)
	predict netbmig_resid01 if year==1960, resid
	reg netbmig  if year==1970, absorb(stateicp)
	predict netbmig_resid02 if year==1970, resid

	replace netbmig_resid=netbmig_resid01 if year==1960
	replace netbmig_resid=netbmig_resid02 if year==1970
	drop netbmig_resid01 netbmig_resid02
	
	/* Boustan (2016): Total number leaving/coming to county: actual and predicted. Note 
	that netbmig is a migration rate (per 100 residents). So, the range is -100 to 
	+whatever. -100 because it is impossible for more than all of the residents to 
	leave. But, on the positive side, the rate is unrestricted, because the growth 
	could be quite high (for a county with 100 blacks in 1940, could have 100,000 
	blacks in 1950 which would be a rate of 1000. */ 
	
	gen totbmig=((bpop_l/100)*netbmig)
	gen totbmig_resid=((bpop_l/100)*netbmig_resid)
	gen weight=netbmig_resid*bpop_l
	
	/* One observation per county, year. */
	
	drop if year==year[_n-1]
	sort countyfips year
	drop if countyfips==.
	rename totbmig actoutmig
	rename totbmig_resid residoutmig
	label var residoutmig "migration, by county-year, south, residualized on state"
	drop _merge
	
	/* Merge with 1940 crosswalks data file. */
	
	cd "$xwalks"
	merge m:1 stateicp countyicp using county1940_crosswalks.dta, keepusing(fips state_name county_name)
	drop if _merge==2 
	g origin_fips=fips
	rename state_name origin_state_name
	rename county_name origin_county_name 
	
	/* Hand correct counties that didn't match using crosswalk file and internet search. */
	
	replace origin_fips = 51067 if countyfips==51620 & _merge==1
	replace origin_fips = 48203 if countyfips==48203 & _merge==1
	replace origin_fips = 51037 if countyfips==54039 & _merge==1
	replace origin_fips = 54041 if countyfips==54041 & _merge==1
	replace origin_fips = 51189 if countyfips==189 & _merge==1
	drop _merge
	
	tostring origin_fips, replace
	keep origin_fips year residoutmig actoutmig netbmig_resid
	
	drop if netbmig_resid==. | residoutmig==.
	
	bysort origin_fips year: gen dup= cond(_N==1,0,_n)
	tab dup
	drop if dup>1
	
	drop dup
	
	cd "$instrument"
	save 3_residstate_act_mig.dta, replace
	
	collapse (sum) netbmig_resid actoutmig residoutmig, by(origin_fips year)
	save "$instrument/3_residstate_act_mig_collapsed.dta", replace	


***2.6. Dropping urban counties.


	/* Load full clean data. */
	cd "$migdata"
	use clean_south_county.dta, clear
	
	/* Predict county-level net migration rate, decade by decade with southern 
	variables chosen by LASSO. Predict net migration rate ("netbmig_pred") based on 
	these vars alone. */
	reg netbmig percot perten perag peragtob tob warfac_pc permin perminot ot if year==1950
	predict netbmig_pred if year==1950
	reg netbmig percot perten perag peragtob tob warfac_pc permin perminot ot if year==1960
	predict netbmig_pred01 if year==1960
	reg netbmig percot perten perag peragtob tob warfac_pc permin perminot ot if year==1970
	predict netbmig_pred02 if year==1970	
	
	replace netbmig_pred=netbmig_pred01 if year==1960
	replace netbmig_pred=netbmig_pred02 if year==1970
	drop netbmig_pred01 netbmig_pred02
	
	/* Boustan (2016): Total number leaving/coming to county: actual and predicted. Note 
	that netbmig is a migration rate (per 100 residents). So, the range is -100 to 
	+whatever. -100 because it is impossible for more than all of the residents to 
	leave. But, on the positive side, the rate is unrestricted, because the growth 
	could be quite high (for a county with 100 blacks in 1940, could have 100,000 
	blacks in 1950 which would be a rate of 1000. */ 
	
	gen totbmig=((bpop_l/100)*netbmig)
	gen totbmig_pred=((bpop_l/100)*netbmig_pred)
	gen weight=netbmig_pred*bpop_l
	
	/* One observation per county, year. */
	
	drop if year==year[_n-1]
	sort countyfips year
	drop if countyfips==.
	rename totbmig actoutmig
	rename totbmig_pred proutmig
	label var proutmig "predicted out migration, by county-year, south"
	drop _merge
	
	/* Merge with 1940 crosswalks data file. */
	
	cd "$xwalks"
	merge m:1 stateicp countyicp using county1940_crosswalks.dta, keepusing(fips state_name county_name ur_code_1990)
	drop if _merge==2 
	
	/* Alternative method dropping top 1% percent urban counties */
	*qui bys state: sum perurb, d
	*drop if perurb > `r(p99)'
	
	/* Drops counties that are NCHS-defined as "central" counties of MSAs of 1 million or more population as of 1990. 
	See replication/data/crosswalks/documentation/urban_rural_county_classification/NCHS_Urbrural_File_Documentation.pdf */
	drop if ur_code_1990==1 
	g origin_fips=fips
	rename state_name origin_state_name
	rename county_name origin_county_name 
	
	/* Hand correct counties that didn't match using crosswalk file and internet search. */
	
	replace origin_fips = 51067 if countyfips==51620 & _merge==1
	replace origin_fips = 48203 if countyfips==48203 & _merge==1
	replace origin_fips = 51037 if countyfips==54039 & _merge==1
	replace origin_fips = 54041 if countyfips==54041 & _merge==1
	replace origin_fips = 51189 if countyfips==189 & _merge==1
	drop _merge
	
	tostring origin_fips, replace
	keep origin_fips year proutmig actoutmig netbmig_pred
	
	drop if netbmig_pred==. | proutmig==.
	
	bysort origin_fips year: gen dup= cond(_N==1,0,_n)
	tab dup
	
	drop dup
	
	cd "$instrument"
	save rur_boustan_predict_mig.dta, replace
	
	collapse (sum) netbmig_pred actoutmig proutmig, by(origin_fips year)
	save "$instrument/rur_boustan_predict_mig_collapsed.dta", replace		
	

***2.7. White southern migration.

    /* Data no disponible
	/* Load full clean data. */
	
	cd "$migdata"
	use clean_south_county_white_nonwhite_mig_1940_1970.dta, clear
	
	/* One observation per county, year. */
	
	drop if year==year[_n-1]
	sort fips year
	drop if fips==.
	rename whitemig actoutmig
	
	/* Merge with 1940 crosswalks data file. */
	keep fips state_name county_name *mig* year
	g origin_fips=fips
	rename state_name origin_state_name
	rename county_name origin_county_name 
	
	/* Hand correct counties that didn't match using crosswalk file and internet search. */
	
	tostring origin_fips, replace
	keep origin_fips year actoutmig
	
	drop if actoutmig==.
	
	bysort origin_fips year: gen dup= cond(_N==1,0,_n)
	tab dup
	
	drop dup
	
	cd "$instrument"
	save 5_white_mig.dta, replace
	
	collapse (sum) actoutmig, by(origin_fips year)
	save "$instrument/5_white_mig_collapsed.dta", replace */

	
/* De aqui en adelante los procedimientos dependen de las bases de datos no disponibles
aun asi se hicieron algunos cambios y comentarios
	
*3. Construct versions of shift-share instrument.	

***3.1. Observed southern county net-migration.


	clear all
	set maxvar 32767
		
	global groups black // took out white
	global origin_id origin_fips
	global origin_id_code origin_fips_code
	global origin_sample origin_sample
	global destination_id city
	global destination_id_code city_code
	global dest_sample dest_sample
	global weights_data "$data/instrument/1_boustan_predict_mig.dta"
	global version 0
	global weight_types act // took out act
	global weight_var outmig
	global start_year 1940
	global panel_length 3
	global shares_dir "$data/instrument/shares" 
	global sharesXweights_dir "$data/instrument" 
	
	cd 
	
	use "$migshares/clean_IPUMS_1935_1940_extract_to_construct_migration_weights.dta", clear
		
	do "$bartik/bartik_generic.do"

***3.2. Original Boustan (2010) predicted southern county net-migration.
	
	clear all
	set maxvar 32000
		
	global groups black // took out white
	global origin_id origin_fips
	global origin_id_code origin_fips_code
	global origin_sample origin_sample
	global destination_id city
	global destination_id_code city_code
	global dest_sample dest_sample
	global weights_data "$data/instrument/1_boustan_predict_mig.dta"
	global version 1
	global weight_types pr // took out act
	global weight_var outmig
	global start_year 1940
	global panel_length 3
	global shares_dir "$data/instrument/shares" 
	global sharesXweights_dir "$data/instrument" 
	
	use "$migshares/clean_IPUMS_1935_1940_extract_to_construct_migration_weights.dta", clear
		
	do "$bartik/bartik_generic.do"	
	
***3.3. Post-LASSO prediction of southern county net-migration.


	clear all
	set maxvar 32000
		
	global groups black // took out white
	global origin_id origin_fips
	global origin_id_code origin_fips_code
	global origin_sample origin_sample
	global destination_id city
	global destination_id_code city_code
	global dest_sample dest_sample
	global weights_data "$data/instrument/2_lasso_boustan_predict_mig.dta"
	global version 2
	global weight_types pr // took out act
	global weight_var outmig
	global start_year 1940
	global panel_length 3
	global shares_dir "$data/instrument/shares" 
	global sharesXweights_dir "$data/instrument" 
	
	use "$migshares/clean_IPUMS_1935_1940_extract_to_construct_migration_weights.dta", clear
		
	do "$bartik/bartik_generic.do"
		
***3.4. 1940 southern state of birth shares.
*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%

* Note: Needs to be run on server for space issues
	cd "$stata"
	local years =  "1940"
	foreach year in `years'{
	clear all
	set maxvar 32000
		
	global groups black // took out white
	global origin_id origin_state_fips
	global origin_sample origin_sample
	global destination_id city
	global destination_id_code city_code
	global dest_sample dest_sample
	global weights_data "$data/instrument/3_lasso_boustan_predict_mig_state.dta"
	global version `year'
	global weight_types pr // took out act
	global weight_var outmig
	global start_year 1940
	global panel_length 3
	global shares_dir "$data/instrument/shares" 
	global sharesXweights_dir "$data/instrument" 
	
	use "${migshares}/clean_IPUMS_`year'_extract_to_construct_migration_weights.dta", clear
		
	do "$bartik/bartik_generic.do"
	}		
		
***3.5. Within state southern mig variation, observed


	clear all
	set maxvar 32000
		
	global groups black // took out white
	global origin_id origin_fips
	global origin_id_code origin_fips_code
	global origin_sample origin_sample
	global destination_id city
	global destination_id_code city_code
	global dest_sample dest_sample
	global weights_data "$data/instrument/3_residstate_act_mig.dta"
	global version 7r
	global weight_types resid // took out act
	global weight_var outmig
	global start_year 1940
	global panel_length 3
	global shares_dir "$data/instrument/shares" 
	global sharesXweights_dir "$data/instrument" 
	
	use "$migshares/clean_IPUMS_1935_1940_extract_to_construct_migration_weights.dta", clear
		
	do "$bartik/bartik_generic.do"
	
***3.6. White southern migration.

	clear all
	set maxvar 32000
		
	global groups white // took out white
	global origin_id origin_fips
	global origin_id_code origin_fips_code
	global origin_sample origin_sample
	global destination_id city
	global destination_id_code city_code
	global dest_sample dest_sample
	global weights_data "$data/instrument/5_white_mig.dta"
	global version 8
	global weight_types act // took out act
	global weight_var outmig
	global start_year 1940
	global panel_length 3
	global shares_dir "$data/instrument/shares" 
	global sharesXweights_dir "$data/instrument" 
	
	use "$migshares/clean_IPUMS_1935_1940_extract_to_construct_migration_weights.dta", clear
		
	do "$bartik/bartik_generic.do"
		
***3.7. Non-urban migration.
*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%

	clear all
	set maxvar 32000
		
	global groups black // took out white
	global origin_id origin_fips
	global origin_id_code origin_fips_code
	global origin_sample origin_sample
	global destination_id city
	global destination_id_code city_code
	global dest_sample dest_sample
	global weights_data "$data/instrument/rur_boustan_predict_mig.dta"
	global version r
	global weight_types pr // took out act
	global weight_var outmig
	global start_year 1940
	global panel_length 3
	global shares_dir "$data/instrument/shares" 
	global sharesXweights_dir "$data/instrument" 
	
	use "$migshares/clean_IPUMS_1935_1940_extract_to_construct_migration_weights.dta", clear
		
	do "$bartik/bartik_generic.do"

***3.8. Placebo migration shocks.
*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%

	forval i=1(1)1000{      	
	use	$xwalks/county1940_crosswalks.dta, clear
	g rndmig`i'=rnormal(0,sqrt(5))
	
	rename fips_str origin_fips
	rename rndmig`i' proutmig
	g year=1940
	
	keep origin_fips year proutmig
	
	tempfile rndmig`i'1940
	
	save `rndmig`i'1940'
	
	cap mkdir $data/instrument/rndmig
	cap mkdir $data/instrument/shares/rndmig

	clear all
	set maxvar 120000
		
	global groups black // took out white
	global origin_id origin_fips
	global origin_id_code origin_fips_code
	global origin_sample origin_sample
	global destination_id city
	global destination_id_code city_code
	global dest_sample dest_sample
	global weights_data "`rndmig`i'1940'"
	global version r`i'
	global weight_types pr // took out act
	global weight_var outmig
	global start_year 1940
	global panel_length 0
	global shares_dir "$data/instrument/shares/rndmig" 
	global sharesXweights_dir "$data/instrument/rndmig" 
	
	use $migshares/clean_IPUMS_1935_1940_extract_to_construct_migration_weights.dta, clear
		
	do $bartik/bartik_generic.do
	}


***3.9. Northern measure of 1940 southern county upward mobility.
*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
	/* Mobility
	use	"$mobdata/clean_county_edu_upm_1940.dta", clear
	
**Base de datos no aparece en el paquete de replicación
**ni puede ser creado en otra parte del proceso
	
	g frac_black_upm=black_upm/black_n
	rename fips_str origin_fips
	rename frac_black_upm proutmig
	g year=1940
	
	keep origin_fips year proutmig
	
	save "$data/instrument/mobility1940.dta", replace

	clear all
	set maxvar 32000
		
	global groups black // took out white
	global origin_id origin_fips
	global origin_id_code origin_fips_code
	global origin_sample origin_sample
	global destination_id city
	global destination_id_code city_code
	global dest_sample dest_sample
	global weights_data "$data/instrument/mobility1940.dta"
	global version m
	global weight_types pr // took out act
	global weight_var outmig
	global start_year 1940
	global panel_length 0
	global shares_dir "$data/instrument/shares" 
	global sharesXweights_dir "$data/instrument" 
	
	use "$migshares/clean_IPUMS_1935_1940_extract_to_construct_migration_weights.dta", clear
		
	do "$bartik/bartik_sum_to_one.do" */

***3.10. Clean and standardize city names and output final instrument measures at the city-level
*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%

	* Version 0
	foreach v in "0"{
	
		local origincode origin_fips
		if `v'>1000{
		local origincode origin_state_fips
		}
		
		use ${migshares}/`v'_black`origincode'1940.dta, clear
		
		order black* total*
		egen totblackmigcity3539=rowmean(total_blackcity*)
		sum totblackmigcity3539
		
		drop total* black*
		save ${instrument}/`v'_city_blackmigshare3539.dta, replace
		
		use ${instrument}/`v'_black_actoutmig`origincode'19401970_collapsed_wide.dta, clear
		merge 1:1 city using ${instrument}/`v'_city_blackmigshare3539.dta, keep(3) nogenerate
		
		save ${instrument}/`v'_black_actmig_1940_1970_wide.dta, replace
	
		use ${instrument}/`v'_black_actmig_1940_1970_wide.dta, clear
		decode city, gen(city_str)
		drop city 
		rename city_str city
		
		*Standardize City Names
		//A - fix spelling and formatting variations
		split city, p(,) g(part)
		replace city = proper(part1) + "," + upper(part2) 
		drop part1 part2
			
		*** Initial cleaning done. Save at this point.
		save ${instrument}/city_crosswalked/`v'_black_actmig_1940_1970_wide_preprocessed.dta, replace
		
		use ${instrument}/city_crosswalked/`v'_black_actmig_1940_1970_wide_preprocessed.dta, clear
		g city_original=city
			
		replace city = "St. Joseph, MO" if city == "Saint Joseph, MO" 
		replace city = "St. Louis, MO" if city == "Saint Louis, MO" 
		replace city = "St. Paul, MN" if city == "Saint Paul, MN" 
		replace city = "McKeesport, PA" if city == "Mckeesport, PA" 
		replace city = "Norristown, PA" if city == "Norristown Borough, PA"
		replace city = "Shenandoah, PA" if city == "Shenandoah Borough, PA"
		replace city = "Jamestown, NY" if city == "Jamestown , NY"
		replace city = "Kensington, PA" if city == "Kensington,"
		replace city = "Oak Park Village, IL" if city == "Oak Park Village,"
		replace city = "Fond du Lac, WI" if city == "Fond Du Lac, WI"
		replace city = "DuBois, PA" if city == "Du Bois, PA"
		replace city = "McKees Rocks, PA" if city == "Mckees Rocks, PA"
		replace city = "McKeesport, PA" if city == "Mckeesport, PA"
		replace city = "Hamtramck, MI" if city == "Hamtramck Village, MI"
		replace city = "Lafayette, IN" if city == "La Fayette, IN"
		replace city = "Schenectady, NY" if city == "Schenectedy, NY"
		replace city = "Wallingford Center, CT" if city == "Wallingford, CT"
		replace city = "Oak Park, IL" if city == "Oak Park Village, IL"
		replace city = "New Kensington, PA" if city == "Kensington, PA"
	
		//B - Replace city names with substitutes in the crosswalk when perfect match with crosswalk impossible
		//B1 - the following cities overlap with their subsitutes
		*replace city = "Silver Lake, NJ" if city == "Belleville, NJ"
		replace city = "Brookdale, NJ" if city == "Bloomfield, NJ" 
		replace city = "Upper Montclair, NJ" if city == "Montclair, NJ"
	
		//B2 - the following cities just share a border with their subsitutes but do not overlap
		replace city = "Glen Ridge, NJ" if city == "Orange, NJ"
		replace city = "Essex Fells, NJ" if city == "West Orange, NJ" 
		replace city = "Bogota, NJ" if city == "Teaneck, NJ" 
	
		//B3 - the following cities do not share a border with their substitutes but are within a few miles
		replace city = "Kenilworth, NJ" if city == "Irvington, NJ"  
		replace city = "Wallington, NJ" if city == "Nutley, NJ" 
		replace city = "Short Hills, NJ" if city == "South Orange, NJ"
		replace city = "Lafayette, IN" if city == "Lafayette, IL"
		   
		*Merge with State Crosswalks
		cd "$xwalks"
		merge 1:1 city using US_place_point_2010_crosswalks.dta, keepusing(cz cz_name)
		replace cz = 19600 if city=="Belleville, NJ"
		replace cz_name = "Newark, NJ" if city=="Belleville, NJ"    
		*Resolve Unmerged Cities
		tab _merge
		
		*Save
		drop if _merge==2
		drop _merge
		cd "$instrument"
		save ${instrument}/city_crosswalked/`v'_black_actmig_1940_1970_wide_xw.dta, replace
	}

	* Versions 1, 2, 1940	
	foreach v in "1" "2" "1940" "r"{
	
		local origincode origin_fips
	
		if "`v'"=="1940" {
		local origincode origin_state_fips
		}
		

		use ${migshares}/`v'_black`origincode'1940.dta, clear
		
		order black* total*
		egen totblackmigcity3539=rowmean(total_blackcity*)
		sum totblackmigcity3539
		
		drop total* black*
		save ${instrument}/`v'_city_blackmigshare3539.dta, replace
		
		use ${instrument}/`v'_black_proutmig`origincode'19401970_collapsed_wide.dta, clear
		merge 1:1 city using ${instrument}/`v'_city_blackmigshare3539.dta, keep(3) nogenerate
		
		save ${instrument}/`v'_black_prmig_1940_1970_wide.dta, replace
	
		use ${instrument}/`v'_black_prmig_1940_1970_wide.dta, clear
		decode city, gen(city_str)
		drop city 
		rename city_str city
		
		*Standardize City Names
			//A - fix spelling and formatting variations
			split city, p(,) g(part)
			replace city = proper(part1) + "," + upper(part2) 
			drop part1 part2
		
		*** Initial cleaning done. Save at this point.
		save ${instrument}/city_crosswalked/`v'_black_prmig_1940_1970_wide_preprocessed.dta, replace
		
		use ${instrument}/city_crosswalked/`v'_black_prmig_1940_1970_wide_preprocessed.dta, clear
		g city_original=city
		
		replace city = "St. Joseph, MO" if city == "Saint Joseph, MO" 
		replace city = "St. Louis, MO" if city == "Saint Louis, MO" 
		replace city = "St. Paul, MN" if city == "Saint Paul, MN" 
		replace city = "McKeesport, PA" if city == "Mckeesport, PA" 
		replace city = "Norristown, PA" if city == "Norristown Borough, PA"
		replace city = "Shenandoah, PA" if city == "Shenandoah Borough, PA"
		replace city = "Jamestown, NY" if city == "Jamestown , NY"
		replace city = "Kensington, PA" if city == "Kensington,"
		replace city = "Oak Park Village, IL" if city == "Oak Park Village,"
		replace city = "Fond du Lac, WI" if city == "Fond Du Lac, WI"
		replace city = "DuBois, PA" if city == "Du Bois, PA"
		replace city = "McKees Rocks, PA" if city == "Mckees Rocks, PA"
		replace city = "McKeesport, PA" if city == "Mckeesport, PA"
		replace city = "Hamtramck, MI" if city == "Hamtramck Village, MI"
		replace city = "Lafayette, IN" if city == "La Fayette, IN"
		replace city = "Schenectady, NY" if city == "Schenectedy, NY"
		replace city = "Wallingford Center, CT" if city == "Wallingford, CT"
		replace city = "Oak Park, IL" if city == "Oak Park Village, IL"
		replace city = "New Kensington, PA" if city == "Kensington, PA"
	
		//B - Replace city names with substitutes in the crosswalk when perfect match with crosswalk impossible
		//B1 - the following cities overlap with their subsitutes
		*replace city = "Silver Lake, NJ" if city == "Belleville, NJ"
		replace city = "Brookdale, NJ" if city == "Bloomfield, NJ" 
		replace city = "Upper Montclair, NJ" if city == "Montclair, NJ"
	
		//B2 - the following cities just share a border with their subsitutes but do not overlap
		replace city = "Glen Ridge, NJ" if city == "Orange, NJ"
		replace city = "Essex Fells, NJ" if city == "West Orange, NJ" 
		replace city = "Bogota, NJ" if city == "Teaneck, NJ" 
	
		//B3 - the following cities do not share a border with their substitutes but are within a few miles
		replace city = "Kenilworth, NJ" if city == "Irvington, NJ"  
		replace city = "Wallington, NJ" if city == "Nutley, NJ" 
		replace city = "Short Hills, NJ" if city == "South Orange, NJ"
		replace city = "Lafayette, IN" if city == "Lafayette, IL"
	   
		*Merge with State Crosswalks
		cd "$xwalks"
		merge 1:1 city using US_place_point_2010_crosswalks.dta, keepusing(cz cz_name)
		replace cz = 19600 if city=="Belleville, NJ"
		replace cz_name = "Newark, NJ" if city=="Belleville, NJ"    
		*Resolve Unmerged Cities
		tab _merge
		
		*Save
		drop if _merge==2
		drop _merge
		cd "$instrument"
		save ${instrument}/city_crosswalked/`v'_black_prmig_1940_1970_wide_xw.dta, replace
	}

	foreach v in "7r" {
		
		local origincode origin_fips
		
		use ${migshares}/`v'_black`origincode'1940.dta, clear
		
		order black* total*
		egen totblackmigcity3539=rowmean(total_blackcity*)
		sum totblackmigcity3539
		
		*egen totblackmig3539=sum(totblackmigcity3539)
		drop total* black*
		save ${instrument}/`v'_city_blackmigshare3539.dta, replace
		
		use ${instrument}/`v'_black_residoutmig`origincode'19401970_collapsed_wide.dta, clear
		merge 1:1 city using ${instrument}/`v'_city_blackmigshare3539.dta, keep(3) nogenerate
		
		save ${instrument}/`v'_black_residmig_1940_1970_wide.dta, replace
		
		use ${instrument}/`v'_black_residmig_1940_1970_wide.dta, clear
		decode city, gen(city_str)
		drop city 
		rename city_str city
		
		*Standardize City Names
		//A - fix spelling and formatting variations
		split city, p(,) g(part)
		replace city = proper(part1) + "," + upper(part2) 
		drop part1 part2
		
		*** Initial cleaning done. Save at this point.
		save ${instrument}/city_crosswalked/`v'_black_residmig_1940_1970_wide_preprocessed.dta, replace
		
		use ${instrument}/city_crosswalked/`v'_black_residmig_1940_1970_wide_preprocessed.dta, clear
		g city_original=city
		
		replace city = "St. Joseph, MO" if city == "Saint Joseph, MO" 
		replace city = "St. Louis, MO" if city == "Saint Louis, MO" 
		replace city = "St. Paul, MN" if city == "Saint Paul, MN" 
		replace city = "McKeesport, PA" if city == "Mckeesport, PA" 
		replace city = "Norristown, PA" if city == "Norristown Borough, PA"
		replace city = "Shenandoah, PA" if city == "Shenandoah Borough, PA"
		replace city = "Jamestown, NY" if city == "Jamestown , NY"
		replace city = "Kensington, PA" if city == "Kensington,"
		replace city = "Oak Park Village, IL" if city == "Oak Park Village,"
		replace city = "Fond du Lac, WI" if city == "Fond Du Lac, WI"
		replace city = "DuBois, PA" if city == "Du Bois, PA"
		replace city = "McKees Rocks, PA" if city == "Mckees Rocks, PA"
		replace city = "McKeesport, PA" if city == "Mckeesport, PA"
		replace city = "Hamtramck, MI" if city == "Hamtramck Village, MI"
		replace city = "Lafayette, IN" if city == "La Fayette, IN"
		replace city = "Schenectady, NY" if city == "Schenectedy, NY"
		replace city = "Wallingford Center, CT" if city == "Wallingford, CT"
		replace city = "Oak Park, IL" if city == "Oak Park Village, IL"
		replace city = "New Kensington, PA" if city == "Kensington, PA"
	
		//B - Replace city names with substitutes in the crosswalk when perfect match with crosswalk impossible
		//B1 - the following cities overlap with their subsitutes
		*replace city = "Silver Lake, NJ" if city == "Belleville, NJ"
		replace city = "Brookdale, NJ" if city == "Bloomfield, NJ" 
		replace city = "Upper Montclair, NJ" if city == "Montclair, NJ"
	
		//B2 - the following cities just share a border with their subsitutes but do not overlap
		replace city = "Glen Ridge, NJ" if city == "Orange, NJ"
		replace city = "Essex Fells, NJ" if city == "West Orange, NJ" 
		replace city = "Bogota, NJ" if city == "Teaneck, NJ" 
	
		//B3 - the following cities do not share a border with their substitutes but are within a few miles
		replace city = "Kenilworth, NJ" if city == "Irvington, NJ"  
		replace city = "Wallington, NJ" if city == "Nutley, NJ" 
		replace city = "Short Hills, NJ" if city == "South Orange, NJ"
		replace city = "Lafayette, IN" if city == "Lafayette, IL"
		
		*Merge with State Crosswalks
		cd "$xwalks"
		merge 1:1 city using US_place_point_2010_crosswalks.dta, keepusing(cz cz_name)
		replace cz = 19600 if city=="Belleville, NJ"
		replace cz_name = "Newark, NJ" if city=="Belleville, NJ"    
		*Resolve Unmerged Cities
		tab _merge
		
		*Save
		drop if _merge==2
		drop _merge
		cd "$instrument"
		save ${instrument}/city_crosswalked/`v'_black_residmig_1940_1970_wide_xw.dta, replace
	}
	
	* Version 8	
	foreach v in  "8" {
		use ${migshares}/`v'_whiteorigin_fips1940.dta, clear
		order white* total*
		egen totwhitemigcity3539=rowmean(total_whitecity*)
		sum totwhitemigcity3539
		drop total* white*
		save ${instrument}/`v'_city_whitemigshare3539.dta, replace
		
		use ${instrument}/`v'_white_actoutmigorigin_fips19401970_collapsed_wide.dta, clear
		merge 1:1 city using ${instrument}/`v'_city_whitemigshare3539.dta, keep(3) nogenerate
		
		save ${instrument}/`v'_white_actmig_1940_1970_wide.dta, replace
	
		use ${instrument}/`v'_white_actmig_1940_1970_wide.dta, clear
		decode city, gen(city_str)
		drop city 
		rename city_str city
		
		*Standardize City Names
			//A - fix spelling and formatting variations
			split city, p(,) g(part)
			replace city = proper(part1) + "," + upper(part2) 
			drop part1 part2
		
	*** Initial cleaning done. Save at this point.
	save ${instrument}/city_crosswalked/`v'_white_actmig_1940_1970_wide_preprocessed.dta, replace
	
	use ${instrument}/city_crosswalked/`v'_white_actmig_1940_1970_wide_preprocessed.dta, clear
	g city_original=city
		
			replace city = "St. Joseph, MO" if city == "Saint Joseph, MO" 
			replace city = "St. Louis, MO" if city == "Saint Louis, MO" 
			replace city = "St. Paul, MN" if city == "Saint Paul, MN" 
			replace city = "McKeesport, PA" if city == "Mckeesport, PA" 
			replace city = "Norristown, PA" if city == "Norristown Borough, PA"
			replace city = "Shenandoah, PA" if city == "Shenandoah Borough, PA"
			replace city = "Jamestown, NY" if city == "Jamestown , NY"
			replace city = "Kensington, PA" if city == "Kensington,"
			replace city = "Oak Park Village, IL" if city == "Oak Park Village,"
			replace city = "Fond du Lac, WI" if city == "Fond Du Lac, WI"
			replace city = "DuBois, PA" if city == "Du Bois, PA"
			replace city = "McKees Rocks, PA" if city == "Mckees Rocks, PA"
			replace city = "McKeesport, PA" if city == "Mckeesport, PA"
			replace city = "Hamtramck, MI" if city == "Hamtramck Village, MI"
			replace city = "Lafayette, IN" if city == "La Fayette, IN"
			replace city = "Schenectady, NY" if city == "Schenectedy, NY"
			replace city = "Wallingford Center, CT" if city == "Wallingford, CT"
			replace city = "Oak Park, IL" if city == "Oak Park Village, IL"
			replace city = "New Kensington, PA" if city == "Kensington, PA"
		
			//B - Replace city names with substitutes in the crosswalk when perfect match with crosswalk impossible
			//B1 - the following cities overlap with their subsitutes
			replace city = "Brookdale, NJ" if city == "Bloomfield, NJ" 
			replace city = "Upper Montclair, NJ" if city == "Montclair, NJ"
		
			//B2 - the following cities just share a border with their subsitutes but do not overlap
			replace city = "Glen Ridge, NJ" if city == "Orange, NJ"
			replace city = "Essex Fells, NJ" if city == "West Orange, NJ" 
			replace city = "Bogota, NJ" if city == "Teaneck, NJ" 
		
			//B3 - the following cities do not share a border with their substitutes but are within a few miles
			replace city = "Kenilworth, NJ" if city == "Irvington, NJ"  
			replace city = "Wallington, NJ" if city == "Nutley, NJ" 
			replace city = "Short Hills, NJ" if city == "South Orange, NJ"
			replace city = "Lafayette, IN" if city == "Lafayette, IL"
	   
		*Merge with State Crosswalks
		cd "$xwalks"
		merge 1:1 city using US_place_point_2010_crosswalks.dta, keepusing(cz cz_name)
		replace cz = 19600 if city=="Belleville, NJ"
		replace cz_name = "Newark, NJ" if city=="Belleville, NJ"    
		*Resolve Unmerged Cities
		tab _merge
		
		*Save
		drop if _merge==2
		drop _merge
		cd "$instrument"
		save ${instrument}/city_crosswalked/`v'_white_actmig_1940_1970_wide_xw.dta, replace
}
		
	* Version m (southern upward mobility in the North)
	foreach v in "m"{
	
	if "`v'"=="m"{
	local svar smob
	}
		local destinationcode city
		
		use ${migshares}/`v'_black`destinationcode'1940.dta, clear
		
		order black* total*
		egen totblackmigcity3539=rowmean(total_blackcity*)
		sum totblackmigcity3539
		
		drop total* black*
		save ${instrument}/`v'_city_blackmigshare3539.dta, replace
		
		use ${instrument}/`v'_black_proutmig`destinationcode'19401940_collapsed_wide.dta, clear
		merge 1:1 city using ${instrument}/`v'_city_blackmigshare3539.dta, keep(3) nogenerate
		
		save ${instrument}/`v'_black_`svar'_1940_1940_wide.dta, replace
	
		use ${instrument}/`v'_black_`svar'_1940_1940_wide.dta, clear
		decode city, gen(city_str)
		drop city 
		rename city_str city
		
		*Standardize City Names
		//A - fix spelling and formatting variations
		split city, p(,) g(part)
		replace city = proper(part1) + "," + upper(part2) 
		drop part1 part2
		
		*** Initial cleaning done. Save at this point.
		save ${instrument}/city_crosswalked/`v'_black_`svar'_1940_1940_wide_preprocessed.dta, replace
		
		use ${instrument}/city_crosswalked/`v'_black_`svar'_1940_1940_wide_preprocessed.dta, clear
		g city_original=city
		
		replace city = "St. Joseph, MO" if city == "Saint Joseph, MO" 
		replace city = "St. Louis, MO" if city == "Saint Louis, MO" 
		replace city = "St. Paul, MN" if city == "Saint Paul, MN" 
		replace city = "McKeesport, PA" if city == "Mckeesport, PA" 
		replace city = "Norristown, PA" if city == "Norristown Borough, PA"
		replace city = "Shenandoah, PA" if city == "Shenandoah Borough, PA"
		replace city = "Jamestown, NY" if city == "Jamestown , NY"
		replace city = "Kensington, PA" if city == "Kensington,"
		replace city = "Oak Park Village, IL" if city == "Oak Park Village,"
		replace city = "Fond du Lac, WI" if city == "Fond Du Lac, WI"
		replace city = "DuBois, PA" if city == "Du Bois, PA"
		replace city = "McKees Rocks, PA" if city == "Mckees Rocks, PA"
		replace city = "McKeesport, PA" if city == "Mckeesport, PA"
		replace city = "Hamtramck, MI" if city == "Hamtramck Village, MI"
		replace city = "Lafayette, IN" if city == "La Fayette, IN"
		replace city = "Schenectady, NY" if city == "Schenectedy, NY"
		replace city = "Wallingford Center, CT" if city == "Wallingford, CT"
		replace city = "Oak Park, IL" if city == "Oak Park Village, IL"
		replace city = "New Kensington, PA" if city == "Kensington, PA"
	
		//B - Replace city names with substitutes in the crosswalk when perfect match with crosswalk impossible
		//B1 - the following cities overlap with their subsitutes
		*replace city = "Silver Lake, NJ" if city == "Belleville, NJ"
		replace city = "Brookdale, NJ" if city == "Bloomfield, NJ" 
		replace city = "Upper Montclair, NJ" if city == "Montclair, NJ"
	
		//B2 - the following cities just share a border with their subsitutes but do not overlap
		replace city = "Glen Ridge, NJ" if city == "Orange, NJ"
		replace city = "Essex Fells, NJ" if city == "West Orange, NJ" 
		replace city = "Bogota, NJ" if city == "Teaneck, NJ" 
	
		//B3 - the following cities do not share a border with their substitutes but are within a few miles
		replace city = "Kenilworth, NJ" if city == "Irvington, NJ"  
		replace city = "Wallington, NJ" if city == "Nutley, NJ" 
		replace city = "Short Hills, NJ" if city == "South Orange, NJ"
		replace city = "Lafayette, IN" if city == "Lafayette, IL"
	   
		*Merge with State Crosswalks
		cd "$xwalks"
		merge 1:1 city using US_place_point_2010_crosswalks.dta, keepusing(cz cz_name)
		replace cz = 19600 if city=="Belleville, NJ"
		replace cz_name = "Newark, NJ" if city=="Belleville, NJ"    
		*Resolve Unmerged Cities
		tab _merge
		
		*Save
		drop if _merge==2
		drop _merge
		cd "$instrument"
		save ${instrument}/city_crosswalked/`v'_black_`svar'_1940_1940_wide_xw.dta, replace
	}	

	
	* Version r`i' (placebo shocks)
	forval i=1(1)1000{
	
		local origincode origin_fips
		
		use ${migshares}/rndmig/r`i'_black`origincode'1940.dta, clear
		
		order black* total*
		egen totblackmigcity3539=rowmean(total_blackcity*)
		sum totblackmigcity3539
		
		drop total* black*
		save ${instrument}/rndmig/r`i'_city_blackmigshare3539.dta, replace
		
		use ${instrument}/rndmig/r`i'_black_proutmig`origincode'19401940_collapsed_wide.dta, clear
		merge 1:1 city using ${instrument}/rndmig/r`i'_city_blackmigshare3539.dta, keep(3) nogenerate
		
		save ${instrument}/rndmig/r`i'_black_prmig_1940_1940_wide.dta, replace
	
		use ${instrument}/rndmig/r`i'_black_prmig_1940_1940_wide.dta, clear
		decode city, gen(city_str)
		drop city 
		rename city_str city
		
		*Standardize City Names
			//A - fix spelling and formatting variations
			split city, p(,) g(part)
			replace city = proper(part1) + "," + upper(part2) 
			drop part1 part2
		
		*** Initial cleaning done. Save at this point.
		cap mkdir $instrument/city_crosswalked/rndmig
		save ${instrument}/city_crosswalked/rndmig/r_`i'_black_prmig_1940_1940_wide_preprocessed.dta, replace
		
		use ${instrument}/city_crosswalked/rndmig/r_`i'_black_prmig_1940_1940_wide_preprocessed.dta, clear
		g city_original=city
		
		replace city = "St. Joseph, MO" if city == "Saint Joseph, MO" 
		replace city = "St. Louis, MO" if city == "Saint Louis, MO" 
		replace city = "St. Paul, MN" if city == "Saint Paul, MN" 
		replace city = "McKeesport, PA" if city == "Mckeesport, PA" 
		replace city = "Norristown, PA" if city == "Norristown Borough, PA"
		replace city = "Shenandoah, PA" if city == "Shenandoah Borough, PA"
		replace city = "Jamestown, NY" if city == "Jamestown , NY"
		replace city = "Kensington, PA" if city == "Kensington,"
		replace city = "Oak Park Village, IL" if city == "Oak Park Village,"
		replace city = "Fond du Lac, WI" if city == "Fond Du Lac, WI"
		replace city = "DuBois, PA" if city == "Du Bois, PA"
		replace city = "McKees Rocks, PA" if city == "Mckees Rocks, PA"
		replace city = "McKeesport, PA" if city == "Mckeesport, PA"
		replace city = "Hamtramck, MI" if city == "Hamtramck Village, MI"
		replace city = "Lafayette, IN" if city == "La Fayette, IN"
		replace city = "Schenectady, NY" if city == "Schenectedy, NY"
		replace city = "Wallingford Center, CT" if city == "Wallingford, CT"
		replace city = "Oak Park, IL" if city == "Oak Park Village, IL"
		replace city = "New Kensington, PA" if city == "Kensington, PA"
	
		//B - Replace city names with substitutes in the crosswalk when perfect match with crosswalk impossible
		//B1 - the following cities overlap with their subsitutes
		*replace city = "Silver Lake, NJ" if city == "Belleville, NJ"
		replace city = "Brookdale, NJ" if city == "Bloomfield, NJ" 
		replace city = "Upper Montclair, NJ" if city == "Montclair, NJ"
	
		//B2 - the following cities just share a border with their subsitutes but do not overlap
		replace city = "Glen Ridge, NJ" if city == "Orange, NJ"
		replace city = "Essex Fells, NJ" if city == "West Orange, NJ" 
		replace city = "Bogota, NJ" if city == "Teaneck, NJ" 
	
		//B3 - the following cities do not share a border with their substitutes but are within a few miles
		replace city = "Kenilworth, NJ" if city == "Irvington, NJ"  
		replace city = "Wallington, NJ" if city == "Nutley, NJ" 
		replace city = "Short Hills, NJ" if city == "South Orange, NJ"
		replace city = "Lafayette, IN" if city == "Lafayette, IL"
	   
		*Merge with State Crosswalks
		cd "$xwalks"
		merge 1:1 city using US_place_point_2010_crosswalks.dta, keepusing(cz cz_name)
		replace cz = 19600 if city=="Belleville, NJ"
		replace cz_name = "Newark, NJ" if city=="Belleville, NJ"    
		*Resolve Unmerged Cities
		tab _merge
		
		*Save
		drop if _merge==2
		drop _merge
		cd "$instrument"
		cap mkdir $instrument/city_crosswalked/rndmig
		save ${instrument}/city_crosswalked/rndmig/r`i'_black_prmig_1940_1940_wide_xw.dta, replace
	}


*4. Assemble final dataset

***4.1. Select sample of cities using complete count 1940 census and CCDB 1944-1977.

	/* Load city population data constructed from complete count 1940 census */
	use "$population/clean_city_population_census_1940.dta", clear // 711 cities in non-South
*	merge 1:1 city using "$population/clean_city_population_ccdb_1944_1977.dta", keepusing(bpop1970 pop1940 pop1970)
	
	/*
	* Analysis of non-matches
	not matched                           789
        from master                       273  (_merge==1) // 273 cities from 1940 census city file do not match
        from using                        516  (_merge==2) // 516 cities from CCDB file do not match because they are Southern or they are non-Southern but do not appear in 1940 Census
	
	Here are the cities that do not appear in 1940 census, are non-southern, and have non-missing data for black pop in 1970: Boise city, ID; East Providence, RI: Huntington Park CA; West
	Haven CT; and Warwick, RI 
	
	Here are the cities that do not appear in 1940 census, are non-southern, and are missing data for black pop in 1970:
	Ardmore, PA
	Arlington, MA
	Arlington, VA
	Belmont, MA
	Belvedere, CA
	Bogota, NJ
	Brookline, MA
	Clarksburg, WV
	Drexel Hill, PA
	Haverford College, PA
	Newport, KY
	Secaucus, NJ
	Watertown, MA
	West Hartford, CT
	Woodbridge, NJ

    matched                               438  (_merge==3)
	
	*/
	
	/* Keep cities large enough (25k+) to appear in CCDB in 1940 and 1970. Results are 
	robust to changing this criterion.*/
	rename bpop1970 bpopc1970 // rename so it is clear these numbers correspond to city populations
	rename pop1970 popc1970 // rename so it is clear these numbers correspond to city populations
	
	/* Butte, MT and Amsterdam, NY received southern black migrants between 1935 and 1940, but are just below pop cutoff for CCDB. 
	Keep them in sample by retrieving 1970 black pop info from Census for these cities */
	replace bpopc1970=38 if city=="Butte, MT" // see Table 27 of published 1970 Census: https://www.census.gov/content/dam/Census/library/working-papers/2005/demo/POP-twps0076.pdf
	replace popc1970=23368 if city=="Butte, MT" // see Table 27 of published 1970 Census: https://www.census.gov/content/dam/Census/library/working-papers/2005/demo/POP-twps0076.pdf
	replace bpopc1970=140 if city=="Amsterdam, NY" // see Table 27 of published 1970 Census: https://www2.census.gov/prod2/decennial/documents/1970a_ny1-02.pdf
	replace popc1970=25524 if city=="Amsterdam, NY" // see Table 27 of published 1970 Census: https://www2.census.gov/prod2/decennial/documents/1970a_ny1-02.pdf
	keep if  bpopc1970!=. & pop1940!=.
	
	/* The following non-southern cities are missing Black population data in 1970 though they have total population data for that year
	city
	Bolingbrook, IL
	Burbank, IL
	Burton, MI
	Farmington Hills, MI
	Grosse Pointe Woods, MI
	Irvine, CA
	Rancho Palos Verdes, CA
	Romulus, MI
	*/	
	
	drop if _merge==2 // Dropping cities in CCDB that do not appear in the 1940 Census list of non-southern cities, see analysis of non-matches above. 
	drop _merge 
	
	
***4.2. Merge in data for instrument.
*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
	
	* Version 0 of the instrument: 1935-1940 black southern migrant location choice X observed total 1940-1970 net-migration for southern counties
	foreach v in "0"{
	merge 1:1 city using ${instrument}/city_crosswalked/`v'_black_actmig_1940_1970_wide_xw.dta
	
	/* Drop cities for which there's no hope of getting predictions for black pop in 
	1970 data for these cities. This set of cities will change depending on the 
	migration matrix used.*/
	drop if _merge==2 
	drop _merge
	
	/* Assume zero change in black pop for cities that black migrants did not move 
	to between 1935 and 1940. Results are robust to changing this criterion. 
	Uncomment "keep if _merge==3" and run again. */
	foreach var of varlist black_actoutmigact*{
	replace `var'=0 if `var'==.
	rename `var' v`v'_`var'
	}
	rename totblackmigcity3539 v`v'_totblackmigcity3539
	}

	* Version 1 of the instrument: 
	*	1935-1940 black southern migrant location choice X total 1940-1970 predicted net-migration for southern counties
	*	Original Boustan (2010) variables for prediction. 
	*	See Boustan (2016) replication files for more details: https://scholar.princeton.edu/lboustan/data-books#ch4.
	
	* Version 2 of the instrument: 
	*	1935-1940 black southern migrant location choice X total 1940-1970 Post-LASSO predicted net-migration for southern counties
	*	See Derenoncourt (2019) Appendix B.2 for more details: https://www.dropbox.com/s/58cv5fv1hsofau8/derenoncourt_2019_appendix.pdf?dl=0
	
	* Version 1940 of the instrument: 
	*	1940 black southern-born state of birth X total 1940-1970 Post-LASSO predicted net-migration for southern states

	* Version r of the instrument: 
	*	1935-1940 black southern migrant location choice X total 1940-1970 Post-LASSO predicted net-migration for southern counties
	* 	excluding the 42 major urban southern counties (NCHS-defined "central" counties of MSAs of 1 million or more population)
	* 	See more here: https://www.cdc.gov/nchs/data/data_acces_files/NCHSUrbruralFileDocumentationInternet2.pdf

/*	foreach v in "1" "2" "1940" "r"{
	merge 1:1 city using  ${instrument}/city_crosswalked/`v'_black_prmig_1940_1970_wide_xw.dta
	
	/* Drop cities for which there's no hope of getting predictions for black pop in 
	1970 data for these cities. This set of cities will change depending on the 
	migration matrix used.*/
	drop if _merge==2 
	drop _merge
	
	/* Assume zero change in black pop for cities that black migrants did not move 
	to between 1935 and 1940. Results are robust to changing this criterion. 
	Uncomment "keep if _merge==3" and run again. */
	foreach var of varlist black_proutmigpr*{
	replace `var'=0 if `var'==.
	rename `var' v`v'_`var'
	}
	rename totblackmigcity3539 v`v'_totblackmigcity3539
	} */
	
	* Version 7r of the instrument: 
	*	1935-1940 black southern migrant location choice X total observed 1940-1970 net-migration for southern counties,
	*	residualized on southern state fixed effects.
	foreach v in "7r" {
	merge 1:1 city using  ${instrument}/city_crosswalked/`v'_black_residmig_1940_1970_wide_xw.dta, keepusing(totblackmigcity3539 black_residoutmigresid*)
	*keep if _merge==3
	
	/* Drop cities for which there's no hope of getting predictions for black pop in 
	1970 data for these cities. This set of cities will change depending on the 
	migration matrix used.*/
	drop if _merge==2 
	drop _merge
	
	/* Assume zero change in black pop for cities that black migrants did not move 
	to between 1935 and 1940. Results are robust to changing this criterion. 
	Uncomment "keep if _merge==3" and run again. */
	
	foreach var of varlist black_residoutmigresid*{
	replace `var'=0 if `var'==.
	rename `var' v`v'_`var'
	}
	rename totblackmigcity3539 v`v'_totblackmigcity3539
	}
	
	* Version 8 of the instrument: 
	*	1935-1940 white southern migrant location choice X total observed 1940-1970 white net-migration for southern counties,
	foreach v in "8" {
	merge 1:1 city using  ${instrument}/city_crosswalked/`v'_white_actmig_1940_1970_wide_xw.dta, keepusing(totwhitemigcity3539 white_actoutmigact*)
	
	/* Drop cities for which there's no hope of getting predictions for black pop in 
	1970 data for these cities. This set of cities will change depending on the 
	migration matrix used.*/
	drop if _merge==2 
	drop _merge

	/* Assume zero change in black pop for cities that black migrants did not move 
	to between 1935 and 1940. Results are robust to changing this criterion. 
	Uncomment "keep if _merge==3" and run again. */
	
	foreach var of varlist white_actoutmigact*{
	replace `var'=0 if `var'==.
	rename `var' v`v'_`var'
	}
	rename totwhitemigcity3539 v`v'_totwhitemigcity3539
	}
	
/*	* Placebo versions of the instrument: 
	*	1935-1940 white southern migrant location choice X normally distributed random shocks,
	*	with mean 0 and variance 5, iterated 1000 times.


	forval i=1(1)1000{
	merge 1:1 city using  ${instrument}/city_crosswalked/rndmig/r`i'_black_prmig_1940_1940_wide_xw.dta 
	*keep if _merge==3
	
	/* Drop cities for which there's no hope of getting predictions for black pop in 
	1970 data for these cities. This set of cities will change depending on the 
	migration matrix used.*/
	drop if _merge==2 
	drop _merge
	
	/* Assume zero change in black pop for cities that black migrants did not move 
	to between 1935 and 1940. Results are robust to changing this criterion. 
	Uncomment "keep if _merge==3" and run again. */
	foreach var of varlist black_proutmigpr*{
	replace `var'=0 if `var'==.
	rename `var' vr`i'_`var'
	}
	rename totblackmigcity3539 vr`i'_totblackmigcity3539
	} */


/*	* Northern CZ measure of 1940 southern county upward mobility: 
	*	1935-1940 black southern migrant location choice X total observed 1940-1970 net-migration for southern counties,
	*	residualized on southern state fixed effects.	
	foreach v in "m" {
	
	if "`v'"=="m"{
	local svar smob
	}
		
	local group "black"
	
	merge 1:1 city using ${instrument}/city_crosswalked/`v'_black_`svar'_1940_1940_wide_xw.dta
	* keep if _merge==3
	
	/* Drop cities for which there's no hope of getting predictions for black southern mob in 1970
	for these cities. This set of cities will change depending on the 
	migration matrix used.*/
	drop if _merge==2
	drop _merge

	/* Assume zero change in black pop for cities that black migrants did not move 
	to between 1935 and 1940. Results are robust to changing this criterion. 
	Uncomment "keep if _merge==3" and run again. */
	foreach var of varlist `group'_proutmigpr*{
	egen mean`svar'_`var'=mean(`var')
	replace `var'=mean`svar'_`var' if `var'==.
	replace `var'=popc1940*`var'
	rename `var' v`v'_`var'
	}
	}	*/
	
/*	keep *_proutmigpr* *_actoutmigact* *_residoutmigresid* popc1940 bpopc1940 popc1970 bpopc1970 *migcity3539 statefip citycode city city_original cz cz_name
	drop if popc1970==. */
	save "$city_sample/GM_city_final_dataset.dta", replace	
	

***4.3. Construct measure of black urban pop change and instrument for black urban in-migration at CZ level.
	/*
	/* Generate measure of black urban in-migration at the CZ level. */
	collapse (sum) *_proutmigpr* *_actoutmigact* *_residoutmigresid* popc1940 bpopc1940 popc1970 bpopc1970 *migcity3539 , by(cz cz_name)
		
	* Actual black pop change in city
	g bpopchange1940_1970=100*(bpopc1970-bpopc1940)/popc1940 */

	* Instrument by version
	* Version 0
	foreach v in "0"{
	g v`v'_bpopchange_pred1940_1970=100*v`v'_black_actoutmigact/popc1940
	
	g v`v'_blackmig3539_share1940=100*v`v'_totblackmigcity3539/popc1940
	}
	
/*	* Versions 1, 2, 1940
	foreach v in "1" "2" "1940"{
	g v`v'_bpopchange_pred1940_1970=100*v`v'_black_proutmigpr/popc1940
	
	g v`v'_blackmig3539_share1940=100*v`v'_totblackmigcity3539/popc1940
	} */
	
/*	* Version r
	foreach v in "r" {
	g v`v'_bpopchange_pred1940_1970=100*v`v'_black_proutmigpr/popc1940
	
	g v`v'_blackmig3539_share1940=100*v`v'_totblackmigcity3539/popc1940
	} */
	
	* Versions 7r
	foreach v in "7r"{
	g v`v'_bpopchange_resid1940_1970=100*v`v'_black_residoutmigresid/popc1940
	
	g v`v'_blackmig3539_share1940=100*v`v'_totblackmigcity3539/popc1940
	}
	
	* Versions 8
	foreach v in "8"{

	g v`v'_wpopchange_pred1940_1970=100*v`v'_white_actoutmigact/popc1940
	
	g v`v'_whitemig3539_share1940=100*v`v'_totwhitemigcity3539/popc1940
	}
		

/*	* Placebo shocks
	forval i=1(1)1000{
	g vr`i'_bpopchange_pred1940_1970=100*vr`i'_black_proutmigpr/popc1940
	
	g vr`i'_blackmig3539_share1940=100*vr`i'_totblackmigcity3539/popc1940
	} */


/*	* Northern CZ measure of 1940 southern county upward mobility
	foreach v in "m"{
	
	if "`v'"=="m"{
	local svar smob
	}
		
	local group "black"
	
	g v`v'_black`svar'1940=v`v'_`group'_proutmigpr/popc1940
	} */
		
***4.4. Merge in all datasets.
	
	global datasets ///
	"$mobdata/clean_cz_mobility_1900_2015.dta" /// <- HISTORICAL & CONTEMPORARY MOBILITY OUTCOMES
	"$jobs/clean_cz_industry_employment_1940_1970.dta" "$schools/clean_cz_med_educd_25plus_1940.dta" /// <- CONTROLS (INDUSTRY MIX & EDUCATION IN DESTINATIONS) 
	"$pf/clean_cz_city_police_per_capita_1920_2007.dta" "$pf/clean_cz_public_finance_1932_2012.dta"  /// <- PUBLIC FINANCE 
	"$political/clean_cz_wallace_share_1968.dta" "$political/clean_cz_weighted_racial_animus.dta" /// <- POLITICAL ECONOMY & RACIAL ANIMUS
	"$population/clean_cz_population_1940_1970.dta" "$population/clean_cz_snq_european_immigration_instrument.dta" "$population/clean_bpopshare_1870_2000_cz.dta" /// <- POPULATION
	"$population/clean_cz_population_density_1940.dta" /// <- POPULATION CONT'D
	"$incarceration/clean_cz_city_murder_rates_1931_1969.dta" /// <- CRIME 
	"$incarceration/clean_cz_jail_rates_1920_1960.dta"  "$incarceration/clean_cz_iob_crime_incarceration.dta"  /// <- CRIME & INCARCERATION CONT'D 
	"$nbhds/clean_cz_neighborhoods.dta" "$nbhds/clean_cz_marriage_income_occscore_1940.dta" /// <- NEIGHBORHOODS
	"$schools/clean_cz_prvschl_share_1920_2010.dta" /// <- PRIVATE SCHOOLS

/*	foreach dataset in "$datasets"{
	merge 1:1 cz using `dataset'
	drop if _merge==2
	drop _merge
	} */

	* Carter (1986) 1960s riots data shared only with author by Robert Margo and William Collins
	cap merge 1:1 cz using "$incarceration/clean_cz_riots_1964_1971.dta" 
	cap drop if _merge==2
	cap drop _merge
	
	/* Get state and region info from cz-to-state_id-to-region crosswalk. 
	merge 1:1 cz using $xwalks/cz_state_region_crosswalk.dta, keepusing(state_id region) keep (3) nogenerate
	replace cz_name="Louisville, KY" if cz==13101 // Fill in Louisville, KY name, which was missing. */
	
	
***4.5. Create rank measure of shock. 

	* OLS
*	xtile GM = bpopchange1940_1970, nq(100) 

	* Instrument by version
	* Version 0
	foreach v in "0"{	
	xtile GM_hat`v' = v`v'_bpopchange_pred1940_1970, nq(100) 
	}
	
/*	* Versions 1, 2, 1940
	foreach v in "1" "2" "1940"{	
	xtile GM_hat`v' = v`v'_bpopchange_pred1940_1970, nq(100) 
	} */
	
/*	* Versions r
	foreach v in "r" {	
	xtile GM_hat`v' = v`v'_bpopchange_pred1940_1970, nq(100) 
	} */
	
	* Versions 7r
	foreach v in "7r"{	
	xtile GM_hat`v' = v`v'_bpopchange_resid1940_1970, nq(100) 
	}	
	
	* Versions 8
	foreach v in "8" {	
	xtile GM_hat`v' = v`v'_wpopchange_pred1940_1970, nq(100) 
	}
	

/*	* Placebo shocks
	forval i=1(1)1000{	
	xtile GM_hatr`i' = vr`i'_bpopchange_pred1940_1970, nq(100) 
	}	*/

		
***4.6. Finalize mechanism variables 	
	
	* Construct Bartik instrument for employment
	local sectorlist ag const fire gov man min nr rtl svc tcu wh

	* Construct Bartik shares: share of each industry located in CZ 
	foreach year in 1940 1970 {
		foreach sector in `sectorlist' {
			egen tot_emp_`sector'`year' = sum(emp_`sector'`year')
			gen empshare_`sector'`year' = emp_`sector'`year' / tot_emp_`sector'`year'
		}
	}
	
/*	* Construct Bartik shock using national leave-one-out growth rates for each industry
	foreach sector in `sectorlist' {
		gen tot_LOO_`sector'1970 =  tot_emp_`sector'1970-emp_`sector'1970[_n]
		gen tot_LOO_`sector'1940 =  tot_emp_`sector'1940-emp_`sector'1940[_n]
		gen gr_LOO_`sector'1940_1970 = tot_LOO_`sector'1970-tot_LOO_`sector'1940
	}

	g empchange_1940_1970=emp_tot1970-emp_tot1940
	g empchange1940_1970_pred = 0
	foreach sector in `sectorlist' {
		replace empchange1940_1970_pred = empchange1940_1970_pred + empshare_`sector'1940*gr_LOO_`sector'1940_1970
	} */
	
/*	* Construct actual and predicted change in employment-to-population (epop) ratio (growth)
	g epopchange_1940_1970=empchange_1940_1970/pop1940
	g epopchange_1940_1970_pred=empchange1940_1970_pred/pop1940
	
	* Convert to percentiles to match functional form in rest of analysis
	xtile emp_hat = epopchange_1940_1970_pred, nq(100)
	xtile emp = epopchange_1940_1970, nq(100) */
	
	* Remove unnecessary vars
	foreach sector in `sectorlist' {
		drop emp_`sector'*
	}	

/*	* Construct European migration shock rank
	xtile eur_mig = wt_instmig_avg, nq(100) */

/*	* Construct quartiles of black population
	xtile bpopquartile = bpopshare1940, nq(4) */
	
/*	* Standardize all mechanism vars	
	qui ds *murder_rate*
	foreach var in `r(varlist)'{
	egen `var'_st = std(`var')
	}
	
	qui ds *polpht* *polshare* *polexppc* *fireexppc* *fireshare* *edushare* *eduexpps* 
	foreach var in `r(varlist)'{
	egen `var'_st = std(`var')
	}

	qui ds *prv*
	foreach var in `r(varlist)'{
	egen `var'_st = std(`var')
	}
	
	qui ds c_whtpop_share*
	foreach var in `r(varlist)'{
	egen `var'_st = std(`var')
	}		

	qui ds *jail_rate* *prison_rate*
	foreach var in `r(varlist)'{
	egen `var'_st = std(`var')
	}

	egen wt_racial_animus_st = std(wt_racial_animus) */
	
	* Generate standardized mech vars for post 1970 averages
	
/*	* Govt spending
	foreach cat in "pol" "fire" "hlthhosp" "sani" "rec" "edu"{
	* Exp share
	egen `cat'share_mean1972_2002=rowmean(`cat'share1972 `cat'share1977 `cat'share1982 `cat'share1987 `cat'share1992 `cat'share1997 `cat'share2002 )
	egen `cat'share_mean1972_2002_st=std(`cat'share_mean1972_2002)	
	} */
	
/*	* Per cap
	foreach cat in "pol" "fire" "hlthhosp" "sani" "rec"{
	egen `cat'exppc_mean1972_2002=rowmean(`cat'exppc1972 `cat'exppc1977 `cat'exppc1982 `cat'exppc1987 `cat'exppc1992 `cat'exppc1997 `cat'exppc2002 )
	egen `cat'exppc_mean1972_2002_st=std(`cat'exppc_mean1972_2002)
	} */
	
/*	* Per pupil edu expenditures
	egen eduexpps_mean1972_2002=rowmean(eduexpps1972 eduexpps1977 eduexpps1982 eduexpps1987 eduexpps1992 eduexpps1997 eduexpps2002 )
	egen eduexpps_mean1972_2002_st=std(eduexpps_mean1972_2002)*/
	
/*	* Murder
	egen murder_mean1931_1943= rowmean(murder_rate1931 murder_rate1943 )
	egen murder_mean1931_1943_st=std(murder_mean1931_1943)
	egen murder_mean1977_2002= rowmean(murder_rate1977 murder_rate1982 murder_rate1987 murder_rate1992  murder_rate1997 murder_rate2002 )
	egen murder_mean1977_2002_st=std(murder_mean1977_2002)*/
	
/*	* Incarceration
	egen total_prison_mean1983_2000=rowmean(total_prison_rate1983 total_prison_rate1984 total_prison_rate1985 total_prison_rate1986 total_prison_rate1987 total_prison_rate1988 total_prison_rate1989 total_prison_rate1990 total_prison_rate1991 total_prison_rate1992 total_prison_rate1993 total_prison_rate1994 total_prison_rate1995 total_prison_rate1996 total_prison_rate1997 total_prison_rate1998 total_prison_rate1999 total_prison_rate2000)
	egen total_prison_mean1983_2000_st=std(total_prison_mean1983_2000)*/

/*	* White private school rates
	egen w_prv_mean1970_2000 = rowmean(w_prv_elemhs_share1970 w_prv_elemhs_share1980 w_prv_elemhs_share1990 w_prv_elemhs_share2000)
	egen w_prv_mean1970_2000_st = std(w_prv_mean1970_2000)*/

/*	* Black private school rates
	egen b_prv_mean1970_2000 = rowmean(b_prv_elemhs_share1970 b_prv_elemhs_share1980 b_prv_elemhs_share1990 b_prv_elemhs_share2000)
	egen b_prv_mean1970_2000_st = std(b_prv_mean1970_2000)*/
	
/*	* Private school rates
	egen prv_mean1970_2000 = rowmean(prv_elemhs_share1970 prv_elemhs_share1980 prv_elemhs_share1990 prv_elemhs_share2000)
	egen prv_mean1970_2000_st = std(prv_mean1970_2000)*/
	
	* Standardize remaining mechanism variables
/*	* Racial segregation
	egen cs_race_theil2000_st= std(cs_race_theil_2000)*/
	
/*	* Income segregation
	egen cs00_seg_inc_st=std(cs00_seg_inc)*/
	
/*	* Commute times
	egen frac_traveltime_lt15_st = std(frac_traveltime_lt15)*/

/*	* Wallace votes
	drop wallace_per_white_vote */
	

***4.7. Clean mobility outcome data


	* Construct change in black men's upward mobility variable
/*	* Create z-score of 1940 measure of black boys' educational upward mobility
	egen mean_bmedu1940=mean(frac_blackm_upm1940)
	egen sd_bmedu1940=sd(frac_blackm_upm1940)
	g bmedu1940_zscore=(frac_blackm_upm1940-mean_bmedu1940)/sd_bmedu1940 */
	
/*	* Create z-score of contemporary measure of black men's income upward mobility
	egen mean_bminc1940=mean(kir_black_male_p50)
	egen sd_bminc1940=sd(kir_black_male_p50)
	g bminc1940_zscore=(kir_black_male_p50-mean_bminc1940)/sd_bminc1940 */

/*	* Construct change in z-score
	g mobchangeb=bminc1940_zscore-bmedu1940_zscore */
	
/*	* Standardize change in z-score
	egen mobchangeb_st=std(mobchangeb) */
	
/*	* Construct racial gap in income upward mobility outcomes by CZ 
	foreach p in "25" "50" "75"{
		g racegap2015_p`p'_cz=kfr_white_pooled_p`p'2015*100-kfr_black_pooled_p`p'2015*100
		} */
		

***4.8. Create regional dummies. 

	tabulate region, gen(reg)	
	
***4.9. Create additional 1940 controls. 

	gen urban_share1940 = popc1940/pop1940
	gen ln_pop_dens1940= log(pop_density1940)
	gen ln_mn_occscore1940=log(mn_occscore1940)
	gen ln_mn_incwage1940 = log(mn_incwage1940)
	
***4.10. Label key variables and save final dataset. 

	la var causal_p25_czkr26 "Expos effect low inc (hh)"
	la var causal_p75_czkr26 "Expos effect high inc (hh)"
	la var causal_p25_czkr26_f "Expos effect low inc girls (hh)"
	la var causal_p75_czkr26_f "Expos effect high inc girls (hh)"
	la var causal_p25_czkr26_m "Expos effect low inc boys (hh)"
	la var causal_p75_czkr26_m "Expos effect high inc boys (hh)"
	
	la var causal_p25_czkir26 "Expos effect low inc (ind)"
	la var causal_p75_czkir26 "Expos effect high inc (ind)"
	la var causal_p25_czkir26_f "Expos effect low inc girls (ind)"
	la var causal_p75_czkir26_f "Expos effect high inc girls (ind)"
	la var causal_p25_czkir26_m "Expos effect low inc boys (ind)"
	la var causal_p75_czkir26_m "Expos effect high inc boys (ind)"
	
	la var mobchangeb_st "Change in black men's mobility standardized"
	la var racegap2015_p50_cz "Racial gap 2015 (CZ)"
	la var racegap2015_p50_ct "Racial gap 2015 (CT)"
	
	la var frac_all_upm1940 "Edu. Upward Mobility 1940"
	la var v2_blackmig3539_share1940 "Black Southern Mig 1935-1940"
	la var reg2 "Midwest"
	la var reg3 "South"
	la var reg4 "West"	
	
	la var GM_hat2 "$\hat{GM}$"
	la var GM "GM"

	save $data/GM_cz_final_dataset.dta, replace

	
*5. Produce main figures and tables.

*STEPS:
	*I. Set specs for baseline analysis. 
	*II. Figures.
	*III. Tables.
	*IV. Appendix figures and tables.
	*V. Estimates cited in text.
		
***5.I. Set specs for baseline analysis.

	global x_ols GM
	global x_iv GM_hat2	
	global baseline_controls frac_all_upm1940 mfg_lfshare1940  v2_blackmig3539_share1940 reg2 reg3 reg4

	run $code/programs/PrintEst.do
		
***5.II. Figures.


	* Figure 1: Black upward mobility in 1940 and 2015

	* (a) Percentage of black teens from median-educated families with 9-plus years of schooling, 1940
	use ${mobdata}/clean_cz_mobility_1900_2015.dta, clear
	replace kfr_black_pooled_p50 = kfr_black_pooled_p50*100
	maptile frac_black_upm1940 if black_n>=10 & kfr_black_pooled_n>=10 & frac_black_upm1940!=. & kfr_black_pooled_p50!=., geo(cz1990) rangecolor(jmpgreen*.15 jmpgreen*1.75) ndfcolor(myslate*1.25) conus 
	cd "$figtab"
	graph export black_edu_mobility_1940_map.png, replace

	* (b) Household income rank of black men and women from below-median-income families, 2015	
	maptile kfr_black_pooled_p50 if black_n>=10 & kfr_black_pooled_n>=10 & frac_black_upm1940!=. & kfr_black_pooled_p50!=., geo(cz1990) rangecolor(jmpgreen*.15 jmpgreen*1.75) ndfcolor(myslate*1.25) conus
	cd "$figtab"
	graph export black_p50_mobility_2015_map.png, replace
	
	
	* Figure 2: Quantiles of urban black share increases, 1940-1970

	use ${data}/GM_cz_final_dataset.dta, clear
	gen mylabel=cz_name+"                  " if regexm(cz_name, "Steubenville") | regexm(cz_name, "Milwaukee") ///
	| regexm(cz_name, "Washington") | regexm(cz_name, "Gary") | regexm(cz_name, "Detroit") 
	replace mylabel="Washington, D.C." if regexm(cz_name, "Washington")
	replace mylabel="Detroit, MI" if regexm(cz_name, "Detroit")
	replace mylabel="Gary, IN" if regexm(cz_name, "Gary")
	replace mylabel="Steubenville, OH" if regexm(cz_name, "Steubenville")
	replace mylabel="Milwaukee, WI" if regexm(cz_name, "Milwaukee")
	graph twoway (scatter bpopchange1940_1970 GM if mylabel!="", legend(off) mcolor(jmporange) msymbol(circle_hollow) ///
	msize(vlarge) mlabel(mylabel) mlabcolor(black) mlabposition(11) graphregion(color(white)) plotregion(ilcolor(white)) ///
	ylabel(,nogrid))  || (scatter bpopchange1940_1970 GM, xtitle("Percentile of urban Black pop increase 40-70") ///
	ytitle("Incr. in urban Black pop '40-70 as ppt of 1940 urban pop") legend(off) mcolor(jmpgreen) graphregion(color(white)) plotregion(ilcolor(white))) 
	cd "$figtab"
	graph export bpopchange_percentiles.png, replace
	
	* Point estimates cited in text:
	
		* median
		summ bpopchange1940_1970, d
		local p50_bpopchng4070 = `r(p50)'
		PrintEst `p50_bpopchng4070' "p50_bpopchng4070" "" " percentage points%" "3.1"
	
		xtile pctbpopchange1940_1970 = bpopchange1940_1970, nq(100)
		
		*Pittsburgh
		summ bpopchange1940_1970 if cz_name=="Pittsburgh, PA"
		local pitt_bpopchng4070 = `r(mean)'
		PrintEst `pitt_bpopchng4070' "pitt_bpopchng4070" "" " percentage points%" "3.1"

		summ pctbpopchange1940_1970 if cz_name=="Pittsburgh, PA"
		local pitt_pctbpopchng4070 = `r(mean)'
		PrintEst `pitt_pctbpopchng4070' "pitt_pctbpopchng4070" "" "rd percentile%" "2.0"
		
		*Detroit
		summ bpopchange1940_1970 if cz_name=="Detroit, MI"
		local detr_bpopchng4070 = `r(mean)'
		PrintEst `detr_bpopchng4070' "detr_bpopchng4070" "" " percentage points%" "3.1"

		summ pctbpopchange1940_1970 if cz_name=="Detroit, MI"
		local detr_pctbpopchng4070 = `r(mean)'
		PrintEst `detr_pctbpopchng4070' "detr_pctbpopchng4070" "" "th percentile%" "2.0"
		
		*Salt Lake City
		summ bpopchange1940_1970 if cz_name=="Salt Lake City, UT"
		local slc_bpopchng4070 = `r(mean)'
		PrintEst `slc_bpopchng4070' "slc_bpopchng4070" "" " percentage points%" "3.1"

		summ pctbpopchange1940_1970 if cz_name=="Salt Lake City, UT"
		local slc_pctbpopchng4070 = `r(mean)'
		PrintEst `slc_pctbpopchng4070' "slc_pctbpopchng4070" "" "th percentile%" "2.0"
		
		*Washington, DC
		summ bpopchange1940_1970 if cz_name=="Washington DC, DC"
		local wadc_bpopchng4070 = `r(mean)'
		PrintEst `wadc_bpopchng4070' "wadc_bpopchng4070" "" " percentage points%" "3.1"

		summ pctbpopchange1940_1970 if cz_name=="Washington DC, DC"
		local wadc_pctbpopchng4070 = `r(mean)'
		PrintEst `wadc_pctbpopchng4070' "wadc_pctbpopchng4070" "" "th percentile%" "2.0"
		

	* Figure 3: Relationship between 1940-1970 Black population change and upward mobility in 2012

	eststo clear
	use ${data}/GM_cz_final_dataset.dta, clear
	local x_ols GM 
	local x_iv GM_hat2
	local y perm_res_p25_kr26
	qui: reg `y' `x_ols' 
	local coeff : di %4.3f _b[`x_ols']
	local coeff_se : di %4.3f _se[`x_ols']
	binscatter `y' `x_ols',  ///
	reportreg lcolor(myslate*1.5) ylabel(,nogrid) mcolor(jmpgreen) xtitle("Percentile of Black pop change 40-70") ///
	ytitle("Avg. Adult HH Inc. Rank in 2012, Parents 25p") caption("Slope = `coeff' (`coeff_se')", ring(0) pos(8))
	PrintEst `coeff' "permres_GM_nocontrols" "" " percentiles" "4.2"
	graph export $figtab/permres_GM_nocontrols.png, replace 
	

	* Figure 4: Shift-share instrument for Great Migration

* Uses data shared only with the author and not contained in the replication files.
/*
	* (a) Composition of recent 1935-1940 black southern migrants in Detroit vs. Baltimore 
	use ${migshares}/2_blackorigin_fips1940.dta, clear
	decode(city), gen(city_str)
	keep if (regexm(city_str, "Detroit")==1 | regexm(city_str, "Baltimore")==1)
	keeporder city_str city blackorigin_fips51760 blackorigin_fips1073 blackorigin_fips5093 blackorigin_fips12031 blackorigin_fips13051 ///
	blackorigin_fips22071 blackorigin_fips28027 blackorigin_fips37119 blackorigin_fips45091 blackorigin_fips48113 ///
	blackorigin_fips21111 blackorigin_fips40109 blackorigin_fips47157 blackorigin_fips54047
	g id=_n
	g citytag=lower(substr(city_str,1,strpos(city_str,",")-1))
	reshape long blackorigin_fips, i(id) j(code)
	g space=0
	gen id_order = .
	replace id_order = 1 if code == 51760 & citytag == "detroit"
	replace id_order = 2 if code == 1073 & citytag == "detroit"
	replace id_order = 3 if code == 5093 & citytag == "detroit"
	replace id_order = 4 if code == 12031 & citytag == "detroit"
	replace id_order = 5 if code == 13051 & citytag == "detroit"
	replace id_order = 6 if code == 22071 & citytag == "detroit"
	replace id_order = 7 if code == 28027 & citytag == "detroit"
	replace id_order = 8 if code == 37119 & citytag == "detroit"
	replace id_order = 9 if code == 45091 & citytag == "detroit"
	replace id_order = 10 if code == 48113 & citytag == "detroit"
	replace id_order = 11 if code == 21111 & citytag == "detroit"
	replace id_order = 12 if code == 40109 & citytag == "detroit"
	replace id_order = 13 if code == 47157 & citytag == "detroit"
	replace id_order = 14 if code == 54047 & citytag == "detroit"
	replace id_order = 15 if code == 51760 & citytag == "baltimore"
	replace id_order = 16 if code == 1073 & citytag == "baltimore"
	replace id_order = 17 if code == 5093 & citytag == "baltimore"
	replace id_order = 18 if code == 12031 & citytag == "baltimore"
	replace id_order = 19 if code == 13051 & citytag == "baltimore"
	replace id_order = 20 if code == 22071 & citytag == "baltimore"
	replace id_order = 21 if code == 28027 & citytag == "baltimore"
	replace id_order = 22 if code == 37119 & citytag == "baltimore"
	replace id_order = 23 if code == 45091 & citytag == "baltimore"
	replace id_order = 24 if code == 48113 & citytag == "baltimore"
	replace id_order = 25 if code == 21111 & citytag == "baltimore"
	replace id_order = 26 if code == 40109 & citytag == "baltimore"
	replace id_order = 27 if code == 47157 & citytag == "baltimore"
	replace id_order = 28 if code == 54047 & citytag == "baltimore"
	twoway (bar blackorigin_fips id_order if id_order == 1,  bcolor(navy) ) ///
	(bar blackorigin_fips id_order if id_order == 2,  bcolor(maroon) text(.12 2 "Alabama", color(black))) ///
	(bar blackorigin_fips id_order if id_order > 2 & id_order < 16,  bcolor(dimgray)) ///
	(bar blackorigin_fips id_order if id_order == 15,  bcolor(navy) text(.115 15 "Virginia", color(black))) ///
	(bar blackorigin_fips id_order if id_order == 16,  bcolor(maroon)) ///
	(bar blackorigin_fips id_order if id_order > 16,  bcolor(dimgray)) ///
	(pcarrowi 0.05 6 .00367647 1 (3), lcolor(gray) mcolor(gray) text(0.053 6 "Virginia", color(black)) ) ///
	(pcarrowi 0.03 19 .00093041 16 (3), lcolor(gray) mcolor(gray) text(0.033 19 "Alabama", color(black)) ) ///
	,legend(off) ysc(r(0(.02).14)) xtitle("") xlabel("") ytitle("") graphregion(color(white)) ///
	ylabel(0 "0%" .02 "2%" .04 "4%" .06 "6%" .08 "8%" .10 "10%" .12 "12%" .14 "14%", angle(horizontal) nogrid) ///
	caption("      Detroit                                    Baltimore          ", ring(0) position(12) size(large)) xsize(8)
	cd "$figtab"
	graph export "detroit_baltimore_migrant_composition.png", replace
	
	* (b) Southern state net-migration, 1940-1970 
	use ${migdata}/raw/south_county.dta, clear
	drop _merge 
	keep  netbmig bpop_l year state
	g mig = (netbmig/100)*bpop_l
	cd "$xwalks"
	statastates, fips(state)
	replace state_name=strlower(state_name)
	keep if (state_name=="alabama" | state_name=="arkansas" | state_name=="florida" | state_name=="georgia" | state_name=="kentucky" ///
	| state_name=="louisiana" | state_name=="mississippi" | state_name=="north carolina" | state_name=="oklahoma" | state_name=="south carolina" ///
	| state_name=="tennessee" | state_name=="texas" | state_name=="virginia" | state_name=="west virginia")
	collapse (sum)  mig , by(state_name year)
	g origin_state_name = upper(state_name)
	
	merge 1:1 origin_state_name year using ${instrument}/3_lasso_boustan_predict_mig_state.dta, keepusing(proutmig) nogenerate
	
	replace origin_state_name = lower(origin_state_name)
	drop state_name
	rename origin_state_name state_name
	
	replace proutmig=mig if proutmig==.
	replace proutmig=proutmig/1000

	append using ${migdata}/raw/fmt_net_migration_south.dta
	
	replace state_name=subinstr(state_name, " ", "", .)

	preserve
	keep if year==1930
	sort mig
	restore
	
	preserve
	keep if year==1970
	sort mig
	restore
		
	reshape wide mig proutmig, i(year) j(state_name) string
	tsset year, delta(10)
	
	keep if year>=1940
	
	foreach var in "migalabama" "migarkansas" "migflorida" "miggeorgia" "migkentucky" "miglouisiana" "migmississippi" "migoklahoma" "migtennessee" "migtexas" "migvirginia"{
	replace `var'=`var'/1000
	gen label_`var' = proper(subinstr("`var'", "mig","",.)) if year==1970
	}
	
	gen label_proutmigalabama = proper(subinstr("proutmigalabama", "proutmig", "Pred Mig ",.)) if year==1970
	gen label_proutmigvirginia = proper(subinstr("proutmigvirginia", "proutmig", "Pred Mig ",.)) if year==1970

	twoway (tsline migalabama , recast(connected) msymbol(none) lcolor(maroon%50) mlabel(label_migalabama) mlabcolor(black)) ///
	(tsline proutmigalabama , recast(connected) msymbol(none) lcolor(maroon)  lwidth(thick) mlabel(label_proutmigalabama) mlabcolor(black))   ///
	(tsline migarkansas , recast(connected) msymbol(none) lcolor(myslate))   ///
	(tsline migflorida , recast(connected) msymbol(none) lcolor(myslate))  ///
	(tsline miggeorgia , recast(connected) msymbol(none)  lcolor(myslate))  ///
	(tsline migkentucky , recast(connected) msymbol(none) lcolor(myslate))  ///
	(tsline miglouisiana , recast(connected) msymbol(none) lcolor(myslate)) ///
	(tsline migmississippi , recast(connected) msymbol(none) lcolor(myslate)) ///
	(tsline migoklahoma , recast(connected) msymbol(none) lcolor(myslate)) ///
	(tsline migtennessee , recast(connected) msymbol(none) lcolor(myslate)) ///
	(tsline migtexas , recast(connected) msymbol(none) lcolor(myslate)) ///
	(tsline migvirginia , recast(connected) msymbol(none)  lcolor(jmpdarkblue%50)   mlabel(label_migvirginia) mlabcolor(black)) /// 	 
	(tsline proutmigvirginia , recast(connected) msymbol(none)  lcolor(jmpdarkblue) lwidth(thick) mlabel(label_proutmigvirginia) mlabcolor(black)), /// 
	xlabel(1940(10)1980) graphregion(color(white)) plotregion(ilcolor(white)) ysc(range(-350 100)) ///
	ytitle("Thousands") xtitle("") ylabel( ,nogrid) yline(0, lcolor(black)) legend(off) 
	cd "$figtab"
	graph export "southern_netmig_1940_1970_2.png", replace	

*/
	

	* Figure 5: First stage on Black population change

	eststo clear
	use ${data}/GM_cz_final_dataset.dta, clear
	local x_ols GM 
	local x_iv GM_hat2

	qui: reg `x_ols' `x_iv' ${baseline_controls} 
	test `x_iv' = 0
	local first_stage_fstat : di %4.2f `r(F)'
	PrintEst `first_stage_fstat' "first_stage_fstat" "" "%" "4.1"

	local first_stage : di %4.3f _b[`x_iv']
	local first_stage_se : di %4.3f _se[`x_iv']
	
	PrintEst `first_stage' "first_stage" "" " percentile" "4.2"

	binscatter `x_ols' `x_iv', controls( ${baseline_controls}) ///
	reportreg lcolor(myslate*1.5) ylabel(,nogrid) mcolor(jmpgreen) ///
	xtitle("Percentile of predicted Black pop change 40-70") ytitle("Percentile of Black pop change 40-70") ///
	caption("F Stat = `first_stage_fstat'", ring(0) pos(7))
	graph export $figtab/first_stage_GM_GM_hat2.png, replace 


	* Figure 6: Great Migration reduced average upward mobility in northern commuting zones

	eststo clear
	use ${data}/GM_cz_final_dataset.dta, clear
	local x_ols GM 
	local x_iv GM_hat2
	local y perm_res_p25_kr26
	qui: reg `y' `x_iv' ${baseline_controls} 
	local coeff : di %4.3f _b[`x_iv']
	local coeff_se : di %4.3f _se[`x_iv']
	binscatter `y' `x_iv' , controls(${baseline_controls}) ///
	reportreg lcolor(myslate*1.5) ylabel(,nogrid) mcolor(jmpgreen) ///
	xtitle("Percentile of predicted Black pop change 40-70") ytitle("Expected Mean Adult HH Inc. Rank, Parents 25p") ///
	caption("Slope = `coeff' (`coeff_se')", ring(0) pos(7)) 
	graph export $figtab/permres_GM_hat2.png, replace 


	* Figure 7: Childhood in Great Migration CZ lowers adult income of children from low income families

	eststo clear
	use ${data}/GM_cz_final_dataset.dta, clear
	local x_ols GM 
	local x_iv GM_hat2
	local y causal_p25_czkr26
	qui: reg `y' `x_iv' ${baseline_controls} [aw=1/(`y'_se^2)]
	local coeff : di %12.4f _b[`x_iv']
	local coeff_se : di %5.4f _se[`x_iv']
	binscatter `y' `x_iv' [aw= 1/(`y'_se^2)], controls(${baseline_controls}) ///
	reportreg lcolor(myslate*1.5) ylabel(,nogrid) mcolor(jmpgreen) ///
	xtitle("Percentile of predicted Black pop change 40-70") ytitle("CZ Expos Effect on HH Inc. Rank, Parents 25p") ///
	caption( "Slope = `coeff' (`coeff_se')", ring(0) pos(7))
	graph export $figtab/causal_GM_hat2.png, replace 


	* Figure 8: Race and gender heterogeneity in impact of Great Migration on Upward Mobility

	use ${data}/GM_cz_final_dataset.dta, clear
		
	local outcomes kfr_black_pooled_p252015 kfr_black_pooled_p752015 kir_black_male_p252015 kir_black_male_p752015 kir_black_female_p252015 kir_black_female_p752015 kfr_white_pooled_p252015 kfr_white_pooled_p752015 kir_white_male_p252015 kir_white_male_p752015 kir_white_female_p252015 kir_white_female_p752015
	foreach y in `outcomes'{
	replace `y'=`y'*100
	}
	
	* Rescale treatment in terms of standard deviations
	qui sum GM_hat2
	replace GM_hat2=GM_hat2/`r(sd)'	
	qui sum GM
	replace GM=GM/`r(sd)'

	foreach y of varlist `outcomes'{
	ivreg2 `y' ($x_ols =$x_iv ) ${baseline_controls} if `y'!=., first 	
	eststo `y'_iv
	}	

	local xvar "GM"
	local type "_iv"
	local scale "xsc(r(-8 8))"
	local xlabel "xla(-8(2)8)"	

	coefplot (kir_black_male_p252015`type', mlabels(`xvar' = 12 "Black men, low inc ") keep(`xvar') ciopts(lcolor(jmpgreen)) mcolor(jmpgreen) mlabcolor(jmpgreen)) ///
	(kir_black_male_p752015`type', mlabels(`xvar' = 12 "Black men, high inc ") keep(`xvar') ciopts(lcolor(jmpgreen)) mcolor(jmpgreen) mlabcolor(jmpgreen)) ///
	(kir_black_female_p252015`type', mlabels(`xvar' = 12 "Black women,  low inc  ") keep(`xvar') xline(0, lcolor(gs8)) ciopts(lcolor(myslate*1.5)) mcolor(myslate*1.5) mlabcolor(myslate*1.5)) ///
	(kir_black_female_p752015`type',mlabels(`xvar' = 12 "Black women,  high inc  ") keep(`xvar') ciopts(lcolor(myslate*1.5)) mcolor(myslate*1.5) mlabcolor(myslate*1.5)) ///
	(kir_white_male_p252015`type', mlabels(`xvar' = 12 "White men, low inc ")  keep(`xvar') ciopts(lcolor(myslate*1.5)) mcolor(myslate*1.5) mlabcolor(myslate*1.5)) ///
	(kir_white_male_p752015`type', mlabels(`xvar' = 12 "White men, high inc ") keep(`xvar') ciopts(lcolor(myslate*1.5)) mcolor(myslate*1.5) mlabcolor(myslate*1.5))  ///
	(kir_white_female_p252015`type', mlabels(`xvar' = 12 "White women,  low inc  ")keep(`xvar') ciopts(lcolor(myslate*1.5)) mcolor(myslate*1.5) mlabcolor(myslate*1.5)) ///
	(kir_white_female_p752015`type', mlabels(`xvar' = 12 "White women,  high inc  ") keep(`xvar') ciopts(lcolor(myslate*1.5)) mcolor(myslate*1.5) mlabcolor(myslate*1.5)), ///
	`scale' `xlabel' graphregion(color(white)) plotregion(ilcolor(white)) ylabel(none,nogrid)  legend(off) ytitle("") ///
	ylabel(none) xtitle("Percentile Change in Average Adult Income Rank in CZ")	
	cd "$figtab"
	graph export GM_race_kir_mobility`type'_coefplot.png, replace	


	* Figure 9: Great Migration CZs have higher segregation, crime, and policing

	cap mkdir "$figtab/tempgraphs"

	use ${data}/GM_cz_final_dataset.dta, clear
	eststo clear
	
	* Rescale treatment in terms of standard deviations
	qui sum GM_hat2
	replace GM_hat2=GM_hat2/`r(sd)'	
	qui sum GM
	replace GM=GM/`r(sd)'

	* Government spending
	foreach cat in "pol" "edu"{
	reg `cat'share1932_st GM_hat2 ${baseline_controls}
	eststo `cat'
	}
	
	foreach cat in "pol"{		
	reg `cat'exppc1932_st GM_hat2 ${baseline_controls}	
	eststo pc`cat'
	}	

	foreach cat in "edu"{		
	reg `cat'expps1932_st GM_hat2 ${baseline_controls}	
	eststo pc`cat'
	}	
	
	* Murder
	reg murder_mean1931_1943_st GM_hat2 ${baseline_controls}
	eststo murder

	* Incarceration
	reg jail_rate1940_st GM_hat2 ${baseline_controls}	
	eststo prison

	* Private school
	reg prv_elemhs_share1920_st GM_hat2 ${baseline_controls}
	eststo prvschl

	local xvar GM_hat2
	
	coefplot (prvschl, mlabels(`xvar' = 12 "Private School")  keep(`xvar')) ///
	(murder, mlabels(`xvar' = 12 "Murder")  keep(`xvar')) ///
	(pol, mlabels(`xvar' = 12 "Police Exp Share") keep(`xvar')) ///
	(pcpol, mlabels(`xvar' = 12 "Police Exp Per Cap")  keep(`xvar')) ///
	(prison, mlabels(`xvar' = 12 "Incarceration") keep(`xvar')) ///
	(edu,mlabels(`xvar' = 12 "Education Exp Share") keep(`xvar')) ///
	(pcedu,mlabels(`xvar' = 12 "Education Exp Per Pupil")  keep(`xvar')), ///
	graphregion(color(white)) plotregion(ilcolor(white)) ylabel(none,nogrid)  legend(off) ytitle("")  ///
	ylabel(none) xtitle("(a) Effects on pre-1940 mechanisms") caption ("Units are standard deviations.", ring(0) pos(4) size(small)) xline(0, lcolor(gs8)) xsc(range(-.5 1)) xlabel(-.5 .5 1) name(p1,replace)
	cd "$figtab"
	graph export GM_locgov_coefplot_pretrends.png, replace	
	
	use ${data}/GM_cz_final_dataset.dta, clear
	eststo clear
	
	* Rescale treatment in terms of standard deviations
	qui sum GM_hat2
	replace GM_hat2=GM_hat2/`r(sd)'	
	qui sum GM
	replace GM=GM/`r(sd)'
	
	* Government spending
	foreach cat in "pol" "fire" "hlthhosp" "sani" "rec" "edu"{	
	reg `cat'share_mean1972_2002_st GM_hat2 ${baseline_controls}
	eststo `cat'
	}
	
	* Murder
	reg murder_mean1977_2002_st GM_hat2 ${baseline_controls} murder_mean1931_1943_st
	eststo murder
	
	* Incarceration
	reg total_prison_mean1983_2000_st GM_hat2 ${baseline_controls} murder_mean1931_1943_st
	eststo prison

	* Private school
	reg w_prv_mean1970_2000_st GM_hat2 ${baseline_controls} murder_mean1931_1943_st
	eststo prvschl
	
	* Black Private school
	reg b_prv_mean1970_2000_st GM_hat2 ${baseline_controls} murder_mean1931_1943_st
	eststo prvschl2

	* Racial segregation
	reg cs_race_theil2000_st GM_hat2 ${baseline_controls} murder_mean1931_1943_st
	eststo rseg
	
	* Income segregation
	reg cs00_seg_inc_st GM_hat2 ${baseline_controls} murder_mean1931_1943_st
	eststo iseg
	
	* Commute times
	reg frac_traveltime_lt15_st GM_hat2 ${baseline_controls} murder_mean1931_1943_st
	eststo ct	

	local xvar GM_hat2

	coefplot (prvschl, mlabels(`xvar' = 12 "   White Private School")  keep(`xvar')) ///
	(prvschl2, mlabels(`xvar' = 12 "Black Private School   ") keep(`xvar')) ///	
	(rseg, mlabels(`xvar' = 12 "Residential Racial Segregation") keep(`xvar')) ///
	(iseg, mlabels(`xvar' = 12 "  Residential Income Segregation") keep(`xvar')) ///
	(ct, mlabels(`xvar' = 12 "Frac w/ short commutes    ") keep(`xvar') ) ///		
	(murder, mlabels(`xvar' = 12 "Murder")  keep(`xvar')) ///	
	(pol, mlabels(`xvar' = 12 "   Police")  keep(`xvar')) ///			
	(prison, mlabels(`xvar' = 12 "Incarceration") keep(`xvar')) ///		
	(edu,mlabels(`xvar' = 12 "Education") keep(`xvar')) ///
	(fire, mlabels(`xvar' = 12 "Fire") keep(`xvar')) ///
	(hlthhosp, mlabels(`xvar' = 12 "Health & Hospitals") keep(`xvar')) ///
	(sani, mlabels(`xvar' = 12 "Sanitation") keep(`xvar')) ///
	(rec, mlabels(`xvar' = 12 "Recreation") keep(`xvar') xline(0, lcolor(gs8))), ///
	graphregion(color(white)) plotregion(ilcolor(white)) ylabel(none,nogrid)  legend(off) ytitle("")  ///
	ylabel(none) xtitle("(b) Effects on post-1970 mechanisms") caption("Units are standard deviations." "Controls for pre-1940 murder rates.", ring(0) pos(4) size(small))  xsc(range(-.5 1)) xlabel(-.5 .5 1) name(p2,replace)
	cd "$figtab"
	graph export GM_locgov_coefplot.png, replace	

	cd $figtab/tempgraphs
	graph combine p1 p2 , rows(1)  graphregion(color(white)) iscale(.75)
	
	graph export $figtab/GM_locgov_coefplot_prepost.png, replace 

	cap rmdir $figtab/tempgraphs
	
*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%	
*2. Tables.
*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%

*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
	* Table 1: Contribution of location versus selection in Great Migration effects
*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
	use ${data}/GM_cz_final_dataset.dta, clear
	
	* Rescale treatment in terms of standard deviations
	qui sum GM_hat2 
	replace GM_hat2=GM_hat2/`r(sd)'	
	qui sum GM 
	replace GM=GM/`r(sd)'
	local ols_SD_wt= `r(sd)'
	di %4.0f `ols_SD_wt'
	PrintEst `ols_SD_wt' "ols_SD_wt" "" "%" "4.0"
	
	qui: ivreg2 perm_res_p25_kr26 ($x_ols = $x_iv )  ${baseline_controls} , first
	local GM_imp_upm = _b[$x_ols]
	di %4.1f `GM_imp_upm'
	PrintEst `GM_imp_upm' "GM_imp_upm" "" "" "4.1"
	
	use ${data}/GM_cz_final_dataset.dta, clear
	
	* Rescale treatment in terms of standard deviations
	qui sum GM_hat2 [aw=1/(causal_p25_czkr26_se^2)]
	replace GM_hat2=GM_hat2/`r(sd)'	
	qui sum GM [aw=1/(causal_p25_czkr26_se^2)]
	replace GM=GM/`r(sd)'
	local causal_SD= `r(sd)'
	di %4.0f `causal_SD'
	PrintEst `causal_SD' "causal_SD" "" "%" "4.0"

	qui: ivreg2 causal_p25_czkr26 ($x_ols = $x_iv )  ${baseline_controls} [aw=1/(causal_p25_czkr26_se^2)], first
	local causal = _b[$x_ols]
	
	local GM_imp_loc_20 = `causal' *20
	local GM_imp_loc_15_53 = `causal'  *15.53
	local GM_imp_loc_14_52 = `causal'  *14.52
	
	PrintEst `GM_imp_loc_20' "GM_imp_loc_20" "" "" "4.1"
	PrintEst `GM_imp_loc_15_53' "GM_imp_loc_15_53" "" "" "4.1"
	PrintEst `GM_imp_loc_14_52' "GM_imp_loc_14_52" "" "" "4.1"

	* Ratio	
	local loc_upm_ratio_20 = (`causal'*20/`GM_imp_upm')*100 
	local loc_upm_ratio_15_53 = (`causal'*15.53/`GM_imp_upm')*100 
	local loc_upm_ratio_14_52 = (`causal'*14.52/`GM_imp_upm')*100 
	
	PrintEst `loc_upm_ratio_20' "loc_upm_ratio_20" "" "%" "4.0"
	PrintEst `loc_upm_ratio_15_53' "loc_upm_ratio_15_53" "" "%" "4.0"
	PrintEst `loc_upm_ratio_14_52' "loc_upm_ratio_14_52" "" "%" "4.0"

	* Percent income (for p25, 1 percentile corresponds to 3.14% of income, see Chetty & Hendren (2018b): https://opportunityinsights.org/wp-content/uploads/2018/03/movers_paper2.pdf)
	local pct_income_permres = abs(`GM_imp_upm')*3.14
	PrintEst `pct_income_permres' "pct_income_permres" "" "%" "4.1"	
	
	* Absolute value of effect for in-text citation
	local GM_imp_upm_abs = abs(`GM_imp_upm')
	PrintEst `GM_imp_upm_abs' "GM_imp_upm_abs" "" "%" "4.1"
	
	local GM_imp_loc_14_52 = abs(`causal'  *14.52)
	PrintEst `GM_imp_loc_14_52' "GM_imp_loc_14_52_abs" "" "%" "4.1"

*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
	* Table 2: Great Migration contribution to northern racial upward mobility gap
*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
	use ${data}/GM_cz_final_dataset.dta, clear

	* Multiply outcomes by 100 to be in percentile units
	foreach var of varlist kfr_white_pooled* kfr_black_pooled*{
	replace `var'=100*`var'
	}

	* Calculate mean of Black upward mobility in the sample
	foreach p in "p25" "p50" "p75"{
	sum kfr_black_pooled_`p', d
	local black_mean_`p'=r(mean)	
	}

	* Create a counterfactual dataset where GM is lowest percentile across all CZs in sample
	replace GM_hat2=1
	replace GM=1
	g id=_n+130
	
	* Replace Black upward mobility as missing in the counterfactual dataset
	foreach r in "black"{
	foreach p in "p25" "p50" "p75"{
	replace kfr_`r'_pooled_`p'=.
	}
	}

	* Calculate racial gap (which will be missing for all values) in the counterfactual dataset
	foreach p in "p25" "p50" "p75"{
	g racegap2015_`p' = kfr_white_pooled_`p'-kfr_black_pooled_`p'
	}
	
	tempfile GM_cf
	
	save `GM_cf'

	* Open a new copy of the original dataset for generating the fitted values
	use ${data}/GM_cz_final_dataset.dta, clear
	
	* Multiply outcomes by 100 to be in percentile units	
	foreach var of varlist kfr_white_pooled* kfr_black_pooled*{
	replace `var'=100*`var'
	}

	* Calculate racial gap in upward mobility in each CZ
	foreach p in "p25" "p50" "p75"{
	g racegap2015_`p' = kfr_white_pooled_`p'-kfr_black_pooled_`p'
	}
	
	g id=_n
	append using `GM_cf'
	
	* 2SLS regression of Black upward mobility on GM; generate fitted vaues
	foreach r in "black" {
	foreach p in "p25" "p50" "p75"{
	ivreg2 kfr_`r'_pooled_`p' ($x_ols = $x_iv) ${baseline_controls}	if kfr_`r'_pooled_`p' !=.  // Only uses original sample in the 2SLS regression
	predict kfr_`r'_pooled_`p'_pred, xb // Generates fitted values of Black upward mobility for all obs using X vars & regression coefficients from step above
	predict kfr_`r'_pooled_`p'_pred_se, stdp	 // Generates standard errors for above fitted values 
	}
	}
	
	* Calculate actual average racial gap in sample CZs
	foreach p in "p25" "p50" "p75"{
	sum racegap2015_`p' if id<131, d           // Only calculate on original sample (not counterfactual sample)
	local rgp_mean_`p'=r(mean)					// Store mean in original sample

	* Compute counterfactual gap by predicting Black outcomes if lowest percentile of GM
	g racegap2015_`p'_pred = kfr_white_pooled_`p'-kfr_black_pooled_`p'_pred  // Calculate racial gap using fitted values	
	sum racegap2015_`p'_pred if id>130, d
	g rgpcf_mean_`p'=r(mean)
	local rgpcf_mean_`p'=r(mean)
	sum kfr_black_pooled_`p'_pred, d
	local blackcf_north_`p'=r(mean)
	}

	* Generate standard errors for counterfactual racial gap
	keep if id>130

	* Variance of counterfactual racial gaps
	foreach p in "p25" "p50" "p75"{	
	egen var_rg_cf_`p'=sum((racegap2015_`p'_pred - rgpcf_mean_`p')^2)
	replace var_rg_cf_`p'=(var_rg_cf_`p')/(_N-1)
	g sqrt_var_rg_cf_`p' = sqrt(var_rg_cf_`p')
	local sqrt_var_rg_cf_`p' = sqrt_var_rg_cf_`p'
	}

	tempname sample
	file open `sample' using "$figtab/text/samplesize.txt", text write replace
	file write `sample' `"`r(N)'"' 
    file close `sample'
	
	foreach p in "p25" "p50" "p75"{
	local text_`p' : di%4.2f `rgp_mean_`p''
	tempname number
	file open `number' using "$figtab/text/rg`p'.txt", text write replace
	file write `number' `"`text_`p''"' 
    file close `number'	
	local text_`p' : di%4.2f `rgpcf_mean_`p''	
	file open `number' using "$figtab/text/rgcf`p'.txt", text write replace
	file write `number' `"`text_`p''"' 
    file close `number'	
	local text_`p' : di%4.2f `sqrt_var_rg_cf_`p''
	file open `number' using "$figtab/text/rg_se`p'.txt", text write replace
	file write `number' `"`text_`p''"' 
    file close `number'	
	local change`p' = (`rgp_mean_`p''-`rgpcf_mean_`p'')*100/`rgp_mean_`p''
	local text_`p' : di %4.0f `change`p''
	file open `number' using "$figtab/text/change`p'.txt", text write replace
	file write `number' `"`text_`p''%"' 
    file close `number'		
	}
	
	foreach p in "p25" "p50" "p75"{
	di ""
	di %4.2f `rgp_mean_`p''
	di %4.2f `rgpcf_mean_`p''	
	di %4.2f `sqrt_var_rg_cf_`p''
	di  %4.0f (`rgp_mean_`p''-`rgpcf_mean_`p'')*100/`rgp_mean_`p''
	di ""
	di %4.2f `black_mean_`p''
	di %4.2f `blackcf_mean_`p''
	}

	* Alternative calculation not assuming zero effect of GM on white upward mobility

	use ${data}/GM_cz_final_dataset.dta, clear

	* Multiply outcomes by 100 to be in percentile units
	foreach var of varlist kfr_white_pooled* kfr_black_pooled*{
	replace `var'=100*`var'
	}

	* Calculate mean of Black upward mobility in the sample
	foreach p in "p25" "p50" "p75"{
	sum kfr_black_pooled_`p', d
	local black_mean_`p'=r(mean)	
	}

	* Create a counterfactual dataset where GM is lowest percentile across all CZs in sample
	replace GM_hat2=1
	replace GM=1
	g id=_n+130
	
	* Replace Black upward mobility as missing in the counterfactual dataset
	foreach r in "black" "white"{
	foreach p in "p25" "p50" "p75"{
	replace kfr_`r'_pooled_`p'=.
	}
	}

	* Calculate racial gap (which will be missing for all values) in the counterfactual dataset
	foreach p in "p25" "p50" "p75"{
	g racegap2015_`p' = kfr_white_pooled_`p'-kfr_black_pooled_`p'
	}
	
	tempfile GM_cf
	
	save `GM_cf'

	* Open a new copy of the original dataset for generating the fitted values
	use ${data}/GM_cz_final_dataset.dta, clear
	
	* Multiply outcomes by 100 to be in percentile units	
	foreach var of varlist kfr_white_pooled* kfr_black_pooled*{
	replace `var'=100*`var'
	}

	* Calculate racial gap in upward mobility in each CZ
	foreach p in "p25" "p50" "p75"{
	g racegap2015_`p' = kfr_white_pooled_`p'-kfr_black_pooled_`p'
	}
	
	g id=_n
	append using `GM_cf'
	
	* 2SLS regression of Black upward mobility on GM; generate fitted vaues
	foreach r in "black" "white"{
	foreach p in "p25" "p50" "p75"{
	ivreg2 kfr_`r'_pooled_`p' ($x_ols = $x_iv) ${baseline_controls}	if kfr_`r'_pooled_`p' !=.  // Only uses original sample in the 2SLS regression
	predict kfr_`r'_pooled_`p'_pred, xb // Generates fitted values of Black upward mobility for all obs using X vars & regression coefficients from step above
	predict kfr_`r'_pooled_`p'_pred_se, stdp	 // Generates standard errors for above fitted values 
	}
	}
	
	* Calculate actual average racial gap in sample CZs
	foreach p in "p25" "p50" "p75"{
	sum racegap2015_`p' if id<131, d           // Only calculate on original sample (not counterfactual sample)
	local rgp_mean_`p'=r(mean)					// Store mean in original sample

	* Compute counterfactual gap by predicting Black outcomes if lowest percentile of GM
	g racegap2015_`p'_pred = kfr_white_pooled_`p'_pred-kfr_black_pooled_`p'_pred  // Calculate racial gap using fitted values	
	sum racegap2015_`p'_pred if id>130, d
	g rgpcf_mean_`p'=r(mean)
	local rgpcf_mean_`p'=r(mean)
	sum kfr_black_pooled_`p'_pred, d
	local blackcf_north_`p'=r(mean)
	}

	* Generate standard errors for counterfactual racial gap
	keep if id>130

	* Variance of counterfactual racial gaps
	foreach p in "p25" "p50" "p75"{	
	egen var_rg_cf_`p'=sum((racegap2015_`p'_pred - rgpcf_mean_`p')^2)
	replace var_rg_cf_`p'=(var_rg_cf_`p')/(_N)
	g sqrt_var_rg_cf_`p' = sqrt(var_rg_cf_`p')
	}

	foreach p in "p25" "p50" "p75"{
	local text_`p' : di%4.2f `rgp_mean_`p''
	tempname number
	file open `number' using "$figtab/text/rg`p'_alt.txt", text write replace
	file write `number' `"`text_`p''"' 
    file close `number'	
	local text_`p' : di%4.2f `rgpcf_mean_`p''	
	file open `number' using "$figtab/text/rgcf`p'_alt.txt", text write replace
	file write `number' `"`text_`p''"' 
    file close `number'	
	local text_`p' : di%4.2f `sqrt_var_rg_cf_`p''
	file open `number' using "$figtab/text/rg_se`p'_alt.txt", text write replace
	file write `number' `"`text_`p''"' 
    file close `number'	
	local change`p' = (`rgp_mean_`p''-`rgpcf_mean_`p'')*100/`rgp_mean_`p''
	local text_`p' : di %4.0f `change`p''
	file open `number' using "$figtab/text/change`p'_alt.txt", text write replace
	file write `number' `"`text_`p''%"' 
    file close `number'		
	}
	
	foreach p in "p25" "p50" "p75"{
	di ""
	di %4.2f `rgp_mean_`p''
	di %4.2f `rgpcf_mean_`p''	
	di %4.2f sqrt_var_rg_cf_`p'
	di  %4.0f (`rgp_mean_`p''-`rgpcf_mean_`p'')*100/`rgp_mean_`p''
	di ""
	di %4.2f `black_mean_`p''
	di %4.2f `blackcf_mean_`p''
	}
	
	* Remove any tempfiles and temp folders 
	cd "$figtab"
	
	local datafiles: dir "$figtab" files "*.dta"
	foreach datafile of local datafiles {
			rm `datafile'
	}

	shell rmdir $figtab/temp
	
*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
	* Table 3: Placebo test of identification strategy using pre-1940 upward mobility and educational attainment
*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
	local outcomes "enrolled_occs1900 enrolled_occs1910 enrolled_occs1920 enrolled_occs1930 enrolled_occs1940  avg_wt_med_educ1940"
	eststo clear
	foreach outcome in `outcomes'{
		use ${data}/GM_cz_final_dataset.dta, clear
		qui sum GM if `outcome'!=.
		local ols_SD=`r(sd)'
		qui: reg `outcome' $x_iv ${baseline_controls}
		sum `outcome' if e(sample)
		estadd scalar basemean=r(mean)
		estadd scalar sd=r(sd)
		eststo `outcome'           
		estadd local hascontrols 	"Y"
		estadd scalar gm_sd=`ols_SD'
		}
	esttab `outcomes' using $figtab/table1.tex, tex replace label nonote nocons nonum  ///
	posthead("&\multicolumn{5}{c}{}&\multicolumn{1}{c}{Median}\\" ///
	"&\multicolumn{5}{c}{Percentage of teens with low}&\multicolumn{1}{c}{adult}\\" ///
	"&\multicolumn{5}{c}{occ. score fathers attending school}&\multicolumn{1}{c}{education}\\" ///
	 "&\multicolumn{1}{c}{1900}&\multicolumn{1}{c}{1910}&\multicolumn{1}{c}{1920}&\multicolumn{1}{c}{1930}&\multicolumn{1}{c}{1940}&\multicolumn{1}{c}{1940}\\" )  ///
	se(%8.3f) b(%8.3f) nostar nofloat nonum nomtitles ///
	keep(GM_hat2) stats(basemean sd gm_sd N hascontrols,fmt(%8.3f %8.3f %8.3f %15.0gc  ) label("Baseline mean" "SD Dep Var" "SD GM" "Observations" "Baseline Controls")) 

*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
	* Table 4: Lower average upward mobility today for low income families in Great Migration CZs 
*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%	
	use ${data}/GM_cz_final_dataset.dta, clear
	
	foreach p in "25" {
	local inc_level "Lower Inc"	
	if "`p'"=="75"{
	local inc_level "Higher Inc"
	}
	
	local outcomes perm_res_p`p'_kr26 perm_res_p`p'_kr26_f perm_res_p`p'_kr26_m perm_res_p`p'_kir26 perm_res_p`p'_kir26_f perm_res_p`p'_kir26_m

	* First Stage	
	eststo clear
	foreach y in `outcomes'{
	use ${data}/GM_cz_final_dataset.dta, clear
	qui sum GM if `y'!=.
	local ols_SD=`r(sd)'
	local ols_SD= `r(sd)'
	di %4.0f `ols_SD'
	PrintEst `ols_SD' "ols_SD" "" "%" "4.0"
	

	reg $x_ols $x_iv  ${baseline_controls}
	eststo `y'
	test $x_iv = 0
	estadd scalar fstat=`r(F)'
	}
	cd "$figtab"
	esttab `outcomes' using "permres_table_p`p'.tex", frag replace  varwidth(25) label se ///
	stats(fstat, labels(  F-Stat)) keep($x_iv) mgroups("\textit{First Stage on GM}", pattern(1 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) nonumber ///
	nostar nomtitle nonotes nolines nogaps prefoot(\cmidrule(lr){2-7}) postfoot(\cmidrule(lr){2-7}) substitute({table} {threeparttable}) 
	
	* OLS 	
	eststo clear
	foreach y in `outcomes'{
	use ${data}/GM_cz_final_dataset.dta, clear
	qui sum GM if `y'!=.
	local ols_SD=`r(sd)'
	reg `y' $x_ols  ${baseline_controls}
	local olsLB`y'=_b[$x_ols] -1.96*_se[$x_ols]
	local olsUB`y'=_b[$x_ols] +1.96*_se[$x_ols]		
	PrintEst `olsLB`y'' "olsLB`y'" "" "%" "4.3"
	PrintEst `olsUB`y'' "olsUB`y'" "" "%" "4.3"		
	eststo `y'
	}
	cd "$figtab"
	esttab `outcomes' using "permres_table_p`p'.tex", frag append  varwidth(25) label se ///
	prehead("\\" "&\multicolumn{3}{c}{Household Income Rank}&\multicolumn{3}{c}{Individual Income Rank}\\" ///
	"&\multicolumn{1}{c}{Pooled}&\multicolumn{1}{c}{Women}&\multicolumn{1}{c}{Men}&\multicolumn{1}{c}{Pooled}&\multicolumn{1}{c}{Women}&\multicolumn{1}{c}{Men} \\\cmidrule(lr){2-7}")  ///
	stats( r2, labels( R-squared)) keep($x_ols) mgroups("\textit{Ordinary Least Squares}", pattern(1 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) nonumber ///
	nostar nomtitle nonotes nolines nogaps prefoot(\cmidrule(lr){2-7}) postfoot(\cmidrule(lr){2-7}) substitute({table} {threeparttable}) 

	* RF 	
	eststo clear
	foreach y in `outcomes'{
	use ${data}/GM_cz_final_dataset.dta, clear
	qui sum GM if `y'!=.
	local ols_SD=`r(sd)'
	reg `y' $x_iv  ${baseline_controls}
	eststo `y'
	}
	cd "$figtab"
	esttab `outcomes' using "permres_table_p`p'.tex", frag append  varwidth(25) label se ///
	stats( r2, labels( R-squared)) keep($x_iv) mgroups("\textit{Reduced Form}", pattern(1 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) nonumber  ///
	nostar nomtitle nonotes nolines nogaps prefoot(\cmidrule(lr){2-7}) postfoot(\cmidrule(lr){2-7}) substitute({table} {threeparttable}) 
	
	* 2SLS 
	eststo clear
	foreach y in `outcomes'{
	use ${data}/GM_cz_final_dataset.dta, clear
	qui sum GM if `y'!=.
	local ols_SD=`r(sd)'
	ivreg2 `y' ($x_ols = $x_iv )  ${baseline_controls}, first
	local GM_`y' = _b[$x_ols]
	local GM_`y'_abs = abs(_b[$x_ols])
	local GM_`y'_se : di %4.3f _se[$x_ols]
	local ivLB`y'=_b[$x_ols] -1.96*_se[$x_ols]
	local ivUB`y'=_b[$x_ols] +1.96*_se[$x_ols]		
	PrintEst `ivLB`y'' "ivLB`y'" "" "%" "4.3"
	PrintEst `ivUB`y'' "ivUB`y'" "" "%" "4.3"	
	PrintEst `GM_`y'' "GM_`y'" "" " percentile points (s.e. = `GM_`y'_se')%" "4.3"
	PrintEst `GM_`y'_abs' "GM_`y'_abs" "" " percentile points (s.e. = `GM_`y'_se')%" "4.3"
	eststo `y'
	sum `y' if e(sample) 
	estadd scalar basemean=r(mean)
	estadd scalar sd=r(sd)	
	estadd scalar gm_sd=`ols_SD'
	}
	
	cd "$figtab"
	esttab `outcomes' using "permres_table_p`p'.tex", frag append  varwidth(25) label se ///
	stats(none) keep($x_ols) mgroups("\textit{Two-stage least squares}", pattern(1 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) nonumber ///
	nostar nomtitle  nonotes nolines nogaps prefoot(\cmidrule(lr){2-7}) substitute({table} {threeparttable}) 
	
	* Footer
	cd "$figtab"
	esttab `outcomes' using "permres_table_p`p'.tex", frag append  varwidth(25) label se ///
	stats( N  basemean sd gm_sd, labels(N "Mean Rank" "SD Rank" "SD GM")) drop(*) nonumber ///
	nostar nomtitle  nonotes nolines nogaps prefoot(\cmidrule(lr){2-7}) substitute({table} {threeparttable}) 
	}

*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
	* Table 5: Childhood exposure to Great Migration CZs lowers upward mobility for low income families 
*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%	
	use ${data}/GM_cz_final_dataset.dta, clear
	
	foreach p in "25" {
	local inc_level "Lower Inc"	
	if "`p'"=="75"{
	local inc_level "Higher Inc"
	}
	
	local outcomes causal_p`p'_czkr26 causal_p`p'_czkr26_f causal_p`p'_czkr26_m causal_p`p'_czkir26 causal_p`p'_czkir26_f causal_p`p'_czkir26_m

	* First Stage	
	eststo clear
	foreach y in `outcomes'{
	use ${data}/GM_cz_final_dataset.dta, clear
	g `y'_wt=1/(`y'_se^2)
	qui sum GM if `y'!=. [aw=`y'_wt]
	local ols_SD=`r(sd)'
	reg $x_ols $x_iv   ${baseline_controls} [aw=`y'_wt]
	eststo `y'
	test $x_iv = 0
	estadd scalar fstat=`r(F)'
	}
	cd "$figtab"
	esttab `outcomes' using "main_table_p`p'.tex", frag replace  varwidth(25) label se ///
	stats(fstat, labels(  F-Stat)) keep($x_iv) mgroups("\textit{First Stage on GM}", pattern(1 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) nonumber ///
	nostar nomtitle nonotes nolines nogaps prefoot(\cmidrule(lr){2-7}) postfoot(\cmidrule(lr){2-7}) substitute({table} {threeparttable}) 
	
	* OLS 	
	eststo clear
	foreach y in `outcomes'{
	use ${data}/GM_cz_final_dataset.dta, clear
	g `y'_wt=1/(`y'_se^2)
	qui sum GM if `y'!=. [aw=`y'_wt]
	local ols_SD=`r(sd)'
	reg `y' $x_ols  ${baseline_controls} [aw=`y'_wt]
	local olsLB`y'=_b[$x_ols] -1.96*_se[$x_ols]
	local olsUB`y'=_b[$x_ols] +1.96*_se[$x_ols]		
	PrintEst `olsLB`y'' "olsLB`y'" "" "%" "12.4"
	PrintEst `olsUB`y'' "olsUB`y'" "" "%" "12.4"		
	eststo `y'
	}
	cd "$figtab"
	esttab `outcomes' using "main_table_p`p'.tex", frag append  varwidth(25) label se ///
	prehead("\\" "&\multicolumn{3}{c}{Household Income Rank}&\multicolumn{3}{c}{Individual Income Rank}\\" ///
	"&\multicolumn{1}{c}{Pooled}&\multicolumn{1}{c}{Women}&\multicolumn{1}{c}{Men}&\multicolumn{1}{c}{Pooled}&\multicolumn{1}{c}{Women}&\multicolumn{1}{c}{Men} \\\cmidrule(lr){2-7}")  ///
	stats( r2, labels( R-squared)) keep($x_ols) mgroups("\textit{Ordinary Least Squares}", pattern(1 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) nonumber ///
	nostar nomtitle nonotes nolines nogaps prefoot(\cmidrule(lr){2-7}) postfoot(\cmidrule(lr){2-7}) substitute({table} {threeparttable}) 

	* RF 	
	eststo clear
	foreach y in `outcomes'{
	use ${data}/GM_cz_final_dataset.dta, clear
	g `y'_wt=1/(`y'_se^2)
	qui sum GM if `y'!=. [aw=`y'_wt]
	local ols_SD=`r(sd)'
	reg `y' $x_iv  ${baseline_controls} [aw=`y'_wt]
	eststo `y'
	}
	cd "$figtab"
	esttab `outcomes' using "main_table_p`p'.tex", frag append  varwidth(25) label se ///
	stats( r2, labels( R-squared)) keep($x_iv) mgroups("\textit{Reduced Form}", pattern(1 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) nonumber  ///
	nostar nomtitle nonotes nolines nogaps prefoot(\cmidrule(lr){2-7}) postfoot(\cmidrule(lr){2-7}) substitute({table} {threeparttable}) 
	
	* 2SLS 	
	eststo clear
	foreach y in `outcomes'{
	use ${data}/GM_cz_final_dataset.dta, clear
	g `y'_wt=1/(`y'_se^2)
	qui sum GM if `y'!=. [aw=`y'_wt]
	local ols_SD=`r(sd)'
	ivreg2 `y' ($x_ols = $x_iv )  ${baseline_controls} [aw=`y'_wt], first
	local GM_`y' = _b[$x_ols]
	local GM_`y'_abs = abs(_b[$x_ols])
	local GM_`y'_se : di %6.4f _se[$x_ols]
	local ivLB`y'=_b[$x_ols] -1.96*_se[$x_ols]
	local ivUB`y'=_b[$x_ols] +1.96*_se[$x_ols]		
	PrintEst `ivLB`y'' "ivLB`y'" "" "%" "12.4"
	PrintEst `ivUB`y'' "ivUB`y'" "" "%" "12.4"	
	PrintEst `GM_`y'' "GM_`y'" "" " percentile points (s.e. = `GM_`y'_se')%" "6.4"
	PrintEst `GM_`y'_abs' "GM_`y'_abs" "" " percentile points (s.e. = `GM_`y'_se')%" "6.4"
	eststo `y'
	sum `y' if e(sample) [aw=`y'_wt]
	estadd scalar basemean=r(mean)
	estadd scalar sd=r(sd)	
	estadd scalar gm_sd=`ols_SD'
	estadd local precisionwt "Y" 
	}
	
	cd "$figtab"
	esttab `outcomes' using "main_table_p`p'.tex", frag append  varwidth(25) label se ///
	stats(none) keep($x_ols) mgroups("\textit{Two-stage least squares}", pattern(1 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) nonumber ///
	nostar nomtitle  nonotes nolines nogaps prefoot(\cmidrule(lr){2-7}) substitute({table} {threeparttable}) 
	
	* Footer
	cd "$figtab"
	esttab `outcomes' using "main_table_p`p'.tex", frag append  varwidth(25) label se ///
	stats( N precisionwt basemean sd gm_sd, labels(N "Precision Wt" "Mean Expos FX" "SD Expos FX" "SD GM")) drop(*) nonumber ///
	nostar nomtitle  nonotes nolines nogaps prefoot(\cmidrule(lr){2-7}) substitute({table} {threeparttable}) 
	}


*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
	* Table 6: Lower average upward mobility today for black households in Great Migration CZs 
	* &
	* Table 7: No Great Migration impact on average upward mobility of white households today
*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%	
	foreach r in "black" "white"{
	local race "Black"	
	if "`r'"=="black"{
	local race "White"
	}
	
	local outcomes kfr_`r'_pooled_p252015 kir_`r'_female_p252015 kir_`r'_male_p252015 kfr_`r'_pooled_p752015 kir_`r'_female_p752015 kir_`r'_male_p752015

	* First Stage	
	eststo clear
	foreach y in `outcomes'{
	use ${data}/GM_cz_final_dataset.dta, clear
	keep if `y'!=.
	* Rescale treatment in terms of standard deviations
	qui sum GM
	local ols_SD=`r(sd)'
	* Rescale outcome in percentile ranks
	replace `y'=100*`y'
	reg $x_ols $x_iv   ${baseline_controls} if `y'!=.
	eststo `y'
	test $x_iv = 0
	estadd scalar fstat=`r(F)'
	}
	cd "${figtab}"
	esttab `outcomes' using "`r'_hh_table.tex", frag replace  varwidth(25) label se ///
	stats(fstat, labels(  F-Stat)) keep($x_iv) mgroups("\textit{First Stage on GM}", pattern(1 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) nonumber ///
	nostar nomtitle nonotes nolines nogaps prefoot(\cmidrule(lr){2-7}) postfoot(\cmidrule(lr){2-7}) substitute({table} {threeparttable}) 
	
	* OLS 	
	eststo clear
	foreach y in `outcomes'{
	use ${data}/GM_cz_final_dataset.dta, clear
	keep if `y'!=.
	* Rescale treatment in terms of standard deviations
	qui sum GM
	local ols_SD=`r(sd)'
	* Rescale outcome in percentile ranks
	replace `y'=100*`y'
	reg `y' $x_ols  ${baseline_controls}  if `y'!=.
	eststo `y'
	}
	cd "${figtab}"
	esttab `outcomes' using "`r'_hh_table.tex", frag append  varwidth(25) label se ///
	prehead("\\" "&\multicolumn{3}{c}{Low Income}&\multicolumn{3}{c}{High Income}\\" ///
	"&\multicolumn{1}{c}{Pooled}&\multicolumn{1}{c}{Women}&\multicolumn{1}{c}{Men}&\multicolumn{1}{c}{Pooled}&\multicolumn{1}{c}{Women}&\multicolumn{1}{c}{Men} \\\cmidrule(lr){2-7}")  ///
	stats( r2, labels( R-squared)) keep($x_ols) mgroups("\textit{Ordinary Least Squares}", pattern(1 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) nonumber ///
	nostar nomtitle nonotes nolines nogaps prefoot(\cmidrule(lr){2-7}) postfoot(\cmidrule(lr){2-7}) substitute({table} {threeparttable}) 

	* RF 	
	eststo clear
	foreach y in `outcomes'{
	use ${data}/GM_cz_final_dataset.dta, clear
	keep if `y'!=.
	* Rescale treatment in terms of standard deviations
	qui sum GM
	local ols_SD=`r(sd)'
	* Rescale outcome in percentile ranks
	replace `y'=100*`y'
	reg `y' $x_iv  ${baseline_controls} if `y'!=.
	eststo `y'
	}
	cd "${figtab}"
	esttab `outcomes' using "`r'_hh_table.tex", frag append  varwidth(25) label se ///
	stats( r2, labels( R-squared)) keep($x_iv) mgroups("\textit{Reduced Form}", pattern(1 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) nonumber  ///
	nostar nomtitle nonotes nolines nogaps prefoot(\cmidrule(lr){2-7}) postfoot(\cmidrule(lr){2-7}) substitute({table} {threeparttable}) 
	
	* 2SLS 	
	eststo clear
	foreach y in `outcomes'{
	use ${data}/GM_cz_final_dataset.dta, clear
	keep if `y'!=.
	* Rescale treatment in terms of standard deviations
	qui sum GM
	local ols_SD=`r(sd)'
	* Rescale outcome in percentile ranks
	replace `y'=100*`y'
	ivreg2 `y' ($x_ols = $x_iv )  ${baseline_controls}  if `y'!=., first
	local GM_`y' = _b[$x_ols]
	local GM_`y'_abs = abs(_b[$x_ols])
	local GM_`y'_se : di %4.3f _se[$x_ols]
	PrintEst `GM_`y'' "GM_`y'" "" " percentile points (s.e. = `GM_`y'_se')%" "4.3"
	PrintEst `GM_`y'_abs' "GM_`y'_abs" "" " percentile points (s.e. = `GM_`y'_se')%" "4.3"
	eststo `y'
	use ${data}/GM_cz_final_dataset.dta, clear
	keep if `y'!=.
	sum `y' 
	estadd scalar basemean=r(mean)
	estadd scalar sd=r(sd)	
	estadd scalar gm_sd=`ols_SD'
	}
	
	cd "${figtab}"
	esttab `outcomes' using "`r'_hh_table.tex", frag append  varwidth(25) label se ///
	stats(none) keep($x_ols) mgroups("\textit{Two-stage least squares}", pattern(1 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) nonumber ///
	nostar nomtitle  nonotes nolines nogaps prefoot(\cmidrule(lr){2-7}) substitute({table} {threeparttable}) 
	
	* Footer
	cd "${figtab}"
	esttab `outcomes' using "`r'_hh_table.tex", frag append  varwidth(25) label se ///
	stats( N  basemean sd gm_sd, labels(N  "Mean Rank" "SD Rank" "SD GM")) drop(*) nonumber ///
	nostar nomtitle  nonotes nolines nogaps prefoot(\cmidrule(lr){2-7}) substitute({table} {threeparttable}) 
	}
	

*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
	* Table 8: Robustness of effects of childhood exposure to Great Migration CZs
*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%	
	global nocon "v2_blackmig3539_share1940"
	global divfe "v2_blackmig3539_share1940 reg2 reg3 reg4"	
	global baseline "frac_all_upm1940 mfg_lfshare1940 v2_blackmig3539_share1940 reg2 reg3 reg4"
	global emp "frac_all_upm1940 mfg_lfshare1940 v2_blackmig3539_share1940 reg2 reg3 reg4 emp_hat"
	global flexbpop40 "frac_all_upm1940 mfg_lfshare1940 v2_blackmig3539_share1940 reg2 reg3 reg4 i.bpopquartile"		
	global swmig "frac_all_upm1940 mfg_lfshare1940 v2_blackmig3539_share1940 reg2 reg3 reg4  GM_hat8"
	global eurmig "frac_all_upm1940 mfg_lfshare1940 v2_blackmig3539_share1940 reg2 reg3 reg4  eur_mig"
	global supmob "frac_all_upm1940 mfg_lfshare1940 v2_blackmig3539_share1940 reg2 reg3 reg4  vm_blacksmob1940"

	local y  causal_p25_czkr26 	

	eststo clear
	use ${data}/GM_cz_final_dataset.dta, clear
	g `y'_wt=1/(`y'_se^2)	

	reg $x_ols $x_iv $nocon [aw=`y'_wt]
	test $x_iv = 0
	local fstat=`r(F)'
	
	reg `y' $x_ols $nocon  [aw=`y'_wt]
	eststo nocon_ols
	estadd local hasdivfe 		"N"
	estadd local hasbaseline 	"N"
	estadd local hasemp		 	"N"
	estadd local hasflexbpop40	"N"
	estadd local hasswmig		"N"
	estadd local haseurmig		"N"
	estadd local hassupmob		"N"	
	estadd local precisionwt 	"Y" 

	ivreg2 `y' ($x_ols = $x_iv ) $nocon  [aw=`y'_wt]
	eststo nocon
	estadd scalar fstat=`fstat'
		
	reg $x_ols $x_iv $divfe [aw=`y'_wt]
	test $x_iv = 0
	local fstat=`r(F)'

	reg `y' $x_ols $divfe    [aw=`y'_wt]
	eststo divfe_ols
	estadd local hasdivfe 		"Y"
	estadd local hasbaseline 	"N"
	estadd local hasemp 		"N"
	estadd local hasflexbpop40	"N"
	estadd local hasswmig		"N"
	estadd local haseurmig		"N"
	estadd local hassupmob		"N"	
	estadd local precisionwt 	"Y" 
	
	ivreg2 `y' ($x_ols = $x_iv )  $divfe    [aw=`y'_wt]
	eststo divfe
	estadd scalar fstat=`fstat'
	
	reg $x_ols $x_iv $baseline [aw=`y'_wt]
	test $x_iv = 0
	local fstat=`r(F)'

	reg `y' $x_ols $baseline [aw=`y'_wt]
	eststo baseline_ols
	estadd local hasdivfe 		"Y"
	estadd local hasbaseline 	"Y"
	estadd local hasemp 		"N"
	estadd local hasflexbpop40	"N"
	estadd local hasswmig		"N"
	estadd local haseurmig		"N"
	estadd local hassupmob		"N"	
	estadd local precisionwt 	"Y" 

	ivreg2 `y' ($x_ols = $x_iv ) $baseline [aw=`y'_wt]
	eststo baseline
	estadd scalar fstat=`fstat'
	
	reg $x_ols $x_iv $emp [aw=`y'_wt]
	test $x_iv = 0
	local fstat=`r(F)'

	reg `y' $x_ols $emp  [aw=`y'_wt]
	eststo emp_ols
	estadd local hasdivfe 		"Y"
	estadd local hasbaseline 	"Y"
	estadd local hasemp 		"Y"
	estadd local hasflexbpop40	"N"
	estadd local hasswmig		"N"
	estadd local haseurmig		"N"
	estadd local hassupmob		"N"	
	estadd local precisionwt 	"Y" 

	ivreg2 `y' ($x_ols = $x_iv ) $emp  [aw=`y'_wt]
	eststo emp
	estadd scalar fstat=`fstat'

	reg $x_ols $x_iv $flexbpop40 [aw=`y'_wt]
	test $x_iv = 0
	local fstat=`r(F)'

	reg `y' $x_ols $flexbpop40  [aw=`y'_wt]
	eststo flexbpop40_ols
	estadd local hasdivfe 		"Y"
	estadd local hasbaseline 	"Y"
	estadd local hasemp 		"N"
	estadd local hasflexbpop40	"Y"
	estadd local hasswmig		"N"
	estadd local haseurmig		"N"
	estadd local hassupmob		"N"	
	estadd local precisionwt 	"Y"

	ivreg2 `y' ($x_ols = $x_iv ) $flexbpop40  [aw=`y'_wt]
	eststo flexbpop40
	estadd scalar fstat=`fstat'
	
	reg $x_ols $x_iv $swmig [aw=`y'_wt]
	test $x_iv = 0
	local fstat=`r(F)'

	reg `y' $x_ols $swmig  [aw=`y'_wt]
	eststo swmig_ols
	estadd local hasdivfe 		"Y"
	estadd local hasbaseline 	"Y"
	estadd local hasemp 		"N"
	estadd local hasflexbpop40	"N"
	estadd local hasswmig		"Y"
	estadd local haseurmig		"N"
	estadd local hassupmob		"N"	
	estadd local precisionwt 	"Y" 

	ivreg2 `y' ($x_ols = $x_iv ) $swmig  [aw=`y'_wt]
	eststo swmig
	estadd scalar fstat=`fstat'

	reg $x_ols $x_iv $eurmig [aw=`y'_wt]
	test $x_iv = 0
	local fstat=`r(F)'

	reg `y' $x_ols $eurmig  [aw=`y'_wt]
	eststo eurmig_ols
	estadd local hasdivfe 		"Y"
	estadd local hasbaseline 	"Y"
	estadd local hasemp 		"N"
	estadd local hasflexbpop40	"N"
	estadd local hasswmig	 	"N"
	estadd local haseurmig		"Y"
	estadd local hassupmob		"N"	
	estadd local precisionwt 	"Y" 
	
	ivreg2 `y' ($x_ols = $x_iv )  $eurmig  [aw=`y'_wt]
	eststo eurmig
	estadd scalar fstat=`fstat'

	reg $x_ols $x_iv $supmob [aw=`y'_wt]
	test $x_iv = 0
	local fstat=`r(F)'

	reg `y' $x_ols $supmob  [aw=`y'_wt]
	eststo supmob_ols
	estadd local hasdivfe 		"Y"
	estadd local hasbaseline 	"Y"
	estadd local hasemp			"N"
	estadd local hasflexbpop40	"N"
	estadd local hasswmig		"N"
	estadd local haseurmig		"N"
	estadd local hassupmob		"Y"	
	estadd local precisionwt 	"Y" 
	
	ivreg2 `y' ($x_ols = $x_iv ) $supmob  [aw=`y'_wt]
	eststo supmob
	estadd scalar fstat= `fstat'
	
	cd "$figtab"
	esttab nocon divfe baseline flexbpop40 supmob  swmig eurmig emp using "main_robust_table_p25.tex", frag replace varwidth(25) label se ///
	stats( fstat, labels("First Stage F-Stat")) keep($x_ols) coeflabel(GM "GM (2SLS)")  nonumber ///
	nostar nomtitle nonotes nolines nogaps  substitute({table} {threeparttable}) prefoot(\cmidrule(lr){2-9})
	
	esttab nocon_ols divfe_ols baseline_ols flexbpop40_ols supmob_ols  swmig_ols eurmig_ols emp_ols  using "main_robust_table_p25.tex", frag append  varwidth(25) label se ///
	prehead("\\") coeflabel(GM "GM (OLS)") ///
	stats( r2  N  precisionwt hasdivfe hasbaseline hasflexbpop40 hassupmob hasswmig haseurmig hasemp , ///
	labels( "R-squared (OLS)" N "Precision Wt" "Census Div FE" "Baseline Controls"  "1940 Black Share Quartile FEs" "Southern Mob" ///
	"White South Mig" "Eur Mig"  "Emp Bartik" )) keep($x_ols)  nonumber  ///
	nostar nomtitle nonotes nolines prefoot(\cmidrule(lr){2-9}) postfoot(\cmidrule(lr){2-9})  substitute({table} {threeparttable}) 

	
*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
	* Table 9: Robustness of Great Migration's effects on black men's upward mobility
*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%	
	global nocon "v2_blackmig3539_share1940"
	global divfe "v2_blackmig3539_share1940 reg2 reg3 reg4"	
	global baseline "frac_all_upm1940 mfg_lfshare1940 v2_blackmig3539_share1940 reg2 reg3 reg4"
	global emp "frac_all_upm1940 mfg_lfshare1940 v2_blackmig3539_share1940 reg2 reg3 reg4 emp_hat"
	global flexbpop40 "frac_all_upm1940 mfg_lfshare1940 v2_blackmig3539_share1940 reg2 reg3 reg4 i.bpopquartile"		
	global swmig "frac_all_upm1940 mfg_lfshare1940 v2_blackmig3539_share1940 reg2 reg3 reg4  GM_hat8"
	global eurmig "frac_all_upm1940 mfg_lfshare1940 v2_blackmig3539_share1940 reg2 reg3 reg4  eur_mig"
	global supmob "frac_all_upm1940 mfg_lfshare1940 v2_blackmig3539_share1940 reg2 reg3 reg4  vm_blacksmob1940"

	local y  kir_black_male_p252015 	
	
	eststo clear
	use ${data}/GM_cz_final_dataset.dta, clear
	replace `y'=100*`y'
	
	reg $x_ols $x_iv $nocon 
	test $x_iv = 0
	local fstat=`r(F)'
	
	reg `y' $x_ols $nocon  
	eststo nocon_ols
	estadd local hasdivfe 		"N"
	estadd local hasbaseline 	"N"
	estadd local hasemp		 	"N"
	estadd local hasflexbpop40	"N"
	estadd local hasswmig		"N"
	estadd local haseurmig		"N"
	estadd local hassupmob		"N"	
	estadd local precisionwt 	"Y" 

	ivreg2 `y' ($x_ols = $x_iv ) $nocon  
	eststo nocon
	estadd scalar fstat=`fstat'
		
	reg $x_ols $x_iv $divfe 
	test $x_iv = 0
	local fstat=`r(F)'

	reg `y' $x_ols $divfe    
	eststo divfe_ols
	estadd local hasdivfe 		"Y"
	estadd local hasbaseline 	"N"
	estadd local hasemp 		"N"
	estadd local hasflexbpop40	"N"
	estadd local hasswmig		"N"
	estadd local haseurmig		"N"
	estadd local hassupmob		"N"	
	estadd local precisionwt 	"Y" 
	
	ivreg2 `y' ($x_ols = $x_iv )  $divfe    
	eststo divfe
	estadd scalar fstat=`fstat'
	
	reg $x_ols $x_iv $baseline 
	test $x_iv = 0
	local fstat=`r(F)'

	reg `y' $x_ols $baseline 
	eststo baseline_ols
	estadd local hasdivfe 		"Y"
	estadd local hasbaseline 	"Y"
	estadd local hasemp 		"N"
	estadd local hasflexbpop40	"N"
	estadd local hasswmig		"N"
	estadd local haseurmig		"N"
	estadd local hassupmob		"N"	
	estadd local precisionwt 	"Y" 

	ivreg2 `y' ($x_ols = $x_iv ) $baseline 
	eststo baseline
	estadd scalar fstat=`fstat'
	
	reg $x_ols $x_iv $emp 
	test $x_iv = 0
	local fstat=`r(F)'

	reg `y' $x_ols $emp  
	eststo emp_ols
	estadd local hasdivfe 		"Y"
	estadd local hasbaseline 	"Y"
	estadd local hasemp 		"Y"
	estadd local hasflexbpop40	"N"
	estadd local hasswmig		"N"
	estadd local haseurmig		"N"
	estadd local hassupmob		"N"	
	estadd local precisionwt 	"Y" 

	ivreg2 `y' ($x_ols = $x_iv ) $emp  
	eststo emp
	estadd scalar fstat=`fstat'

	reg $x_ols $x_iv $flexbpop40 
	test $x_iv = 0
	local fstat=`r(F)'

	reg `y' $x_ols $flexbpop40  
	eststo flexbpop40_ols
	estadd local hasdivfe 		"Y"
	estadd local hasbaseline 	"Y"
	estadd local hasemp 		"N"
	estadd local hasflexbpop40	"Y"
	estadd local hasswmig		"N"
	estadd local haseurmig		"N"
	estadd local hassupmob		"N"	
	estadd local precisionwt 	"Y"

	ivreg2 `y' ($x_ols = $x_iv ) $flexbpop40  
	eststo flexbpop40
	estadd scalar fstat=`fstat'
	
	reg $x_ols $x_iv $swmig 
	test $x_iv = 0
	local fstat=`r(F)'

	reg `y' $x_ols $swmig  
	eststo swmig_ols
	estadd local hasdivfe 		"Y"
	estadd local hasbaseline 	"Y"
	estadd local hasemp 		"N"
	estadd local hasflexbpop40	"N"
	estadd local hasswmig		"Y"
	estadd local haseurmig		"N"
	estadd local hassupmob		"N"	
	estadd local precisionwt 	"Y" 

	ivreg2 `y' ($x_ols = $x_iv ) $swmig  
	eststo swmig
	estadd scalar fstat=`fstat'

	reg $x_ols $x_iv $eurmig 
	test $x_iv = 0
	local fstat=`r(F)'

	reg `y' $x_ols $eurmig  
	eststo eurmig_ols
	estadd local hasdivfe 		"Y"
	estadd local hasbaseline 	"Y"
	estadd local hasemp 		"N"
	estadd local hasflexbpop40	"N"
	estadd local hasswmig	 	"N"
	estadd local haseurmig		"Y"
	estadd local hassupmob		"N"	
	estadd local precisionwt 	"Y" 
	
	ivreg2 `y' ($x_ols = $x_iv )  $eurmig  
	eststo eurmig
	estadd scalar fstat=`fstat'

	reg $x_ols $x_iv $supmob 
	test $x_iv = 0
	local fstat=`r(F)'

	reg `y' $x_ols $supmob  
	eststo supmob_ols
	estadd local hasdivfe 		"Y"
	estadd local hasbaseline 	"Y"
	estadd local hasemp			"N"
	estadd local hasflexbpop40	"N"
	estadd local hasswmig		"N"
	estadd local haseurmig		"N"
	estadd local hassupmob		"Y"	
	estadd local precisionwt 	"Y" 
	
	ivreg2 `y' ($x_ols = $x_iv ) $supmob  
	eststo supmob
	estadd scalar fstat= `fstat'
	
	cd "$figtab"
	esttab nocon divfe baseline flexbpop40 supmob  swmig eurmig emp using "main_robust_table_bmp25.tex", frag replace varwidth(25) label se ///
	stats( fstat, labels("First Stage F-Stat")) keep($x_ols) coeflabel(GM "GM (2SLS)")  nonumber ///
	nostar nomtitle nonotes nolines nogaps  substitute({table} {threeparttable}) prefoot(\cmidrule(lr){2-9})
	
	esttab nocon_ols divfe_ols baseline_ols flexbpop40_ols supmob_ols  swmig_ols eurmig_ols emp_ols  using "main_robust_table_bmp25.tex", frag append  varwidth(25) label se ///
	prehead("\\") coeflabel(GM "GM (OLS)") ///
	stats( r2  N  precisionwt hasdivfe hasbaseline hasflexbpop40 hassupmob hasswmig haseurmig hasemp , ///
	labels( "R-squared (OLS)" N "Precision Wt" "Census Div FE" "Baseline Controls"  "1940 Black Share Quartile FEs" "Southern Mob" ///
	"White South Mig" "Eur Mig"  "Emp Bartik" )) keep($x_ols)  nonumber  ///
	nostar nomtitle nonotes nolines prefoot(\cmidrule(lr){2-9}) postfoot(\cmidrule(lr){2-9})  substitute({table} {threeparttable}) 
			
*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%	
*V. Estimates cited in text.
*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------% 
	
	/* Share of non-southern continental US total and Black population in sample CZs
	use $data/GM_cz_final_dataset.dta, clear
	merge 1:1 cz using $population/clean_bpopshare_1870_2000_cz.dta, keepusing(cz totpop2000 bpop2000 stateabbrv)
	g original_sample=(_merge==3)
	drop _merge
	g south= (stateabbrv == "AL" |stateabbrv == "AR" |stateabbrv == "FL" |stateabbrv == "GA" |stateabbrv == "KY" |stateabbrv == "LA" |stateabbrv == "MS" |stateabbrv == "NC" | ///
	stateabbrv == "OK" |stateabbrv == "SC" |stateabbrv == "TN" |stateabbrv == "TX" |stateabbrv == "WV" |stateabbrv == "VA" )
	drop if south==1
	
	g no_cz=1
	
	collapse (sum) totpop2000 bpop2000 no_cz, by(original_sample)
	g id=1
	reshape wide totpop2000 bpop2000 no_cz, i(id) j(original_sample)
	
	g share_pop_original_sample = totpop20001/(totpop20001+totpop20000)
	g share_bpop_original_sample = bpop20001/(bpop20001+bpop20000)
	
	local share_pop_original_sample = share_pop_original_sample *100
	local share_bpop_original_sample = share_bpop_original_sample *100
	PrintEst `share_pop_original_sample' "share_pop_original_sample" "" "%" "4.0"	
	PrintEst `share_bpop_original_sample' "share_bpop_original_sample" "" "%" "4.0"	
	
	* Correlation between 1940 upward mobility and 2015 upward mobility for continental US (where both measures are available--721 CZs)
	use $mobdata/clean_cz_mobility_1900_2015.dta, clear
	corr frac_all_upm1940 kfr_pooled_pooled_p252015
	local mobmeasurecorr : di r(rho)
	PrintEst `mobmeasurecorr' "mobmeasurecorr" "" "%" "4.2"
	
	* Correlation between 2015 income and educational upward mobility for continental US
	use $mobdata/clean_cz_mobility_1900_2015.dta, clear
	corr kfr_pooled_pooled_p252015 hs_pooled_pooled_p25
	local incedumobcorr : di r(rho)
	PrintEst `incedumobcorr' "incedumobcorr" "" "%" "4.2"

	* Correlation between Great Migration and baseline 1940 covariates
	use ${data}/GM_cz_final_dataset.dta, clear
	corr GM frac_all_upm1940
	local gmeduupmcorr: di r(rho)
	PrintEst `gmeduupmcorr' "gmeduupmcorr" "" "%" "4.2"
	
	corr GM mfg_lfshare1940
	local gmmfgsharcorr: di r(rho)
	PrintEst `gmmfgsharcorr' "gmmfgsharcorr" "" "%" "4.2"
	
	corr GM v2_blackmig3539_share1940
	local gmrecblkmig: di r(rho)
	PrintEst `gmrecblkmig' "gmrecblkmig" "" "%" "4.2"	
	
	corr bpopchange1940_1970 bpopshare1940
	local bpopchangebpopshare1940: di r(rho)
	PrintEst `bpopchangebpopshare1940' "bpopchangebpopshare1940" "" "%" "4.2"	
	
	
	*Correlation between years of manufacturing share
	use ${data}/GM_cz_final_dataset.dta, clear
	
	corr mfg_lfshare1950 mfg_lfshare1940
	local mfgshar4050corr: di r(rho)
	PrintEst `mfgshar4050corr' "mfgshar4050corr" "" "%" "4.2"
	
	corr mfg_lfshare1970 mfg_lfshare1940
	local mfgshar4070corr: di r(rho)
	PrintEst `mfgshar4070corr' "mfgshar4070corr" "" "%" "4.2" */
	
	
	
	
	*/
	
	