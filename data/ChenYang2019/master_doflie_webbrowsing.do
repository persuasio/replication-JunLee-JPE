*** Do-file: replication of Chen and Yang (2018)
*** Part 2: based on web browsing data
*** December 2018


clear all
set more off


*** Define root name for data folder
	local directory_path 	`""'
	cd "`directory_path'/Outregs"

	

*** Prepare NYTimes wedge

	preserve
	use "`directory_path'/nytimes_newslog.dta", clear
	drop if date == .
	drop if date < 20443 | date > 20918
	gen censored_prop_day = censored_topic / total_article
	// code weekends and holidays as 0
	replace censored_prop_day = 0  	if censored_prop_day == .
	egen week = cut(date), at(20443(7)20919)
	bysort week: egen censored_prop_week = mean(censored_prop_day)
	keep week censored_prop_week
	duplicates drop
	drop if week == .
	save "`directory_path'/Outregs/nytimes_newslog_wedge.dta", replace
	


*** Prepare panel survey completion indicator

	preserve
	use "`directory_path'/panelsurvey_raw.dta", clear
	keep responseID_wave1 panelmerged_wave2 panelmerged_wave3
	save "`directory_path'/Outregs/panelmerged_indicator.dta", replace
	restore
	


*** Figure 1

	// load raw data
	use "`directory_path'/browsingtime_nytimes.dta", clear
	
	// reshape data (standardized variable name)
	drop if visit_day == .
	keep responseID_wave1 visit_day time_spent_min_nytimes
	reshape wide time_spent_min_nytimes, i(responseID_wave1) j(visit_day)

	// generate summary stats: average browsing time after 1st quiz
	gen time_spent_min_nytimes_sum = 0
	gen time_spent_ext_nytimes_sum = 0
	gen total_days = 0
	forvalues i = 20443/20918 {
		gen time_spent_ext_nytimes`i' = (time_spent_min_nytimes`i' > 0) & time_spent_min_nytimes`i' != . 
		replace time_spent_min_nytimes_sum = time_spent_min_nytimes_sum + time_spent_min_nytimes`i'
		replace time_spent_ext_nytimes_sum = time_spent_ext_nytimes_sum + time_spent_ext_nytimes`i'
		replace total_days = total_days + 1
		}
	gen time_spent_min_nytimes_avg = time_spent_min_nytimes_sum / total_days
	gen time_spent_ext_nytimes_avg = time_spent_ext_nytimes_sum / total_days
	
	// generate summary stats: percentage of weeks who have at least one visit to NYTimes
	forvalues i = 20443(7)20918 {
		gen time_week_min_nytimes`i' = 0
		gen time_week_ext_nytimes`i' = 0
		gen time_week_ett_nytimes`i' = 0
		forvalues j = 0/6 {
			local k = `i' + `j'
			replace time_week_min_nytimes`i' = time_week_min_nytimes`i' + time_spent_min_nytimes`k'
			replace time_week_ext_nytimes`i' = 1 	if time_spent_ext_nytimes`k' == 1
			replace time_week_ett_nytimes`i' = time_week_ett_nytimes`i' + time_spent_ext_nytimes`k'
			}
		}

	// generate indicator for active browser (based on ett)
	forvalues i = 20443(7)20918 {
		gen time_week_ett2_nytimes`i' = (time_week_ett_nytimes`i' >= 2) & time_week_ett_nytimes`i' != .
		}
		
	gen time_week_min_nytimes_sum = 0
	gen time_week_ext_nytimes_sum = 0
	gen time_week_ett_nytimes_sum = 0
	gen time_week_day_nytimes_sum = 0
	
	// generate total after 1st quiz
	gen total_weeks = 0
	forvalues i = 20485(7)20918 {
		replace time_week_min_nytimes_sum = time_week_min_nytimes_sum + time_week_min_nytimes`i'
		replace time_week_ext_nytimes_sum = time_week_ext_nytimes_sum + time_week_ext_nytimes`i'
		replace time_week_ett_nytimes_sum = time_week_ett_nytimes_sum + time_week_ett2_nytimes`i'
		replace time_week_day_nytimes_sum = time_week_day_nytimes_sum + time_week_ett_nytimes`i'
		replace total_weeks = total_weeks + 1
		}
	
	// deal with missing days
	foreach i in 20709 20793 20828 {
		replace time_week_min_nytimes`i' = time_week_min_nytimes`i' * 7 / 6
		}
		
	gen time_week_min_nytimes_avg = time_week_min_nytimes_sum / total_weeks
	gen time_week_ext_nytimes_avg = time_week_ext_nytimes_sum / total_weeks		
	gen time_week_ett_nytimes_avg = time_week_ett_nytimes_sum / total_weeks			
	gen time_week_day_nytimes_avg = time_week_day_nytimes_sum / total_weeks			
	
	// merge in adoption date
	merge 1:1 responseID_wave1 using "`directory_path'/vpn_date_adoption.dta"
	drop if _merge == 2
	drop _merge
	
	// merge in student information
	merge 1:1 responseID using "`directory_path'/vpn_treatment_roster"
	keep if _merge == 3
	drop _merge
	
	// reshape data
	keep if vpn_date_adoption != .
	keep treatment_newsletter time_week_min_nytimes*
	drop time_week_min_nytimes_sum time_week_min_nytimes_avg
	collapse time_week_min_nytimes*, by(treatment_newsletter)
	reshape long time_week_min_nytimes, i(treatment_newsletter) j(week)
	reshape wide time_week_min_nytimes, i(week) j(treatment_newsletter)

	// merge nytimes_wedge
	merge 1:1 week using "`directory_path'/Outregs/nytimes_newslog_wedge.dta"
	keep if _merge == 3
	drop _merge

	// intensive margin plot: full
	keep if week < 20912
	sum time_week_min_nytimes1
	local max = `r(max)'
	gen c1 = 12.5  if week >= 20443 & week <= 20485
	gen c2 = 12.5  if week >= 20485 & week <= 20527
	twoway 	(area c1 week, color(gs15) lcolor(gs15) fintensity(inten100)) ///
			(area c2 week, color(gs13) lcolor(gs13) fintensity(inten100)) ///
			(connected time_week_min_nytimes0 week, yaxis(1) lcolor(navy) msymbol(O) mcolor(navy)) ///
			(connected time_week_min_nytimes1 week, yaxis(1) lwidth(thick) lcolor(cranberry) msymbol(O) mcolor(cranberry)) ///
			(connected censored_prop_week week, yaxis(2) lpattern(dash) lwidth(medthick) lcolor(gs6) msymbol(th) mcolor(gs6)), ///
			xtitle("") ///
			ytitle("Total browsing time on NYTimes per week (min)", axis(1)) ///
			ytitle("% pol. sensitive articles on NYTimes", axis(2)) ///
			yscale(r(0 12) axis(1)) ///
			ylabel(0(4)12, axis(1)) ///
			yscale(r(0 0.30) axis(2)) ///
			ylabel(0 "0" 0.10 "10" 0.20 "20" 0.30 "30", axis(2)) ///
			xlabel(20454 "2016-01" 20485 "2016-02" 20514 "2016-03" 20545 "2016-04" 20575 "2016-05" 20606 "2016-06" 20636 "2016-07" 20667 "2016-08" 20698 "2016-09" 20728 "2016-10" 20759 "2016-11" 20789 "2016-12" 20820 "2017-01" 20851 "2017-02" 20879 "2017-03" 20910 "2017-04") ///
			legend(order(3 "Access" 4 "Access + Encour." 5 "% sensitive articles on NYTimes") col(3)) ///
			xsize(16) ///
			graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
			saving(figure_browsing_nyt_intensive_full, replace)
	graph export "figure_browsing_nyt_intensive_full.pdf", replace

	


*** Table 2 and Table A.8

	preserve
	clear
	gen i = ""
	gen mu_a_all = .
	gen sd_a_all = .
	gen mu_ae_all = .
	gen sd_ae_all = .
	gen pvalue_all = .
	gen mu_a_endline = .
	gen sd_a_endline = .
	gen mu_ae_endline = .
	gen sd_ae_endline = .
	gen pvalue_endline = .
	save "`directory_path'/Outregs/browsing_activities.dta", replace
	restore

* 	Panel A: adopter

	use "`directory_path'/browsingtime_totaltime.dta", clear
	
	// reshape data
	drop if visit_day == .
	keep responseID_wave1 visit_day time_spent_min_total_active2min
	rename time_spent_min_total_active2min 	time_spent_min_total_active
	reshape wide time_spent_min_total_active, i(responseID_wave1) j(visit_day)
	
	// merge in student information
	merge 1:1 responseID using "`directory_path'/vpn_treatment_roster"
	keep if treatment_vpn == 1
	drop _merge

	// fill in missing days
	foreach i in 20713 20793 20830 {
		replace time_spent_min_total_active`i' = .
		}
		
	// impute zeros
	forvalues i = 20439/20918 {
		replace time_spent_min_total_active`i' = 0 		if time_spent_min_total_active`i' == .
		}
		
	// merge in survey participation info
	merge 1:1 responseID using "`directory_path'/Outregs/panelmerged_indicator"
	keep if _merge == 3
	drop _merge
	
	// indicator for adopters
	forvalues i = 20439/20918 {
		gen time_spent_pos`i' = (time_spent_min_total_active`i' > 0)
		}
	gen time_spent_pos_sum = 0
	forvalues i = 20439/20918 {
		replace time_spent_pos_sum = time_spent_pos_sum + time_spent_pos`i'
		}
	gen adopter = (time_spent_pos_sum > 0)

	preserve
	gen i = "adopter"

	// test extensive margin of adoption: among all students
	ttest adopter, by(treatment_newsletter)
	gen mu_a_all = `r(mu_1)'
	gen sd_a_all = `r(sd_1)'
	gen mu_ae_all = `r(mu_2)'
	gen sd_ae_all = `r(sd_2)'
	gen pvalue_all = `r(p)'
	
	// test extensive margin of adoption: among endline participants
	ttest adopter if panelmerged_wave3 == 1, by(treatment_newsletter)
	gen mu_a_endline = `r(mu_1)'
	gen sd_a_endline = `r(sd_1)'
	gen mu_ae_endline = `r(mu_2)'
	gen sd_ae_endline = `r(sd_2)'
	gen pvalue_endline = `r(p)'
	
	keep i mu_a_all sd_a_all mu_ae_all sd_ae_all pvalue_all mu_a_endline sd_a_endline mu_ae_endline sd_ae_endline pvalue_endline
	duplicates drop
	append using "`directory_path'/Outregs/browsing_activities.dta"
	save "`directory_path'/Outregs/browsing_activities.dta", replace
	restore
	
	

* 	Panel A: active_user

	use "`directory_path'/browsingtime_totaltime.dta", clear
	
	// keep post 2nd quiz data (March 2016) only
	keep if visit_day > 20516

	// reshape data
	drop if visit_day == .
	keep responseID_wave1 visit_day time_spent_min_total_active2min
	rename time_spent_min_total_active2min 	time_spent_min_total_active
	reshape wide time_spent_min_total_active, i(responseID_wave1) j(visit_day)
	
	// merge in student information
	merge 1:1 responseID using "`directory_path'/vpn_treatment_roster"
	keep if treatment_vpn == 1
	drop _merge

	// fill in missing days
	foreach i in 20713 20793 20830 {
		replace time_spent_min_total_active`i' = .
		}
		
	// impute zeros
	forvalues i = 20517/20918 {
		replace time_spent_min_total_active`i' = 0 		if time_spent_min_total_active`i' == .
		}
		
	// merge in survey participation info
	merge 1:1 responseID using "`directory_path'/Outregs/panelmerged_indicator"
	keep if _merge == 3
	drop _merge
	
	// indicator for active users
	forvalues i = 20517/20918 {
		gen time_spent_pos`i' = (time_spent_min_total_active`i' > 0)
		}
	gen time_spent_pos_sum = 0
	forvalues i = 20517/20918 {
		replace time_spent_pos_sum = time_spent_pos_sum + time_spent_pos`i'
		}
	gen active_user = (time_spent_pos_sum > 40)

	// genearate user level averages
	gen time_spent_min_total_sum = 0
	forvalues i = 20517/20918 {
		replace time_spent_min_total_sum = time_spent_min_total_sum + time_spent_min_total_active`i'
		}
	gen time_spent_min_total_avg = time_spent_min_total_sum / 402

	
	preserve
	gen i = "activeuser"

	// test extensive margin of active_user: among all students
	ttest active_user, by(treatment_newsletter)
	gen mu_a_all = `r(mu_1)'
	gen sd_a_all = `r(sd_1)'
	gen mu_ae_all = `r(mu_2)'
	gen sd_ae_all = `r(sd_2)'
	gen pvalue_all = `r(p)'
	
	// test extensive margin of active_user: among endline participants
	ttest active_user if panelmerged_wave3 == 1, by(treatment_newsletter)
	gen mu_a_endline = `r(mu_1)'
	gen sd_a_endline = `r(sd_1)'
	gen mu_ae_endline = `r(mu_2)'
	gen sd_ae_endline = `r(sd_2)'
	gen pvalue_endline = `r(p)'
	
	keep i mu_a_all sd_a_all mu_ae_all sd_ae_all pvalue_all mu_a_endline sd_a_endline mu_ae_endline sd_ae_endline pvalue_endline
	duplicates drop
	append using "`directory_path'/Outregs/browsing_activities.dta"
	save "`directory_path'/Outregs/browsing_activities.dta", replace
	restore

	

* 	Panel B/C: total browsing time

	preserve
	gen i = "totaltime_active"
	
	// t-test by treatment: among active users
	ttest time_spent_min_total_avg if active_user == 1, by(treatment_newsletter)
	gen mu_a_all = `r(mu_1)'
	gen sd_a_all = `r(sd_1)'
	gen mu_ae_all = `r(mu_2)'
	gen sd_ae_all = `r(sd_2)'
	gen pvalue_all = `r(p)'

	// t-test by treatment: among active users & endline participants
	ttest time_spent_min_total_avg if active_user == 1 & panelmerged_wave3 == 1, by(treatment_newsletter)
	gen mu_a_endline = `r(mu_1)'
	gen sd_a_endline = `r(sd_1)'
	gen mu_ae_endline = `r(mu_2)'
	gen sd_ae_endline = `r(sd_2)'
	gen pvalue_endline = `r(p)'
	
	keep i mu_a_all sd_a_all mu_ae_all sd_ae_all pvalue_all mu_a_endline sd_a_endline mu_ae_endline sd_ae_endline pvalue_endline
	duplicates drop
	append using "`directory_path'/Outregs/browsing_activities.dta"
	save "`directory_path'/Outregs/browsing_activities.dta", replace
	restore

	
	preserve
	gen i = "totaltime_all"
	
	// t-test by treatment: among active users
	ttest time_spent_min_total_avg, by(treatment_newsletter)
	gen mu_a_all = `r(mu_1)'
	gen sd_a_all = `r(sd_1)'
	gen mu_ae_all = `r(mu_2)'
	gen sd_ae_all = `r(sd_2)'
	gen pvalue_all = `r(p)'

	// t-test by treatment: among active users & endline participants
	ttest time_spent_min_total_avg if panelmerged_wave3 == 1, by(treatment_newsletter)
	gen mu_a_endline = `r(mu_1)'
	gen sd_a_endline = `r(sd_1)'
	gen mu_ae_endline = `r(mu_2)'
	gen sd_ae_endline = `r(sd_2)'
	gen pvalue_endline = `r(p)'
	
	keep i mu_a_all sd_a_all mu_ae_all sd_ae_all pvalue_all mu_a_endline sd_a_endline mu_ae_endline sd_ae_endline pvalue_endline
	duplicates drop
	append using "`directory_path'/Outregs/browsing_activities.dta"
	save "`directory_path'/Outregs/browsing_activities.dta", replace
	restore
	


* 	Panel B/C: all categories except for top foreign news

	foreach w in google youtube facebook twitter nytimes information wikipedia entertainment porn {

		use "`directory_path'/browsingtime_`w'.dta", clear
		
		// keep post 2nd quiz data (March 2016) only
		keep if visit_day > 20516

		// reshape data (standardized variable name)
		drop if visit_day == .
		rename time_spent_min_`w' 	time_spent_min_total_active
		keep responseID_wave1 visit_day time_spent_min_total_active
		reshape wide time_spent_min_total_active, i(responseID_wave1) j(visit_day)

		// merge in student information
		merge 1:1 responseID using "`directory_path'/vpn_treatment_roster"
		keep if treatment_vpn == 1
		drop _merge

		// fill in missing days
		foreach i in 20713 20793 20830 {
			replace time_spent_min_total_active`i' = .
			}
			
		// impute zeros
		forvalues i = 20517/20918 {
			replace time_spent_min_total_active`i' = 0 		if time_spent_min_total_active`i' == .
			}
			
		// merge in survey participation info
		merge 1:1 responseID using "`directory_path'/Outregs/panelmerged_indicator"
		keep if _merge == 3
		drop _merge

		// merge indicator for active user
		merge 1:1 responseID using "`directory_path'/vpn_browsing_active_user.dta"
		drop _merge
		
		// genearate user level averages
		gen time_spent_min_total_sum = 0
		forvalues i = 20517/20918 {
			replace time_spent_min_total_sum = time_spent_min_total_sum + time_spent_min_total_active`i'
			}
		gen time_spent_min_total_avg = time_spent_min_total_sum / 402
		
		preserve
		gen i = "`w'_active"
		
		// t-test by treatment: among active users
		ttest time_spent_min_total_avg if active_user == 1, by(treatment_newsletter)
		gen mu_a_all = `r(mu_1)'
		gen sd_a_all = `r(sd_1)'
		gen mu_ae_all = `r(mu_2)'
		gen sd_ae_all = `r(sd_2)'
		gen pvalue_all = `r(p)'

		// t-test by treatment: among active users & endline participants
		ttest time_spent_min_total_avg if active_user == 1 & panelmerged_wave3 == 1, by(treatment_newsletter)
		gen mu_a_endline = `r(mu_1)'
		gen sd_a_endline = `r(sd_1)'
		gen mu_ae_endline = `r(mu_2)'
		gen sd_ae_endline = `r(sd_2)'
		gen pvalue_endline = `r(p)'
		
		keep i mu_a_all sd_a_all mu_ae_all sd_ae_all pvalue_all mu_a_endline sd_a_endline mu_ae_endline sd_ae_endline pvalue_endline
		duplicates drop
		append using "`directory_path'/Outregs/browsing_activities.dta"
		save "`directory_path'/Outregs/browsing_activities.dta", replace
		restore

		
		preserve
		gen i = "`w'_all"
		
		// t-test by treatment: among active users
		ttest time_spent_min_total_avg, by(treatment_newsletter)
		gen mu_a_all = `r(mu_1)'
		gen sd_a_all = `r(sd_1)'
		gen mu_ae_all = `r(mu_2)'
		gen sd_ae_all = `r(sd_2)'
		gen pvalue_all = `r(p)'

		// t-test by treatment: among active users & endline participants
		ttest time_spent_min_total_avg if panelmerged_wave3 == 1, by(treatment_newsletter)
		gen mu_a_endline = `r(mu_1)'
		gen sd_a_endline = `r(sd_1)'
		gen mu_ae_endline = `r(mu_2)'
		gen sd_ae_endline = `r(sd_2)'
		gen pvalue_endline = `r(p)'
		
		keep i mu_a_all sd_a_all mu_ae_all sd_ae_all pvalue_all mu_a_endline sd_a_endline mu_ae_endline sd_ae_endline pvalue_endline
		duplicates drop
		append using "`directory_path'/Outregs/browsing_activities.dta"
		save "`directory_path'/Outregs/browsing_activities.dta", replace
		restore
	
	}
	

	
* 	Panel B/C: top foreign news

	// load NYTimes data
	use "`directory_path'/browsingtime_nytimes.dta", clear
	
	// merge with other top foreign news sites
	foreach var in economist reddit cnn theguardian huffingtonpost foxnews bbc bloomberg wsj usatoday reuters nbcnews ft {
		merge 1:1 responseID_wave1 visit_day using "`directory_path'/browsingtime_`var'.dta"
		drop _merge
		}
	
	// keep post 2nd quiz data (March 2016) only
	keep if visit_day > 20516

	// generate total browsing time
	gen time_spent_min_foreignnews = time_spent_min_nytimes
	foreach var in economist reddit cnn theguardian huffingtonpost foxnews bbc bloomberg wsj usatoday reuters nbcnews ft {
		replace time_spent_min_foreignnews = time_spent_min_foreignnews + time_spent_min_`var'
		}

	// reshape data (standardized variable name)
	drop if visit_day == .
	rename time_spent_min_foreignnews 	time_spent_min_total_active
	keep responseID_wave1 visit_day time_spent_min_total_active
	reshape wide time_spent_min_total_active, i(responseID_wave1) j(visit_day)

	// merge in student information
	merge 1:1 responseID using "`directory_path'/vpn_treatment_roster"
	keep if treatment_vpn == 1
	drop _merge

	// fill in missing days
	foreach i in 20713 20793 20830 {
		replace time_spent_min_total_active`i' = .
		}
		
	// impute zeros
	forvalues i = 20517/20918 {
		replace time_spent_min_total_active`i' = 0 		if time_spent_min_total_active`i' == .
		}
		
	// merge in survey participation info
	merge 1:1 responseID using "`directory_path'/Outregs/panelmerged_indicator"
	keep if _merge == 3
	drop _merge

	// merge indicator for active user
	merge 1:1 responseID using "`directory_path'/vpn_browsing_active_user.dta"
	drop _merge
	
	// genearate user level averages
	gen time_spent_min_total_sum = 0
	forvalues i = 20517/20918 {
		replace time_spent_min_total_sum = time_spent_min_total_sum + time_spent_min_total_active`i'
		}
	gen time_spent_min_total_avg = time_spent_min_total_sum / 402
	
	preserve
	gen i = "foreignnews_active"
	
	// t-test by treatment: among active users
	ttest time_spent_min_total_avg if active_user == 1, by(treatment_newsletter)
	gen mu_a_all = `r(mu_1)'
	gen sd_a_all = `r(sd_1)'
	gen mu_ae_all = `r(mu_2)'
	gen sd_ae_all = `r(sd_2)'
	gen pvalue_all = `r(p)'

	// t-test by treatment: among active users & endline participants
	ttest time_spent_min_total_avg if active_user == 1 & panelmerged_wave3 == 1, by(treatment_newsletter)
	gen mu_a_endline = `r(mu_1)'
	gen sd_a_endline = `r(sd_1)'
	gen mu_ae_endline = `r(mu_2)'
	gen sd_ae_endline = `r(sd_2)'
	gen pvalue_endline = `r(p)'
	
	keep i mu_a_all sd_a_all mu_ae_all sd_ae_all pvalue_all mu_a_endline sd_a_endline mu_ae_endline sd_ae_endline pvalue_endline
	duplicates drop
	append using "`directory_path'/Outregs/browsing_activities.dta"
	save "`directory_path'/Outregs/browsing_activities.dta", replace
	restore

	
	preserve
	gen i = "foreignnews_all"
	
	// t-test by treatment: among active users
	ttest time_spent_min_total_avg, by(treatment_newsletter)
	gen mu_a_all = `r(mu_1)'
	gen sd_a_all = `r(sd_1)'
	gen mu_ae_all = `r(mu_2)'
	gen sd_ae_all = `r(sd_2)'
	gen pvalue_all = `r(p)'

	// t-test by treatment: among active users & endline participants
	ttest time_spent_min_total_avg if panelmerged_wave3 == 1, by(treatment_newsletter)
	gen mu_a_endline = `r(mu_1)'
	gen sd_a_endline = `r(sd_1)'
	gen mu_ae_endline = `r(mu_2)'
	gen sd_ae_endline = `r(sd_2)'
	gen pvalue_endline = `r(p)'
	
	keep i mu_a_all sd_a_all mu_ae_all sd_ae_all pvalue_all mu_a_endline sd_a_endline mu_ae_endline sd_ae_endline pvalue_endline
	duplicates drop
	append using "`directory_path'/Outregs/browsing_activities.dta"
	save "`directory_path'/Outregs/browsing_activities.dta", replace
	restore

	

	
*** Figure A.4

	// load raw data
	use "`directory_path'/browsingtime_nytimes.dta", clear
	
	// reshape data (standardized variable name)
	drop if visit_day == .
	keep responseID_wave1 visit_day time_spent_min_nytimes
	reshape wide time_spent_min_nytimes, i(responseID_wave1) j(visit_day)

	// generate summary stats: average browsing time after 1st quiz
	gen time_spent_min_nytimes_sum = 0
	gen time_spent_ext_nytimes_sum = 0
	gen total_days = 0
	forvalues i = 20443/20918 {
		gen time_spent_ext_nytimes`i' = (time_spent_min_nytimes`i' > 0) & time_spent_min_nytimes`i' != . 
		replace time_spent_min_nytimes_sum = time_spent_min_nytimes_sum + time_spent_min_nytimes`i'
		replace time_spent_ext_nytimes_sum = time_spent_ext_nytimes_sum + time_spent_ext_nytimes`i'
		replace total_days = total_days + 1
		}
	gen time_spent_min_nytimes_avg = time_spent_min_nytimes_sum / total_days
	gen time_spent_ext_nytimes_avg = time_spent_ext_nytimes_sum / total_days
	
	// generate summary stats: percentage of weeks who have at least one visit to NYTimes
	forvalues i = 20443(7)20918 {
		gen time_week_min_nytimes`i' = 0
		gen time_week_ext_nytimes`i' = 0
		gen time_week_ett_nytimes`i' = 0
		forvalues j = 0/6 {
			local k = `i' + `j'
			replace time_week_min_nytimes`i' = time_week_min_nytimes`i' + time_spent_min_nytimes`k'
			replace time_week_ext_nytimes`i' = 1 	if time_spent_ext_nytimes`k' == 1
			replace time_week_ett_nytimes`i' = time_week_ett_nytimes`i' + time_spent_ext_nytimes`k'
			}
		}

	// generate indicator for active browser (based on ett)
	forvalues i = 20443(7)20918 {
		gen time_week_ett2_nytimes`i' = (time_week_ett_nytimes`i' >= 2) & time_week_ett_nytimes`i' != .
		}
		
	gen time_week_min_nytimes_sum = 0
	gen time_week_ext_nytimes_sum = 0
	gen time_week_ett_nytimes_sum = 0
	gen time_week_day_nytimes_sum = 0
	
	// generate total after 1st quiz
	gen total_weeks = 0
	forvalues i = 20485(7)20918 {
		replace time_week_min_nytimes_sum = time_week_min_nytimes_sum + time_week_min_nytimes`i'
		replace time_week_ext_nytimes_sum = time_week_ext_nytimes_sum + time_week_ext_nytimes`i'
		replace time_week_ett_nytimes_sum = time_week_ett_nytimes_sum + time_week_ett2_nytimes`i'
		replace time_week_day_nytimes_sum = time_week_day_nytimes_sum + time_week_ett_nytimes`i'
		replace total_weeks = total_weeks + 1
		}
	
	// deal with missing days
	foreach i in 20709 20793 20828 {
		replace time_week_min_nytimes`i' = time_week_min_nytimes`i' * 7 / 6
		}
		
	gen time_week_min_nytimes_avg = time_week_min_nytimes_sum / total_weeks
	gen time_week_ext_nytimes_avg = time_week_ext_nytimes_sum / total_weeks		
	gen time_week_ett_nytimes_avg = time_week_ett_nytimes_sum / total_weeks			
	gen time_week_day_nytimes_avg = time_week_day_nytimes_sum / total_weeks			
	
	// merge in adoption date
	merge 1:1 responseID_wave1 using "`directory_path'/vpn_date_adoption.dta"
	drop if _merge == 2
	drop _merge
	
	// merge in student information
	merge 1:1 responseID using "`directory_path'/vpn_treatment_roster"
	keep if treatment_vpn == 1
	drop _merge
	
	// generate indicator for ever adopted
	gen adopted = (vpn_date_adoption != .)

	// adoption count by day: newsletter == 1
	preserve
	keep if treatment_newsletter == 1
	
	gen c = 1
	egen count_nl_yes = sum(c)
	bysort vpn_date_adoption: egen adoption_nl_yes = sum(adopted)
	keep vpn_date_adoption adoption_nl_yes count_nl_yes
	duplicates drop
	drop if vpn_date_adoption == .
	save "`directory_path'/Outregs/adoption_count_nl_yes.dta", replace
	restore
	
	// adoption count by day: newsletter == 0
	preserve
	keep if treatment_newsletter == 0
	
	gen c = 1
	egen count_nl_no = sum(c)
	bysort vpn_date_adoption: egen adoption_nl_no = sum(adopted)
	keep vpn_date_adoption adoption_nl_no count_nl_no
	duplicates drop
	drop if vpn_date_adoption == .
	save "`directory_path'/Outregs/adoption_count_nl_no.dta", replace
	restore
		
	// merge in data
	preserve
	use "`directory_path'/vpn_log_dates.dta", clear
	drop if dates == ""
	keep if dates_code <= 20603
	keep dates dates_code
	// resolve 2016-02-19
	replace dates = "2016-02-19" 	if dates_code == 20503
	duplicates drop
	
	rename dates_code vpn_date_adoption
	
	merge 1:1 vpn_date_adoption using "`directory_path'/Outregs/adoption_count_nl_yes.dta"
	drop _merge
	merge 1:1 vpn_date_adoption using "`directory_path'/Outregs/adoption_count_nl_no.dta"
	drop _merge
	erase "`directory_path'/Outregs/adoption_count_nl_yes.dta"
	erase "`directory_path'/Outregs/adoption_count_nl_no.dta"
	
	replace adoption_nl_yes = 0 	if adoption_nl_yes == .
	replace adoption_nl_no = 0 		if adoption_nl_no == .
	
	// switch to percentage
	egen count_nl_yes1 = mean(count_nl_yes) 
	drop count_nl_yes
	rename count_nl_yes1 count_nl_yes
	egen count_nl_no1 = mean(count_nl_no) 
	drop count_nl_no
	rename count_nl_no1 count_nl_no	
	
	// cut data by weeks
	egen vpn_date_interval = cut(vpn_date_adoption), at(20438(5)20603)
	bysort vpn_date_interval: egen adoption_total_nl_yes = sum(adoption_nl_yes)
	bysort vpn_date_interval: egen adoption_total_nl_no = sum(adoption_nl_no)
	
	gen adoption_pc_nl_yes = adoption_total_nl_yes / count_nl_yes
	gen adoption_pc_nl_no = adoption_total_nl_no / count_nl_no
	
	keep vpn_date_interval adoption_pc_nl_yes adoption_pc_nl_no 
	duplicates drop
	
	// generate cumulative adoption
	sort vpn_date_interval
	gen adoption_cum_nl_yes = adoption_pc_nl_yes[1]
	replace adoption_cum_nl_yes = adoption_pc_nl_yes[_n] + adoption_cum_nl_yes[_n-1] if _n > 1
	gen adoption_cum_nl_no = adoption_pc_nl_no[1]
	replace adoption_cum_nl_no = adoption_pc_nl_no[_n] + adoption_cum_nl_no[_n-1] if _n > 1

	qui su adoption_cum_nl_yes
	local max = `r(max)'
	gen c1 = `max' + 0.03  if vpn_date_interval >= 20438 & vpn_date_interval <= 20488
	gen c2 = `max' + 0.03  if vpn_date_interval >= 20488 & vpn_date_interval <= 20537

	twoway 	(area c1 vpn_date_interval, color(gs15) lcolor(gs15) fintensity(inten100)) ///
			(area c2 vpn_date_interval, color(gs13) lcolor(gs13) fintensity(inten100)) ///
			(connected adoption_cum_nl_no vpn_date_interval, lcolor(navy) msymbol(O) mcolor(navy)) ///
			(connected adoption_cum_nl_yes vpn_date_interval, lwidth(thick) lcolor(cranberry) msymbol(O) mcolor(cranberry)), ///
			text(0.69 20463 "Unincentivized" "encouragement" "period", size(small)) ///
			text(0.69 20510 "Incentivized" "encouragement" "period", size(small)) ///
		    ytitle("Cumulative activation rate") ///
			xtitle("") ///
			xlabel(20454 "2016-01" 20485 "2016-02" 20514 "2016-03" 20545 "2016-04" 20575 "2016-05") ///
			legend(on order(3 4) label(3 "Access") label(4 "Access + Encour.")) ///
			graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
			saving(figure_adoption_cumulative_all, replace)
	
	graph export "figure_adoption_cumulative_all.pdf", replace
	
	restore	
	
	

	
*** Figure A.5

	// load raw data
	use "`directory_path'/browsingtime_nytimes.dta", clear
	
	// reshape data (standardized variable name)
	drop if visit_day == .
	keep responseID_wave1 visit_day time_spent_min_nytimes
	reshape wide time_spent_min_nytimes, i(responseID_wave1) j(visit_day)

	// generate summary stats: average browsing time after 1st quiz
	gen time_spent_min_nytimes_sum = 0
	gen time_spent_ext_nytimes_sum = 0
	gen total_days = 0
	forvalues i = 20443/20918 {
		gen time_spent_ext_nytimes`i' = (time_spent_min_nytimes`i' > 0) & time_spent_min_nytimes`i' != . 
		replace time_spent_min_nytimes_sum = time_spent_min_nytimes_sum + time_spent_min_nytimes`i'
		replace time_spent_ext_nytimes_sum = time_spent_ext_nytimes_sum + time_spent_ext_nytimes`i'
		replace total_days = total_days + 1
		}
	gen time_spent_min_nytimes_avg = time_spent_min_nytimes_sum / total_days
	gen time_spent_ext_nytimes_avg = time_spent_ext_nytimes_sum / total_days
	
	// generate summary stats: percentage of weeks who have at least one visit to NYTimes
	forvalues i = 20443(7)20918 {
		gen time_week_min_nytimes`i' = 0
		gen time_week_ext_nytimes`i' = 0
		gen time_week_ett_nytimes`i' = 0
		forvalues j = 0/6 {
			local k = `i' + `j'
			replace time_week_min_nytimes`i' = time_week_min_nytimes`i' + time_spent_min_nytimes`k'
			replace time_week_ext_nytimes`i' = 1 	if time_spent_ext_nytimes`k' == 1
			replace time_week_ett_nytimes`i' = time_week_ett_nytimes`i' + time_spent_ext_nytimes`k'
			}
		}

	// generate indicator for active browser (based on ett)
	forvalues i = 20443(7)20918 {
		gen time_week_ett2_nytimes`i' = (time_week_ett_nytimes`i' >= 2) & time_week_ett_nytimes`i' != .
		}
		
	gen time_week_min_nytimes_sum = 0
	gen time_week_ext_nytimes_sum = 0
	gen time_week_ett_nytimes_sum = 0
	gen time_week_day_nytimes_sum = 0
	
	// generate total after 1st quiz
	gen total_weeks = 0
	forvalues i = 20485(7)20918 {
		replace time_week_min_nytimes_sum = time_week_min_nytimes_sum + time_week_min_nytimes`i'
		replace time_week_ext_nytimes_sum = time_week_ext_nytimes_sum + time_week_ext_nytimes`i'
		replace time_week_ett_nytimes_sum = time_week_ett_nytimes_sum + time_week_ett2_nytimes`i'
		replace time_week_day_nytimes_sum = time_week_day_nytimes_sum + time_week_ett_nytimes`i'
		replace total_weeks = total_weeks + 1
		}
	
	// deal with missing days
	foreach i in 20709 20793 20828 {
		replace time_week_min_nytimes`i' = time_week_min_nytimes`i' * 7 / 6
		}
		
	gen time_week_min_nytimes_avg = time_week_min_nytimes_sum / total_weeks
	gen time_week_ext_nytimes_avg = time_week_ext_nytimes_sum / total_weeks		
	gen time_week_ett_nytimes_avg = time_week_ett_nytimes_sum / total_weeks			
	gen time_week_day_nytimes_avg = time_week_day_nytimes_sum / total_weeks			
	
	// merge in adoption date
	merge 1:1 responseID_wave1 using "`directory_path'/vpn_date_adoption.dta"
	drop if _merge == 2
	drop _merge
	
	// merge in student information
	merge 1:1 responseID using "`directory_path'/vpn_treatment_roster"
	keep if treatment_vpn == 1
	drop _merge
	
	replace time_week_day_nytimes_avg = 0 if vpn_date_adoption == .
	replace time_week_min_nytimes_avg = 0 if vpn_date_adoption == .

	cumul time_week_day_nytimes_avg if treatment_newsletter == 1, gen(time_week_day_nytimes_avg_cum1)
	cumul time_week_day_nytimes_avg if treatment_newsletter == 0, gen(time_week_day_nytimes_avg_cum0)
	sort time_week_day_nytimes_avg_cum1 time_week_day_nytimes_avg_cum0
	twoway 	(line time_week_day_nytimes_avg_cum1 time_week_day_nytimes_avg, mcolor(cranberry) lpattern(solid) lwidth(thick) lcolor(cranberry)) /// 
			(line time_week_day_nytimes_avg_cum0 time_week_day_nytimes_avg, mcolor(navy) lcolor(navy)), ///
			yline(0.454, lpattern(dash) lcolor(navy)) ///
			yline(0.305, lpattern(dash) lwidth(thick) lcolor(cranberry)) ///
			ylabel(, grid) ytitle("") ///
			xtitle("Average # days visit on NYT each week") ///
			legend(on order(2 1) label(1 "Access + Encour.") label(2 "Access")) ///
			graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
			saving(figure_cdf_time_week_day_nytimes_avg_all, replace)			
	graph export "figure_cdf_time_week_day_nytimes_avg_all.pdf", replace
					
	cumul time_week_min_nytimes_avg if treatment_newsletter == 1, gen(time_week_min_nytimes_avg_cum1)
	cumul time_week_min_nytimes_avg if treatment_newsletter == 0, gen(time_week_min_nytimes_avg_cum0)
	sort time_week_min_nytimes_avg_cum1 time_week_min_nytimes_avg_cum0
	twoway 	(line time_week_min_nytimes_avg_cum1 time_week_min_nytimes_avg, mcolor(cranberry) lpattern(solid) lwidth(thick) lcolor(cranberry)) /// 
			(line time_week_min_nytimes_avg_cum0 time_week_min_nytimes_avg, mcolor(navy) lcolor(navy)), ///
			yline(0.454, lpattern(dash) lcolor(navy)) ///
			yline(0.305, lpattern(dash) lwidth(thick) lcolor(cranberry)) ///
			ylabel(, grid) ytitle("") ///
			xtitle("Average mins spent on NYT each week") ///
			legend(on order(2 1) label(1 "Access + Encour.") label(2 "Access")) ///
			graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
			saving(figure_cdf_time_week_min_nytimes_avg_all, replace)			
	graph export "figure_cdf_time_week_min_nytimes_avg_all.pdf", replace
	
	// combine graphs
	grc1leg figure_cdf_time_week_day_nytimes_avg_all.gph figure_cdf_time_week_min_nytimes_avg_all.gph, ///
		cols(2) ///
		scale(1.1) ///
		ysize(3.5) xsize(15) ///
		legendfrom(figure_cdf_time_week_day_nytimes_avg_all.gph) ///
		graphregion(fcolor(white) ilcolor(white) lcolor(white))
	graph export "figure_cdf_time_week_nytimes_all.pdf", replace




*** Figure A.6

	// load data: NYTimes
	use "`directory_path'/browsingtime_nytimes.dta", clear

	// merge with other foreign news sites
	foreach var in economist reddit cnn theguardian huffingtonpost foxnews bbc bloomberg wsj usatoday reuters nbcnews ft {
		merge 1:1 responseID_wave1 visit_day using "`directory_path'/browsingtime_`var'.dta"
		drop _merge
		}

	// generate total browsing time
	gen time_spent_min_foreignnews = time_spent_min_nytimes
	foreach var in economist reddit cnn theguardian huffingtonpost foxnews bbc bloomberg wsj usatoday reuters nbcnews ft {
		replace time_spent_min_foreignnews = time_spent_min_foreignnews + time_spent_min_`var'
		}	

	// reshape data (standardized variable name)
	drop if visit_day == .
	keep responseID_wave1 visit_day time_spent_min_foreignnews
	reshape wide time_spent_min_foreignnews, i(responseID_wave1) j(visit_day)

	// generate summary stats: average browsing time after 1st quiz
	gen time_spent_min_foreignnews_sum = 0
	gen time_spent_ext_foreignnews_sum = 0
	gen total_days = 0
	forvalues i = 20443/20918 {
		gen time_spent_ext_foreignnews`i' = (time_spent_min_foreignnews`i' > 0) & time_spent_min_foreignnews`i' != . 
		replace time_spent_min_foreignnews_sum = time_spent_min_foreignnews_sum + time_spent_min_foreignnews`i'
		replace time_spent_ext_foreignnews_sum = time_spent_ext_foreignnews_sum + time_spent_ext_foreignnews`i'
		replace total_days = total_days + 1
		}
	gen time_spent_min_foreignnews_avg = time_spent_min_foreignnews_sum / total_days
	gen time_spent_ext_foreignnews_avg = time_spent_ext_foreignnews_sum / total_days
	
	// generate summary stats: percentage of weeks who have at least one visit to foreignnews
	forvalues i = 20443(7)20918 {
		gen time_week_min_foreignnews`i' = 0
		gen time_week_ext_foreignnews`i' = 0
		gen time_week_ett_foreignnews`i' = 0
		forvalues j = 0/6 {
			local k = `i' + `j'
			replace time_week_min_foreignnews`i' = time_week_min_foreignnews`i' + time_spent_min_foreignnews`k'
			replace time_week_ext_foreignnews`i' = 1 	if time_spent_ext_foreignnews`k' == 1
			replace time_week_ett_foreignnews`i' = time_week_ett_foreignnews`i' + time_spent_ext_foreignnews`k'
			}
		}

	// generate indicator for active browser (based on ett)
	forvalues i = 20443(7)20918 {
		gen time_week_ett2_foreignnews`i' = (time_week_ett_foreignnews`i' >= 2) & time_week_ett_foreignnews`i' != .
		}
		
	gen time_week_min_foreignnews_sum = 0
	gen time_week_ext_foreignnews_sum = 0
	gen time_week_ett_foreignnews_sum = 0
	gen time_week_day_foreignnews_sum = 0
	
	// generate total after 1st quiz
	gen total_weeks = 0
	forvalues i = 20485(7)20918 {
		replace time_week_min_foreignnews_sum = time_week_min_foreignnews_sum + time_week_min_foreignnews`i'
		replace time_week_ext_foreignnews_sum = time_week_ext_foreignnews_sum + time_week_ext_foreignnews`i'
		replace time_week_ett_foreignnews_sum = time_week_ett_foreignnews_sum + time_week_ett2_foreignnews`i'
		replace time_week_day_foreignnews_sum = time_week_day_foreignnews_sum + time_week_ett_foreignnews`i'
		replace total_weeks = total_weeks + 1
		}
	
	// deal with missing days
	foreach i in 20709 20793 20828 {
		replace time_week_min_foreignnews`i' = time_week_min_foreignnews`i' * 7 / 6
		}
		
	gen time_week_min_foreignnews_avg = time_week_min_foreignnews_sum / total_weeks
	gen time_week_ext_foreignnews_avg = time_week_ext_foreignnews_sum / total_weeks		
	gen time_week_ett_foreignnews_avg = time_week_ett_foreignnews_sum / total_weeks			
	gen time_week_day_foreignnews_avg = time_week_day_foreignnews_sum / total_weeks			
	
	// merge in adoption date
	merge 1:1 responseID_wave1 using "`directory_path'/vpn_date_adoption.dta"
	drop if _merge == 2
	drop _merge
	
	// merge in student information
	merge 1:1 responseID using "`directory_path'/vpn_treatment_roster"
	keep if treatment_vpn == 1
	drop _merge
	

* 	extract and merge data
	preserve
	keep if vpn_date_adoption != .
	keep treatment_newsletter time_week_ett2_foreignnews*
	collapse time_week_ett2_foreignnews*, by(treatment_newsletter)
	reshape long time_week_ett2_foreignnews, i(treatment_newsletter) j(week)
	reshape wide time_week_ett2_foreignnews, i(week) j(treatment_newsletter)
	save "`directory_path'/Outregs/vpn_log_master_dynamic_foreignnews.dta", replace
	restore

	preserve
	keep if vpn_date_adoption != .
	keep treatment_newsletter time_week_min_foreignnews*
	drop time_week_min_foreignnews_sum time_week_min_foreignnews_avg
	collapse time_week_min_foreignnews*, by(treatment_newsletter)
	reshape long time_week_min_foreignnews, i(treatment_newsletter) j(week)
	reshape wide time_week_min_foreignnews, i(week) j(treatment_newsletter)
	
	merge 1:1 week using "`directory_path'/Outregs/vpn_log_master_dynamic_foreignnews.dta"
	keep if _merge == 3
	drop _merge
	save "`directory_path'/Outregs/vpn_log_master_dynamic_foreignnews.dta", replace
	restore

	
* 	merge in nytimes data

	// load raw data
	use "`directory_path'/browsingtime_nytimes.dta", clear
	
	// reshape data (standardized variable name)
	drop if visit_day == .
	keep responseID_wave1 visit_day time_spent_min_nytimes
	reshape wide time_spent_min_nytimes, i(responseID_wave1) j(visit_day)

	// generate summary stats: average browsing time after 1st quiz
	gen time_spent_min_nytimes_sum = 0
	gen time_spent_ext_nytimes_sum = 0
	gen total_days = 0
	forvalues i = 20443/20918 {
		gen time_spent_ext_nytimes`i' = (time_spent_min_nytimes`i' > 0) & time_spent_min_nytimes`i' != . 
		replace time_spent_min_nytimes_sum = time_spent_min_nytimes_sum + time_spent_min_nytimes`i'
		replace time_spent_ext_nytimes_sum = time_spent_ext_nytimes_sum + time_spent_ext_nytimes`i'
		replace total_days = total_days + 1
		}
	gen time_spent_min_nytimes_avg = time_spent_min_nytimes_sum / total_days
	gen time_spent_ext_nytimes_avg = time_spent_ext_nytimes_sum / total_days
	
	// generate summary stats: percentage of weeks who have at least one visit to NYTimes
	forvalues i = 20443(7)20918 {
		gen time_week_min_nytimes`i' = 0
		gen time_week_ext_nytimes`i' = 0
		gen time_week_ett_nytimes`i' = 0
		forvalues j = 0/6 {
			local k = `i' + `j'
			replace time_week_min_nytimes`i' = time_week_min_nytimes`i' + time_spent_min_nytimes`k'
			replace time_week_ext_nytimes`i' = 1 	if time_spent_ext_nytimes`k' == 1
			replace time_week_ett_nytimes`i' = time_week_ett_nytimes`i' + time_spent_ext_nytimes`k'
			}
		}

	// generate indicator for active browser (based on ett)
	forvalues i = 20443(7)20918 {
		gen time_week_ett2_nytimes`i' = (time_week_ett_nytimes`i' >= 2) & time_week_ett_nytimes`i' != .
		}
		
	gen time_week_min_nytimes_sum = 0
	gen time_week_ext_nytimes_sum = 0
	gen time_week_ett_nytimes_sum = 0
	gen time_week_day_nytimes_sum = 0
	
	// generate total after 1st quiz
	gen total_weeks = 0
	forvalues i = 20485(7)20918 {
		replace time_week_min_nytimes_sum = time_week_min_nytimes_sum + time_week_min_nytimes`i'
		replace time_week_ext_nytimes_sum = time_week_ext_nytimes_sum + time_week_ext_nytimes`i'
		replace time_week_ett_nytimes_sum = time_week_ett_nytimes_sum + time_week_ett2_nytimes`i'
		replace time_week_day_nytimes_sum = time_week_day_nytimes_sum + time_week_ett_nytimes`i'
		replace total_weeks = total_weeks + 1
		}
	
	// deal with missing days
	foreach i in 20709 20793 20828 {
		replace time_week_min_nytimes`i' = time_week_min_nytimes`i' * 7 / 6
		}
		
	gen time_week_min_nytimes_avg = time_week_min_nytimes_sum / total_weeks
	gen time_week_ext_nytimes_avg = time_week_ext_nytimes_sum / total_weeks		
	gen time_week_ett_nytimes_avg = time_week_ett_nytimes_sum / total_weeks			
	gen time_week_day_nytimes_avg = time_week_day_nytimes_sum / total_weeks			
	
	// merge in adoption date
	merge 1:1 responseID_wave1 using "`directory_path'/vpn_date_adoption.dta"
	drop if _merge == 2
	drop _merge
	
	// merge in student information
	merge 1:1 responseID using "`directory_path'/vpn_treatment_roster"
	keep if _merge == 3
	drop _merge
	
	// reshape data
	preserve
	keep if vpn_date_adoption != .
	keep treatment_newsletter time_week_ett2_nytimes*
	collapse time_week_ett2_nytimes*, by(treatment_newsletter)
	reshape long time_week_ett2_nytimes, i(treatment_newsletter) j(week)
	reshape wide time_week_ett2_nytimes, i(week) j(treatment_newsletter)
	save "`directory_path'/Outregs/vpn_log_master_dynamic_nytimes.dta", replace
	restore

	preserve
	keep if vpn_date_adoption != .
	keep treatment_newsletter time_week_min_nytimes*
	drop time_week_min_nytimes_sum time_week_min_nytimes_avg
	collapse time_week_min_nytimes*, by(treatment_newsletter)
	reshape long time_week_min_nytimes, i(treatment_newsletter) j(week)
	reshape wide time_week_min_nytimes, i(week) j(treatment_newsletter)
	
	merge 1:1 week using "`directory_path'/Outregs/vpn_log_master_dynamic_nytimes.dta"
	keep if _merge == 3
	drop _merge
	save "`directory_path'/Outregs/vpn_log_master_dynamic_nytimes.dta", replace
	restore

	
* 	merge
	use "`directory_path'/Outregs/vpn_log_master_dynamic_foreignnews.dta", clear
	merge 1:1 week using "`directory_path'/Outregs/vpn_log_master_dynamic_nytimes.dta"
	keep if _merge == 3
	drop _merge
	save "`directory_path'/Outregs/vpn_log_master_dynamic_foreignnews.dta", replace

	
* 	intensive margin plot: full
	preserve
	use "`directory_path'/Outregs/vpn_log_master_dynamic_foreignnews.dta", clear
	keep if week < 20912
	sum time_week_min_nytimes1
	local max = `r(max)'
	gen c1 = 12.5  if week >= 20443 & week <= 20485
	gen c2 = 12.5  if week >= 20485 & week <= 20527
	twoway 	(area c1 week, color(gs15) lcolor(gs15) fintensity(inten100)) ///
			(area c2 week, color(gs13) lcolor(gs13) fintensity(inten100)) ///
			(connected time_week_min_nytimes0 week, yaxis(1) lcolor(navy) msymbol(O) mcolor(navy)) ///
			(connected time_week_min_nytimes1 week, yaxis(1) lcolor(cranberry) msymbol(O) mcolor(cranberry)) ///
			(connected time_week_min_foreignnews0 week, yaxis(1) lcolor(navy) lpattern(dash) msymbol(O) mcolor(navy)) ///
			(connected time_week_min_foreignnews1 week, yaxis(1) lpattern(dash) lcolor(cranberry) msymbol(O) mcolor(cranberry)), ///
			xtitle("") ///
			ytitle("Total browsing time per week (min)") ///
			yscale(r(0 12)) ///
			ylabel(0(4)12) ///
			xlabel(20454 "2016-01" 20485 "2016-02" 20514 "2016-03" 20545 "2016-04" 20575 "2016-05" 20606 "2016-06" 20636 "2016-07" 20667 "2016-08" 20698 "2016-09" 20728 "2016-10" 20759 "2016-11" 20789 "2016-12" 20820 "2017-01" 20851 "2017-02" 20879 "2017-03" 20910 "2017-04") ///
			legend(order(3 "Access: NYTimes" 4 "Access + Encour.: NYTimes" 5 "Access: all foreign news" 6 "Access + Encour.: all foreign news") col(4)) ///
			xsize(16) ///
			graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
			saving(figure_browsing_foreignnews_intensive_full, replace)
	graph export "figure_browsing_foreignnews_intensive_full.pdf", replace
	restore	
	
	
	
	
*** Figure A.7

	// load raw data
	use "`directory_path'/browsingtime_nytimes.dta", clear
	
	// reshape data (standardized variable name)
	drop if visit_day == .
	keep responseID_wave1 visit_day time_spent_min_nytimes
	reshape wide time_spent_min_nytimes, i(responseID_wave1) j(visit_day)

	// generate summary stats: average browsing time after 1st quiz
	gen time_spent_min_nytimes_sum = 0
	gen time_spent_ext_nytimes_sum = 0
	gen total_days = 0
	forvalues i = 20443/20918 {
		gen time_spent_ext_nytimes`i' = (time_spent_min_nytimes`i' > 0) & time_spent_min_nytimes`i' != . 
		replace time_spent_min_nytimes_sum = time_spent_min_nytimes_sum + time_spent_min_nytimes`i'
		replace time_spent_ext_nytimes_sum = time_spent_ext_nytimes_sum + time_spent_ext_nytimes`i'
		replace total_days = total_days + 1
		}
	gen time_spent_min_nytimes_avg = time_spent_min_nytimes_sum / total_days
	gen time_spent_ext_nytimes_avg = time_spent_ext_nytimes_sum / total_days
	
	// generate summary stats: percentage of weeks who have at least one visit to NYTimes
	forvalues i = 20443(7)20918 {
		gen time_week_min_nytimes`i' = 0
		gen time_week_ext_nytimes`i' = 0
		gen time_week_ett_nytimes`i' = 0
		forvalues j = 0/6 {
			local k = `i' + `j'
			replace time_week_min_nytimes`i' = time_week_min_nytimes`i' + time_spent_min_nytimes`k'
			replace time_week_ext_nytimes`i' = 1 	if time_spent_ext_nytimes`k' == 1
			replace time_week_ett_nytimes`i' = time_week_ett_nytimes`i' + time_spent_ext_nytimes`k'
			}
		}

	// generate indicator for active browser (based on ett)
	forvalues i = 20443(7)20918 {
		gen time_week_ett2_nytimes`i' = (time_week_ett_nytimes`i' >= 2) & time_week_ett_nytimes`i' != .
		}
		
	gen time_week_min_nytimes_sum = 0
	gen time_week_ext_nytimes_sum = 0
	gen time_week_ett_nytimes_sum = 0
	gen time_week_day_nytimes_sum = 0
	
	// generate total after 1st quiz
	gen total_weeks = 0
	forvalues i = 20485(7)20918 {
		replace time_week_min_nytimes_sum = time_week_min_nytimes_sum + time_week_min_nytimes`i'
		replace time_week_ext_nytimes_sum = time_week_ext_nytimes_sum + time_week_ext_nytimes`i'
		replace time_week_ett_nytimes_sum = time_week_ett_nytimes_sum + time_week_ett2_nytimes`i'
		replace time_week_day_nytimes_sum = time_week_day_nytimes_sum + time_week_ett_nytimes`i'
		replace total_weeks = total_weeks + 1
		}
	
	// deal with missing days
	foreach i in 20709 20793 20828 {
		replace time_week_min_nytimes`i' = time_week_min_nytimes`i' * 7 / 6
		}
		
	gen time_week_min_nytimes_avg = time_week_min_nytimes_sum / total_weeks
	gen time_week_ext_nytimes_avg = time_week_ext_nytimes_sum / total_weeks		
	gen time_week_ett_nytimes_avg = time_week_ett_nytimes_sum / total_weeks			
	gen time_week_day_nytimes_avg = time_week_day_nytimes_sum / total_weeks			
	
	// merge in adoption date
	merge 1:1 responseID_wave1 using "`directory_path'/vpn_date_adoption.dta"
	drop if _merge == 2
	drop _merge
	
	// merge in student information
	merge 1:1 responseID using "`directory_path'/vpn_treatment_roster"
	keep if _merge == 3
	drop _merge
	
	// reshape data
	keep if vpn_date_adoption != .
	keep treatment_newsletter time_week_ett2_nytimes*
	collapse time_week_ett2_nytimes*, by(treatment_newsletter)
	reshape long time_week_ett2_nytimes, i(treatment_newsletter) j(week)
	reshape wide time_week_ett2_nytimes, i(week) j(treatment_newsletter)

	// merge nytimes_wedge
	merge 1:1 week using "`directory_path'/Outregs/nytimes_newslog_wedge.dta"
	keep if _merge == 3
	drop _merge

	// extensive margin plot: full
	keep if week < 20912
	sum time_week_ett2_nytimes1
	local max = `r(max)'
	gen c1 = 0.805  if week >= 20443 & week <= 20485
	gen c2 = 0.805  if week >= 20485 & week <= 20527
	twoway 	(area c1 week, color(gs15) lcolor(gs15) fintensity(inten100)) ///
			(area c2 week, color(gs13) lcolor(gs13) fintensity(inten100)) ///
			(connected time_week_ett2_nytimes0 week, yaxis(1) lcolor(navy) msymbol(O) mcolor(navy)) ///
			(connected time_week_ett2_nytimes1 week, yaxis(1) lwidth(thick) lcolor(cranberry) msymbol(O) mcolor(cranberry)) ///
			(connected censored_prop_week week, yaxis(2) lpattern(dash) lwidth(medthick) lcolor(gs6) msymbol(th) mcolor(gs6)), ///
			xtitle("") ///
			ytitle("% students actively browsing NYTimes", axis(1)) ///
			ytitle("% pol. sensitive articles on NYTimes", axis(2)) ///
			yscale(r(0 0.30) axis(2)) ///
			ylabel(0 "0" 0.10 "10" 0.20 "20" 0.30 "30", axis(2)) ///
			xlabel(20454 "2016-01" 20485 "2016-02" 20514 "2016-03" 20545 "2016-04" 20575 "2016-05" 20606 "2016-06" 20636 "2016-07" 20667 "2016-08" 20698 "2016-09" 20728 "2016-10" 20759 "2016-11" 20789 "2016-12" 20820 "2017-01" 20851 "2017-02" 20879 "2017-03" 20910 "2017-04") ///
			legend(order(3 "Access" 4 "Access + Encour." 5 "% sensitive articles on NYTimes") col(3)) ///
			xsize(16) ///
			graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
			saving(figure_browsing_nyt_extensive_full, replace)
	graph export "figure_browsing_nyt_extensive_full.pdf", replace
	
	
	

*** Figure A.20

	// load raw data
	use "`directory_path'/browsingtime_nytimes.dta", clear
	
	// reshape data (standardized variable name)
	drop if visit_day == .
	keep responseID_wave1 visit_day time_spent_min_nytimes
	reshape wide time_spent_min_nytimes, i(responseID_wave1) j(visit_day)

	// generate summary stats: average browsing time after 1st quiz
	gen time_spent_min_nytimes_sum = 0
	gen time_spent_ext_nytimes_sum = 0
	gen total_days = 0
	forvalues i = 20443/20918 {
		gen time_spent_ext_nytimes`i' = (time_spent_min_nytimes`i' > 0) & time_spent_min_nytimes`i' != . 
		replace time_spent_min_nytimes_sum = time_spent_min_nytimes_sum + time_spent_min_nytimes`i'
		replace time_spent_ext_nytimes_sum = time_spent_ext_nytimes_sum + time_spent_ext_nytimes`i'
		replace total_days = total_days + 1
		}
	gen time_spent_min_nytimes_avg = time_spent_min_nytimes_sum / total_days
	gen time_spent_ext_nytimes_avg = time_spent_ext_nytimes_sum / total_days
	
	// generate summary stats: percentage of weeks who have at least one visit to NYTimes
	forvalues i = 20443(7)20918 {
		gen time_week_min_nytimes`i' = 0
		gen time_week_ext_nytimes`i' = 0
		gen time_week_ett_nytimes`i' = 0
		forvalues j = 0/6 {
			local k = `i' + `j'
			replace time_week_min_nytimes`i' = time_week_min_nytimes`i' + time_spent_min_nytimes`k'
			replace time_week_ext_nytimes`i' = 1 	if time_spent_ext_nytimes`k' == 1
			replace time_week_ett_nytimes`i' = time_week_ett_nytimes`i' + time_spent_ext_nytimes`k'
			}
		}

	// generate indicator for active browser (based on ett)
	forvalues i = 20443(7)20918 {
		gen time_week_ett2_nytimes`i' = (time_week_ett_nytimes`i' >= 2) & time_week_ett_nytimes`i' != .
		}
		
	gen time_week_min_nytimes_sum = 0
	gen time_week_ext_nytimes_sum = 0
	gen time_week_ett_nytimes_sum = 0
	gen time_week_day_nytimes_sum = 0
	
	// generate total after 1st quiz
	gen total_weeks = 0
	forvalues i = 20485(7)20918 {
		replace time_week_min_nytimes_sum = time_week_min_nytimes_sum + time_week_min_nytimes`i'
		replace time_week_ext_nytimes_sum = time_week_ext_nytimes_sum + time_week_ext_nytimes`i'
		replace time_week_ett_nytimes_sum = time_week_ett_nytimes_sum + time_week_ett2_nytimes`i'
		replace time_week_day_nytimes_sum = time_week_day_nytimes_sum + time_week_ett_nytimes`i'
		replace total_weeks = total_weeks + 1
		}
	
	// deal with missing days
	foreach i in 20709 20793 20828 {
		replace time_week_min_nytimes`i' = time_week_min_nytimes`i' * 7 / 6
		}
		
	gen time_week_min_nytimes_avg = time_week_min_nytimes_sum / total_weeks
	gen time_week_ext_nytimes_avg = time_week_ext_nytimes_sum / total_weeks		
	gen time_week_ett_nytimes_avg = time_week_ett_nytimes_sum / total_weeks			
	gen time_week_day_nytimes_avg = time_week_day_nytimes_sum / total_weeks			
	
	// merge in adoption date
	merge 1:1 responseID_wave1 using "`directory_path'/vpn_date_adoption.dta"
	drop if _merge == 2
	drop _merge
	
	// merge in student information
	merge 1:1 responseID using "`directory_path'/vpn_treatment_roster"
	keep if treatment_vpn == 1
	drop _merge
	
	// generate indicator for ever adopted
	gen adopted = (vpn_date_adoption != .)

	// merge in treatment_expiration info
	merge 1:1 responseID_wave1 using "`directory_path'/vpn_treatment_expiration.dta"
	replace treatment_vpnexpiration = 0 	if treatment_vpnexpiration == .
	keep if _merge == 3
	keep if treatment_vpnexpiration != .
	
	// adoption count by day: treatment_vpnexpiration == 1
	preserve
	keep if treatment_vpnexpiration == 1
	
	gen c = 1
	egen count_nl_yes = sum(c)
	bysort vpn_date_adoption: egen adoption_nl_yes = sum(adopted)
	keep vpn_date_adoption adoption_nl_yes count_nl_yes
	duplicates drop
	drop if vpn_date_adoption == .
	save "`directory_path'/Outregs/adoption_count_treatment_vpnexpiration_yes.dta", replace
	restore
	
	// adoption count by day: treatment_vpnexpiration == 0
	preserve
	keep if treatment_vpnexpiration == 0
	
	gen c = 1
	egen count_nl_no = sum(c)
	bysort vpn_date_adoption: egen adoption_nl_no = sum(adopted)
	keep vpn_date_adoption adoption_nl_no count_nl_no
	duplicates drop
	drop if vpn_date_adoption == .
	save "`directory_path'/Outregs/adoption_count_treatment_vpnexpiration_no.dta", replace
	restore
		
	// merge in data
	preserve
	use "`directory_path'/vpn_log_dates.dta", clear
	drop if dates == ""
	keep if dates_code <= 20603
	keep dates dates_code
	// resolve 2016-02-19
	replace dates = "2016-02-19" 	if dates_code == 20503
	duplicates drop
	
	rename dates_code vpn_date_adoption
	
	merge 1:1 vpn_date_adoption using "`directory_path'/Outregs/adoption_count_treatment_vpnexpiration_yes.dta"
	drop _merge
	merge 1:1 vpn_date_adoption using "`directory_path'/Outregs/adoption_count_treatment_vpnexpiration_no.dta"
	drop _merge
	erase "`directory_path'/Outregs/adoption_count_treatment_vpnexpiration_yes.dta"
	erase "`directory_path'/Outregs/adoption_count_treatment_vpnexpiration_no.dta"
	
	replace adoption_nl_yes = 0 	if adoption_nl_yes == .
	replace adoption_nl_no = 0 		if adoption_nl_no == .
	
	// switch to percentage
	egen count_nl_yes1 = mean(count_nl_yes) 
	drop count_nl_yes
	rename count_nl_yes1 count_nl_yes
	egen count_nl_no1 = mean(count_nl_no) 
	drop count_nl_no
	rename count_nl_no1 count_nl_no	
	
	// cut data by weeks
	egen vpn_date_interval = cut(vpn_date_adoption), at(20438(5)20603)
	bysort vpn_date_interval: egen adoption_total_nl_yes = sum(adoption_nl_yes)
	bysort vpn_date_interval: egen adoption_total_nl_no = sum(adoption_nl_no)
	
	gen adoption_pc_nl_yes = adoption_total_nl_yes / count_nl_yes
	gen adoption_pc_nl_no = adoption_total_nl_no / count_nl_no
	
	keep vpn_date_interval adoption_pc_nl_yes adoption_pc_nl_no 
	duplicates drop
	
	// generate cumulative adoption
	sort vpn_date_interval
	gen adoption_cum_nl_yes = adoption_pc_nl_yes[1]
	replace adoption_cum_nl_yes = adoption_pc_nl_yes[_n] + adoption_cum_nl_yes[_n-1] if _n > 1
	gen adoption_cum_nl_no = adoption_pc_nl_no[1]
	replace adoption_cum_nl_no = adoption_pc_nl_no[_n] + adoption_cum_nl_no[_n-1] if _n > 1

	qui su adoption_cum_nl_yes
	local max = `r(max)'
	gen c1 = `max' + 0.03  if vpn_date_interval >= 20438 & vpn_date_interval <= 20488
	gen c2 = `max' + 0.03  if vpn_date_interval >= 20488 & vpn_date_interval <= 20537

	twoway 	(area c1 vpn_date_interval, color(gs15) lcolor(gs15) fintensity(inten100)) ///
			(area c2 vpn_date_interval, color(gs13) lcolor(gs13) fintensity(inten100)) ///
			(connected adoption_cum_nl_no vpn_date_interval, msymbol(O) mcolor(navy) lwidth(medthick) lcolor(navy)) ///
			(connected adoption_cum_nl_yes vpn_date_interval, msymbol(O) mcolor(gray) lwidth(thick) lcolor(gray)), ///
			text(0.75 20463 "Unincentivized" "encouragement" "period", size(small)) ///
			text(0.75 20510 "Incentivized" "encouragement" "period", size(small)) ///
		    ytitle("Cumulative activation rate") ///
			xtitle("") ///
			xlabel(20454 "2016-01" 20485 "2016-02" 20514 "2016-03" 20545 "2016-04" 20575 "2016-05") ///
			legend(on order(3 4) label(3 "Expiration date uninformed") label(4 "Expiration date informed")) ///
			graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
			saving(figure_adoption_cumulative_expirationtreatment, replace)
	
	graph export "figure_adoption_cumulative_expirationtreatment.pdf", replace
	
	restore	

	

*** Table A.10

* 	Panel A

	// load raw data
	use "`directory_path'/browsingtime_nytimes.dta", clear
	
	// reshape data (standardized variable name)
	drop if visit_day == .
	keep responseID_wave1 visit_day time_spent_min_nytimes
	reshape wide time_spent_min_nytimes, i(responseID_wave1) j(visit_day)

	// generate summary stats: average browsing time after 1st quiz
	gen time_spent_min_nytimes_sum = 0
	gen time_spent_ext_nytimes_sum = 0
	gen total_days = 0
	forvalues i = 20443/20918 {
		gen time_spent_ext_nytimes`i' = (time_spent_min_nytimes`i' > 0) & time_spent_min_nytimes`i' != . 
		replace time_spent_min_nytimes_sum = time_spent_min_nytimes_sum + time_spent_min_nytimes`i'
		replace time_spent_ext_nytimes_sum = time_spent_ext_nytimes_sum + time_spent_ext_nytimes`i'
		replace total_days = total_days + 1
		}
	gen time_spent_min_nytimes_avg = time_spent_min_nytimes_sum / total_days
	gen time_spent_ext_nytimes_avg = time_spent_ext_nytimes_sum / total_days
	
	// generate summary stats: percentage of weeks who have at least one visit to NYTimes
	forvalues i = 20443(7)20918 {
		gen time_week_min_nytimes`i' = 0
		gen time_week_ext_nytimes`i' = 0
		gen time_week_ett_nytimes`i' = 0
		forvalues j = 0/6 {
			local k = `i' + `j'
			replace time_week_min_nytimes`i' = time_week_min_nytimes`i' + time_spent_min_nytimes`k'
			replace time_week_ext_nytimes`i' = 1 	if time_spent_ext_nytimes`k' == 1
			replace time_week_ett_nytimes`i' = time_week_ett_nytimes`i' + time_spent_ext_nytimes`k'
			}
		}

	// generate indicator for active browser (based on ett)
	forvalues i = 20443(7)20918 {
		gen time_week_ett2_nytimes`i' = (time_week_ett_nytimes`i' >= 2) & time_week_ett_nytimes`i' != .
		}
		
	gen time_week_min_nytimes_sum = 0
	gen time_week_ext_nytimes_sum = 0
	gen time_week_ett_nytimes_sum = 0
	gen time_week_day_nytimes_sum = 0
	
	// generate total after 1st quiz
	gen total_weeks = 0
	forvalues i = 20485(7)20918 {
		replace time_week_min_nytimes_sum = time_week_min_nytimes_sum + time_week_min_nytimes`i'
		replace time_week_ext_nytimes_sum = time_week_ext_nytimes_sum + time_week_ext_nytimes`i'
		replace time_week_ett_nytimes_sum = time_week_ett_nytimes_sum + time_week_ett2_nytimes`i'
		replace time_week_day_nytimes_sum = time_week_day_nytimes_sum + time_week_ett_nytimes`i'
		replace total_weeks = total_weeks + 1
		}
	
	// deal with missing days
	foreach i in 20709 20793 20828 {
		replace time_week_min_nytimes`i' = time_week_min_nytimes`i' * 7 / 6
		}
		
	gen time_week_min_nytimes_avg = time_week_min_nytimes_sum / total_weeks
	gen time_week_ext_nytimes_avg = time_week_ext_nytimes_sum / total_weeks		
	gen time_week_ett_nytimes_avg = time_week_ett_nytimes_sum / total_weeks			
	gen time_week_day_nytimes_avg = time_week_day_nytimes_sum / total_weeks			
	
	// merge in adoption date
	merge 1:1 responseID_wave1 using "`directory_path'/vpn_date_adoption.dta"
	drop if _merge == 2
	drop _merge
	
	// merge in student information
	merge 1:1 responseID using "`directory_path'/vpn_treatment_roster"
	keep if _merge == 3
	drop _merge
	
	// keep relevant data
	keep if vpn_date_adoption != .
	keep time_week_min_nytimes* responseID_wave1 treatment_newsletter
	drop time_week_min_nytimes_sum time_week_min_nytimes_avg

	// reshape data
	reshape long time_week_min_nytimes, i(responseID_wave1) j(week)
	
	// merge nytimes_wedge
	merge m:1 week using "`directory_path'/Outregs/nytimes_newslog_wedge.dta"
	keep if _merge == 3
	drop _merge

	// regression: group-AE main
	reg time_week_min_nytimes censored_prop_week if treatment_newsletter == 1 & week > 20527, absorb(responseID_wave1) r
	preserve
	keep censored_prop_week week
	duplicates drop
	su censored_prop_week if week > 20527
	restore
	su time_week_min_nytimes if treatment_newsletter == 1 & week > 20527
	
	// regression: group-AE robustness - remove Trump
	reg time_week_min_nytimes censored_prop_week if treatment_newsletter == 1 & week > 20527 & week != 20765, absorb(responseID_wave1) r
	preserve
	keep censored_prop_week week
	duplicates drop
	su censored_prop_week if week > 20527 & week != 20765
	restore
	su time_week_min_nytimes if treatment_newsletter == 1 & week > 20527 & week != 20765
	
	// regression: rgroup-AE obustness - remove week long holidays
	reg time_week_min_nytimes censored_prop_week if treatment_newsletter == 1 & week > 20527 & week != 20730 & week != 20849, absorb(responseID_wave1) r	
	preserve
	keep censored_prop_week week
	duplicates drop
	su censored_prop_week if week > 20527 & week != 20730 & week != 20849
	restore
	su time_week_min_nytimes if treatment_newsletter == 1 & week > 20527 & week != 20730 & week != 20849

	// regression: group-AE robustness - remove week Trump & long holidays
	reg time_week_min_nytimes censored_prop_week if treatment_newsletter == 1 & week != 20765 & week > 20527 & week != 20730 & week != 20849, absorb(responseID_wave1) r	
	preserve
	keep censored_prop_week week
	duplicates drop
	su censored_prop_week if week > 20527 & week != 20765 & week != 20730 & week != 20849
	restore
	su time_week_min_nytimes if treatment_newsletter == 1 & week > 20527 & week != 20765 & week != 20730 & week != 20849

	// regression: group-A main
	reg time_week_min_nytimes censored_prop_week if treatment_newsletter == 0 & week > 20527, absorb(responseID_wave1) r
	preserve
	keep censored_prop_week week
	duplicates drop
	su censored_prop_week if week > 20527
	restore
	su time_week_min_nytimes if treatment_newsletter == 0 & week > 20527
	
	// regression: group-A robustness - remove Trump
	reg time_week_min_nytimes censored_prop_week if treatment_newsletter == 0 & week > 20527 & week != 20765, absorb(responseID_wave1) r
	preserve
	keep censored_prop_week week
	duplicates drop
	su censored_prop_week if week > 20527 & week != 20765
	restore
	su time_week_min_nytimes if treatment_newsletter == 0 & week > 20527 & week != 20765
	
	// regression: group-A robustness - remove week long holidays
	reg time_week_min_nytimes censored_prop_week if treatment_newsletter == 0 & week > 20527 & week != 20730 & week != 20849, absorb(responseID_wave1) r	
	preserve
	keep censored_prop_week week
	duplicates drop
	su censored_prop_week if week > 20527 & week != 20730 & week != 20849
	restore
	su time_week_min_nytimes if treatment_newsletter == 0 & week > 20527 & week != 20730 & week != 20849

	// regression: group-A robustness - remove week Trump & long holidays
	reg time_week_min_nytimes censored_prop_week if treatment_newsletter == 0 & week != 20765 & week > 20527 & week != 20730 & week != 20849, absorb(responseID_wave1) r	
	preserve
	keep censored_prop_week week
	duplicates drop
	su censored_prop_week if week > 20527 & week != 20765 & week != 20730 & week != 20849
	restore
	su time_week_min_nytimes if treatment_newsletter == 0 & week > 20527 & week != 20765 & week != 20730 & week != 20849




* 	Panel B

	// load raw data
	use "`directory_path'/browsingtime_nytimes.dta", clear
	
	// merge with other top foreign news sites
	foreach var in economist reddit cnn theguardian huffingtonpost foxnews bbc bloomberg wsj usatoday reuters nbcnews ft {
		merge 1:1 responseID_wave1 visit_day using "`directory_path'/browsingtime_`var'.dta"
		drop _merge
		}
	
	// generate total browsing time
	gen time_spent_min_foreignnews = time_spent_min_nytimes
	foreach var in economist reddit cnn theguardian huffingtonpost foxnews bbc bloomberg wsj usatoday reuters nbcnews ft {
		replace time_spent_min_foreignnews = time_spent_min_foreignnews + time_spent_min_`var'
		}

	// generate non-NYTimes browsing time
	gen time_spent_min_nonnyt = time_spent_min_foreignnews - time_spent_min_nytimes
	
	// reshape data (standardized variable name)
	drop if visit_day == .
	keep responseID_wave1 visit_day time_spent_min_nonnyt
	reshape wide time_spent_min_nonnyt, i(responseID_wave1) j(visit_day)

	// generate summary stats: average browsing time after 1st quiz
	gen time_spent_min_nonnyt_sum = 0
	gen time_spent_ext_nonnyt_sum = 0
	gen total_days = 0
	forvalues i = 20443/20918 {
		gen time_spent_ext_nonnyt`i' = (time_spent_min_nonnyt`i' > 0) & time_spent_min_nonnyt`i' != . 
		replace time_spent_min_nonnyt_sum = time_spent_min_nonnyt_sum + time_spent_min_nonnyt`i'
		replace time_spent_ext_nonnyt_sum = time_spent_ext_nonnyt_sum + time_spent_ext_nonnyt`i'
		replace total_days = total_days + 1
		}
	gen time_spent_min_nonnyt_avg = time_spent_min_nonnyt_sum / total_days
	gen time_spent_ext_nonnyt_avg = time_spent_ext_nonnyt_sum / total_days
	
	// generate summary stats: percentage of weeks who have at least one visit to NYTimes
	forvalues i = 20443(7)20918 {
		gen time_week_min_nonnyt`i' = 0
		gen time_week_ext_nonnyt`i' = 0
		gen time_week_ett_nonnyt`i' = 0
		forvalues j = 0/6 {
			local k = `i' + `j'
			replace time_week_min_nonnyt`i' = time_week_min_nonnyt`i' + time_spent_min_nonnyt`k'
			replace time_week_ext_nonnyt`i' = 1 	if time_spent_ext_nonnyt`k' == 1
			replace time_week_ett_nonnyt`i' = time_week_ett_nonnyt`i' + time_spent_ext_nonnyt`k'
			}
		}

	// generate indicator for active browser (based on ett)
	forvalues i = 20443(7)20918 {
		gen time_week_ett2_nonnyt`i' = (time_week_ett_nonnyt`i' >= 2) & time_week_ett_nonnyt`i' != .
		}
		
	gen time_week_min_nonnyt_sum = 0
	gen time_week_ext_nonnyt_sum = 0
	gen time_week_ett_nonnyt_sum = 0
	gen time_week_day_nonnyt_sum = 0
	
	// generate total after 1st quiz
	gen total_weeks = 0
	forvalues i = 20485(7)20918 {
		replace time_week_min_nonnyt_sum = time_week_min_nonnyt_sum + time_week_min_nonnyt`i'
		replace time_week_ext_nonnyt_sum = time_week_ext_nonnyt_sum + time_week_ext_nonnyt`i'
		replace time_week_ett_nonnyt_sum = time_week_ett_nonnyt_sum + time_week_ett2_nonnyt`i'
		replace time_week_day_nonnyt_sum = time_week_day_nonnyt_sum + time_week_ett_nonnyt`i'
		replace total_weeks = total_weeks + 1
		}
	
	// deal with missing days
	foreach i in 20709 20793 20828 {
		replace time_week_min_nonnyt`i' = time_week_min_nonnyt`i' * 7 / 6
		}
		
	gen time_week_min_nonnyt_avg = time_week_min_nonnyt_sum / total_weeks
	gen time_week_ext_nonnyt_avg = time_week_ext_nonnyt_sum / total_weeks		
	gen time_week_ett_nonnyt_avg = time_week_ett_nonnyt_sum / total_weeks			
	gen time_week_day_nonnyt_avg = time_week_day_nonnyt_sum / total_weeks			
	
	// merge in adoption date
	merge 1:1 responseID_wave1 using "`directory_path'/vpn_date_adoption.dta"
	drop if _merge == 2
	drop _merge
	
	// merge in student information
	merge 1:1 responseID using "`directory_path'/vpn_treatment_roster"
	keep if _merge == 3
	drop _merge
	
	// keep relevant data
	keep if vpn_date_adoption != .
	keep time_week_min_nonnyt* responseID_wave1 treatment_newsletter
	drop time_week_min_nonnyt_sum time_week_min_nonnyt_avg

	// reshape data
	reshape long time_week_min_nonnyt, i(responseID_wave1) j(week)
	
	// merge nytimes_wedge
	merge m:1 week using "`directory_path'/Outregs/nytimes_newslog_wedge.dta"
	keep if _merge == 3
	drop _merge

	// regression: group-AE main
	reg time_week_min_nonnyt censored_prop_week if treatment_newsletter == 1 & week > 20527, absorb(responseID_wave1) r
	preserve
	keep censored_prop_week week
	duplicates drop
	su censored_prop_week if week > 20527
	restore
	su time_week_min_nonnyt if treatment_newsletter == 1 & week > 20527
	
	// regression: group-AE robustness - remove Trump
	reg time_week_min_nonnyt censored_prop_week if treatment_newsletter == 1 & week > 20527 & week != 20765, absorb(responseID_wave1) r
	preserve
	keep censored_prop_week week
	duplicates drop
	su censored_prop_week if week > 20527 & week != 20765
	restore
	su time_week_min_nonnyt if treatment_newsletter == 1 & week > 20527 & week != 20765
	
	// regression: rgroup-AE obustness - remove week long holidays
	reg time_week_min_nonnyt censored_prop_week if treatment_newsletter == 1 & week > 20527 & week != 20730 & week != 20849, absorb(responseID_wave1) r	
	preserve
	keep censored_prop_week week
	duplicates drop
	su censored_prop_week if week > 20527 & week != 20730 & week != 20849
	restore
	su time_week_min_nonnyt if treatment_newsletter == 1 & week > 20527 & week != 20730 & week != 20849

	// regression: group-AE robustness - remove week Trump & long holidays
	reg time_week_min_nonnyt censored_prop_week if treatment_newsletter == 1 & week != 20765 & week > 20527 & week != 20730 & week != 20849, absorb(responseID_wave1) r	
	preserve
	keep censored_prop_week week
	duplicates drop
	su censored_prop_week if week > 20527 & week != 20765 & week != 20730 & week != 20849
	restore
	su time_week_min_nonnyt if treatment_newsletter == 1 & week > 20527 & week != 20765 & week != 20730 & week != 20849

	// regression: group-A main
	reg time_week_min_nonnyt censored_prop_week if treatment_newsletter == 0 & week > 20527, absorb(responseID_wave1) r
	preserve
	keep censored_prop_week week
	duplicates drop
	su censored_prop_week if week > 20527
	restore
	su time_week_min_nonnyt if treatment_newsletter == 0 & week > 20527
	
	// regression: group-A robustness - remove Trump
	reg time_week_min_nonnyt censored_prop_week if treatment_newsletter == 0 & week > 20527 & week != 20765, absorb(responseID_wave1) r
	preserve
	keep censored_prop_week week
	duplicates drop
	su censored_prop_week if week > 20527 & week != 20765
	restore
	su time_week_min_nonnyt if treatment_newsletter == 0 & week > 20527 & week != 20765
	
	// regression: group-A robustness - remove week long holidays
	reg time_week_min_nonnyt censored_prop_week if treatment_newsletter == 0 & week > 20527 & week != 20730 & week != 20849, absorb(responseID_wave1) r	
	preserve
	keep censored_prop_week week
	duplicates drop
	su censored_prop_week if week > 20527 & week != 20730 & week != 20849
	restore
	su time_week_min_nonnyt if treatment_newsletter == 0 & week > 20527 & week != 20730 & week != 20849

	// regression: group-A robustness - remove week Trump & long holidays
	reg time_week_min_nonnyt censored_prop_week if treatment_newsletter == 0 & week != 20765 & week > 20527 & week != 20730 & week != 20849, absorb(responseID_wave1) r	
	preserve
	keep censored_prop_week week
	duplicates drop
	su censored_prop_week if week > 20527 & week != 20765 & week != 20730 & week != 20849
	restore
	su time_week_min_nonnyt if treatment_newsletter == 0 & week > 20527 & week != 20765 & week != 20730 & week != 20849
	
