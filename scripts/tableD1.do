************
* SCRIPT:   tableD1.do
* PURPOSE:  Creates Table D1
*
* ACKNOWLEDGMENT
*	The orginal dataset "FoxNewsFinalDataQJE" is from DellaVigna and Kaplan (QJE, 2007).
*	The dataset is available at http://eml.berkeley.edu/~sdellavi/index.html.
*   We thank the authors of the paper to make their data available online. 
*
* NOTES
* 	An alternative measure of the dependent variable is used. 
*   That is, Y = (Republican voting share as a share of the voting-age population) 
************

use "$Persuasion/data/DellaVignaKaplan2007/FoxNewsFinalDataQJE.dta", clear

*** Generate macros of controls
** Census Controls
foreach tc in "2000" "00m90" {
	local contrcens`tc'= "pop`tc' hs`tc' hsp`tc' college`tc' male`tc' black`tc' hisp`tc' empl`tc' unempl`tc' married`tc' income`tc' urban`tc'"
	local contrcenssh`tc'stat= "pop`tc'stat hsp`tc' college`tc' black`tc' hisp`tc' unempl`tc' urban`tc'"
}
** Cable Controls
foreach tc in "2000" {
	local contrcbl`tc'= "poptot`tc'd2-poptot`tc'd10 noch`tc'd2-noch`tc'd10"
}

** Keep only basic sample
keep if sample12000==1

* Table V -- Presidential Effects -- Robustness and Presidential Delayed Effects
* Lagged vote share

gen new_y2000 = reppresvotes2000/pop18p2000
gen new_y1996 = reppresvotes1996/pop18p1996

sum new_y2000 new_y1996 [aweight=totpresvotes1996]

**************************************************
*** Estimating parameters for persuation rates ***

local fixed_effect_list = "diststate countystate"

foreach fixed_effect in  `fixed_effect_list' {
	
matrix results_`fixed_effect' = J(4,1,.)	

*** regression of Y on Z
areg new_y2000 new_y1996 foxnews2000 `contrcens2000' `contrcens00m90' `contrcbl2000' ///
	if sample12000==1 [aweight=pop18p1996], a(`fixed_effect') robust cluster(account2000)

matrix reg_b = e(b)
predict yhat, xbd 

*** predict Y given Z=1 and Z=0

scalar coeff_b = reg_b[1,2]
scalar list coeff_b
gen yhat1 = yhat + coeff_b - coeff_b*foxnews2000
gen yhat0 = yhat - coeff_b*foxnews2000

replace yhat1 = min(max(yhat1,0),1)
replace yhat0 = min(max(yhat0,0),1)

gen thetahat_num = (yhat1-yhat0)
gen thetahat_den = (1-yhat0)
replace thetahat_den = max(thetahat_den, 1e-8)

*** Lower bound for theta 

sum thetahat_num if auddiaryScar > 0 [aweight=auddiaryScar]
scalar avg_num = r(mean)

sum thetahat_den if auddiaryScar > 0 [aweight=auddiaryScar]
scalar avg_den = r(mean)

scalar avg_lb = avg_num/avg_den
matrix results_`fixed_effect'[1,1] = avg_lb

*** regression of T on Z
areg foxanyScar foxnews2000 `contrcens2000' `contrcens00m90' `contrcbl2000' ///
	if sample12000==1 [aweight=auddiaryScar], a(`fixed_effect') robust cluster(account2000)
matrix reg_b = e(b)
predict ehat, xbd 

*** predict exposure rates given Z=1 and Z=0

scalar coeff_b = reg_b[1,1]
scalar list coeff_b
gen ehat1 = ehat + coeff_b - coeff_b*foxnews2000
gen ehat0 = ehat - coeff_b*foxnews2000

replace ehat1 = min(max(ehat1,0),1)
replace ehat0 = min(max(ehat0,0),1)

*** compute the upper bound for theta

gen ub_num1 = yhat1 + 1 - ehat1
gen ub_num2 = yhat0 - ehat0
gen ub_num = min(1, ub_num1) - max(0,ub_num2)
gen ub_den = 1 - max(0,ub_num2)
replace ub_den = max(ub_den, 1e-8)

sum ub_num if auddiaryScar > 0 [aweight=auddiaryScar]
scalar ub_avg_num = r(mean)

sum ub_den if auddiaryScar > 0 [aweight=auddiaryScar]
scalar ub_avg_den = r(mean)

matrix results_`fixed_effect'[2,1] = ub_avg_num/ub_avg_den

*** compute the lower bound for theta_local

gen late_weight = (ehat1-ehat0)
gen theta_local_den = min(thetahat_den, late_weight)
sum theta_local_den if auddiaryScar > 0 [aweight=auddiaryScar]
scalar local_den = r(mean)
scalar theta_local = avg_num/local_den

matrix results_`fixed_effect'[3,1] = theta_local
matrix results_`fixed_effect'[4,1] = 1

drop yhat* theta* ehat* ub_* late*

}

*******************************************
****  Printing result in a LaTeX file *****
*******************************************
						  					  
matrix results = (results_diststate , results_countystate)

frmttable using "$Persuasion/results/tableD1", statmat(results) tex replace sdec(3) ///
	ctitles("", U.S. House district fixed effects, County fixed effects )  ///
	rtitles("APR (LB)" \ "APR (UB)" \ "LPR (LB)" \ "LPR (UB)" )  ///
	title("TABLE D1. Persuasion Rates: Fox News Effects")
	
