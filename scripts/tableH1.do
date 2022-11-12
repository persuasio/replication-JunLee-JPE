************
* SCRIPT:   tableH1.do
* PURPOSE:  Creates Table H1
*
* ACKNOWLEDGMENT
*	The orginal dataset "NTV_Individual_Data.dta" is from Enikolopov, Petrova, and Zhuravskaya (AER, 2011).
************

use "$Persuasion/data/EnikolopovPetrovaZhuravskaya2011/NTV_Individual_Data.dta", clear

local basic = "logpop98 wage98_ln"
local sociodem = "male age educ1 married consump"

local party_list = "Unity OVR"

*** create outcome variables based on directional treatment ***

gen y_vote_Unity = 0
replace y_vote_Unity = 1 if vote_Unity == 0   

gen y_vote_OVR = 0
replace y_vote_OVR = 1 if vote_OVR == 1   

*** create a binary instrument ***

svyset [pweight=kishweig]

_pctile tvmaxtveloss5050powerA [pweight=kishweig], p(50)	 
scalar IV_p50 = r(r1)	 

	gen IV_NTV = 0
replace IV_NTV = 1 if tvmaxtveloss5050powerA > IV_p50

*** predict Z ***

qui probit IV_NTV `sociodem' `basic' [pweight=kishweig]  
predict phat_Z1

*** predict T given Z=1 and Z=0 ***

qui probit Watches_NTV_1999 `sociodem' `basic' [pweight=kishweig] if IV_NTV==1 
predict phat_T_Z1

qui probit Watches_NTV_1999 `sociodem' `basic' [pweight=kishweig] if IV_NTV==0 
predict phat_T_Z0

*** input for ATE estimation ***

foreach party in  `party_list'    {

gen newvote_`party' = 0
replace newvote_`party' = 1 if y_vote_`party' == 1  

*** predict Y given Z=1 and Z=0 ***

qui probit newvote_`party' `sociodem' `basic' [pweight=kishweig] if IV_NTV==1 
predict phat_Y_Z1_`party'

qui probit newvote_`party' `sociodem' `basic' [pweight=kishweig] if IV_NTV==0 
predict phat_Y_Z0_`party'

*** compute relevant quantities ***

gen      w1_`party' = IV_NTV/phat_Z1
gen      w0_`party' = (1-IV_NTV)/(1-phat_Z1)
gen      t1_`party' = w1_`party'*(newvote_`party'-phat_Y_Z1_`party') 
gen      t0_`party' = w0_`party'*(newvote_`party'-phat_Y_Z0_`party')

*** average out with respect to covariates ***

svyset [pweight=kishweig]

svy: mean phat_Y_Z1_`party'
scalar beta1 = e(b)[1,1] 

svy: mean phat_Y_Z0_`party'
scalar beta0 = e(b)[1,1] 

*** compute the lower bound on the average persuasion rate ***

scalar ate_lb = (beta1 - beta0)/(1-beta0)	

matrix results_`party' = J(3,1,0)
matrix results_`party'[1,1] = ate_lb

*** compute the standard errors and one-sided 80% confidence intervals ***

gen infl_t1_`party' = t1_`party' + (phat_Y_Z1_`party' - beta1)
gen infl_t0_`party' = t0_`party' + (phat_Y_Z0_`party' - beta0)
gen infl_`party' = (1/(1-beta0))*(infl_t1_`party' - ((1-beta1)/(1-beta0))*infl_t0_`party')

svy: mean infl_`party'
matrix result_se = e(V)
matrix results_`party'[2,1] = sqrt(result_se[1,1])				 
matrix results_`party'[3,1] = results_`party'[1,1]-1.645*results_`party'[2,1]	

}

*******************************************
****  Printing result in a LaTeX file *****
*******************************************
						  
matrix results = (results_Unity, results_OVR)

frmttable using "$Persuasion/results/tableH1", statmat(results) tex replace sdec(3) ///
	title("Table H1") ///
	ctitles("", Not Voting for Unity, Voting for Unity )  ///
	rtitles("Point estimate of the lower bound" \ "Standard error of the lower bound" \ ///
	"One-sided 95\% confidence interval" ) 
