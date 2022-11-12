************
* SCRIPT:   tableD2.do
* PURPOSE:  Creates Table D2
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

*** predict the exposure rate, that is, e(X,Z) = P(T=1|X,Z) ***
*** specification: linear in X and cubic in Z, while interacting some of X and Z ***

reg Watches_NTV_1999 `sociodem' `basic' ///
    c.tvmaxtveloss5050powerA c.tvmaxtveloss5050powerA#c.tvmaxtveloss5050powerA /// 
    c.tvmaxtveloss5050powerA#c.tvmaxtveloss5050powerA#c.tvmaxtveloss5050powerA /// 
	c.tvmaxtveloss5050powerA#i.male c.tvmaxtveloss5050powerA#c.age ///
	c.tvmaxtveloss5050powerA#i.educ1 c.tvmaxtveloss5050powerA#i.married ///
	c.tvmaxtveloss5050powerA#c.consump ///
	[pweight=kishweig]
predict phat
*** fitted probability is truncated to be between 0 and 1 ***
replace phat = 0 if phat < 0 & phat != .
replace phat = 1 if phat > 1 & phat != .

*** estimate the marginal persuasion rates ***

foreach party in `party_list' {

	* generate Y*(1-T)
    gen notwatch_vote_`party' = 0
	replace notwatch_vote_`party' = 1 if y_vote_`party' == 1 & Watches_NTV_1999==0
		
	* linear probability model of Y given X and e(X,Z)
	* specification: linear in X and cubic in e(X,Z), while interacting some of X and e(X,Z)
	reg y_vote_`party' `sociodem' `basic' ///
			c.phat c.phat#c.phat c.phat#c.phat#c.phat  ///
			c.phat#i.male c.phat#c.age c.phat#i.educ1 c.phat#i.married c.phat#c.consump ///
			[pweight=kishweig]
	
	margins [pweight=kishweig], dydx(phat) at(phat = (0.4(0.1)0.6))
	matrix num = r(b)

	* linear probability model of Y*(1-T) given X and e(X,Z)
	* specification: linear in X and cubic in e(X,Z), while interacting some of X and e(X,Z)
	reg notwatch_vote_`party' `sociodem' `basic' ///
			c.phat c.phat#c.phat c.phat#c.phat#c.phat  ///
			c.phat#i.male c.phat#c.age c.phat#i.educ1 c.phat#i.married c.phat#c.consump ///
			[pweight=kishweig]
	
	margins [pweight=kishweig], dydx(phat) at(phat = (0.4(0.1)0.6))
	matrix den = J(1,3,1) + r(b)
	
	mata: st_matrix("mte", st_matrix("num") :/ st_matrix("den"))
	matrix list mte

	matrix results_`party' = mte'
	
}

matrix results = (results_Unity, results_OVR)

frmttable using "$Persuasion/results/tableD2", statmat(results) tex replace sdec(3) ///
	ctitles("", "Not Vote for Unity", "Vote for OVR" )  ///
	rtitles("v = 0.4" \ "v = 0.5" \ "v = 0.6") ///
	title("Estimates of Marginal Persuasion Rates")
