************
* SCRIPT:   tableE2.do
* PURPOSE:  Creates Table E2
*
* ACKNOWLEDGMENT
*	The orginal dataset "CharityOutputQJE" is from DellaVigna, List, and Malmendier (QJE, 2012).
*	The data collected in DLM are available at http://eml.berkeley.edu/~sdellavi/index.html. 
************

** Load data with experimental observations
use "$Persuasion/data/DellaVignaListMalmendier2012/CharityOutputQJE.dta", clear
memory 
su

** Code treatment with 2-week warning (2Ww) together with other Warning (flyer) treatments 
replace treatment="W" if treatment=="2Ww"

** Initial Sample
count
count if (charity=="Ecu" | charity=="LaRabida")
count if (charity=="Sv2008")
count if (charity=="Sv2009")

*** Drop two groups of observations which the solicitors did not contact
* I. Households with no solicitor sign
* II. Households where solicitor could not knock on door (big dog barking, house for sale, etc.)
* Important to eliminate these observations because these households are disproportionally
* in No-Warning treatment (since some are excluded after flyering), 
* and have no giving/no survey completion
* nosolsign==1 indicates that these households have a sign that they do not accept solicitors
count if nosol==1
ta treatment nosol,row
drop if nosol==1
drop nosol
count if toel==1
ta treatment toel,row
drop if toeliminate==1
drop toeliminate
** Drop 5 solicitor-days with serious inconsistencies
** treatment is no-flyer, but solicitors report flyer on door throughout
drop if sol=="Angelena" & (date==mdy(7,27,2008)|date==mdy(7,13,2008))
drop if sol=="Shedora" & (date==mdy(7,27,2008)|date==mdy(8,10,2008))
drop if sol=="Tehmur" & date==mdy(6,1,2008)
drop if sol=="Phillip" & date==mdy(8,9,2008)
drop if sol=="Robert" & date==mdy(7,13,2008) & (hour==11|hour==13)
count

** Final Sample
ta treatment,m
count
count if (charity=="Ecu" | charity=="LaRabida")
count if (charity=="Sv2008")
count if (charity=="Sv2009")

**** Information in Experimental Design Section
* Number of solicitors -- 92
codebook solicitor
* Number of households
ta location year

** Distinguish three waves of charity data
gen month=month(date)
gen dwave=1 if ((month==7 | month==8) & charity =="LaRabida") & year==2008
replace dwave=0 if ((month==4 | month==5 | month==6) | (charity=="Ecu" & month==7)) & year==2008
replace dwave=2 if (month==9 | month==10) & year==2008
* Randomization by wave
ta date dwave
bys dwave: ta treatment charity if charityall=="All", column
bys dwave: ta treatment charity if charityall=="All", row

** Generate fixed effects
egen sodate=concat(solicitor date)
egen treatmentby=concat(treatment charity)
egen byte grsol=group(solicitor)
egen byte grdatloc=group(date location)
egen byte grdatlocsol=group(date location solicitor)
egen byte grhour=group(hour)
egen byte grarea=group(area_rank)
drop solicitor area* hour 

char treatment[omit] "Nw"
replace amt_donate=0 if amt==.

* Save final file
*save CharityOutputQJEFinal,replace

foreach x in "W" "Oo" {
	gen d`x'=treatment=="`x'"
	gen d`x'Ecu=treatment=="`x'" & charity=="Ecu"
	gen d`x'Lar=treatment=="`x'" & charity=="LaRabida"
	}	
gen dEcu=charity=="Ecu"
foreach x in "0d5m" "0d10m" "5d5m" "5d10m" "10d10m" "10d5m" {
	foreach y in "Nw" "W" "Oo" {
		replace treatment="`y'`x'" if treatment=="`y'-`x'"
		}
	}
foreach x in "Nw0d10m" "W0d10m" "W0d5m" "W10d10m" {
	gen d`x'08=(treatment=="`x'" & year==2008)
	}
foreach x in "Nw0d5m" "Nw5d5m" "W0d10m" "W0d5m" "W10d5m" "W5d5m" "Oo0d5m" "Oo5d5m" {
	gen d`x'09=treatment=="`x'" & year==2009
	}
order dNw0d10m08 dW0d10m08 dW0d5m08 dW10d10m08
order dNw5d5m09 dW0d10m09 dW0d5m09 dW5d5m09 dW10d5m09 dOo0d5m09 dOo5d5m09
* drop omitted category
drop dNw0d10m08 dNw0d5m09
compress

**************************************************
*** Estimating parameters for persuation rates ***

local charity_list = "LaRabida Ecu"

foreach charity in  `charity_list' {

xi: reg answer dW dOo i.grsol i.grdatloc i.grhour i.grarea if charity=="`charity'", robust cluster(sodate)
matrix reg_b = e(b)
predict ahat if charity=="`charity'"
scalar coeff_b1 = reg_b[1,1]
scalar coeff_b2 = reg_b[1,2]
scalar list coeff_b1 coeff_b2
gen ahat_Nw = ahat - coeff_b1*dW - coeff_b2*dOo
gen ahat_W  = ahat_Nw + coeff_b1
gen ahat_Oo = ahat_Nw + coeff_b2
replace ahat_Nw = min(max(ahat_Nw,0),1)
replace ahat_W = min(max(ahat_W,0),1)
replace ahat_Oo = min(max(ahat_Oo,0),1)

xi: reg saidyes dW dOo i.grsol i.grdatloc i.grhour i.grarea if charity=="`charity'", robust cluster(sodate)
matrix reg_b = e(b)
predict yhat if charity=="`charity'" 
scalar coeff_b1 = reg_b[1,1]
scalar coeff_b2 = reg_b[1,2]
scalar list coeff_b1 coeff_b2
gen yhat_Nw = yhat - coeff_b1*dW - coeff_b2*dOo
gen yhat_W  = yhat_Nw + coeff_b1
gen yhat_Oo = yhat_Nw + coeff_b2
replace yhat_Nw = min(max(yhat_Nw,0),1)
replace yhat_W = min(max(yhat_W,0),1)
replace yhat_Oo = min(max(yhat_Oo,0),1)

gen lb_Nw = yhat_Nw
gen lb_W  = yhat_W
gen lb_Oo = yhat_Oo
  
gen ub_Nw = yhat_Nw + 1 - ahat_Nw
gen ub_W  = yhat_W  + 1 - ahat_W
gen ub_Oo = yhat_Oo + 1 - ahat_Oo

matrix results_`charity' = J(3,5,.)
  
mean yhat_Nw ahat_Nw lb_Nw ub_Nw if charity=="`charity'" 
matrix mean_tmp=e(b)
matrix late_tmp=mean_tmp[1,1]/mean_tmp[1,2] 
matrix results_`charity'[1,1] = mean_tmp[1,1]
matrix results_`charity'[1,2] = mean_tmp[1,2]
matrix results_`charity'[1,3] = mean_tmp[1,3]
matrix results_`charity'[1,4] = mean_tmp[1,4]
matrix results_`charity'[1,5] = late_tmp

mean yhat_W ahat_W lb_W ub_W       if charity=="`charity'"  
matrix mean_tmp=e(b)
matrix late_tmp=mean_tmp[1,1]/mean_tmp[1,2] 
matrix results_`charity'[2,1] = mean_tmp[1,1]
matrix results_`charity'[2,2] = mean_tmp[1,2]
matrix results_`charity'[2,3] = mean_tmp[1,3]
matrix results_`charity'[2,4] = mean_tmp[1,4]
matrix results_`charity'[2,5] = late_tmp

mean yhat_Oo ahat_Oo lb_Oo ub_Oo   if charity=="`charity'"   
matrix mean_tmp=e(b)
matrix late_tmp=mean_tmp[1,1]/mean_tmp[1,2] 
matrix results_`charity'[3,1] = mean_tmp[1,1]
matrix results_`charity'[3,2] = mean_tmp[1,2]
matrix results_`charity'[3,3] = mean_tmp[1,3]
matrix results_`charity'[3,4] = mean_tmp[1,4]
matrix results_`charity'[3,5] = late_tmp

drop ahat ahat_* yhat yhat_* lb_* ub_*

}

*******************************************
****  Printing result in a LaTeX file *****
*******************************************
						  
matrix results = 100*(results_LaRabida \ results_Ecu)

frmttable using "$Persuasion/results/tableE2", statmat(results) tex replace sdec(1) ///
	ctitles("", P(Y=1|Z=1), e(1), APR (LB), APR (UB), LPR )  ///
	rtitles("La Rabida: Baseline" \ "La Rabida: Flyer" \ "La Rabida: Opt-Out" \ ///
			"Ecu: Baseline" \ "Ecu: Flyer" \ "Ecu: Opt-Out")  ///
	title("Persuasive Effect by Treatment in DLM")
	

