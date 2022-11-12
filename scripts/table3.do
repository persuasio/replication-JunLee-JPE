************
* SCRIPT:   table3.do
* PURPOSE:  Creates Table 3
*
* ACKNOWLEDGMENT
*	The orginal dataset "publicdata.dta" is from Gerber, Karlan, and Bergan (2009, AEJ Applied).
************

use "$Persuasion/data/GerberKarlanBergan2009/publicdata.dta", clear

**************************************
****       Setup                 *****
**************************************

drop if times==1	/* drop observations with The Washing Times treatment */


gen data_avail = 0	/* criterion to be included in the estimation sample */
replace data_avail = 1 if survey == 1 & !missing(voteddem_all) & !missing(readsome) 

keep if data_avail==1         /* keep observations in the estimation sample */

*** two way tables by treatment ***
tab voteddem_all readsome if post == 1, matcell(treat)

frmttable using "$Persuasion/results/table3", statmat(treat) tex replace sdec(0) ///
	title("Table 3: The Washington Post (Z = 1)") ///
	ctitles("", "T=0", "T=1" )  ///
	rtitles("Y=0" \ "Y=1" )

tab voteddem_all readsome if post == 0, matcell(control)

frmttable using "$Persuasion/results/table3", statmat(control) tex replace sdec(0) ///
	title("Table 3: Control (Z = 0)") ///
	ctitles("", "T=0", "T=1" )  ///
	rtitles("Y=0" \ "Y=1" )	
