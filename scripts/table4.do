************
* SCRIPT:   table4.do
* PURPOSE:  Creates Table 4
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

gen   case_id = _n           /* generate Case ID */
gen outcome = voteddem_all   /* Y_i */
gen   treat = readsome       /* T_i */
gen   instr = post           /* Z_i */  

*** generate variables used in estimating persuation rates ***

gen a_u = outcome*treat + (1-treat)
gen b_l = outcome*(1-treat)

gen a_u2 = outcome + (1-treat)
gen b_l2 = outcome - treat

*** specify a variety of critical values for different purposes ***

scalar alpha_level = 0.2     /* the (1-alpha_level) confidence interval */  
scalar cv_cns1 = invnormal(1-alpha_level)   /* one-sided critical value */
scalar cv_cns2 = invnormal(1-alpha_level/2) /* two-sided critical value */

scalar pretest_level = 0.001
scalar cv_cns_pretest = invnormal(1-pretest_level/4) 
/* note that pretest_level is divided by 2 since we pretest two hypotheses */

scalar adjusted = alpha_level-pretest_level
scalar cv_cns1_adj = invnormal(1-adjusted)   /* one-sided critical value after pretesting */
scalar cv_cns2_adj = invnormal(1-adjusted/2) /* two-sided critical value after pretesting */

*************************
****  ITT Analysis  *****
*************************

/* ITT estimate on outcome */

reg outcome instr, robust

scalar itt_coef = _b[instr]
scalar itt_se = _se[instr]

scalar itt_lb = max(0,itt_coef - cv_cns2*itt_se)
scalar itt_ub = min(1,itt_coef + cv_cns2*itt_se)
scalar list itt_lb itt_coef itt_ub

*** saving results ***
matrix results =  (itt_lb, itt_coef, itt_coef, itt_ub)

/* ITT estimate on treatment */

reg treat instr, robust

scalar itt_coef_treat = _b[instr]
scalar itt_se_treat = _se[instr]

scalar itt_lb = max(0,itt_coef_treat - cv_cns2*itt_se_treat)
scalar itt_ub = min(1,itt_coef_treat + cv_cns2*itt_se_treat)
scalar list itt_lb itt_coef itt_ub

*** saving results ***
matrix results = results \ (itt_lb, itt_coef_treat, itt_coef_treat, itt_ub)
 

**********************************************************************
****  Average persuation rate when (Y_i, T_i, Z_i) are observed  *****
**********************************************************************

reg outcome instr, robust

/* estimate of the lower bound for average persuation rate */ 
nlcom _b[instr]/(1-_b[_cons])	 

matrix  lower_bound_est = r(b)
matrix lower_bound_avar = r(V)
scalar lower_bound_coef = lower_bound_est[1,1]
scalar   lower_bound_se = sqrt(lower_bound_avar[1,1])

reg a_u instr 
est store a_u_reg
reg b_l instr 
est store b_l_reg
suest a_u_reg b_l_reg, vce(cluster case_id) 

/* estimate of the uppe bound for average persuation rate */  
nlcom ([a_u_reg_mean]_cons + [a_u_reg_mean]instr - [b_l_reg_mean]_cons)/(1-[b_l_reg_mean]_cons)

matrix  upper_bound_est = r(b)
matrix upper_bound_avar = r(V)
scalar upper_bound_coef = upper_bound_est[1,1]
scalar   upper_bound_se = sqrt(upper_bound_avar[1,1])

/* compute the critical value using Stoye (2009) */  
scalar correction_term = (upper_bound_coef-lower_bound_coef)/max(upper_bound_se,lower_bound_se)
egen cvtmp = fill("0 0.01")
replace cvtmp = . if cvtmp > cv_cns2+0.01
gen difftmp = abs(normal(cvtmp + correction_term) - normal(-cvtmp) - (1-alpha_level))
replace cvtmp = . if difftmp > 0.002
sum cvtmp
scalar cv_cns_stoye = r(mean)
drop cvtmp difftmp

scalar lower_bound_end = max(0,lower_bound_coef - cv_cns_stoye*lower_bound_se)
scalar upper_bound_end = min(1,upper_bound_coef + cv_cns_stoye*upper_bound_se)

scalar list lower_bound_end lower_bound_coef upper_bound_coef upper_bound_end


*** saving results ***
matrix results = results \ (lower_bound_end, lower_bound_coef, upper_bound_coef, upper_bound_end) 


****************************************************************************************
****  Average persuation rate when (Y_i,T_i) and (Y_i,Z_i) are separately observed *****
****************************************************************************************

reg a_u2 instr, robust 

lincom _cons + instr
scalar coef_a_u2 = r(estimate)
scalar   se_a_u2 = r(se)

reg b_l2 instr, robust 
scalar coef_b_l2 = _b[_cons]
matrix      avar = e(V) 
scalar   se_b_l2 = sqrt(el(avar,2,2))

scalar upper_bound2_coef = (min(1,coef_a_u2) - max(0,coef_b_l2))/(1-max(0,coef_b_l2))


*** pre-testing stage ***

scalar a_u2_pretest = coef_a_u2 + cv_cns_pretest*se_a_u2
scalar b_l2_pretest = coef_b_l2 + cv_cns_pretest*se_b_l2

*** N.B. It turns out that (a_u2_pretest < 1) & (b_l2_pretest <= 0) ***
 
scalar upper_bound2_se = se_a_u2

/* compute the critical value using Stoye (2009) again */  
scalar correction_term = (upper_bound2_coef-lower_bound_coef)/max(upper_bound2_se,lower_bound_se)
egen cvtmp = fill("0 0.01")
replace cvtmp = . if cvtmp > cv_cns2_adj+0.01
gen difftmp = abs(normal(cvtmp + correction_term) - normal(-cvtmp) - (1-adjusted))
replace cvtmp = . if difftmp > 0.002
sum cvtmp
scalar cv_cns2_stoye = r(mean)

scalar lower_bound2_end = max(0,lower_bound_coef - cv_cns2_stoye*lower_bound_se)
scalar upper_bound2_end = min(1,upper_bound2_coef + cv_cns2_stoye*upper_bound2_se)

scalar list lower_bound2_end lower_bound_coef upper_bound2_coef upper_bound2_end


*** saving results ***
matrix results = results \ (lower_bound2_end, lower_bound_coef, upper_bound2_coef, upper_bound2_end) 


****************************************************************
****  Average persuation rate when (Y_i,Z_i) are observed  *****
****************************************************************

scalar lower_bound3_end = max(0,lower_bound_coef - cv_cns1*lower_bound_se)
scalar list lower_bound3_end lower_bound_coef


*** saving results ***
matrix results = results \ (lower_bound3_end, lower_bound_coef, 1, 1) 


**********************************************************************
****  Local persuation rate when (Y_i, T_i, Z_i) are observed    *****
**********************************************************************

reg outcome instr 
est store num_reg

gen den_lpr = (1-outcome)*(1-treat)
reg den_lpr instr 
est store den_reg

suest num_reg den_reg, vce(cluster case_id) 
nlcom [num_reg_mean]instr/(-[den_reg_mean]instr)

matrix  lpr_est = r(b)
matrix lpr_avar = r(V)
scalar lpr_coef = lpr_est[1,1]
scalar   lpr_se = sqrt(lpr_avar[1,1])


scalar lpr_lb = max(0,lpr_coef - cv_cns2*lpr_se)
scalar lpr_ub = min(1,lpr_coef + cv_cns2*lpr_se)
scalar list lpr_lb lpr_coef lpr_ub

*** saving results ***
matrix results = results \ (lpr_lb, lpr_coef,  lpr_coef, lpr_ub) 

****************************************************************************************
****  Local persuation rate when (Y_i,T_i) and (Y_i,Z_i) are separately observed   *****
****************************************************************************************

** theta^*_L = max(LATE,theta_L) 

ivregress 2sls outcome (treat = instr), robust 
scalar lpr2_coef = _b[treat]
matrix      avar = e(V) 
scalar   lpr2_se = sqrt(el(avar,1,1))


scalar lpr_lb_LATE = max(0,lpr2_coef - cv_cns2*lpr2_se)
scalar lpr_lb_LBND = max(0,lower_bound_coef - cv_cns2*lower_bound_se)

scalar lpr2_max = max(lpr2_coef, lower_bound_coef)
scalar lpr2_lb = max(lpr_lb_LATE, lpr_lb_LBND)

scalar list lpr2_lb lpr2_max 

*** saving results ***
matrix results = results \ (lpr2_lb, lpr2_max, 1, 1) 

****************************************************************
****  Local persuation rate when (Y_i,Z_i) are observed    *****
****************************************************************

scalar lpr3_lb = lower_bound3_end 
scalar list lpr3_lb lower_bound_coef

*** saving results ***
matrix results = results \ (lpr3_lb, lower_bound_coef, 1, 1) 

*******************
****  LATE    *****
*******************

scalar lpr_ub_LATE = min(1,lpr2_coef + cv_cns2*lpr2_se)
scalar list lpr_lb_LATE lpr2_coef lpr_ub_LATE

*** saving results ***
matrix results = results \ (lpr_lb_LATE, lpr2_coef,  lpr2_coef, lpr_ub_LATE) 

************************************************
****  LATE  when (Y_i,Z_i) are observed    *****
************************************************

scalar lpr3_lb_LATE = max(0,itt_coef - cv_cns1*itt_se)
scalar list lpr3_lb_LATE itt_coef

*** saving results ***
matrix results = results \ (lpr3_lb_LATE, itt_coef,  1, 1) 

**********************************************************************
****  Average persuation rate                                    *****
****	conditional on voting without persuasive treatement      *****
****	when (Y_i, T_i, Z_i) are observed                        *****
**********************************************************************


gen yt00 = (1-voted)*(1-treat)

reg outcome instr
est store multi_num_reg

reg yt00 instr
est store multi_den_reg

suest multi_num_reg multi_den_reg, vce(cluster case_id) 
*** estimating p_1(1|1) + p_0(1,0|0)
lincom [multi_num_reg_mean]instr+[multi_num_reg_mean]_cons+[multi_den_reg_mean]_cons
	
scalar multi_pretest_coef = r(estimate)
scalar multi_pretest_se = r(se)	
	
*** pre-testing stage ***
scalar multi_pretest = multi_pretest_coef + cv_cns_pretest*multi_pretest_se


*** estimating the lower bound	
suest multi_num_reg multi_den_reg, vce(cluster case_id) 
nlcom [multi_num_reg_mean]instr/(1-[multi_num_reg_mean]_cons-[multi_den_reg_mean]_cons)		

matrix  multi_lower_bound_est = r(b)
matrix multi_lower_bound_avar = r(V)
scalar multi_lower_bound_coef = multi_lower_bound_est[1,1]
scalar   multi_lower_bound_se = sqrt(lower_bound_avar[1,1])


scalar  multi_lower_bound_end = max(0,multi_lower_bound_coef - cv_cns1_adj*multi_lower_bound_se)
scalar list multi_lower_bound_end multi_lower_bound_coef


*** saving results ***
matrix results = results \ (multi_lower_bound_end, multi_lower_bound_coef, 1, 1) 


*******************************************
****  Printing result in a LaTeX file *****
*******************************************
						  
matrix results_est = (results[1..11,2], results[1..11,3], results[1..11,1], results[1..11,4])

frmttable using "$Persuasion/results/table4", statmat(results_est) tex replace sdec(4) ///
	title("Table 4") ///
	ctitles("", Estimate (LB), Estimate (UB), CI (LB), CI (UB) )  ///
	rtitles("Intent-to-Treat (Voting)" \ "Intent-to-Treat (Reading newspaper)" \ ///
	"Average Persuasion Rate (Case 1)" \ "Average Persuasion Rate (Case 2)" \ "Average Persuasion Rate (Case 3)" \    ///
	"Local Persuasion Rate (Case 1)" \ "Local Persuasion Rate (Case 2)" \ "Local Persuasion Rate (Case 3)" \ /// 
	"LATE (Cases 1 and 2)" \ "LATE (Case 3)" \ "Non-binary case") 
