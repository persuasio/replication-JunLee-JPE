************
* SCRIPT:   tableE1.do
* PURPOSE:  Creates Table E1
*
* ACKNOWLEDGMENT
*	The summary statistic used here is from Landry, Lange, List, Price, and Rupp (QJE, 2006).
*   Specifically, the first three rows of Table II in their paper provide information to compute
*   the following quantities:
*   - P(Y=1|Z=1) is obtained by (# of households that contributed)/(Total households approached); 
*   - e(1) = P(Z=1) is obtained by (Total households home)/(Total households approached).  
************

* input for table E1
matrix landry_et_al = ///
  ( 1186, 446, 113 \ ///
	1282, 453, 67 \ ///
	963, 363, 165 \ ///
	1402, 493, 177)

matrix U = J(rowsof(landry_et_al),1,1)
matrix All = U'*landry_et_al
matrix landry_et_al = ( landry_et_al \ All )	
	
matrix results = J(1,5,.)
	
foreach j of numlist 1/5 {	

	display `j'

	scalar pr_y1_z1 = landry_et_al[`j',3]/landry_et_al[`j',1]
	scalar pr_z1 = landry_et_al[`j',2]/landry_et_al[`j',1]
	scalar theta_lb = pr_y1_z1
	scalar theta_ub = pr_y1_z1 + 1 - pr_z1
	scalar theta_local = pr_y1_z1/pr_z1

	matrix results = results \ (pr_y1_z1, pr_z1, theta_lb, theta_ub, theta_local) 

}

matrix results = 100*results[2..6,1..5]

frmttable using "$Persuasion/results/tableE1", statmat(results) tex replace sdec(1) ///
	ctitles("", P(Y=1|Z=1), e(1), APR (LB), APR (UB), LPR )  ///
	rtitles("VCM" \ "VCM with seed money" \ ///
	"Single-prize lottery" \ "Multiple-prize lottery" \    ///
	"All")  ///
	title("Persuasive Effect by Treatment in Landry et al. (2006)")
	
