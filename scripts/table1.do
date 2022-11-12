************
* SCRIPT:   table1.do
* PURPOSE:  Creates Table 1
************

* input for table 1
matrix turnout = ///
  ( 0.472, 0.448, 0.279, 0 \ ///
	0.310, 0.286, 0.293, 0 \ ///
	0.711, 0.660, 0.737, 0 \ ///
	0.416, 0.405, 0.414, 0 \ ///
	0.455, 0.435, 0.800, 0 \ ///
	0.700, 0.690, 0.250, 0 )
	
matrix results = J(1,5,.)
	
foreach j of numlist 1/6 {	

	display `j'

	scalar y1 = turnout[`j',1]
	scalar y0 = turnout[`j',2]
	scalar e1 = turnout[`j',3]
	scalar e0 = turnout[`j',4]

	calc4persuasio y1 y0 e1 e0
	
	scalar persuasion_DK = (y1-y0)/((e1-e0)*(1-y0))

	matrix results = results \ (persuasion_DK, r(apr_lb), r(apr_ub), r(lpr_lb), r(lpr_ub)) 

}

matrix results = 100*results[2..7,1..5]

frmttable using "$Persuasion/results/table1", statmat(results) tex replace sdec(1) ///
	ctitles("", DK, APR (LB), APR (UB), LPR (LB), LPR (UB) )  ///
	rtitles("Green and Gerber (2000)" \ "Green, Gerber, and Nickerson (2003)" \ ///
	"Green and Gerber (2001)" \ "Green and Gerber (2001)" \    ///
	"Gentzkow (2006)" \ "Gentzkow, Shapiro, and Sinkinson (2011)")  ///
	title("Persuasion Rates: Papers on Voter Turnout")
