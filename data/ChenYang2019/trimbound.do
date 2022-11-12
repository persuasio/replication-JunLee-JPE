# delimit ;

set more 1;


capture program drop trimbound;
program define trimbound; /* syntax: depvar treatmentvar wgt sortvar */;
			  /* depvar is the outcome variable */;
			  /* treatmentvar is the binary treatment variable */;
			  /* wgt is the weight */;
			  /* sortvar is the variable that will break "ties" if the outcome variable is discrete */;
preserve;     /* preserve all the data; will restore later */;
	local y="`1'"; /* local macros */;
	local d="`2'";
	local wgt="`3'";
	local sv="`4'";
	tempvar s swgt upper lower oldd; /* selection variable, and sum of weights, upper and lower*/;

	keep if `d'~=.; /* delete all observations that don't have a treatment status */;

	gen `s'=0;  /* generate a variable for selection */;
	replace `s'=1 if `y'~=.; /* s=1 for nonmissing outcomes */;
	
	sum `y' [aw=`wgt']; /* summarize depvar */;
	sum `s' [aw=`wgt']; /* summarize selection variable */;
	sum `d' [aw=`wgt']; /* summarize treatment variable */;
	local ED=r(mean);

	quietly reg `s' `d' [aw=`wgt']; /* first determine which group has a higher fraction non-missing */;

	if _b[`d']>=0 {;  /* if the treatment has a higher fraction of non-missing */;

		local p=_b[`d']/(_b[`d']+_b[_cons]); /* trimming fraction p */;
		local alpha=_b[_cons]; /* alpha coefficient */;

		quietly sum `y' if `d'==1 [aw=`wgt']; /* summarize depvar if d equals 1 */;
		local sumw=r(sum_w); /* macro that has sum of weights */;

		sort `d' `y' `sv'; /*sort the data.... nonmissing lhrwage is first */;
		quietly by `d': gen `swgt'=sum(`wgt'); /* generate running sum of weights */;
		quietly by `d': gen `upper'=1 if `swgt'>`p'*`sumw'; /* generate variable to indicate upper distribution */;
		quietly by `d': gen `lower'=1 if `swgt'<(1-`p')*`sumw'; /* generate variable to indicate lower distribution */;

		
		quietly su `s' if `d'==0 [aw=`wgt']; /* for control */;
		local controlN=r(N);
		quietly su `s' if `d'==1 [aw=`wgt']; /* and for treatments */;
		local treatN=r(N); /* local macro for number of treatment observations */;

		quietly su `y' if `d'==0 [aw=`wgt']; local yc=r(mean); local ycvar=r(Var)/r(N); local UcontrolN=r(N); /* sum for the untrimmed control */;
		quietly su `y' if `d'==1 [aw=`wgt']; local yt=r(mean); local ytvar=r(Var)/r(N); local UtreatN=r(N); /* sum for the untrimmed treatment */;
		quietly su `y' if `d'==1 & `upper'==1 [aw=`wgt']; local uyt=r(mean); local uytN=r(N) /* sum for the UB treatment */;
			local  ucomp1=r(Var)/r(N); /* first component of variance */;
			local  ucomp2=((r(min)-r(mean))^2)*(`p')/((1-`p')*`UtreatN');  /* second component of variance */;
			local  alphat=`alpha'/(1-`p');
			local  pvar=((1-`p')^2)*(((1-`alphat')/((`treatN')*(`alphat')))+((1-`alpha')/((`controlN')*(`alpha'))));	
			local  ucomp3=(((r(min)-r(mean))^2)/((1-`p')^2))*`pvar'; /* third component*/;
			local ubse=((`ucomp1')+(`ucomp2')+(`ucomp3'))^.5; /* UB std. error */;
			local lcutoff=r(min); /* lower cutoff */;
		quietly su `y' if `d'==1 & `lower'==1 [aw=`wgt']; local lyt=r(mean); local lytN=r(N) /* sum for the LB treatment */;
			local  lcomp1=r(Var)/r(N); /* first component of variance */;
			local  lcomp2=((r(max)-r(mean))^2)*`p'/((1-`p')*(`UtreatN'));  /* second component of variance */;
			local  lcomp3=(((r(max)-r(mean))^2)/((1-`p')^2))*`pvar'; /* third component*/;
			local lbse=((`lcomp1')+(`lcomp2')+(`lcomp3'))^.5; /* LB std. error */;
			local ucutoff=r(max); /* upper cutoff */;

		display "Treatment group must be trimmed";
		display "  Total Number of Control Observations: "; count if `d'==0;
		display "Total Number of Treatment Observations: "; count if `d'==1;
		display "Total Number of non-missing Treatment: " `UtreatN';
		display "         Fraction Non-missing, Control: " `alpha';
		display "       Fraction Non-missing, Treatment: " `alpha'/(1-`p');
		display "        Trimming Proportion (Std. Err): " `p' " (" `pvar'^.5 ")";
		display "                          Lower cutoff: " `lcutoff';
		display "                          Upper cutoff: " `ucutoff';
		display "     Untrimmed Control Mean (Std. Err): " `yc' " (" `ycvar'^.5 ")";
		display "   Untrimmed Treatment Mean (Std. Err): " `yt' " (" `ytvar'^.5 ")";
		display "           Untrimmed Effect (Std. Err): " `yt'-`yc' " (" sqrt(`ytvar'+`ycvar') ")";
		display "Upper bound Treatment (Std. Err) [Obs]: " `uyt' " (" `ubse' ") [" `uytN' "]";
		display "                           Component 1: " `ucomp1'^.5;
		display "                           Component 2: " `ucomp2'^.5;
		display "                           Component 3: " `ucomp3'^.5;
		display "Lower bound Treatment (Std. Err) [Obs]: " `lyt' " (" `lbse' ") [" `lytN' "]";
		display "                           Component 1: " `lcomp1'^.5;
		display "                           Component 2: " `lcomp2'^.5;
		display "                           Component 3: " `lcomp3'^.5;
		display "         Upper bound Effect (Std. Err): " `uyt'-`yc' " (" ((`ycvar')+(`ubse')^2)^.5 ")";
		display "         Lower bound Effect (Std. Err): " `lyt'-`yc' " (" ((`ycvar')+(`lbse')^2)^.5 ")";


	};


	if _b[`d']<0 {;  /* if the control has a higher fraction of non-missing */;

		replace `d'=1-`d';		

		quietly reg `s' `d' [aw=`wgt']; /* first determine which group has a higher fraction non-missing */;

		local p=_b[`d']/(_b[`d']+_b[_cons]); /* trimming fraction p */;
		local alpha=_b[_cons]; /* alpha coefficient */;

		quietly sum `y' if `d'==1 [aw=`wgt']; /* summarize depvar if d equals 1 */;
		local sumw=r(sum_w); /* macro that has sum of weights */;

		sort `d' `y' `sv'; /*sort the data.... nonmissing lhrwage is first */;
		quietly by `d': gen `swgt'=sum(`wgt'); /* generate running sum of weights */;
		quietly by `d': gen `upper'=1 if `swgt'>`p'*`sumw'; /* generate variable to indicate upper distribution */;
		quietly by `d': gen `lower'=1 if `swgt'<(1-`p')*`sumw'; /* generate variable to indicate lower distribution */;

		
		quietly su `s' if `d'==0 [aw=`wgt']; /* for control */;
		local controlN=r(N);
		quietly su `s' if `d'==1 [aw=`wgt']; /* and for treatments */;
		local treatN=r(N); /* local macro for number of treatment observations */;

		quietly su `y' if `d'==0 [aw=`wgt']; local yc=r(mean); local ycvar=r(Var)/r(N); local UcontrolN=r(N); /* sum for the untrimmed control */;
		quietly su `y' if `d'==1 [aw=`wgt']; local yt=r(mean); local ytvar=r(Var)/r(N); local UtreatN=r(N); /* sum for the untrimmed treatment */;
		quietly su `y' if `d'==1 & `upper'==1 [aw=`wgt']; local uyt=r(mean); local uytN=r(N) /* sum for the UB treatment */;
			local  ucomp1=r(Var)/r(N); /* first component of variance */;
			local  ucomp2=((r(min)-r(mean))^2)*(`p')/((1-`p')*`UtreatN');  /* second component of variance */;
			local  alphat=`alpha'/(1-`p');
			local  pvar=((1-`p')^2)*(((1-`alphat')/((`treatN')*(`alphat')))+((1-`alpha')/((`controlN')*(`alpha'))));	
			local  ucomp3=(((r(min)-r(mean))^2)/((1-`p')^2))*`pvar'; /* third component*/;
			local ubse=((`ucomp1')+(`ucomp2')+(`ucomp3'))^.5; /* UB std. error */;
			local lcutoff=r(min); /* lower cutoff */;
		quietly su `y' if `d'==1 & `lower'==1 [aw=`wgt']; local lyt=r(mean); local lytN=r(N) /* sum for the LB treatment */;
			local  lcomp1=r(Var)/r(N); /* first component of variance */;
			local  lcomp2=((r(max)-r(mean))^2)*`p'/((1-`p')*(`UtreatN'));  /* second component of variance */;
			local  lcomp3=(((r(max)-r(mean))^2)/((1-`p')^2))*`pvar'; /* third component*/;
			local lbse=((`lcomp1')+(`lcomp2')+(`lcomp3'))^.5; /* LB std. error */;
			local ucutoff=r(max); /* upper cutoff */;



		display "Control group must be trimmed";
		display "  Total Number of Control Observations: "; count if `d'==1;
		display "Total Number of Treatment Observations: "; count if `d'==0;
		display "Total Number of non-missing Treatment: " `UtreatN';
		display "         Fraction Non-missing, Control: " `alpha'/(1-`p');
		display "       Fraction Non-missing, Treatment: " `alpha';
		display "        Trimming Proportion (Std. Err): " `p' " (" `pvar'^.5 ")";
		display "                          Lower cutoff: " `lcutoff';
		display "                          Upper cutoff: " `ucutoff';
		display "     Untrimmed Control Mean (Std. Err): " `yt' " (" `ytvar'^.5 ")";
		display "   Untrimmed Treatment Mean (Std. Err): " `yc' " (" `ycvar'^.5 ")";
		display "           Untrimmed Effect (Std. Err): " `yc'-`yt' " (" sqrt(`ytvar'+`ycvar') ")";
		display "  Upper bound Control (Std. Err) [Obs]: " `uyt' " (" `ubse' ") [" `uytN' "]";
		display "                           Component 1: " `ucomp1'^.5;
		display "                           Component 2: " `ucomp2'^.5;
		display "                           Component 3: " `ucomp3'^.5;
		display "  Lower bound Control (Std. Err) [Obs]: " `lyt' " (" `lbse' ") [" `lytN' "]";
		display "                           Component 1: " `lcomp1'^.5;
		display "                           Component 2: " `lcomp2'^.5;
		display "                           Component 3: " `lcomp3'^.5;
		display "         Upper bound Effect (Std. Err): " `yc' - `lyt' " (" ((`ycvar')+(`lbse')^2)^.5 ")";
		display "         Lower bound Effect (Std. Err): " `yc' - `uyt' " (" ((`ycvar')+(`ubse')^2)^.5 ")";


	};


end;





