************
* SCRIPT:   table2.do
* PURPOSE:  Creates Table 2
************

	use "$Persuasion/data/ChenYang2019/ChenYang2019.dta", clear

	matrix results = J(1, 5, .)
	qui reg active_user treatment_vpnonly treatment_nlonly treatment_vpnnl, r
	local share_compliers = _b[treatment_vpnnl]
			
	// persuasion rates for panel variables
	foreach Y in info_foreign_website ///	A.1.2 Ranked high: foreign websites
				 info_freq_website_for ///	A.1.6 Freq. of visiting foreign websites for info.
				 az_belief_media_value ///	A.3 Valuation of access to foreign media outlets
				 az_belief_media_trust ///	A.4 Trust in non-domestic media outlets
				 bias_domestic ///			A.5.1 Degree of censorship on domestic news outlets	
				 bias_foreign ///	A.5.2 Degree of censorship on foreign news outlets	
				 az_belief_media_justif /// A.6 Censorship unjustified
				 bias_dom_govt_policy_t1 ///  A.7.1 Domestic cens. driven by govt. policies
				 bias_for_govt_policy_t1 /// A.7.2 Foreign cens. driven by govt. policies
		{
				
		qui reg `Y'_w3_p treatment_vpnonly treatment_nlonly treatment_vpnnl, r
		qui summ `Y'_w1_p if treatment_vpnnl == 1
		local persuasion_DK = _b[treatment_vpnnl]/(`share_compliers'*(1-`r(mean)'))
				
		persuasio apr `Y'_w3_p active_user treatment_vpnnl ///
					if (treatment_main == 1) | (treatment_main == 3), ///
					method("bootstrap") nboot($nbt)

		matrix pers_res = (e(apr_est) \ e(apr_ci))
		
		persuasio lpr `Y'_w3_p active_user treatment_vpnnl ///
					if (treatment_main == 1) | (treatment_main == 3), ///
					method("bootstrap") nboot($nbt)
					
		matrix CY = (`persuasion_DK' \ .)			
		matrix pers_res = (pers_res, ( (e(lpr_est), .) \ e(lpr_ci)))
		matrix results = results \ ( CY , pers_res)
				
		}

	// persuasion rates for non-panel variables
	foreach Y in ///	           
	vpn_purchase_wmt_record ///	  A.2.1 Purchase discounted tool we offered
	vpn_purchase_yes ///		  A.2.2 Purchase any tool
	{
				
		qui reg `Y'_p treatment_vpnonly treatment_nlonly treatment_vpnnl, r
		local persuasion_DK = _b[treatment_vpnnl]/(`share_compliers'*(1-_b[_cons]))
				
		persuasio apr `Y'_p active_user treatment_vpnnl ///
					if (treatment_main == 1) | (treatment_main == 3), ///
					method("bootstrap") nboot($nbt)

		matrix pers_res = (e(apr_est) \ e(apr_ci))
		
		persuasio lpr `Y'_p active_user treatment_vpnnl ///
					if (treatment_main == 1) | (treatment_main == 3), ///
					method("bootstrap") nboot($nbt)
					
		matrix CY = (`persuasion_DK' \ .)			
		matrix pers_res = (pers_res, ( (e(lpr_est), .) \ e(lpr_ci)))
		matrix results = results \ ( CY , pers_res)
				
	}
	
	matrix results = results[2..23,1..5] 
	matrix results = results[1..4, 1..5] \ results[19..22, 1..5] \ results[5..18, 1..5]
	
	matrix results = 100*results
	
	frmttable using "$Persuasion/results/table2", statmat(results) tex replace sdec(1) ///
	title("Table 2") ///
	ctitles("CY19", APR (LB), APR (UB), LPR, "")
