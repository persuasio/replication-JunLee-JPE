************
* SCRIPT:   CY19_data.do
* PURPOSE:  Creates a dataset for Table 2
*
* ACKNOWLEDGMENT
*	This script is a slightly modified version of the do file
*   "master_dofile_panelsurvey.do" included in the replication files
*   for Chen and Yang (2019, AER). 
*   Before running this code, it is necessary to download the original dataset from 
*   the AER webpage at https://www.aeaweb.org/articles?id=10.1257/aer.20171765
*   and store them at "$Persuasion/data/ChenYang2019".
************

*** Do-file: replication of Chen and Yang (2018)
*** Part 1: based on student survey data
*** December 2018

clear all
set more off

*** Define root name for data folder
	cd "$Persuasion/data/ChenYang2019"
		
*** Program Anderson z-score generation

	cap program drop andersonz 
	program andersonz
	syntax varlist [if] [in], GENerate(name)

	* 	Convert to effect sizes
		local varliste ""
		foreach X of varlist `varlist' {
			qui sum `X' `if'
			gen `X'e = `X'/r(sd) `if'
			local csd = r(sd)
			qui sum `X'
			replace `X'e = `X'e - r(mean)/`csd'
			local varliste "`varliste' `X'e" 
			}
	
	*   Generate GLS weighting matrix
		matrix accum R = `varliste', nocons dev
		matrix R = R/r(N)
		matrix R = syminv(R)
		local counter1 = 1
		matrix J = J(colsof(R),1,1)
		while `counter1' <= colsof(R) {
			matrix T = R[`counter1',1..colsof(R)]
			matrix A = T*J
			local weight`counter1'f = A[1,1]
			local counter1 = `counter1' + 1
			}
	
	*   Generate total number of variables per obs
	*   As opposed to Anderson's code, do not replace missing with zeros
		cap drop sample
		gen sample = 0
		local counter1 = 1
		foreach X of varlist `varlist' {
			replace sample = sample + `weight`counter1'f' if `X'~=.
	*		replace `X'e = 0 if `X' ==.
			local counter1 = `counter1' + 1
			}
	
	*   Apply weights to outcomes
		quietly gen `generate' = 0 `if'
		
		local counter1 = 1
		foreach X of varlist `varlist' {
			replace `generate' = `X'e * (`weight`counter1'f') + `generate' `if'
			local counter1 = `counter1' + 1
			}
			
		replace `generate' = `generate' / sample `if'
		cap drop sample
		cap drop `varliste'
		
	* 	Re-standardize aggregate z-score
		egen z_`generate' = std(`generate')
		drop `generate'
		rename z_`generate' `generate'
	
	end	
	
	

*** Load merged raw data for student survey

	use panelsurvey_raw, clear
	
	

	
*** Define variables

	local var_identity_w1 					`"responseID_wave1"'
	local var_university 					`"university"'
	local var_treatment 					`"treatment_vpn treatment_newsletter vpn_current_paid_user treatment_vpnexpiration"'
	
	
* 	A: Media related beliefs, attitudes, and behaviors
	local var_info_ranking_reg_w1 			`"info_domestic_website_w1 info_foreign_website_w1 info_social_media_dom_w1 info_word_of_mouth_w1"'
	local var_info_ranking_reg_w2 			`"info_domestic_website_w2 info_foreign_website_w2 info_social_media_dom_w2 info_social_media_for_w2 info_word_of_mouth_w2"'
	local var_info_ranking_reg_w3 			`"info_domestic_website_w3 info_foreign_website_w3 info_social_media_dom_w3 info_social_media_for_w3 info_word_of_mouth_w3"'
	local var_info_ranking_fig_w23 			`"info_domestic_website info_foreign_website info_social_media_dom info_word_of_mouth info_social_media_for"'
	
	local var_info_freq_reg_w1 				`"info_freq_website_for_w1"'
	local var_info_freq_reg_w2 				`"info_freq_website_for_w2"'
	local var_info_freq_reg_w3 				`"info_freq_website_for_w3"'
	local var_info_freq_fig_w123 			`"info_freq_website_for"'
	
	local var_vpnusage_own 					`"vpn_usage vpn_neveruser vpn_payment_now vpn_stop_when vpn_stop_why vpn_stop_why_text"'
	
	local var_vpnusage_rm_all 				`"vpn_roommate_usage"'
	
	local var_vpn_purchase_w3 				`"vpn_purchase vpn_purchase_yes vpn_purchase_wmt vpn_purchase_premium"'
	local var_vpn_purchase_reg_w3 			`"vpn_purchase_wmt_record vpn_purchase_yes"'
		
	local var_media_valuation_reg_w1 		`"wtp_vpn_w1 added_value_foreign_media_w1"'
	local var_media_valuation_reg_w2 		`"wtp_vpn_w2 added_value_foreign_media_w2"'
	local var_media_valuation_reg_w3 		`"wtp_vpn_w3 added_value_foreign_media_w3"'
	local var_media_valuation_fig_w123 		`"wtp_vpn added_value_foreign_media az_belief_media_value"'
	
	local var_media_trust_reg_w1 			`"trust_media_dom_state_w1 trust_media_dom_private_w1 trust_media_foreign_w1"'
	local var_media_trust_reg_w2 			`"trust_media_dom_state_w2 trust_media_dom_private_w2 trust_media_foreign_w2"'
	local var_media_trust_reg_w3 			`"trust_media_dom_state_w3 trust_media_dom_private_w3 trust_media_foreign_w3"'
	local var_media_trust_fig_w123 			`"trust_media_dom_state trust_media_dom_private trust_media_foreign az_belief_media_trust"'
	
	local var_percmediabias_reg_di_cn_w2 	`"distneutral_cn_neg_cn_w2 distneutral_cn_pos_cn_w2 distneutral_us_neg_cn_w2 distneutral_us_pos_cn_w2"'
	local var_percmediabias_reg_di_us_w2 	`"distneutral_cn_neg_us_w2 distneutral_cn_pos_us_w2 distneutral_us_neg_us_w2 distneutral_us_pos_us_w2"'
	local var_percmediabias_reg_ce_cn_w2 	`"bias_cn_neg_cn_cens_w2 bias_cn_pos_cn_cens_w2 bias_us_neg_cn_cens_w2 bias_us_pos_cn_cens_w2"'
	local var_percmediabias_reg_ce_us_w2 	`"bias_cn_neg_us_cens_w2 bias_cn_pos_us_cens_w2 bias_us_neg_us_cens_w2 bias_us_pos_us_cens_w2"'
	local var_percmediabias_fig_w12 		`"distneutral_cn_neg_cn distneutral_cn_pos_cn distneutral_us_neg_cn distneutral_us_pos_cn distneutral_cn_neg_us distneutral_cn_pos_us distneutral_us_neg_us distneutral_us_pos_us bias_cn_neg_cn_cens bias_cn_pos_cn_cens bias_us_neg_cn_cens bias_us_pos_cn_cens bias_cn_neg_us_cens bias_cn_pos_us_cens bias_us_neg_us_cens bias_us_pos_us_cens"'
	
	local var_censor_justif_reg_w1 			`"censor_just_dom_economic_w1 censor_just_dom_political_w1 censor_just_dom_social_w1 censor_just_for_w1"'
	local var_censor_justif_reg_w2 			`"censor_just_dom_economic_w2 censor_just_dom_political_w2 censor_just_dom_social_w2 censor_just_for_w2 censor_just_porn_w2"'
	local var_censor_justif_reg_w3 			`"censor_just_dom_economic_w3 censor_just_dom_political_w3 censor_just_dom_social_w3 censor_just_for_w3 censor_just_porn_w3"'
	local var_censor_justif_reg_fp_w1 		`"censor_just_dom_economic_w1 censor_just_dom_political_w1 censor_just_dom_social_w1 censor_just_for_w1"'
	local var_censor_justif_reg_fp_w2 		`"censor_just_dom_economic_w2 censor_just_dom_political_w2 censor_just_dom_social_w2 censor_just_for_w2"'
	local var_censor_justif_reg_fp_w3 		`"censor_just_dom_economic_w3 censor_just_dom_political_w3 censor_just_dom_social_w3 censor_just_for_w3"'
	local var_censor_justif_fig_w123 		`"censor_just_dom_economic censor_just_dom_political censor_just_dom_social censor_just_for az_belief_media_justif_pf"'
	local var_censor_justif_fig_w23 		`"censor_just_porn"'
	
	local var_censor_level_reg_w1			`"bias_domestic_w1 bias_foreign_w1"'
	local var_censor_level_reg_w2			`"bias_domestic_w2 bias_foreign_w2"'
	local var_censor_level_reg_w3			`"bias_domestic_w3 bias_foreign_w3"'
	local var_censor_driver_reg_dom_w1 		`"bias_dom_govt_policy_t1_w1 bias_dom_firm_interest_t1_w1 bias_dom_media_pref_t1_w1 bias_dom_reader_demand_t1_w1"'
	local var_censor_driver_reg_dom_w2 		`"bias_dom_govt_policy_t1_w2 bias_dom_firm_interest_t1_w2 bias_dom_media_pref_t1_w2 bias_dom_reader_demand_t1_w2"'
	local var_censor_driver_reg_dom_w3 		`"bias_dom_govt_policy_t1_w3 bias_dom_firm_interest_t1_w3 bias_dom_media_pref_t1_w3 bias_dom_reader_demand_t1_w3"'
	local var_censor_driver_reg_for_w1		`"bias_for_govt_policy_t1_w1 bias_for_firm_interest_t1_w1 bias_for_media_pref_t1_w1 bias_for_reader_demand_t1_w1"'
	local var_censor_driver_reg_for_w2		`"bias_for_govt_policy_t1_w2 bias_for_firm_interest_t1_w2 bias_for_media_pref_t1_w2 bias_for_reader_demand_t1_w2"'
	local var_censor_driver_reg_for_w3		`"bias_for_govt_policy_t1_w3 bias_for_firm_interest_t1_w3 bias_for_media_pref_t1_w3 bias_for_reader_demand_t1_w3"'
	local var_censor_level_fig_w123 		`"bias_domestic bias_foreign"'
	local var_censor_driver_fig_w123 		`"bias_dom_govt_policy_t1 bias_dom_firm_interest_t1 bias_dom_media_pref_t1 bias_dom_reader_demand_t1 bias_for_govt_policy_t1 bias_for_firm_interest_t1 bias_for_media_pref_t1 bias_for_reader_demand_t1"'
	
	
* 	B: Knowledge
	local var_knowledge_news_reg_w1 		`"news_c_stock_rise news_c_exchange_depreciate news_c_train_brazil_peru news_c_army_reduce news_c_taiwan_election news_c_china_us_network_coop news_c_lijiacheng_china news_c_nanjing_memory_list"'
	local var_knowledge_news_reg_w2 		`"news_c_panamapapers news_c_stockcrash news_c_topinequality news_c_cpicensorship news_c_laborunrest news_c_tenyearshk news_c_economistcensor news_c_waterpollution news_c_applefbi news_c_taiwanelection news_c_yihehotel news_perccorrect_w2 news_perccor_cen_w2 news_perccor_unc_w2 news_perccor_qui"'
	local var_knowledge_news_reg_quiz 		`"news_c_topinequality news_c_cpicensorship news_c_laborunrest news_c_waterpollution"'
	local var_knowledge_news_reg_cen_w1 	`"news_c_stock_rise news_c_train_brazil_peru news_c_taiwan_election news_c_china_us_network_coop"'
	local var_knowledge_news_reg_cen_w2 	`"news_c_panamapapers news_c_tenyearshk news_c_stockcrash news_c_economistcensor"'
	local var_knowledge_news_reg_cen_w3 	`"news_c_coalprod news_c_trumpchina news_c_xiaojianhua news_c_xijiangcar news_c_chinanorway news_c_womenrights news_c_hkceelection"'
	local var_knowledge_news_reg_unc_w1 	`"news_c_exchange_depreciate news_c_army_reduce news_c_lijiacheng_china news_c_nanjing_memory_list"'
	local var_knowledge_news_reg_unc_w2 	`"news_c_applefbi news_c_taiwanelection news_c_yihehotel"'
	local var_knowledge_news_reg_unc_w3 	`"news_c_northkoreacoal news_c_birdflu news_c_ethiopiatrain news_c_foreignreserve"'
	local var_knowledge_news_fig_w123 		`"news_perccor_cen news_perccor_unc"'
	
	local var_knowledge_people_reg_tocens 	`"people_puzhiqiang_w2 people_renzhiqiang_w2 people_huangzhifeng_w2"'
	local var_knowledge_people_reg_censor 	`"people_lizehou_w2 people_chenguangcheng_w2 people_lixiaolin_w2"'
	local var_knowledge_people_reg_uncens 	`"people_maoyushi_w2 people_honghuang_w2 people_liuqiangdong_w2"'
	local var_knowledge_people_reg_fake 	`"people_jialequn_w2"'
	local var_knowledge_people_censor_w1 	`"people_lizehou_w1 people_huangzhifeng_w1 people_chenguangcheng_w1 people_lixiaolin_w1"'
	local var_knowledge_people_uncens_w1 	`"people_puzhiqiang_w1 people_renzhiqiang_w1 people_maoyushi_w1 people_honghuang_w1 people_liuqiangdong_w1"'
	local var_knowledge_people_fig_w12 		`"people_puzhiqiang people_lizehou people_huangzhifeng people_chenguangcheng people_lixiaolin people_renzhiqiang people_maoyushi people_honghuang people_liuqiangdong people_jialequn"'
	
	local var_knowledge_prot_reg_w2 		`"protest_2012_hk_curriculum_w2 protest_2014_umbrella_w2 protest_2016_mongkok_riot_w2 protest_2014_sun_flower_w2 protest_2014_europe_square_w2 protest_2010_arabic_spring_w2 protest_2014_crimea_vote_w2 protest_2010_catal_indep_w2 protest_2011_tmrw_parade_w2 protest_pcheard_total_w2 protest_pcheard_china_w2 protest_pcheard_foreign_w2"'
	local var_knowledge_prot_reg_w3 		`"protest_2012_hk_curriculum_w3 protest_2014_umbrella_w3 protest_2016_mongkok_riot_w3 protest_2014_sun_flower_w3 protest_2014_europe_square_w3 protest_2010_arabic_spring_w3 protest_2014_crimea_vote_w3 protest_2010_catal_indep_w3 protest_2011_tmrw_parade_w3 protest_2017_women_march_w3 protest_pcheard_total_w3 protest_pcheard_china_w3 protest_pcheard_foreign_w3"'
	local var_knowledge_prot_reg_chi_w1 	`"protest_2012_hk_curriculum_w1 protest_2014_umbrella_w1 protest_2014_sun_flower_w1"'
	local var_knowledge_prot_reg_chi_w2  	`"protest_2012_hk_curriculum_w2 protest_2014_umbrella_w2 protest_2016_mongkok_riot_w2 protest_2014_sun_flower_w2"'
	local var_knowledge_prot_reg_chi_w3  	`"protest_2012_hk_curriculum_w3 protest_2014_umbrella_w3 protest_2016_mongkok_riot_w3 protest_2014_sun_flower_w3"'
	local var_knowledge_prot_reg_for_w1 	`"protest_2014_europe_square_w1 protest_2010_arabic_spring_w1 protest_2014_crimea_vote_w1 protest_2010_catal_indep_w1"'
	local var_knowledge_prot_reg_for_w2 	`"protest_2014_europe_square_w2 protest_2010_arabic_spring_w2 protest_2014_crimea_vote_w2 protest_2010_catal_indep_w2"'
	local var_knowledge_prot_reg_for_w3 	`"protest_2014_europe_square_w3 protest_2010_arabic_spring_w3 protest_2014_crimea_vote_w3 protest_2010_catal_indep_w3 protest_2017_women_march_w3"'
	local var_knowledge_prot_reg_fak_w1 	`"protest_2011_tmrw_parade_w1"'
	local var_knowledge_prot_reg_fak_w2 	`"protest_2011_tmrw_parade_w2"'
	local var_knowledge_prot_reg_fak_w3 	`"protest_2011_tmrw_parade_w3"'
	local var_knowledge_prot_fig_w123 		`"protest_2014_europe_square protest_2014_sun_flower protest_2010_arabic_spring protest_2014_crimea_vote protest_2012_hk_curriculum protest_2010_catal_indep protest_2014_umbrella protest_2011_tmrw_parade protest_pcheard_china protest_pcheard_foreign"' 
	local var_knowledge_prot_fig_w23  		`"protest_2016_mongkok_riot"'
	local var_knowledge_prot_fig_w3  		`"protest_2017_women_march"'
	
	local var_knowledge_time_w2 			`"news_time_firstclick news_time_lastclick news_time_submit news_time_totalclick pppr_time_firstclick pppr_time_lastclick pppr_time_submit pppr_time_totalclick"'
	
	local var_knowledge_past_w1 			`"event_1994_beijing_jianguomen event_1994_kelamayi_fire event_2003_sunzhigang event_2008_guizhou_wengan event_2008_sanlu_milk_powder event_2010_hefei_wetland event_2011_guangdong_wukan event_2013_fudan_poison event_2014_caixin_xuxiao event_2015_sichuan_linshui_road"'
	
	local var_knowledge_meta_reg_w1 		`"familiar_china_issues_self_w1 familiar_china_others_w1"'
	local var_knowledge_meta_reg_w2 		`"familiar_china_issues_self_w2 familiar_china_others_w2"'
	local var_knowledge_meta_reg_w3 		`"familiar_china_issues_self_w3 familiar_china_others_w3"'
	local var_knowledge_meta_fig_w123 		`"familiar_china_issues_self familiar_china_others az_knowledge_meta"'

	
* 	C: Economic beliefs
	local var_econ_guess_reg_cn_perf_w1 	`"guess_gdp_growth_china_w1 guess_stock_index_sh_w1"'
	local var_econ_guess_reg_cn_perf_w2 	`"guess_gdp_growth_china_w2 guess_stock_index_sh_w2"'
	local var_econ_guess_reg_cn_perf_w3 	`"guess_gdp_growth_china_w3 guess_stock_index_sh_w3"'
	local var_econ_guess_reg_cn_conf_w1 	`"guess_gdp_growth_china_con_w1 guess_stock_index_sh_con_w1"'
	local var_econ_guess_reg_cn_conf_w2 	`"guess_gdp_growth_china_con_w2 guess_stock_index_sh_con_w2"'
	local var_econ_guess_reg_cn_conf_w3 	`"guess_gdp_growth_china_con_w3 guess_stock_index_sh_con_w3"'
	local var_econ_guess_reg_us_perf_w2 	`"guess_gdp_growth_us_w2 guess_stock_index_dj_w2"'
	local var_econ_guess_reg_us_perf_w3 	`"guess_gdp_growth_us_w3 guess_stock_index_dj_w3"'
	local var_econ_guess_reg_us_conf_w2 	`"guess_gdp_growth_us_con_w2 guess_stock_index_dj_con_w2"'
	local var_econ_guess_reg_us_conf_w3 	`"guess_gdp_growth_us_con_w3 guess_stock_index_dj_con_w3"'
	local var_econ_guess_fig_w123 			`"guess_gdp_growth_china guess_gdp_growth_china_con guess_stock_index_sh guess_stock_index_sh_con az_belief_econ_perf_cn az_belief_econ_conf_cn"'
	local var_econ_guess_fig_w23 			`"guess_gdp_growth_us guess_gdp_growth_us_con guess_stock_index_dj guess_stock_index_dj_con az_belief_econ_perf_us az_belief_econ_conf_us"'
	

* 	D: Political attitudes
	local var_demand_change_reg_w1 			`"inst_change_econ_w1 inst_change_poli_w1"'
	local var_demand_change_reg_w2 			`"inst_change_econ_w2 inst_change_poli_w2"'
	local var_demand_change_reg_w3 			`"inst_change_econ_w3 inst_change_poli_w3"'
	local var_demand_change_fig_w123 		`"inst_change_econ inst_change_poli az_belief_instchange"'
	
	local var_trust_inst_reg_govt_w1 		`"trust_central_govt_w1 trust_provincial_govt_w1 trust_local_govt_w1"'
	local var_trust_inst_reg_govt_w2 		`"trust_central_govt_w2 trust_provincial_govt_w2 trust_local_govt_w2"'
	local var_trust_inst_reg_govt_w3 		`"trust_central_govt_w3 trust_provincial_govt_w3 trust_local_govt_w3"'
	local var_trust_inst_reg_foreign_w1 	`"trust_japan_govt_w1 trust_us_govt_w1"'
	local var_trust_inst_reg_foreign_w2 	`"trust_japan_govt_w2 trust_us_govt_w2"'
	local var_trust_inst_reg_foreign_w3 	`"trust_japan_govt_w3 trust_us_govt_w3"'
	local var_trust_inst_reg_finance_w2 	`"trust_financial_domestic_w2 trust_financial_foreign_w2"'
	local var_trust_inst_reg_ngo_w2 		`"trust_ngo_w2"'
	local var_trust_inst_reg_copo_w2 		`"trust_court_w2 trust_police_w2"'
	local var_trust_inst_fig_w123  			`"trust_central_govt trust_provincial_govt trust_local_govt trust_japan_govt trust_us_govt az_belief_trust_govt az_belief_trust_foreign"'
	local var_trust_inst_fig_w12 			`"trust_court trust_police trust_ngo trust_financial_domestic trust_financial_foreign"'
	
	local var_eval_govt_reg_w1 				`"eval_govt_economic_w1 eval_govt_dom_politics_w1 eval_govt_for_relations_w1"'
	local var_eval_govt_reg_w2 				`"eval_govt_economic_w2 eval_govt_dom_politics_w2 eval_govt_for_relations_w2"'
	local var_eval_govt_reg_w3 				`"eval_govt_economic_w3 eval_govt_dom_politics_w3 eval_govt_for_relations_w3"'
	local var_eval_govt_fig_w123 			`"eval_govt_economic eval_govt_dom_politics eval_govt_for_relations az_belief_evalgovt"'
	
	local var_eval_criteria_reg_w2 			`"revalgovt_election_w2 revalgovt_economy_w2 revalgovt_equality_w2 revalgovt_ruleoflaw_w2 revalgovt_human_rights_w2 revalgovt_freedom_speech_w2 revalgovt_global_power_w2 revalgovt_fair_history_w2"'
	local var_eval_criteria_fig_w12 		`"revalgovt_election revalgovt_economy revalgovt_equality revalgovt_ruleoflaw revalgovt_human_rights revalgovt_freedom_speech revalgovt_global_power revalgovt_fair_history"'
	
	local var_severity_reg_w2 				`"severity_welfare_w2 severity_employment_w2 severity_pollution_w2 severity_inequality_w2 severity_corruption_w2 severity_dscrm_minority_w2"'
	local var_severity_fig_w12 				`"severity_welfare severity_employment severity_pollution severity_inequality severity_corruption severity_dscrm_minority"'
	
	local var_democracy_reg_w2 				`"china_interest_group_w2 china_rate_democracy_w2 china_rate_humanrights_w2 importance_live_in_demo_w2"'
	local var_democracy_reg_w3 				`"importance_live_in_demo_w3"'
	local var_democracy_reg_fp_w1 			`"importance_live_in_demo_w1"'
	local var_democracy_reg_fp_w2 			`"importance_live_in_demo_w2"'
	local var_democracy_reg_fp_w3 			`"importance_live_in_demo_w3"'
	local var_democracy_fig_w123 			`"importance_live_in_demo az_belief_democracy_fp"'
	local var_democracy_fig_w12 			`"china_interest_group china_rate_democracy china_rate_humanrights"'
	
	local var_contro_justi_reg_policy_w2 	`"justify_minority_policy_w2 justify_hukou_w2 justify_one_child_w2 justify_hongkong_policy_w2 justify_taiwan_policy_w2 justify_violence_stability_w2 justify_receive_refugee_w2 justify_reduce_pollution_w2 justify_gaokao_w2 justify_soe_privatize_w2"'
	local var_contro_justi_reg_liberal_w2 	`"justify_homo_marriage_w2 justify_legal_prostitute_w2 justify_abortion_w2 justify_exmaritalsex justify_transgene_w2 justify_soft_drug_w2"'
	local var_contro_justi_fig_w12 			`"justify_minority_policy justify_reduce_pollution justify_hukou justify_one_child justify_gaokao justify_hongkong_policy justify_taiwan_policy justify_transgene justify_receive_refugee justify_soe_privatize justify_homo_marriage justify_legal_prostitute justify_abortion justify_soft_drug justify_violence_stability"'
	
	local var_willing_fight_reg_w1 			`"willing_against_illi_govt_w1 willing_report_mis_w1 willing_protect_weak_w1"'
	local var_willing_fight_reg_w2 			`"willing_against_illi_govt_w2 willing_report_mis_w2 willing_protect_weak_w2"'
	local var_willing_fight_reg_w3 			`"willing_against_illi_govt_w3 willing_report_mis_w3 willing_protect_weak_w3"'
	local var_willing_fight_fig_123 		`"willing_against_illi_govt willing_report_mis willing_protect_weak az_belief_willing"'
	
	local var_interest_reg_w2 				`"interest_economic_w2 interest_politics_w2"'
	local var_interest_fig_w12 				`"interest_economic interest_politics"'
	
	local var_patriotism_reg_w2 			`"proud_being_chinese_w2"'
	local var_patriotism_fig_w12 			`"proud_being_chinese"'
	
	local var_fear_critgovt_reg_w2 			`"fear_critic_govt_self_w2"'
	local var_fear_critgovt_fig_w12 		`"fear_critic_govt_self"'
	
	
* 	E: Behaviors
	local var_socialinteract_reg_w1 		`"frequency_talk_politic_w1 frequency_persuade_friends_w1"'
	local var_socialinteract_reg_w2 		`"frequency_talk_politic_w2 frequency_persuade_friends_w2"'
	local var_socialinteract_reg_w3 		`"frequency_talk_politic_w3 frequency_persuade_friends_w3"'
	local var_socialinteract_fig_w123 		`"frequency_talk_politic frequency_persuade_friends az_var_socialinteract"'
	
	local var_socialknowledge_w1 			`"know_attitudes_relatives know_attitudes_schoolmates know_attitudes_outschool_friends"'
	
	local var_polparticipation_reg_w2 		`"participate_ngo_w2 participate_social_protest_w2 participate_plan_vote_w2 participate_complain_school_w2"'
	local var_polparticipation_reg_w3 		`"participate_social_protest_w3 participate_plan_vote_w3 participate_complain_school_w3"'
	local var_polparticipation_reg_pf_w1 	`"participate_social_protest_w1 participate_plan_vote_w1 participate_complain_school_w1"'
	local var_polparticipation_reg_pf_w2 	`"participate_social_protest_w2 participate_plan_vote_w2 participate_complain_school_w2"'
	local var_polparticipation_reg_pf_w3 	`"participate_social_protest_w3 participate_plan_vote_w3 participate_complain_school_w3"'
	local var_polparticipation_fig_w123 	`"participate_social_protest participate_plan_vote participate_complain_school az_var_polparticipation"'
	local var_polparticipation_fig_w12 		`"participate_ngo"'
	
	local var_polengagement_w1  			`"ccp_member ccp_join_year ccp_prep_member participate_stuunion participate_stuunion_ever participate_tuanwei participate_tuanwei_ever"'

	local var_planaftergrad_reg_w1 			`"plan_grad_gradschool_dom_w1 plan_grad_foreignmaster_w1 plan_grad_foreignphd_w1 plan_grad_military_w1 plan_grad_work_w1"'
	local var_planaftergrad_reg_w2 			`"plan_grad_gradschool_dom_w2 plan_grad_foreignmaster_w2 plan_grad_foreignphd_w2 plan_grad_military_w2 plan_grad_work_w2"'
	local var_planaftergrad_reg_w3 			`"plan_grad_gradschool_dom_w3 plan_grad_foreignmaster_w3 plan_grad_foreignphd_w3 plan_grad_military_w3 plan_grad_work_w3"'
	local var_planaftergrad_fig_w123 		`"plan_grad_gradschool_dom plan_grad_foreignmaster plan_grad_foreignphd plan_grad_military plan_grad_work"'
	
	local var_career_sector_reg_w1 			`"cp_t3_national_civil_w1 cp_t3_local_civil_w1 cp_t3_military_w1 cp_t3_chinese_private_w1 cp_t3_for_firm_w1 cp_t3_soe_w1 cp_t3_institutional_w1 cp_t3_entrepreneur_w1"'
	local var_career_sector_reg_w2 			`"cp_t3_national_civil_w2 cp_t3_local_civil_w2 cp_t3_military_w2 cp_t3_chinese_private_w2 cp_t3_for_firm_w2 cp_t3_soe_w2 cp_t3_institutional_w2 cp_t3_entrepreneur_w2"'
	local var_career_sector_reg_w3 			`"cp_t3_national_civil_w3 cp_t3_local_civil_w3 cp_t3_military_w3 cp_t3_chinese_private_w3 cp_t3_for_firm_w3 cp_t3_soe_w3 cp_t3_institutional_w3 cp_t3_entrepreneur_w3"'
	local var_career_sector_fig_w123 		`"cp_t1_national_civil cp_t1_local_civil cp_t1_military cp_t1_chinese_private cp_t1_for_firm cp_t1_soe cp_t1_institutional cp_t1_entrepreneur"'
	
	local var_career_loc_reg_w1 			`"cloc_beijing_w1 cloc_shanghai_w1 cloc_gzsz_w1 cloc_tjcq_w1 cloc_hkmc_w1 cloc_taiwan_w1 cloc_dom_w1 cloc_for_w1"'
	local var_career_loc_reg_w2 			`"cloc_beijing_w2 cloc_shanghai_w2 cloc_gzsz_w2 cloc_tjcq_w2 cloc_hkmc_w2 cloc_taiwan_w2 cloc_dom_w2 cloc_for_w2"'
	local var_career_loc_reg_w3 			`"cloc_beijing_w3 cloc_shanghai_w3 cloc_gzsz_w3 cloc_tjcq_w3 cloc_hkmc_w3 cloc_taiwan_w3 cloc_dom_w2 cloc_for_w3"'
	local var_career_loc_reg_w123 			`"cloc_beijing cloc_shanghai cloc_gzsz cloc_tjcq cloc_hkmc cloc_taiwan cloc_dom cloc_for"'
	
	local var_stock_invest_reg_w1  			`"stock_participation_w1"'
	local var_stock_invest_reg_w2  			`"stock_participation_w2"'
	local var_stock_invest_reg_w3  			`"stock_participation_w3"'
	local var_stock_invest_fig_w123 		`"stock_participation"'
	
	
* 	F: Demographics, background characteristics, and fundamental preferences
		
	// demographics
	local var_demog_reg 					`"gender birth_year ethnicity_han birthplace_coastal residence_coastal hukou_urban religion_religious ccp_member university_elite hs_track_science department_ssh siblings_total father_edu_hsabove father_ccp mother_edu_hsabove mother_ccp hh_income domestic_english_atleast4 foreign_english_yes travel_hktaiwan travel_foreign_yes"'
	local var_demog_reg_personal 			`"gender birth_year ethnicity_han birthplace_coastal residence_coastal hukou_urban religion_religious ccp_member"'
	local var_demog_reg_education 			`"university_elite hs_track_science department_ssh"'
	local var_demog_reg_english 			`"domestic_english_atleast4 foreign_english_yes"'
	local var_demog_reg_travel 				`"travel_hktaiwan travel_foreign_yes"'
	local var_demog_reg_household 			`"siblings_total father_edu_hsabove work_father_govt father_ccp mother_edu_hsabove work_mother_govt mother_ccp hh_income"'
	
	local var_demog_reg_imbalance 			`"residence_coastal hs_track_science father_ccp mother_ccp gift_amount risk_preference_ce"'
	
	// fundamental preferences
	local var_preference_risk 				`"willing_risk risk_preference_ce risk_lottery_choice"'
	local var_preference_time 				`"willing_future procrastinate"'
	local var_preference_altruism 			`"willing_goodcauses donate_amount"'
	local var_preference_reciprocity 		`"willing_returnfavor best_intentions gift_amount willing_punish_you willing_punish_others willing_revenge"'

		
* 	G: List experiment
	local var_listexp_all 					`"list_exp_direct list_trust_direct_count list_trust_direct_yes list_trust_veiled_count"'
	local var_listexp_reg 					`"list_count_trust"'
	
	
* 	Y: Overall effects
	local var_overall_w1 					`"info_freq_website_for_w1 az_belief_media_value_w1 az_belief_media_trust_w1 bias_domestic_w1 az_belief_media_justif_w1 bias_foreign_r_w1 news_perccor_cen_w1 news_perccor_unc_w1 protest_pcheard_china_w1 protest_pcheard_foreign_w1 az_knowledge_meta_w1 az_belief_econ_perf_cn_r_w1 az_belief_econ_conf_cn_w1 az_belief_instchange_w1 az_belief_trust_govt_r_w1 az_belief_trust_foreign_w1 az_belief_evalgovt_r_w1 importance_live_in_demo_r_w1 az_belief_willing_w1 az_var_socialinteract_w1 az_var_polparticipation_w1 plan_grad_gradschool_dom_w1 plan_grad_foreignmaster_w1 plan_grad_foreignphd_w1 plan_grad_military_w1 plan_grad_work_w1 cloc_beijing_w1 cloc_shanghai_w1 cloc_gzsz_w1 cloc_tjcq_w1 cloc_hkmc_w1 cloc_taiwan_w1 cloc_dom_w2 cloc_for_w1 stock_participation_w1"'
	local var_overall_w2 					`"info_freq_website_for_w2 az_belief_media_value_w2 az_belief_media_trust_w2 bias_domestic_w2 az_belief_media_justif_w2 bias_foreign_r_w2 news_perccor_cen_w2 news_perccor_unc_w2 protest_pcheard_china_w2 protest_pcheard_foreign_w2 az_knowledge_meta_w2 az_belief_econ_perf_cn_r_w2 az_belief_econ_conf_cn_w2 az_belief_instchange_w2 az_belief_trust_govt_r_w2 az_belief_trust_foreign_w2 az_belief_evalgovt_r_w2 importance_live_in_demo_r_w2 az_belief_willing_w2 az_var_socialinteract_w2 az_var_polparticipation_w2 plan_grad_gradschool_dom_w2 plan_grad_foreignmaster_w2 plan_grad_foreignphd_w2 plan_grad_military_w2 plan_grad_work_w2 cloc_beijing_w2 cloc_shanghai_w2 cloc_gzsz_w2 cloc_tjcq_w2 cloc_hkmc_w2 cloc_taiwan_w2 cloc_dom_w2 cloc_for_w2 stock_participation_w2"'
	local var_overall_w3 					`"info_freq_website_for_w3 az_belief_media_value_w3 az_belief_media_trust_w3 bias_domestic_w3 az_belief_media_justif_w3 bias_foreign_r_w3 news_perccor_cen_w3 news_perccor_unc_w3 protest_pcheard_china_w3 protest_pcheard_foreign_w3 az_knowledge_meta_w3 az_belief_econ_perf_cn_r_w3 az_belief_econ_conf_cn_w3 az_belief_instchange_w3 az_belief_trust_govt_r_w3 az_belief_trust_foreign_w3 az_belief_evalgovt_r_w3 importance_live_in_demo_r_w3 az_belief_willing_w3 az_var_socialinteract_w3 az_var_polparticipation_w3 plan_grad_gradschool_dom_w3 plan_grad_foreignmaster_w3 plan_grad_foreignphd_w3 plan_grad_military_w3 plan_grad_work_w3 cloc_beijing_w3 cloc_shanghai_w3 cloc_gzsz_w3 cloc_tjcq_w3 cloc_hkmc_w3 cloc_taiwan_w3 cloc_dom_w3 cloc_for_w3 stock_participation_w3"'

	
	
	
*** Prepare data: merge additional variables

* 	merge in vpn usage data: adoption dates
	merge 1:1 responseID_wave1 using vpn_date_adoption
	drop _merge
		

* 	merge in indicator of vpn active user
	merge 1:1 responseID_wave1 using vpn_browsing_active_user
	drop _merge
	replace active_user = 0 if active_user == .
	
	
* 	merge in vpn purchase data (from service provider server)
	merge 1:1 responseID_wave1 using vpn_purchase
	drop _merge
	
		


*** Prepare data: generate additional variables
	
* 	wave 1

	// foreach news, indicator correct or not
	gen news_c_stock_rise 				= (news_stock_rise == 1) if news_stock_rise != .
	gen news_c_exchange_depreciate 		= (news_exchange_depreciate == 1) if news_exchange_depreciate != .
	gen news_c_train_brazil_peru  		= (news_train_brazil_peru == 0) if news_train_brazil_peru != .
	gen news_c_army_reduce 				= (news_army_reduce == 1) if news_army_reduce != .
	gen news_c_taiwan_election 			= (news_taiwan_election == 0) if news_taiwan_election != .
	gen news_c_china_us_network_coop 	= (news_china_us_network_coop == 0) if news_china_us_network_coop != .
	gen news_c_lijiacheng_china 		= (news_lijiacheng_china == 0) if news_lijiacheng_china != .
	gen news_c_nanjing_memory_list 		= (news_nanjing_memory_list == 1) if news_nanjing_memory_list != .

	// generate count variable: # news answered correctly
	gen news_totalcorrect = news_c_stock_rise + news_c_exchange_depreciate + news_c_train_brazil_peru + news_c_army_reduce + news_c_taiwan_election + news_c_china_us_network_coop + news_c_lijiacheng_china + news_c_nanjing_memory_list
	
	// generate percentage correct
	gen news_perccorrect = news_totalcorrect / 8
	
	// correct count in percentage: censored news
	gen news_perccor_cen = (news_c_taiwan_election + news_c_lijiacheng_china) / 2
	
	// correct count in percentage: uncensored news
	gen news_perccor_unc = (news_c_stock_rise + news_c_exchange_depreciate + news_c_train_brazil_peru + news_c_army_reduce + news_c_china_us_network_coop + news_c_nanjing_memory_list) / 6

	// generate count variable: # people heard before
	gen people_totalheard = people_puzhiqiang_w1 + people_lizehou_w1 + people_huangzhifeng_w1 + people_chenguangcheng_w1 + people_lixiaolin_w1 + people_renzhiqiang_w1 + people_maoyushi_w1 + people_honghuang_w1 + people_liuqiangdong_w1
	gen people_pcheard = (people_puzhiqiang_w1 + people_lizehou_w1 + people_huangzhifeng_w1 + people_chenguangcheng_w1 + people_lixiaolin_w1 + people_renzhiqiang_w1 + people_maoyushi_w1 + people_honghuang_w1 + people_liuqiangdong_w1) / 9

	// generate count variable: # protests heard before
	gen protest_totalheard = protest_2014_europe_square_w1 + protest_2014_sun_flower_w1 + protest_2010_arabic_spring_w1 + protest_2014_crimea_vote_w1 + protest_2012_hk_curriculum_w1 + protest_2010_catal_indep_w1 + protest_2014_umbrella_w1
	gen protest_pcheard_total = (protest_2014_europe_square_w1 + protest_2014_sun_flower_w1 + protest_2010_arabic_spring_w1 + protest_2014_crimea_vote_w1 + protest_2012_hk_curriculum_w1 + protest_2010_catal_indep_w1 + protest_2014_umbrella_w1) / 9
	
	// generate percentage by category: full panel version
	gen protest_pcheard_china = (protest_2014_sun_flower_w1 + protest_2012_hk_curriculum_w1 + protest_2014_umbrella_w1) / 3
	gen protest_pcheard_foreign = (protest_2014_europe_square_w1 + protest_2010_arabic_spring_w1 + protest_2014_crimea_vote_w1 + protest_2010_catal_indep_w1) / 4

	// political participation: recode into ever participated
	gen participate_stuunion_ever = (participate_stuunion < 3) if participate_stuunion != .
	gen participate_tuanwei_ever = (participate_tuanwei < 3) if participate_tuanwei != .
	
	// generate distance from neutrality
	foreach var in us_pos_us us_pos_cn cn_pos_us cn_pos_cn us_neg_us us_neg_cn cn_neg_us cn_neg_cn {
		gen distneutral_`var' = abs(bias_`var'_media_w1 - 4) if bias_`var'_media_w1 != .
		}	
		
	// separate censorship and bias
	foreach var in us_pos_us us_pos_cn cn_pos_us cn_pos_cn us_neg_us us_neg_cn cn_neg_us cn_neg_cn {
		gen bias_`var'_cens = (bias_`var'_media_w1 == 1) if bias_`var'_media_w1 != .
		}
	
	// generate indicator for top categories: bias_for
	foreach var in bias_for_govt_policy bias_for_firm_interest bias_for_media_pref bias_for_reader_demand {
		gen `var'_t1 = (`var'_w1 == 1) if `var'_w1 != .
		}
		
	// generate indicator for top categories: bias_dom
	foreach var in bias_dom_govt_policy bias_dom_firm_interest bias_dom_media_pref bias_dom_reader_demand {
		gen `var'_t1 = (`var'_w1 == 1) if `var'_w1 != .
		}
		
	// generate relative weights for evaluating government performance
	gen evalgovt_total_w1 = evalgovt_election_w1 + evalgovt_economy_w1 + evalgovt_equality_w1 + evalgovt_ruleoflaw_w1 + evalgovt_human_rights_w1 + evalgovt_freedom_speech_w1 + evalgovt_global_power_w1 + evalgovt_fair_history_w1
	foreach var in evalgovt_election evalgovt_economy evalgovt_equality evalgovt_ruleoflaw evalgovt_human_rights evalgovt_freedom_speech evalgovt_global_power evalgovt_fair_history {
		gen r`var' = `var'_w1 / evalgovt_total_w1
		}

	// generate dummies for plan_graduation
	gen plan_grad_gradschool_dom 	= (plan_graduation_w1 == 1) if plan_graduation_w1 != .
	gen plan_grad_foreignmaster 	= (plan_graduation_w1 == 2) if plan_graduation_w1 != .
	gen plan_grad_foreignphd 		= (plan_graduation_w1 == 3) if plan_graduation_w1 != .
	gen plan_grad_military 			= (plan_graduation_w1 == 4) if plan_graduation_w1 != .
	gen plan_grad_work 				= (plan_graduation_w1 == 5) if plan_graduation_w1 != .
	gen plan_grad_dontknow 			= (plan_graduation_w1 == 6) if plan_graduation_w1 != .
	
	// generate top1 and top3 dummies for career preferences
	gen work_top3_missing = (work_top3_national_civil_w1 == . & work_top3_local_civil_w1 == . & work_top3_military_w1 == . & work_top3_chinese_private_w1 == . & work_top3_for_firm_w1 == . & work_top3_soe_w1 == . & work_top3_institutional_w1 == . & work_top3_entrepreneur_w1 == .)
	
	foreach var in national_civil local_civil military chinese_private for_firm soe institutional entrepreneur {
		gen cp_t1_`var' = (work_top3_`var'_w1 == 1) if work_top3_missing == 0
		gen cp_t3_`var' = (work_top3_`var'_w1 != .) if work_top3_missing == 0
		}
	drop work_top3_missing
	
	// generate dummies for career location preferences
	gen cloc_beijing 		= (place_top_w1 == 1) if place_top_w1 != .
	gen cloc_shanghai 		= (place_top_w1 == 2) if place_top_w1 != .
	gen cloc_gzsz 			= (place_top_w1 == 3 | place_top_w1 == 4) if place_top_w1 != .
	gen cloc_tjcq 			= (place_top_w1 == 5 | place_top_w1 == 6) if place_top_w1 != .
	gen cloc_hkmc 			= (place_top_w1 == 7 | place_top_w1 == 8) if place_top_w1 != .
	gen cloc_taiwan 		= (place_top_w1 == 9) if place_top_w1 != .
	gen cloc_dom 			= (place_top_w1 == 10) if place_top_w1 != .
	gen cloc_for  			= (place_top_w1 == 11) if place_top_w1 != .

	
	// add subscript: wave 1
	foreach var in news_totalcorrect news_perccorrect news_perccor_cen news_perccor_unc people_totalheard people_pcheard protest_totalheard protest_pcheard_total protest_pcheard_china protest_pcheard_foreign distneutral_us_pos_us distneutral_us_pos_cn distneutral_cn_pos_us distneutral_cn_pos_cn distneutral_us_neg_us distneutral_us_neg_cn distneutral_cn_neg_us distneutral_cn_neg_cn bias_us_pos_us_cens bias_us_pos_cn_cens bias_cn_pos_us_cens bias_cn_pos_cn_cens bias_us_neg_us_cens bias_us_neg_cn_cens bias_cn_neg_us_cens bias_cn_neg_cn_cens bias_for_govt_policy_t1 bias_for_firm_interest_t1 bias_for_media_pref_t1 bias_for_reader_demand_t1 bias_dom_govt_policy_t1 bias_dom_firm_interest_t1 bias_dom_media_pref_t1 bias_dom_reader_demand_t1 revalgovt_election revalgovt_economy revalgovt_equality revalgovt_ruleoflaw revalgovt_human_rights revalgovt_freedom_speech revalgovt_global_power revalgovt_fair_history plan_grad_gradschool_dom plan_grad_foreignmaster plan_grad_foreignphd plan_grad_military plan_grad_work plan_grad_dontknow cp_t1_national_civil cp_t3_national_civil cp_t1_local_civil cp_t3_local_civil cp_t1_military cp_t3_military cp_t1_chinese_private cp_t3_chinese_private cp_t1_for_firm cp_t3_for_firm cp_t1_soe cp_t3_soe cp_t1_institutional cp_t3_institutional cp_t1_entrepreneur cp_t3_entrepreneur cloc_beijing cloc_shanghai cloc_gzsz cloc_tjcq cloc_hkmc cloc_taiwan cloc_dom cloc_for {
		rename `var' 	`var'_w1
		}
	

* 	wave 2
	
	// foreach news, indicator correct or not
	gen news_c_topinequality 		= (news_topinequality == 0) if news_topinequality != .
	gen news_c_cpicensorship 	 	= (news_cpicensorship == 1) if news_cpicensorship != .
	gen news_c_laborunrest 			= (news_laborunrest == 1) if news_laborunrest != .
	gen news_c_waterpollution 		= (news_waterpollution == 0) if news_waterpollution != .
	gen news_c_stockcrash 			= (news_stockcrash == 1) if news_stockcrash != .
	gen news_c_taiwanelection 		= (news_taiwanelection == 0) if news_taiwanelection != .
	gen news_c_applefbi 			= (news_applefbi == 0) if news_applefbi != .
	gen news_c_tenyearshk 			= (news_tenyearshk == 1) if news_tenyearshk != .
	gen news_c_panamapapers 		= (news_panamapapers == 1) if news_panamapapers != .
	gen news_c_yihehotel 			= (news_yihehotel == 1) if news_yihehotel != .
	gen news_c_economistcensor 		= (news_economistcensor == 1) if news_economistcensor != .

	// generate count variable: # news answered correctly
	gen news_totalcorrect = news_c_topinequality + news_c_cpicensorship + news_c_laborunrest + news_c_waterpollution + news_c_stockcrash + news_c_taiwanelection + news_c_applefbi + news_c_tenyearshk + news_c_panamapapers + news_c_yihehotel + news_c_economistcensor

	// generate percentage correct
	gen news_perccorrect = news_totalcorrect / 11
	
	// correct count in percentage: quiz questions
	gen news_perccor_qui = (news_c_topinequality + news_c_cpicensorship + news_c_laborunrest + news_c_waterpollution) / 4
	
	// correct count in percentage: censored news
	gen news_perccor_cen = (news_c_stockcrash + news_c_panamapapers + news_c_tenyearshk + news_c_economistcensor) / 4
	
	// correct count in percentage: uncensored news
	gen news_perccor_unc = (news_c_taiwanelection + news_c_applefbi + news_c_yihehotel) / 3
		
	// generate count variable: # people heard before
	gen people_totalheard = people_puzhiqiang_w2 + people_lizehou_w2 + people_huangzhifeng_w2 + people_chenguangcheng_w2 + people_lixiaolin_w2 + people_renzhiqiang_w2 + people_maoyushi_w2 + people_honghuang_w2 + people_liuqiangdong_w2
	gen people_pcheard = (people_puzhiqiang_w2 + people_lizehou_w2 + people_huangzhifeng_w2 + people_chenguangcheng_w2 + people_lixiaolin_w2 + people_renzhiqiang_w2 + people_maoyushi_w2 + people_honghuang_w2 + people_liuqiangdong_w2) / 9
		
	// percentage of names heard by category
	gen people_perchd_censored = (people_lixiaolin_w2 + people_huangzhifeng_w2 + people_chenguangcheng_w2 + people_lizehou_w2) / 4
	gen people_perchd_uncensor = (people_maoyushi_w2 + people_liuqiangdong_w2 + people_honghuang_w2) / 3
	gen people_perchd_unctocen = (people_puzhiqiang_w2 + people_renzhiqiang_w2) / 2

	// generate count variable: # events heard before
	gen protest_totalheard = protest_2014_europe_square_w2 + protest_2014_sun_flower_w2 + protest_2010_arabic_spring_w2 + protest_2014_crimea_vote_w2 + protest_2012_hk_curriculum_w2 + protest_2010_catal_indep_w2 + protest_2014_umbrella_w2 + protest_2016_mongkok_riot_w2
	gen protest_pcheard_total = (protest_2014_europe_square_w2 + protest_2014_sun_flower_w2 + protest_2010_arabic_spring_w2 + protest_2014_crimea_vote_w2 + protest_2012_hk_curriculum_w2 + protest_2010_catal_indep_w2 + protest_2014_umbrella_w2 + protest_2016_mongkok_riot_w2) / 8
	
	// generate percentage by category: full panel version
	gen protest_pcheard_china = (protest_2014_sun_flower_w2 + protest_2012_hk_curriculum_w2 + protest_2014_umbrella_w2) / 3
	gen protest_pcheard_foreign = (protest_2014_europe_square_w2 + protest_2010_arabic_spring_w2 + protest_2014_crimea_vote_w2 + protest_2010_catal_indep_w2) / 4
		
	// generate distance from neutrality
	foreach var in us_pos_us us_pos_cn cn_pos_us cn_pos_cn us_neg_us us_neg_cn cn_neg_us cn_neg_cn {
		gen distneutral_`var' = abs(bias_`var'_media_w2 - 4) if bias_`var'_media_w2 != .
		}	
		
	// separate censorship and bias
	foreach var in us_pos_us us_pos_cn cn_pos_us cn_pos_cn us_neg_us us_neg_cn cn_neg_us cn_neg_cn {
		gen bias_`var'_cens = (bias_`var'_media_w2 == 1) if bias_`var'_media_w2 != .
		}
		
	// generate indicator for top categories: bias_for
	foreach var in bias_for_govt_policy bias_for_firm_interest bias_for_media_pref bias_for_reader_demand {
		gen `var'_t1 = (`var'_w2 == 1) if `var'_w2 != .
		}
		
	// generate indicator for top categories: bias_dom
	foreach var in bias_dom_govt_policy bias_dom_firm_interest bias_dom_media_pref bias_dom_reader_demand {
		gen `var'_t1 = (`var'_w2 == 1) if `var'_w2 != .
		}

	// generate relative weights for evaluating government performance
	gen evalgovt_total_w2 = evalgovt_election_w2 + evalgovt_economy_w2 + evalgovt_equality_w2 + evalgovt_ruleoflaw_w2 + evalgovt_human_rights_w2 + evalgovt_freedom_speech_w2 + evalgovt_global_power_w2 + evalgovt_fair_history_w2
	foreach var in evalgovt_election evalgovt_economy evalgovt_equality evalgovt_ruleoflaw evalgovt_human_rights evalgovt_freedom_speech evalgovt_global_power evalgovt_fair_history {
		gen r`var' = `var'_w2 / evalgovt_total_w2
		}

	// generate dummies for plan_graduation
	gen plan_grad_gradschool_dom 	= (plan_graduation_w2 == 1) if plan_graduation_w2 != .
	gen plan_grad_foreignmaster 	= (plan_graduation_w2 == 2) if plan_graduation_w2 != .
	gen plan_grad_foreignphd 		= (plan_graduation_w2 == 3) if plan_graduation_w2 != .
	gen plan_grad_military 			= (plan_graduation_w2 == 4) if plan_graduation_w2 != .
	gen plan_grad_work 				= (plan_graduation_w2 == 5) if plan_graduation_w2 != .
	gen plan_grad_dontknow 			= (plan_graduation_w2 == 6) if plan_graduation_w2 != .

	// generate top1 and top3 dummies for career preferences
	gen work_top3_missing = (work_top3_national_civil_w2 == . & work_top3_local_civil_w2 == . & work_top3_military_w2 == . & work_top3_chinese_private_w2 == . & work_top3_for_firm_w2 == . & work_top3_soe_w2 == . & work_top3_institutional_w2 == . & work_top3_entrepreneur_w2 == .)

	foreach var in national_civil local_civil military chinese_private for_firm soe institutional entrepreneur {
		gen cp_t1_`var' = (work_top3_`var'_w2 == 1) if work_top3_missing == 0
		gen cp_t3_`var' = (work_top3_`var'_w2 != .) if work_top3_missing == 0
		}
	drop work_top3_missing

	// generate dummies for career location preferences
	gen cloc_beijing 		= (place_top_w2 == 1) if place_top_w2 != .
	gen cloc_shanghai 		= (place_top_w2 == 2) if place_top_w2 != .
	gen cloc_gzsz 			= (place_top_w2 == 3 | place_top_w2 == 4) if place_top_w2 != .
	gen cloc_tjcq 			= (place_top_w2 == 5 | place_top_w2 == 6) if place_top_w2 != .
	gen cloc_hkmc 			= (place_top_w2 == 7 | place_top_w2 == 8) if place_top_w2 != .
	gen cloc_taiwan 		= (place_top_w2 == 9) if place_top_w2 != .
	gen cloc_dom 			= (place_top_w2 == 10) if place_top_w2 != .
	gen cloc_for  			= (place_top_w2 == 11) if place_top_w2 != .

	
	// add subscript: wave 2
	foreach var in news_totalcorrect news_perccorrect news_perccor_cen news_perccor_unc people_totalheard people_pcheard people_perchd_censored people_perchd_uncensor people_perchd_unctocen protest_totalheard protest_pcheard_total protest_pcheard_china protest_pcheard_foreign distneutral_us_pos_us distneutral_us_pos_cn distneutral_cn_pos_us distneutral_cn_pos_cn distneutral_us_neg_us distneutral_us_neg_cn distneutral_cn_neg_us distneutral_cn_neg_cn bias_us_pos_us_cens bias_us_pos_cn_cens bias_cn_pos_us_cens bias_cn_pos_cn_cens bias_us_neg_us_cens bias_us_neg_cn_cens bias_cn_neg_us_cens bias_cn_neg_cn_cens bias_for_govt_policy_t1 bias_for_firm_interest_t1 bias_for_media_pref_t1 bias_for_reader_demand_t1 bias_dom_govt_policy_t1 bias_dom_firm_interest_t1 bias_dom_media_pref_t1 bias_dom_reader_demand_t1 revalgovt_election revalgovt_economy revalgovt_equality revalgovt_ruleoflaw revalgovt_human_rights revalgovt_freedom_speech revalgovt_global_power revalgovt_fair_history plan_grad_gradschool_dom plan_grad_foreignmaster plan_grad_foreignphd plan_grad_military plan_grad_work plan_grad_dontknow cp_t1_national_civil cp_t3_national_civil cp_t1_local_civil cp_t3_local_civil cp_t1_military cp_t3_military cp_t1_chinese_private cp_t3_chinese_private cp_t1_for_firm cp_t3_for_firm cp_t1_soe cp_t3_soe cp_t1_institutional cp_t3_institutional cp_t1_entrepreneur cp_t3_entrepreneur cloc_beijing cloc_shanghai cloc_gzsz cloc_tjcq cloc_hkmc cloc_taiwan cloc_dom cloc_for {
		rename `var' 	`var'_w2
		}
	
	
* 	wave 3

	// foreach news, indicator correct or not
	gen news_c_coalprod  		= (news_coalprod == 0) if news_coalprod != .
	gen news_c_trumpchina 		= (news_trumpchina == 1) if news_trumpchina != .
	gen news_c_xiaojianhua 		= (news_xiaojianhua == 0) if news_xiaojianhua != .
	gen news_c_northkoreacoal  	= (news_northkoreacoal == 1) if news_northkoreacoal != .
	gen news_c_birdflu 			= (news_birdflu == 1) if news_birdflu != .
	gen news_c_ethiopiatrain 	= (news_ethiopiatrain == 0) if news_ethiopiatrain != .
	gen news_c_foreignreserve 	= (news_foreignreserve == 0) if news_foreignreserve != .
	gen news_c_xijiangcar 		= (news_xijiangcar == 0) if news_xijiangcar != .
	gen news_c_chinanorway 		= (news_chinanorway == 1) if news_chinanorway != .
	gen news_c_womenrights 		= (news_womenrights == 0) if news_womenrights != .
	gen news_c_hkceelection 	= (news_hkceelection == 0) if news_hkceelection != .
	
	// generate count variable: # news answered correctly
	gen news_totalcorrect = news_c_coalprod + news_c_trumpchina + news_c_xiaojianhua + news_c_northkoreacoal + news_c_birdflu + news_c_ethiopiatrain + news_c_foreignreserve + news_c_xijiangcar + news_c_chinanorway + news_c_womenrights + news_c_hkceelection

	// generate percentage correct
	gen news_perccorrect = news_totalcorrect / 11
	
	// correct count in percentage: censored news
	gen news_perccor_cen = (news_c_coalprod + news_c_trumpchina + news_c_xiaojianhua + news_c_xijiangcar + news_c_chinanorway + news_c_womenrights + news_c_hkceelection) / 7
	
	// correct count in percentage: uncensored news
	gen news_perccor_unc = (news_c_northkoreacoal + news_c_birdflu + news_c_ethiopiatrain + news_c_foreignreserve) / 4
	
	// generate count variable: # events heard before
	gen protest_totalheard = protest_2014_europe_square_w3 + protest_2014_sun_flower_w3 + protest_2010_arabic_spring_w3 + protest_2014_crimea_vote_w3 + protest_2012_hk_curriculum_w3 + protest_2010_catal_indep_w3 + protest_2014_umbrella_w3 + protest_2016_mongkok_riot_w3 + protest_2017_women_march_w3
	gen protest_pcheard_total = (protest_2014_europe_square_w3 + protest_2014_sun_flower_w3 + protest_2010_arabic_spring_w3 + protest_2014_crimea_vote_w3 + protest_2012_hk_curriculum_w3 + protest_2010_catal_indep_w3 + protest_2014_umbrella_w3 + protest_2016_mongkok_riot_w3 + protest_2017_women_march_w3) / 9
	
	// generate percentage by category: full panel version
	gen protest_pcheard_china = (protest_2014_sun_flower_w3 + protest_2012_hk_curriculum_w3 + protest_2014_umbrella_w3) / 3
	gen protest_pcheard_foreign = (protest_2014_europe_square_w3 + protest_2010_arabic_spring_w3 + protest_2014_crimea_vote_w3 + protest_2010_catal_indep_w3) / 4

	// generate indicator for top categories: bias_for
	foreach var in bias_for_govt_policy bias_for_firm_interest bias_for_media_pref bias_for_reader_demand {
		gen `var'_t1 = (`var'_w3 == 1) if `var'_w3 != .
		}
		
	// generate indicator for top categories: bias_dom
	foreach var in bias_dom_govt_policy bias_dom_firm_interest bias_dom_media_pref bias_dom_reader_demand {
		gen `var'_t1 = (`var'_w3 == 1) if `var'_w3 != .
		}

	// generate indicators for regression
	gen vpn_purchase_yes = (vpn_purchase < 6) if vpn_purchase != .
	
	gen vpn_purchase_wmt = ((vpn_purchase_yes == 1) & (vpn_purchase != 1)) if vpn_purchase != .
	
	gen vpn_purchase_premium = .
	replace vpn_purchase_premium = 0 	if (vpn_purchase == 6)
	replace vpn_purchase_premium = 1 	if (vpn_purchase == 1)
	replace vpn_purchase_premium = 2 	if (vpn_purchase == 4)
	replace vpn_purchase_premium = 3 	if (vpn_purchase == 3)
	replace vpn_purchase_premium = 4 	if (vpn_purchase == 5)	

	// generate dummies for plan_graduation
	gen plan_grad_gradschool_dom 	= (plan_graduation_w3 == 1) if plan_graduation_w3 != .
	gen plan_grad_foreignmaster 	= (plan_graduation_w3 == 2) if plan_graduation_w3 != .
	gen plan_grad_foreignphd 		= (plan_graduation_w3 == 3) if plan_graduation_w3 != .
	gen plan_grad_military 			= (plan_graduation_w3 == 4) if plan_graduation_w3 != .
	gen plan_grad_work 				= (plan_graduation_w3 == 5) if plan_graduation_w3 != .
	gen plan_grad_dontknow 			= (plan_graduation_w3 == 6) if plan_graduation_w3 != .

	// generate top1 and top3 dummies for career preferences
	gen work_top3_missing = (work_top3_national_civil_w3 == . & work_top3_local_civil_w3 == . & work_top3_military_w3 == . & work_top3_chinese_private_w3 == . & work_top3_for_firm_w3 == . & work_top3_soe_w3 == . & work_top3_institutional_w3 == . & work_top3_entrepreneur_w3 == .)

	foreach var in national_civil local_civil military chinese_private for_firm soe institutional entrepreneur {
		gen cp_t1_`var' = (work_top3_`var'_w3 == 1) if work_top3_missing == 0
		gen cp_t3_`var' = (work_top3_`var'_w3 != .) if work_top3_missing == 0
		}
	drop work_top3_missing

	// generate dummies for career location preferences
	gen cloc_beijing 		= (place_top_w3 == 1) if place_top_w3 != .
	gen cloc_shanghai 		= (place_top_w3 == 2) if place_top_w3 != .
	gen cloc_gzsz 			= (place_top_w3 == 3 | place_top_w3 == 4) if place_top_w3 != .
	gen cloc_tjcq 			= (place_top_w3 == 5 | place_top_w3 == 6) if place_top_w3 != .
	gen cloc_hkmc 			= (place_top_w3 == 7 | place_top_w3 == 8) if place_top_w3 != .
	gen cloc_taiwan 		= (place_top_w3 == 9) if place_top_w3 != .
	gen cloc_dom 			= (place_top_w3 == 10) if place_top_w3 != .
	gen cloc_for  			= (place_top_w3 == 11) if place_top_w3 != .

	
	// add subscript: wave 3
	foreach var in news_totalcorrect news_perccorrect news_perccor_cen news_perccor_unc protest_totalheard protest_pcheard_total protest_pcheard_china protest_pcheard_foreign bias_for_govt_policy_t1 bias_for_firm_interest_t1 bias_for_media_pref_t1 bias_for_reader_demand_t1 bias_dom_govt_policy_t1 bias_dom_firm_interest_t1 bias_dom_media_pref_t1 bias_dom_reader_demand_t1 plan_grad_gradschool_dom plan_grad_foreignmaster plan_grad_foreignphd plan_grad_military plan_grad_work plan_grad_dontknow cp_t1_national_civil cp_t3_national_civil cp_t1_local_civil cp_t3_local_civil cp_t1_military cp_t3_military cp_t1_chinese_private cp_t3_chinese_private cp_t1_for_firm cp_t3_for_firm cp_t1_soe cp_t3_soe cp_t1_institutional cp_t3_institutional cp_t1_entrepreneur cp_t3_entrepreneur cloc_beijing cloc_shanghai cloc_gzsz cloc_tjcq cloc_hkmc cloc_taiwan cloc_dom cloc_for {
		rename `var' 	`var'_w3
		}

	// generate percentage of total news quiz
	gen news_perccor_cen_all = (news_perccor_cen_w2 + news_perccor_cen_w3)/2
	

		
* 	generate certainty equivalent from Falk et al. 2015 staircase risk preference

	gen risk_preference_ce = .
	replace risk_preference_ce = 32 	if (risk_preference_1 == 1 & risk_preference_17 == 1 & risk_preference_25 == 1 & risk_preference_29 == 1 & risk_preference_31 == 1)
	replace risk_preference_ce = 31 	if (risk_preference_1 == 1 & risk_preference_17 == 1 & risk_preference_25 == 1 & risk_preference_29 == 1 & risk_preference_31 == 2)
	replace risk_preference_ce = 30 	if (risk_preference_1 == 1 & risk_preference_17 == 1 & risk_preference_25 == 1 & risk_preference_29 == 2 & risk_preference_30 == 1)
	replace risk_preference_ce = 29 	if (risk_preference_1 == 1 & risk_preference_17 == 1 & risk_preference_25 == 1 & risk_preference_29 == 2 & risk_preference_30 == 2)
	replace risk_preference_ce = 28 	if (risk_preference_1 == 1 & risk_preference_17 == 1 & risk_preference_25 == 2 & risk_preference_26 == 1 & risk_preference_27 == 1)
	replace risk_preference_ce = 27 	if (risk_preference_1 == 1 & risk_preference_17 == 1 & risk_preference_25 == 2 & risk_preference_26 == 1 & risk_preference_27 == 2)
	replace risk_preference_ce = 26 	if (risk_preference_1 == 1 & risk_preference_17 == 1 & risk_preference_25 == 2 & risk_preference_26 == 2 & risk_preference_28 == 1)
	replace risk_preference_ce = 25 	if (risk_preference_1 == 1 & risk_preference_17 == 1 & risk_preference_25 == 2 & risk_preference_26 == 2 & risk_preference_28 == 2)
	replace risk_preference_ce = 24 	if (risk_preference_1 == 1 & risk_preference_17 == 2 & risk_preference_18 == 1 & risk_preference_22 == 1 & risk_preference_23 == 1)
	replace risk_preference_ce = 23 	if (risk_preference_1 == 1 & risk_preference_17 == 2 & risk_preference_18 == 1 & risk_preference_22 == 1 & risk_preference_23 == 2)
	replace risk_preference_ce = 22 	if (risk_preference_1 == 1 & risk_preference_17 == 2 & risk_preference_18 == 1 & risk_preference_22 == 2 & risk_preference_24 == 1)
	replace risk_preference_ce = 21 	if (risk_preference_1 == 1 & risk_preference_17 == 2 & risk_preference_18 == 1 & risk_preference_22 == 2 & risk_preference_24 == 2)
	replace risk_preference_ce = 20 	if (risk_preference_1 == 1 & risk_preference_17 == 2 & risk_preference_18 == 2 & risk_preference_19 == 1 & risk_preference_20 == 1)
	replace risk_preference_ce = 19 	if (risk_preference_1 == 1 & risk_preference_17 == 2 & risk_preference_18 == 2 & risk_preference_19 == 1 & risk_preference_20 == 2)
	replace risk_preference_ce = 18 	if (risk_preference_1 == 1 & risk_preference_17 == 2 & risk_preference_18 == 2 & risk_preference_19 == 2 & risk_preference_21 == 1)
	replace risk_preference_ce = 17 	if (risk_preference_1 == 1 & risk_preference_17 == 2 & risk_preference_18 == 2 & risk_preference_19 == 2 & risk_preference_21 == 2)
	replace risk_preference_ce = 16 	if (risk_preference_1 == 2 & risk_preference_2  == 1 & risk_preference_10 == 1 & risk_preference_14 == 1 & risk_preference_15 == 1)
	replace risk_preference_ce = 15 	if (risk_preference_1 == 2 & risk_preference_2  == 1 & risk_preference_10 == 1 & risk_preference_14 == 1 & risk_preference_15 == 2)
	replace risk_preference_ce = 14 	if (risk_preference_1 == 2 & risk_preference_2  == 1 & risk_preference_10 == 1 & risk_preference_14 == 2 & risk_preference_16 == 1)
	replace risk_preference_ce = 13 	if (risk_preference_1 == 2 & risk_preference_2  == 1 & risk_preference_10 == 1 & risk_preference_14 == 2 & risk_preference_16 == 2)
	replace risk_preference_ce = 12 	if (risk_preference_1 == 2 & risk_preference_2  == 1 & risk_preference_10 == 2 & risk_preference_11 == 1 & risk_preference_13 == 1)
	replace risk_preference_ce = 11 	if (risk_preference_1 == 2 & risk_preference_2  == 1 & risk_preference_10 == 2 & risk_preference_11 == 1 & risk_preference_13 == 2)
	replace risk_preference_ce = 10 	if (risk_preference_1 == 2 & risk_preference_2  == 1 & risk_preference_10 == 2 & risk_preference_11 == 2 & risk_preference_12 == 1)
	replace risk_preference_ce = 9 		if (risk_preference_1 == 2 & risk_preference_2  == 1 & risk_preference_10 == 2 & risk_preference_11 == 2 & risk_preference_12 == 2)
	replace risk_preference_ce = 8 		if (risk_preference_1 == 2 & risk_preference_2  == 2 & risk_preference_3  == 1 & risk_preference_4  == 1 & risk_preference_5  == 1)
	replace risk_preference_ce = 7 		if (risk_preference_1 == 2 & risk_preference_2  == 2 & risk_preference_3  == 1 & risk_preference_4  == 1 & risk_preference_5  == 2)
	replace risk_preference_ce = 6 		if (risk_preference_1 == 2 & risk_preference_2  == 2 & risk_preference_3  == 1 & risk_preference_4  == 2 & risk_preference_6  == 1)
	replace risk_preference_ce = 5 		if (risk_preference_1 == 2 & risk_preference_2  == 2 & risk_preference_3  == 1 & risk_preference_4  == 2 & risk_preference_6  == 2)
	replace risk_preference_ce = 4 		if (risk_preference_1 == 2 & risk_preference_2  == 2 & risk_preference_3  == 2 & risk_preference_7  == 1 & risk_preference_8  == 1)
	replace risk_preference_ce = 3 		if (risk_preference_1 == 2 & risk_preference_2  == 2 & risk_preference_3  == 2 & risk_preference_7  == 1 & risk_preference_8  == 2)
	replace risk_preference_ce = 2 		if (risk_preference_1 == 2 & risk_preference_2  == 2 & risk_preference_3  == 2 & risk_preference_7  == 2 & risk_preference_9  == 1)
	replace risk_preference_ce = 1 		if (risk_preference_1 == 2 & risk_preference_2  == 2 & risk_preference_3  == 2 & risk_preference_7  == 2 & risk_preference_9  == 2)

	
* 	generate certainty equivalent from Falk et al. 2015 staircase risk preference

	gen time_preference_fe = .
	replace time_preference_fe = 32 	if (time_preference_1 == 1 & time_preference_17 == 1 & time_preference_18 == 1 & time_preference_22 == 1 & time_preference_23 == 1)
	replace time_preference_fe = 31 	if (time_preference_1 == 1 & time_preference_17 == 1 & time_preference_18 == 1 & time_preference_22 == 1 & time_preference_23 == 2)
	replace time_preference_fe = 30 	if (time_preference_1 == 1 & time_preference_17 == 1 & time_preference_18 == 1 & time_preference_22 == 2 & time_preference_24 == 1)
	replace time_preference_fe = 29 	if (time_preference_1 == 1 & time_preference_17 == 1 & time_preference_18 == 1 & time_preference_22 == 2 & time_preference_24 == 2)
	replace time_preference_fe = 28 	if (time_preference_1 == 1 & time_preference_17 == 1 & time_preference_18 == 2 & time_preference_19 == 1 & time_preference_20 == 1)
	replace time_preference_fe = 27 	if (time_preference_1 == 1 & time_preference_17 == 1 & time_preference_18 == 2 & time_preference_19 == 1 & time_preference_20 == 2)
	replace time_preference_fe = 26 	if (time_preference_1 == 1 & time_preference_17 == 1 & time_preference_18 == 2 & time_preference_19 == 2 & time_preference_21 == 1)
	replace time_preference_fe = 25 	if (time_preference_1 == 1 & time_preference_17 == 1 & time_preference_18 == 2 & time_preference_19 == 2 & time_preference_21 == 2)
	replace time_preference_fe = 24 	if (time_preference_1 == 1 & time_preference_17 == 2 & time_preference_25 == 1 & time_preference_29 == 1 & time_preference_31 == 1)
	replace time_preference_fe = 23 	if (time_preference_1 == 1 & time_preference_17 == 2 & time_preference_25 == 1 & time_preference_29 == 1 & time_preference_31 == 2)
	replace time_preference_fe = 22 	if (time_preference_1 == 1 & time_preference_17 == 2 & time_preference_25 == 1 & time_preference_29 == 2 & time_preference_30 == 1)
	replace time_preference_fe = 21 	if (time_preference_1 == 1 & time_preference_17 == 2 & time_preference_25 == 1 & time_preference_29 == 2 & time_preference_30 == 2)
	replace time_preference_fe = 20 	if (time_preference_1 == 1 & time_preference_17 == 2 & time_preference_25 == 2 & time_preference_26 == 1 & time_preference_28 == 1)
	replace time_preference_fe = 19 	if (time_preference_1 == 1 & time_preference_17 == 2 & time_preference_25 == 2 & time_preference_26 == 1 & time_preference_28 == 2)
	replace time_preference_fe = 18 	if (time_preference_1 == 1 & time_preference_17 == 2 & time_preference_25 == 2 & time_preference_26 == 2 & time_preference_27 == 1)
	replace time_preference_fe = 17 	if (time_preference_1 == 1 & time_preference_17 == 2 & time_preference_25 == 2 & time_preference_26 == 2 & time_preference_27 == 2)
	replace time_preference_fe = 16 	if (time_preference_1 == 2 & time_preference_2  == 1 & time_preference_10 == 1 & time_preference_14 == 1 & time_preference_16 == 1)
	replace time_preference_fe = 15 	if (time_preference_1 == 2 & time_preference_2  == 1 & time_preference_10 == 1 & time_preference_14 == 1 & time_preference_16 == 2)
	replace time_preference_fe = 14 	if (time_preference_1 == 2 & time_preference_2  == 1 & time_preference_10 == 1 & time_preference_14 == 2 & time_preference_15 == 1)
	replace time_preference_fe = 13 	if (time_preference_1 == 2 & time_preference_2  == 1 & time_preference_10 == 1 & time_preference_14 == 2 & time_preference_15 == 2)
	replace time_preference_fe = 12 	if (time_preference_1 == 2 & time_preference_2  == 1 & time_preference_10 == 2 & time_preference_11 == 1 & time_preference_13 == 1)
	replace time_preference_fe = 11 	if (time_preference_1 == 2 & time_preference_2  == 1 & time_preference_10 == 2 & time_preference_11 == 1 & time_preference_13 == 2)
	replace time_preference_fe = 10 	if (time_preference_1 == 2 & time_preference_2  == 1 & time_preference_10 == 2 & time_preference_11 == 2 & time_preference_12 == 1)
	replace time_preference_fe = 9 		if (time_preference_1 == 2 & time_preference_2  == 1 & time_preference_10 == 2 & time_preference_11 == 2 & time_preference_12 == 2)
	replace time_preference_fe = 8 		if (time_preference_1 == 2 & time_preference_2  == 2 & time_preference_3  == 1 & time_preference_7  == 1 & time_preference_8  == 1)
	replace time_preference_fe = 7 		if (time_preference_1 == 2 & time_preference_2  == 2 & time_preference_3  == 1 & time_preference_7  == 1 & time_preference_8  == 2)
	replace time_preference_fe = 6 		if (time_preference_1 == 2 & time_preference_2  == 2 & time_preference_3  == 1 & time_preference_7  == 2 & time_preference_9  == 1)
	replace time_preference_fe = 5 		if (time_preference_1 == 2 & time_preference_2  == 2 & time_preference_3  == 1 & time_preference_7  == 2 & time_preference_9  == 2)
	replace time_preference_fe = 4 		if (time_preference_1 == 2 & time_preference_2  == 2 & time_preference_3  == 2 & time_preference_4  == 1 & time_preference_6  == 1)
	replace time_preference_fe = 3		if (time_preference_1 == 2 & time_preference_2  == 2 & time_preference_3  == 2 & time_preference_4  == 1 & time_preference_6  == 2)
	replace time_preference_fe = 2		if (time_preference_1 == 2 & time_preference_2  == 2 & time_preference_3  == 2 & time_preference_4  == 2 & time_preference_5  == 1)
	replace time_preference_fe = 1		if (time_preference_1 == 2 & time_preference_2  == 2 & time_preference_3  == 2 & time_preference_4  == 2 & time_preference_5  == 2)
		
	
	
* 	demographic and background characteristics

	gen ethnicity_han = (ethnicity == 1) if ethnicity != .
	gen hukou_urban = (hukou == 1) if hukou != .
	gen siblings_total = number_bro_younger + number_bro_older + number_sis_younger + number_sis_older
	gen father_edu_hsabove = (father_education > 3) if father_education != .
	gen mother_edu_hsabove = (mother_education > 3) if mother_education != .
	gen religion_religious = (religion > 1) if religion != .
	gen hs_track_science = (hs_track == 1) 	if hs_track != .

	// indicator of university
	gen university_elite = (university == 1) if university != .
	
	// recode birth_year
	recode birth_year (1=1990) (2=1991) (3=1992) (4=1993) (5=1994) (6=1995) (7=1996) (8=1997) (9=1998) (10=1999) (11=2000) (12=2001) (13=2002) (14=2003) (15=2004) (16=2005)

	// indicator of coastal provinces
	gen birthplace_coastal = (birthplace_province == "Beijing" | birthplace_province == "Fujian" | birthplace_province == "Guangdong" | birthplace_province == "Hainan" | birthplace_province == "Hebei" | birthplace_province == "Jiangsu" | birthplace_province == "Shandong" | birthplace_province == "Shanghai" | birthplace_province == "Tianjin" | birthplace_province == "Zhejiang" | birthplace_province == "Non-mainland")
	gen residence_coastal = (residence_province == "Beijing" | residence_province == "Fujian" | residence_province == "Guangdong" | residence_province == "Hainan" | residence_province == "Hebei" | residence_province == "Jiangsu" | residence_province == "Shandong" | residence_province == "Shanghai" | residence_province == "Tianjin" | residence_province == "Zhejiang" | residence_province == "Non-mainland")
	
	// indicator of english credentials
	gen domestic_english_atleast4 = (english_qual_domestic > 1) if english_qual_domestic != .
	gen foreign_english_yes = (english_qual_foreign < 3) if english_qual_foreign != .
	
	// indicator of parents' work sector
	gen work_father_govt = (father_work <= 4) if father_work != .
	gen work_mother_govt = (mother_work <= 4) if mother_work != .
	
	// indicator of travel experience
	gen travel_foreign_yes = (travel_foreign > 1) if travel_foreign != .
	
	
	
* 	treatment related indicators

	// indicators and dummies
	gen treatment_control 		= (treatment_newsletter == 0 & treatment_vpn == 0 & vpn_current_paid_user == 0)
	gen treatment_vpnonly 		= (treatment_newsletter == 0 & treatment_vpn == 1 & vpn_current_paid_user == 0)
	gen treatment_nlonly 		= (treatment_newsletter == 1 & treatment_vpn == 0 & vpn_current_paid_user == 0)
	gen treatment_vpnnl  		= (treatment_newsletter == 1 & treatment_vpn == 1 & vpn_current_paid_user == 0)
	gen treatment_user 			= (vpn_current_paid_user == 1)
	
	gen treatment_master = .
	replace treatment_master = 1 	if treatment_control == 1
	replace treatment_master = 2 	if treatment_vpnonly == 1
	replace treatment_master = 3 	if treatment_nlonly == 1
	replace treatment_master = 4 	if treatment_vpnnl == 1
	replace treatment_master = 5 	if treatment_user == 1
	label def treatment_master 1 "Group-C" 2 "Group-A" 3 "Group-CE" 4 "Group-AE" 5 "Existing users"
	label value treatment_master treatment_master
	
	gen treatment_main = .
	replace treatment_main = 1 	if (treatment_control == 1 | treatment_nlonly == 1)
	replace treatment_main = 2 	if treatment_vpnonly == 1
	replace treatment_main = 3 	if treatment_vpnnl == 1
	replace treatment_main = 4 	if treatment_user == 1
	label def treatment_main 1 "Control" 2 "Access" 3 "Access + Encour." 4 "Existing users"
	label value treatment_main treatment_main

	// treatment interactions	
	gen treatvpnXnl = treatment_vpn * treatment_newsletter

	
	
* 	vpn adoption status
	gen vpn_adopted = (vpn_date_adoption != .)	if vpn_date_adoption != .
	replace vpn_adopted = 0 					if vpn_adopted != 1 & (treatment_master == 2 | treatment_master == 4)

		
* 	roommates' vpn usage

	// generate social learning primitives
	gen vpn_roommate_existing = vpn_roommate_usage_w1
	gen vpn_roommate_new_w2 = vpn_roommate_usage_w2 - vpn_roommate_usage_w1
	gen vpn_roommate_new_w3 = vpn_roommate_usage_w3 - vpn_roommate_usage_w1
	
	// top code at 2
	replace vpn_roommate_existing = 2 		if vpn_roommate_existing > 2
	replace vpn_roommate_new_w2 = 0 		if vpn_roommate_new_w2 < 0
	replace vpn_roommate_new_w2 = 2 		if vpn_roommate_new_w2 > 2
	replace vpn_roommate_new_w3 = 0 		if vpn_roommate_new_w3 < 0
	replace vpn_roommate_new_w3 = 2 		if vpn_roommate_new_w3 > 2
	
	// finalize primitives
	gen soclearning_ownaccess = (treatment_master >= 4)		if treatment_master != .
	gen soclearning_rm_new_w2 = (vpn_roommate_new_w2) 		if vpn_roommate_new_w2 != .
	gen soclearning_rm_new_w3 = (vpn_roommate_new_w3) 		if vpn_roommate_new_w3 != .
	gen soclearning_rm_existing = (vpn_roommate_existing) 	if vpn_roommate_existing != .
	gen soclearning_ownXnew_w2 = soclearning_ownaccess * soclearning_rm_new_w2
	gen soclearning_ownXnew_w3 = soclearning_ownaccess * soclearning_rm_new_w3
		
		
		
* 	list experiment

	// list experiment: indicator for provision of veil
	gen list_exp_veiled = 1 - list_exp_direct

	// list experiment: calculated overall count for direct group
	foreach c in trust {
		gen list_countd_`c'_direct = list_`c'_direct_count + list_`c'_direct_yes
		gen list_count_`c' = .
		replace list_count_`c' = list_countd_`c'_direct 	if list_exp_direct == 1
		replace list_count_`c' = list_`c'_veiled_count 		if list_exp_direct == 0
		}
	

	
			
*** Prepare data: generate az-scores

* 	A: Beliefs and attitudes regarding media
	
	// A.2: Purchase of censorship circumvention tools
	andersonz `var_vpn_purchase_reg_w3', gen(az_var_vpnpurchase)

	// A.3: Valuation of access to foreign media outlets
	andersonz `var_media_valuation_reg_w1', gen(az_belief_media_value_w1)
	andersonz `var_media_valuation_reg_w2', gen(az_belief_media_value_w2)
	andersonz `var_media_valuation_reg_w3', gen(az_belief_media_value_w3)
	
	// A.4: Trust in media outlets
	andersonz `var_media_trust_reg_w1', gen(az_belief_media_trust_w1)
	andersonz `var_media_trust_reg_w2', gen(az_belief_media_trust_w2)
	andersonz `var_media_trust_reg_w3', gen(az_belief_media_trust_w3)
	
	// A.5: Calibration of news outlets' level of censorship and biases
	andersonz `var_percmediabias_reg_ce_cn_w2', gen(az_belief_media_cens_cn)
	andersonz `var_percmediabias_reg_ce_us_w2', gen(az_belief_media_cens_us)
	andersonz `var_percmediabias_reg_di_cn_w2', gen(az_belief_media_bias_cn)
	andersonz `var_percmediabias_reg_di_us_w2', gen(az_belief_media_bias_us)
	
	// A.6: Justification of media censorship
	andersonz `var_censor_justif_reg_w1', gen(az_belief_media_justif_w1)
	andersonz `var_censor_justif_reg_w2', gen(az_belief_media_justif_w2)
	andersonz `var_censor_justif_reg_w3', gen(az_belief_media_justif_w3)
	andersonz `var_censor_justif_reg_fp_w1', gen(az_belief_media_justif_pf_w1)
	andersonz `var_censor_justif_reg_fp_w2', gen(az_belief_media_justif_pf_w2)
	andersonz `var_censor_justif_reg_fp_w3', gen(az_belief_media_justif_pf_w3)
	
	// A.7: Belief regarding drivers of media censorship
	andersonz `var_censor_driver_reg_dom_w2', gen(az_belief_media_dri_dom_w2)
	andersonz `var_censor_driver_reg_for_w2', gen(az_belief_media_dri_for_w2)
	
	
* 	B: Knowledge

	// B.1: Current news events covered in the demand treatment
	andersonz `var_knowledge_news_reg_quiz', gen(az_knowledge_news_quiz)
	
	// B.2: Current news events not covered in the demand treatment
	andersonz `var_knowledge_news_reg_cen_w2', gen(az_knowledge_news_censored)
	andersonz `var_knowledge_news_reg_unc_w2', gen(az_knowledge_news_uncensor)
	andersonz `var_knowledge_news_reg_cen_w1', gen(az_knowledge_news_censored_w1)
	andersonz `var_knowledge_news_reg_unc_w1', gen(az_knowledge_news_uncensor_w1)
	
	// B.3: Awareness of notable figures
	andersonz `var_knowledge_people_reg_tocens', gen(az_knowledge_people_tocens)
	andersonz `var_knowledge_people_reg_censor', gen(az_knowledge_people_censor)
	andersonz `var_knowledge_people_reg_uncens', gen(az_knowledge_people_uncens)
	andersonz `var_knowledge_people_censor_w1', gen(az_knowledge_people_censor_w1)
	andersonz `var_knowledge_people_uncens_w1', gen(az_knowledge_people_uncens_w1)

	// B.4: Awareness of protest events
	andersonz `var_knowledge_prot_reg_chi_w2', gen(az_knowledge_protest_china)
	andersonz `var_knowledge_prot_reg_for_w2', gen(az_knowledge_protest_forei)
	andersonz `var_knowledge_prot_reg_chi_w1', gen(az_knowledge_protest_china_w1)
	
	// B.5: Meta-knowledge
	andersonz `var_knowledge_meta_reg_w1', gen(az_knowledge_meta_w1)
	andersonz `var_knowledge_meta_reg_w2', gen(az_knowledge_meta_w2)
	andersonz `var_knowledge_meta_reg_w3', gen(az_knowledge_meta_w3)
	

* 	C: Economic beliefs
	
	// C.1: Belief on economic performance in China
	andersonz `var_econ_guess_reg_cn_perf_w1', gen(az_belief_econ_perf_cn_w1)
	andersonz `var_econ_guess_reg_cn_perf_w2', gen(az_belief_econ_perf_cn_w2)
	andersonz `var_econ_guess_reg_cn_perf_w3', gen(az_belief_econ_perf_cn_w3)

	// C.2: Confidence on guesses regarding economic performance in China
	andersonz `var_econ_guess_reg_cn_conf_w1', gen(az_belief_econ_conf_cn_w1)
	andersonz `var_econ_guess_reg_cn_conf_w2', gen(az_belief_econ_conf_cn_w2)
	andersonz `var_econ_guess_reg_cn_conf_w3', gen(az_belief_econ_conf_cn_w3)

	// C.3: Belief on economic performance in the US
	andersonz `var_econ_guess_reg_us_perf_w2', gen(az_belief_econ_perf_us_w2)
	andersonz `var_econ_guess_reg_us_perf_w3', gen(az_belief_econ_perf_us_w3)
	
	// C.4: Confidence on guesses regarding economic performance in the US
	andersonz `var_econ_guess_reg_us_conf_w2', gen(az_belief_econ_conf_us_w2)
	andersonz `var_econ_guess_reg_us_conf_w3', gen(az_belief_econ_conf_us_w3)
	
	
* 	D: Political attitudes

	// D.1: Demand for institutional change
	andersonz `var_demand_change_reg_w1', gen(az_belief_instchange_w1)
	andersonz `var_demand_change_reg_w2', gen(az_belief_instchange_w2)
	andersonz `var_demand_change_reg_w3', gen(az_belief_instchange_w3)
				
	// D.2: Trust in institutions
	andersonz `var_trust_inst_reg_govt_w1', gen(az_belief_trust_govt_w1)
	andersonz `var_trust_inst_reg_govt_w2', gen(az_belief_trust_govt_w2)
	andersonz `var_trust_inst_reg_govt_w3', gen(az_belief_trust_govt_w3)
	andersonz `var_trust_inst_reg_foreign_w1', gen(az_belief_trust_foreign_w1)
	andersonz `var_trust_inst_reg_foreign_w2', gen(az_belief_trust_foreign_w2)
	andersonz `var_trust_inst_reg_foreign_w3', gen(az_belief_trust_foreign_w3)
	andersonz `var_trust_inst_reg_copo_w2', gen(az_belief_trust_copo_w2)
	
	// D.3: Evaluation of government???s performance
	andersonz `var_eval_govt_reg_w1', gen(az_belief_evalgovt_w1)
	andersonz `var_eval_govt_reg_w2', gen(az_belief_evalgovt_w2)
	andersonz `var_eval_govt_reg_w3', gen(az_belief_evalgovt_w3)

	// D.4: Performance evaluation criteria
	andersonz `var_eval_criteria_reg_w2', gen(az_belief_evalcrit_w2)
	
	// D.5: Evaluation of severity of socioeconomic issues
	andersonz `var_severity_reg_w2', gen(az_belief_severity_w2)

	// D.6: Evaluation of democracy and human rights protection in China
	andersonz `var_democracy_reg_fp_w1', gen(az_belief_democracy_fp_w1)
	andersonz `var_democracy_reg_fp_w2', gen(az_belief_democracy_fp_w2)
	andersonz `var_democracy_reg_fp_w3', gen(az_belief_democracy_fp_w3)
	andersonz `var_democracy_reg_w2', gen(az_belief_democracy_w2)
	
	// D.7: Justification of controversial policies and issues
	andersonz `var_contro_justi_reg_policy_w2', gen(az_belief_justify_policy_w2)
	andersonz `var_contro_justi_reg_liberal_w2', gen(az_belief_justify_liberal_w2)
	
	// D.8: Willingness to act
	andersonz `var_willing_fight_reg_w1', gen(az_belief_willing_w1)
	andersonz `var_willing_fight_reg_w2', gen(az_belief_willing_w2)
	andersonz `var_willing_fight_reg_w3', gen(az_belief_willing_w3)

	// D.9: Interest in politics and economics
	andersonz `var_interest_reg_w2', gen(az_belief_interest_w2)
	
	
	
* 	E: behaviors

	// E.1: Social interactions
	andersonz `var_socialinteract_reg_w1', gen(az_var_socialinteract_w1)
	andersonz `var_socialinteract_reg_w2', gen(az_var_socialinteract_w2)
	andersonz `var_socialinteract_reg_w3', gen(az_var_socialinteract_w3)	
	
	// E.2: Participation
	andersonz `var_polparticipation_reg_pf_w1', gen(az_var_polparticipation_w1)
	andersonz `var_polparticipation_reg_pf_w2', gen(az_var_polparticipation_w2)
	andersonz `var_polparticipation_reg_pf_w3', gen(az_var_polparticipation_w3)

	
	
* 	F: Demographics, background characteristics, and fundamental preferences
	
	// F.1: Personal characteristics
	andersonz `var_demog_reg_personal', gen(az_demographics_personal)

	// F.2: Educational background
	andersonz `var_demog_reg_education', gen(az_demographics_education)
	
	// F.3: English ability and oversea travel experiences
	andersonz `var_demog_reg_english', gen(az_demographics_english)
	andersonz `var_demog_reg_travel', gen(az_demographics_travel)
	
	// F.4: Household characteristics
	andersonz `var_demog_reg_household', gen(az_demographics_household)
	
	// F.5: Fundamental preferences
	andersonz `var_preference_risk', gen(az_preference_risk)
	andersonz `var_preference_time', gen(az_preference_time)
	andersonz `var_preference_altruism', gen(az_preference_altruism)
	andersonz `var_preference_reciprocity', gen(az_preference_reciprocity)

	
	
* 	X: Overall effect index

	// flip variables to generate overall impact index
	forvalues i = 1/3 {
		gen bias_foreign_r_w`i' = 10 - bias_foreign_w`i'
		gen az_belief_econ_perf_cn_r_w`i' = 0 - az_belief_econ_perf_cn_w`i'
		gen az_belief_trust_govt_r_w`i' = 0 - az_belief_trust_govt_w`i'
		gen az_belief_evalgovt_r_w`i' = 0 - az_belief_evalgovt_w`i'
		gen importance_live_in_demo_r_w`i' = 10 - importance_live_in_demo_w`i'
		}

	andersonz `var_overall_w1', gen(az_overall_w1)
	andersonz `var_overall_w2', gen(az_overall_w2)
	andersonz `var_overall_w3', gen(az_overall_w3)
	
	// five main categories: wave 1
	andersonz info_freq_website_for_w1 az_belief_media_value_w1 az_belief_media_trust_w1 az_belief_media_justif_pf_w1 bias_domestic_w1 bias_foreign_r_w1, gen(az_overall_a_w1)
	andersonz news_perccor_cen_w1 news_perccor_unc_w1 protest_pcheard_china_w1 protest_pcheard_foreign_w1 az_knowledge_meta_w1, gen(az_overall_b_w1)
	andersonz az_belief_econ_perf_cn_r_w1 az_belief_econ_conf_cn_w1, gen(az_overall_c_w1)
	andersonz az_belief_instchange_w1 az_belief_trust_govt_r_w1 az_belief_trust_foreign_w1 az_belief_evalgovt_r_w1 importance_live_in_demo_r_w1 az_belief_willing_w1, gen(az_overall_d_w1)
	andersonz az_var_socialinteract_w1 az_var_polparticipation_w1 plan_grad_foreignmaster_w1 cloc_for_w1 stock_participation_w1, gen(az_overall_e_w1)
	
	// five main categories: wave 2
	andersonz info_freq_website_for_w2 az_belief_media_value_w2 az_belief_media_trust_w2 az_belief_media_justif_pf_w2 bias_domestic_w2 bias_foreign_r_w2, gen(az_overall_a_w2)
	andersonz news_perccor_cen_w2 news_perccor_unc_w2 protest_pcheard_china_w2 protest_pcheard_foreign_w2 az_knowledge_meta_w2, gen(az_overall_b_w2)
	andersonz az_belief_econ_perf_cn_r_w2 az_belief_econ_conf_cn_w2, gen(az_overall_c_w2)
	andersonz az_belief_instchange_w2 az_belief_trust_govt_r_w2 az_belief_trust_foreign_w2 az_belief_evalgovt_r_w2 importance_live_in_demo_r_w2 az_belief_willing_w2, gen(az_overall_d_w2)
	andersonz az_var_socialinteract_w2 az_var_polparticipation_w2 plan_grad_foreignmaster_w2 cloc_for_w2 stock_participation_w2, gen(az_overall_e_w2)

	// five main categories: wave 3
	andersonz info_freq_website_for_w3 az_belief_media_value_w3 az_belief_media_trust_w3 az_belief_media_justif_pf_w3 bias_domestic_w3 bias_foreign_r_w3, gen(az_overall_a_w3)
	andersonz news_perccor_cen_w3 news_perccor_unc_w3 protest_pcheard_china_w3 protest_pcheard_foreign_w3 az_knowledge_meta_w3, gen(az_overall_b_w3)
	andersonz az_belief_econ_perf_cn_r_w3 az_belief_econ_conf_cn_w3, gen(az_overall_c_w3)
	andersonz az_belief_instchange_w3 az_belief_trust_govt_r_w3 az_belief_trust_foreign_w3 az_belief_evalgovt_r_w3 importance_live_in_demo_r_w3 az_belief_willing_w3, gen(az_overall_d_w3)
	andersonz az_var_socialinteract_w3 az_var_polparticipation_w3 plan_grad_foreignmaster_w3 cloc_for_w3 stock_participation_w3, gen(az_overall_e_w3)
	
	// mega z-scores: wave 3
	andersonz az_overall_a_w3 az_overall_b_w3 az_overall_c_w3 az_overall_d_w3 az_overall_e_w3, gen(az_overall_all_w3)
	
		
	
	
*** Prepare data: generate heterogeneity cuts

	// demographics (category F)
	gen h_gender = (gender == 1) 							if gender != .
	gen h_birth_year = (birth_year < 1996) 					if birth_year != .
	gen h_ethnicity_han = (ethnicity_han == 1) 				if ethnicity_han != .
	gen h_birthplace_coastal = (birthplace_coastal == 1) 	if birthplace_coastal != .
	gen h_residence_coastal = (residence_coastal == 1) 		if residence_coastal != .
	gen h_hukou_urban = (hukou_urban == 1) 					if hukou_urban != .
	gen h_religion_religious = (religion_religious == 1) 	if religion_religious != .
	gen h_ccp_member = (ccp_member == 1) 					if ccp_member != .
	gen h_university_elite = (university_elite == 1) 		if university_elite != .
	gen h_hs_track_science = (hs_track_science == 1) 		if hs_track_science != .
	gen h_department_ssh = (department_ssh == 1) 			if department_ssh != .
	gen h_domestic_english_atleast4 = (domestic_english_atleast4 == 1) 	if domestic_english_atleast4 != .
	gen h_foreign_english_yes = (foreign_english_yes == 1) 	if foreign_english_yes != .
	gen h_travel_hktaiwan = (travel_hktaiwan == 1) 			if travel_hktaiwan != .
	gen h_travel_foreign_yes = (travel_foreign_yes == 1) 	if travel_foreign_yes != .
	gen h_siblings_total = (siblings_total > 0) 			if siblings_total != .
	gen h_father_edu_hsabove = (father_edu_hsabove == 1) 	if father_edu_hsabove != .
	gen h_work_father_govt = (work_father_govt == 1) 		if work_father_govt != .
	gen h_father_ccp = (father_ccp == 1) 					if father_ccp != .
	gen h_mother_edu_hsabove = (mother_edu_hsabove == 1) 	if mother_edu_hsabove != .
	gen h_work_mother_govt = (work_mother_govt == 1) 		if work_mother_govt != .
	gen h_mother_ccp = (mother_ccp == 1) 					if mother_ccp != .
	gen h_hh_income = (hh_income >= 75000) 					if hh_income != .
	gen h_az_preference_risk = (az_preference_risk > 0) 	if az_preference_risk != .
	gen h_az_preference_time = (az_preference_time > 0) 	if az_preference_risk != .
	gen h_az_preference_altruism = (az_preference_altruism > 0) 		if az_preference_altruism != .
	gen h_az_preference_reciprocity = (az_preference_reciprocity > 0) 	if az_preference_reciprocity != .
	
	// baseline outcomes
	gen h_az_overall_a_w1 = (az_overall_a_w1 > 0) 			if az_overall_a_w1 != .
	gen h_az_overall_b_w1 = (az_overall_b_w1 > 0) 			if az_overall_b_w1 != .
	gen h_az_overall_c_w1 = (az_overall_c_w1 > 0) 			if az_overall_c_w1 != .
	gen h_az_overall_d_w1 = (az_overall_d_w1 > 0) 			if az_overall_d_w1 != .
	gen h_az_overall_e_w1 = (az_overall_e_w1 > 0) 			if az_overall_e_w1 != .
	gen h_az_belief_media_value_w1 = (az_belief_media_value_w1 > 0) 		if az_belief_media_value_w1 != .
	gen h_az_belief_media_trust_w1 = (az_belief_media_trust_w1 > 0) 		if az_belief_media_trust_w1 != .
	gen h_az_knowledge_news_cens_w1 = (az_knowledge_news_censored_w1 > 0) 	if az_knowledge_news_censored_w1 != .
	gen h_az_knowledge_news_unce_w1 = (az_knowledge_news_uncensor_w1 > 0) 	if az_knowledge_news_uncensor_w1 != .
	gen h_az_knowledge_pp_censor_w1 = (az_knowledge_people_censor_w1 > 0) 	if az_knowledge_people_censor_w1 != .
	gen h_az_knowledge_pp_uncens_w1 = (az_knowledge_people_uncens_w1 > 0) 	if az_knowledge_people_uncens_w1 != .
	gen h_az_knowledge_pr_china_w1 = (az_knowledge_protest_china_w1 > 0) 	if az_knowledge_protest_china_w1 != .
	gen h_az_belief_trust_govt_w1 = (az_belief_trust_govt_w1 > 0) 			if az_belief_trust_govt_w1 != .
		
		
	
	
*** Prepare data: generate variable labels

	// A. Beleifs and attitudes regarding media
	local l_wtp_vpn 						`"WTP for uncensored Internet access ($/month)"'
	local l_added_value_foreign_media 		`"Value added of foreign media access"'
	local l_az_belief_media_value 			`"Valuation of access to foreign media outlets"'
	
	local l_trust_media_dom_state 			`"Distrust in domestic state-owned media"'
	local l_trust_media_dom_private 		`"Distrust in domestic privately-owned media"'
	local l_trust_media_foreign 			`"Trust in foreign media"'
	local l_az_belief_media_trust 			`"Trust in non-domestic media outlets"'
	
	local l_bias_domestic 					`"Degree of censorship on domestic news outlets"'
	local l_bias_foreign 					`"Degree of censorship on foreign news outlets"'
	local l_bias_dom_govt_policy_t1 		`"Domestic cens. driven by govt. policies"'
	local l_bias_dom_firm_interest_t1 		`"Domestic cens. driven by corp. interest"'
	local l_bias_dom_media_pref_t1 			`"Domestic cens. driven by media???s ideology"'
	local l_bias_dom_reader_demand_t1		`"Domestic cens. driven by readers??? demand"'
	local l_bias_for_govt_policy_t1 		`"Foreign cens. driven by govt. policies"'
	local l_bias_for_firm_interest_t1 		`"Foreign cens. driven by corp. interest"'
	local l_bias_for_media_pref_t1 			`"Foreign cens. driven by media???s ideology"'
	local l_bias_for_reader_demand_t1 		`"Foreign cens. driven by readers??? demand"'

	local l_distneutral_cn_neg_cn 			`"Bias: Chinese media on neg. news in China"'
	local l_distneutral_cn_pos_cn 			`"Bias: Chinese media on pos. news in China"'
	local l_distneutral_us_neg_cn 			`"Bias: Chinese media on neg. news in US"'
	local l_distneutral_us_pos_cn 			`"Bias: Chinese media on pos. news in US"'
	local l_distneutral_cn_neg_us 			`"Bias: US media on neg. news in China"'
	local l_distneutral_cn_pos_us 			`"Bias: US media on pos. news in China"'
	local l_distneutral_us_neg_us 			`"Bias: US media on neg. news in US"'
	local l_distneutral_us_pos_us 			`"Bias: US media on pos. news in US"'

	local l_bias_cn_neg_cn_cens 			`"Censorship: Chinese media on neg. news in China"'
	local l_bias_cn_pos_cn_cens 			`"Censorship: Chinese media on pos. news in China"'
	local l_bias_us_neg_cn_cens 			`"Censorship: Chinese media on neg. news in US"'
	local l_bias_us_pos_cn_cens 			`"Censorship: Chinese media on pos. news in US"'
	local l_bias_cn_neg_us_cens 			`"Censorship: US media on neg. news in China"'
	local l_bias_cn_pos_us_cens 			`"Censorship: US media on pos. news in China"'
	local l_bias_us_neg_us_cens 			`"Censorship: US media on neg. news in US"'
	local l_bias_us_pos_us_cens 			`"Censorship: US media on pos. news in US"'
	
	local l_censor_just_dom_economic 		`"Unjustified: censoring economic news"'
	local l_censor_just_dom_political 		`"Unjustified: censoring political news"'
	local l_censor_just_dom_social 			`"Unjustified: censoring social news"'
	local l_censor_just_for 				`"Unjustified: censoring foreign news"'
	local l_censor_just_porn 				`"Unjustified: censoring pornography"'
	local l_az_belief_media_justif_pf 		`"Censorship is unjustified"'
	
	
	// B. Knowledge
	local l_panamapapers 					`"Panama Papers"'
	local l_tenyearshk  					`"HK independence"'
	local l_stockcrash  					`"2016 stock mkt crash"'
	local l_economistcensor 				`"Censoring Economist"'
	local l_coalprod 						`"Steel prod. & pollution"'
	local l_trumpchina 						`"Trump trademark in China"'
	local l_xiaojianhua 					`"Jianhua Xiao kidnap"'
	local l_xijiangcar 						`"Tracking Xinjiang cars"'
	local l_chinanorway 					`"China-Norway relations"'
	local l_womenrights 					`"Women rights activisits"'
	local l_hkceelection 					`"HK CE election"'
	local l_news_perccor_cen 				`"% quizzes answered correctly: sensitive"'
	local l_news_perccor_unc 				`"% quizzes answered correctly: non-sensitive"'

	local l_people_puzhiqiang 				`"Aware of Zhiqiang Pu"'
	local l_people_lizehou 					`"Aware of Zehou Li"'
	local l_people_huangzhifeng 			`"Aware of Joshua Wong"'
	local l_people_chenguangcheng 			`"Aware of Guangcheng Cheng"'
	local l_people_lixiaolin 				`"Aware of Xiaolin Li"'
	local l_people_renzhiqiang 				`"Aware of Zhiqiang Ren"'
	local l_people_maoyushi 				`"Aware of Yushi Mao"'
	local l_people_honghuang 				`"Aware of Huang Hong"'
	local l_people_liuqiangdong 			`"Aware of Qiangdong Liu"'
	local l_people_jialequn 				`"Aware of Lequn Jia"'
	
	local l_protest_2014_europe_square 		`"Aware of 2014 Ukrainian Euromaidan Revolution"'
	local l_protest_2014_sun_flower 		`"Aware of 2014 Taiwan Sunflower Stud. Movement"'
	local l_protest_2010_arabic_spring 		`"Aware of 2010 Arab Spring"'
	local l_protest_2014_crimea_vote 		`"Aware of 2014 Crimean Status Referendum"'
	local l_protest_2012_hk_curriculum 		`"Aware of 2012 HK Anti-National Curr. Movement"'
	local l_protest_2010_catal_indep 		`"Aware of 2010 Catalonian Indep. Movement"'
	local l_protest_2014_umbrella 			`"Aware of 2014 HK Umbrella Revolution"'
	local l_protest_2011_tmrw_parade 		`"Aware of 2011 Tomorrow Revolution"'
	local l_protest_2016_mongkok_riot 		`"Aware of 2016 HK Mong Kok Revolution"'
	local l_protest_pcheard_china 			`"Awareness of protests in Greater China"'
	local l_protest_pcheard_foreign			`"Awareness of foreign protests"'
	
	local l_familiar_china_issues_self 		`"Informedness of issues in China"'
	local l_familiar_china_others 			`"Greater informedness than peers"'
	local l_az_knowledge_meta 				`"Self-assessment of knowledge level"'
	
	
	// C. Economic beliefs
	local l_guess_gdp_growth_china 			`"Guess on GDP growth rate in China"'
	local l_guess_stock_index_sh 			`"Guess on year-end SSCI"'
	local l_az_belief_econ_perf_cn 			`"Optimistic belief of Chinese economy"'
	local l_az_belief_econ_conf_cn 			`"Confidence of guesses on Chinese economy"'
	
	local l_guess_gdp_growth_china_con 		`"Confidence of China GDP guess"'
	local l_guess_stock_index_sh_con 		`"Confidence of SSCI guess"'

	local l_guess_gdp_growth_us 			`"Guess on GDP growth rate in US"'
	local l_guess_stock_index_dj			`"Guess on year-end DJI"'

	local l_guess_gdp_growth_us_con 		`"Confidence of US GDP guess"'
	local l_guess_stock_index_dj_con 		`"Confidence of DJI guess"'

	
	// D. Political attitudes
	local l_inst_change_econ 				`"Economic institution needs changes"'
	local l_inst_change_poli 				`"Political institution needs changes"'
	local l_az_belief_instchange 			`"Demand for institutional change"'
	
	local l_trust_central_govt 				`"Trust in central govt. of China"'
	local l_trust_provincial_govt 			`"Trust in provincial govt. of China"'
	local l_trust_local_govt 				`"Trust in local govt. of China"'
	local l_az_belief_trust_govt 			`"Trust in Chinese govt."'
	local l_trust_japan_govt 				`"Trust in central govt. of Japan"'
	local l_trust_us_govt 					`"Trust in federal govt. of US"'
	local l_az_belief_trust_foreign 		`"Trust in foreign govt."'
	local l_trust_court 					`"Trust in court"'
	local l_trust_police 					`"Trust in police"'
	local l_trust_ngo 						`"Trust in NGOs"'
	local l_trust_financial_domestic 		`"Trust in domestic financial inst."'
	local l_trust_financial_foreign 		`"Trust in foreign financial inst."'
	
	local l_eval_govt_economic 				`"Satisfaction of economic dev."'
	local l_eval_govt_dom_politics 			`"Satisfaction of domestic politics"'
	local l_eval_govt_for_relations 		`"Satisfaction of diplomatic affairs"'
	local l_az_belief_evalgovt 				`"Satisfaction of govt. performance"'
	
	local l_revalgovt_election 				`"Eval. importance: universal suffrage"'
	local l_revalgovt_economy 				`"Eval. importance: economic dev."'
	local l_revalgovt_equality 				`"Eval. importance: income and wealth equality"'
	local l_revalgovt_ruleoflaw 			`"Eval. importance: rule of law"'
	local l_revalgovt_human_rights 			`"Eval. importance: civil and human rights"'
	local l_revalgovt_freedom_speech 		`"Eval. importance: freedom of speech"'
	local l_revalgovt_global_power 			`"Eval. importance: intl. affairs"'
	local l_revalgovt_fair_history 			`"Eval. importance: handle history fairly"'

	local l_severity_welfare 				`"Severity: social security and welfare"'
	local l_severity_employment 			`"Severity: employments"'
	local l_severity_pollution 				`"Severity: environmental pollution"'
	local l_severity_inequality 			`"Severity: wealth inequality"'
	local l_severity_corruption 			`"Severity: govt. corruption"'
	local l_severity_dscrm_minority 		`"Severity: minority discrimination"'

	local l_importance_live_in_demo 		`"Living in democracy is not important"'
	local l_china_interest_group 			`"China cares interest for masses"'
	local l_china_rate_democracy 			`"Level of democracy in China"'
	local l_china_rate_humanrights 			`"Level of human rights protection"'
	local l_az_belief_democracy_fp 			`"Living in democracy is not important"'
	
	local l_justify_minority_policy 		`"Justified: minority policies"'
	local l_justify_reduce_pollution 		`"Justified: prod. cut to reduce pollution"'
	local l_justify_hukou 					`"Justified: migration restrictions"'
	local l_justify_one_child 				`"Justified: one-child policy"'
	local l_justify_gaokao 					`"Justified: college admission policies"'
	local l_justify_hongkong_policy 		`"Justified: policy towards HK"'
	local l_justify_taiwan_policy 			`"Justified: policy towards Taiwan"'
	local l_justify_transgene 				`"Justified: transgenetic food"'
	local l_justify_receive_refugee 		`"Justified: refusal of DPRK refugees"'
	local l_justify_soe_privatize 			`"Justified: privatization of SOEs"'
	local l_justify_homo_marriage 			`"Justified: legal. of homosexual marriages"'
	local l_justify_legal_prostitute 		`"Justified: legal. of prostitution"'
	local l_justify_abortion 				`"Justified: abortion"'
	local l_justify_soft_drug 				`"Justified: soft drugs usage"'
	local l_justify_violence_stability 		`"Justified: govt. use of violence"'

	local l_willing_against_illi_govt 		`"Willing to battle illegal govt. acts"'
	local l_willing_report_mis 				`"Willing to report govt. misconduct"'
	local l_willing_protect_weak 			`"Willing to stand up for the weak"'
	local l_az_belief_willing 				`"Willingness to act"'
	
	local l_interest_economic 				`"Interest in economics"'
	local l_interest_politics 				`"Interest in politics"'

	local l_proud_being_chinese 			`"Proud of being Chinese"'

	local l_fear_critic_govt_self 			`"Fear to criticize govt. in public"'
	
	
	// E. Behaviors
	local l_info_domestic_website 			`"Important info source: domestic websites"'
	local l_info_foreign_website 			`"Important info source: foreign websites"'
	local l_info_social_media_dom 			`"Important info source: domestic social media"'
	local l_info_word_of_mouth 				`"Important info source: word of mouth"'
	local l_info_social_media_for 			`"Important info source: foreign social media"'

	local l_info_freq_website_for 			`"Frequently visit foreign websites for info"'

	local l_participate_social_protest 		`"Participated in social protests"'
	local l_participate_plan_vote 			`"Plan to vote in PCR election"'
	local l_participate_complain_school 	`"Filed complaints to school"'
	local l_participate_ngo 				`"Participated in NGO activities"'
	local l_az_var_polparticipation 		`"Participation behaviors"'
	
	local l_frequency_talk_politic 			`"Frequency of discussing poli. with friends"'
	local l_frequency_persuade_friends 		`"Frequency of persuading others"'
	local l_az_var_socialinteract 			`"Social interaction in politics"'
	
	local l_plan_grad_gradschool_dom 		`"Plan: grad. school in China"'
	local l_plan_grad_foreignmaster 		`"Plan: master degree abroad"'
	local l_plan_grad_foreignphd 			`"Plan: PhD degree abroad"'
	local l_plan_grad_military 				`"Plan: military in China"'
	local l_plan_grad_work 					`"Plan: work right away"'

	local l_cp_t1_national_civil 			`"Sector pref.: national civil service"'
	local l_cp_t1_local_civil 				`"Sector pref.: local civil service"'
	local l_cp_t1_military 					`"Sector pref.: military"'
	local l_cp_t1_chinese_private 			`"Sector pref.: private firm in China"'
	local l_cp_t1_for_firm 					`"Sector pref.: private firm in China"'
	local l_cp_t1_soe 						`"Sector pref.: SOEs"'
	local l_cp_t1_institutional 			`"Sector pref.: inst. organizations"'
	local l_cp_t1_entrepreneur 				`"Sector pref.: entrepreneurship"'

	local l_cloc_beijing 					`"Location pref.: Beijing"'
	local l_cloc_shanghai 					`"Location pref.: Shanghai"'
	local l_cloc_gzsz 						`"Location pref.: Shenzhen/Guangzhou"'
	local l_cloc_tjcq 						`"Location pref.: tier 2 cities in central"'
	local l_cloc_hkmc 						`"Location pref.: HK and Macau"'
	local l_cloc_taiwan 					`"Location pref.: Taiwan"'
	local l_cloc_dom 						`"Location pref.: other cities in China"'
	local l_cloc_for 						`"Location pref.: foreign cities"'
	local l_stock_participation 			`"Currently invested in Chinese stock mkt."'

/* The changes start */	
	
****************************
**  Code for Table A.13   ** 		
****************************

	// keep wave 3 non-existing users only
	keep if panelmerged_wave3 == 1
	keep if treatment_user != 1
	
	// transform into dummy indicators: uncensored belief - above median
	foreach Y in info_foreign_website_w3 info_freq_website_for_w3 az_belief_media_value_w3 az_belief_media_trust_w3 bias_domestic_w3 az_belief_media_justif_w3 news_perccor_cen_w3 news_perccor_unc_w3 protest_pcheard_china_w3 protest_pcheard_foreign_w3 az_knowledge_meta_w3 az_belief_econ_conf_cn_w3 az_belief_econ_perf_us_w3 az_belief_econ_conf_us_w3 az_belief_instchange_w3 az_belief_trust_foreign_w3 az_belief_willing_w3 frequency_talk_politic_w3 frequency_persuade_friends_w3 {
		egen `Y'_m = median(`Y')
		gen `Y'_p = (`Y' > `Y'_m)
		}

	foreach Y in info_foreign_website_w1 info_freq_website_for_w1 az_belief_media_value_w1 az_belief_media_trust_w1 bias_domestic_w1 az_belief_media_justif_w1 news_perccor_cen_w1 news_perccor_unc_w1 protest_pcheard_china_w1 protest_pcheard_foreign_w1 az_knowledge_meta_w1 az_belief_econ_conf_cn_w1 az_belief_instchange_w1 az_belief_trust_foreign_w1 az_belief_willing_w1 frequency_talk_politic_w1 frequency_persuade_friends_w1 {
		egen `Y'_m = median(`Y')
		gen `Y'_p = (`Y' > `Y'_m)
		}

	// transform into dummy indicators: uncensored belief - below median
	foreach Y in bias_foreign_w3 importance_live_in_demo_w3 az_belief_econ_perf_cn_w3 az_belief_trust_govt_w3 az_belief_evalgovt_w3 {
		egen `Y'_m = median(`Y')
		gen `Y'_p = (`Y' < `Y'_m)
		}

	foreach Y in bias_foreign_w1 importance_live_in_demo_w1 az_belief_econ_perf_cn_w1 az_belief_trust_govt_w1 az_belief_evalgovt_w1 {
		egen `Y'_m = median(`Y')
		gen `Y'_p = (`Y' < `Y'_m)
		}

	// already dummy, keep as it is
	rename participate_complain_school_w1 par_complain_school_w1
	rename participate_complain_school_w3 par_complain_school_w3

	foreach Y in vpn_purchase_wmt_record vpn_purchase_yes bias_dom_govt_policy_t1_w3 protest_2011_tmrw_parade_w3 participate_social_protest_w3 participate_plan_vote_w3 par_complain_school_w3 plan_grad_foreignmaster_w3 cp_t3_for_firm_w3 cloc_for_w3 {
		gen `Y'_p = `Y'
		}

	foreach Y in bias_dom_govt_policy_t1_w1 protest_2011_tmrw_parade_w1 participate_social_protest_w1 participate_plan_vote_w1 par_complain_school_w1 plan_grad_foreignmaster_w1 cp_t3_for_firm_w1 cloc_for_w1 {
		gen `Y'_p = `Y'
		}

	// already dummy, flip sign
	foreach Y in bias_for_govt_policy_t1_w3 stock_participation_w3 {
		gen `Y'_p = 1- `Y'
		}

	foreach Y in bias_for_govt_policy_t1_w1 stock_participation_w1 {
		gen `Y'_p = 1- `Y'
		}

keep treatment_main active_user treatment_control ///
	 treatment_vpnonly treatment_nlonly treatment_vpnnl /// 
	 info_foreign_website_w3_p ///	A.1.2 Ranked high: foreign websites
	 info_freq_website_for_w3_p ///	A.1.6 Freq. of visiting foreign websites for info.
     az_belief_media_value_w3_p ///	A.3 Valuation of access to foreign media outlets
	 az_belief_media_trust_w3_p ///	A.4 Trust in non-domestic media outlets
	 bias_domestic_w3_p ///			A.5.1 Degree of censorship on domestic news outlets	
	 bias_foreign_w3_p ///	A.5.2 Degree of censorship on foreign news outlets	
	 az_belief_media_justif_w3_p /// A.6 Censorship unjustified
	 bias_dom_govt_policy_t1_w3_p ///  A.7.1 Domestic cens. driven by govt. policies
	 bias_for_govt_policy_t1_w3_p /// A.7.2 Foreign cens. driven by govt. policies
	 info_foreign_website_w1_p ///	A.1.2 Ranked high: foreign websites
	 info_freq_website_for_w1_p ///	A.1.6 Freq. of visiting foreign websites for info.
	 az_belief_media_value_w1_p ///	A.3 Valuation of access to foreign media outlets
	 az_belief_media_trust_w1_p ///	A.4 Trust in non-domestic media outlets
	 bias_domestic_w1_p ///			A.5.1 Degree of censorship on domestic news outlets	
	 bias_foreign_w1_p ///	A.5.2 Degree of censorship on foreign news outlets	
	 az_belief_media_justif_w1_p /// A.6 Censorship unjustified
	 bias_dom_govt_policy_t1_w1_p ///  A.7.1 Domestic cens. driven by govt. policies
	 bias_for_govt_policy_t1_w1_p /// A.7.2 Foreign cens. driven by govt. policies				 
	 vpn_purchase_wmt_record_p ///	  A.2.1 Purchase discounted tool we offered
	 vpn_purchase_yes_p ///		  A.2.2 Purchase any tool
		
save ChenYang2019, replace
		
**************
** ANALYSES ** 		
**************


*** Set directory for output

*	cd "`directory_path'/Outregs"
	
/* The changes end */	

/*

*** Figure 2

	foreach Y in wtp_vpn {
	
		preserve
		keep if panelmerged_wave2 == 1 & panelmerged_wave3 == 1
		collapse (mean) `Y'_w3 `Y'_w2 `Y'_w1 (sd) s_`Y'_w3 = `Y'_w3 s_`Y'_w2 = `Y'_w2 s_`Y'_w1 = `Y'_w1 (count) n_`Y'_w3 = `Y'_w3 n_`Y'_w2 = `Y'_w2 n_`Y'_w1 = `Y'_w1, by(treatment_main)
		gen h_`Y'_w3 = `Y'_w3 + invttail(n_`Y'_w3-1,0.05)*(s_`Y'_w3 / sqrt(n_`Y'_w3))
		gen l_`Y'_w3 = `Y'_w3 - invttail(n_`Y'_w3-1,0.05)*(s_`Y'_w3 / sqrt(n_`Y'_w3))	
		gen h_`Y'_w2 = `Y'_w2 + invttail(n_`Y'_w2-1,0.05)*(s_`Y'_w2 / sqrt(n_`Y'_w2))
		gen l_`Y'_w2 = `Y'_w2 - invttail(n_`Y'_w2-1,0.05)*(s_`Y'_w2 / sqrt(n_`Y'_w2))	
		gen h_`Y'_w1 = `Y'_w1 + invttail(n_`Y'_w1-1,0.05)*(s_`Y'_w1 / sqrt(n_`Y'_w1))
		gen l_`Y'_w1 = `Y'_w1 - invttail(n_`Y'_w1-1,0.05)*(s_`Y'_w1 / sqrt(n_`Y'_w1))	
		drop s_* n_*
		drop if treatment_main == .
		
		reshape long `Y'_w h_`Y'_w l_`Y'_w, i(treatment_main) j(wave)
		replace wave = 4 if wave == 3
		
		twoway 	(connected `Y'_w wave if treatment_main == 1, msymbol(o) mcolor(gs9) lpattern(dash) lcolor(gs9) lwidth(medthick)) /// 
				(connected `Y'_w wave if treatment_main == 2, msymbol(d) mcolor(gs3) lpattern(dash) lcolor(gs3) lwidth(medthick)) /// 
				(connected `Y'_w wave if treatment_main == 3, msymbol(o) mcolor(cranberry) lpattern(line) lcolor(cranberry) lwidth(vthick)) /// 
				(connected `Y'_w wave if treatment_main == 4, msymbol(o) mcolor(emidblue) lpattern(line) lcolor(emidblue) lwidth(thick)) /// 
				(rcap h_`Y'_w l_`Y'_w wave), /// 
				xline(1.2) ///
				title("`l_`Y''", size(medium)) ///
				ytitle("") ///
				xtitle("") ///
				xscale(r(0.8 4.2)) ///
				xlabel(1 "Nov. 2015" 2 "Apr. 2016" 4 "May 2017", noticks) ///
				legend(order(1 "Control" 2 "Access" 3 "Access + Encour." 4 "Existing users") col(2) size(small)) ///
				graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
				saving(figure_panel_`Y'_s, replace)
				
		graph export "figure_panel_`Y'_s.pdf", replace
		restore
		
	}
	
		
		
		
*** Figure 3

	preserve
	collapse (mean) m_vpn_purchase_yes = vpn_purchase_yes m_vpn_purchase_wmt = vpn_purchase_wmt_record (sd) s_vpn_purchase_yes = vpn_purchase_yes s_vpn_purchase_wmt = vpn_purchase_wmt_record (count) n_vpn_purchase_yes = vpn_purchase_yes n_vpn_purchase_wmt = vpn_purchase_wmt_record, by(treatment_main)
	gen h_vpn_purchase_yes = m_vpn_purchase_yes + invttail(n_vpn_purchase_yes-1,0.05)*(s_vpn_purchase_yes/sqrt(n_vpn_purchase_yes))
	gen l_vpn_purchase_yes = m_vpn_purchase_yes - invttail(n_vpn_purchase_yes-1,0.05)*(s_vpn_purchase_yes/sqrt(n_vpn_purchase_yes))
	gen h_vpn_purchase_wmt = m_vpn_purchase_wmt + invttail(n_vpn_purchase_wmt-1,0.05)*(s_vpn_purchase_wmt/sqrt(n_vpn_purchase_wmt))
	gen l_vpn_purchase_wmt = m_vpn_purchase_wmt - invttail(n_vpn_purchase_wmt-1,0.05)*(s_vpn_purchase_wmt/sqrt(n_vpn_purchase_wmt))
	
	gen benchmark_1 = 0.65
	gen benchmark_2 = 0.47
	gen benchmark_x = treatment_main - 0.5
	
	twoway 	(bar m_vpn_purchase_yes treatment_main if treatment_main == 1, barwidth(0.75) color(gs11) fintensity(40)) ///
			(bar m_vpn_purchase_yes treatment_main if treatment_main == 2, barwidth(0.75) color(gs11) fintensity(40)) ///
			(bar m_vpn_purchase_yes treatment_main if treatment_main == 3, barwidth(0.75) color(cranberry) fintensity(45)) ///
			(bar m_vpn_purchase_yes treatment_main if treatment_main == 4, barwidth(0.75) color(emidblue) fintensity(35)) ///
			(bar m_vpn_purchase_wmt treatment_main if treatment_main == 1, barwidth(0.75) color(gs11) fintensity(80))  ///
			(bar m_vpn_purchase_wmt treatment_main if treatment_main == 2, barwidth(0.75) color(gs11) fintensity(80))  ///
			(bar m_vpn_purchase_wmt treatment_main if treatment_main == 3, barwidth(0.75) color(cranberry) fintensity(90))  ///
			(bar m_vpn_purchase_wmt treatment_main if treatment_main == 4, barwidth(0.75) color(emidblue) fintensity(70))  ///
			(line benchmark_1 benchmark_x if benchmark_x >= 2, lcolor(erose) lpattern(dash) lwidth(thick)) ///
			(line benchmark_2 benchmark_x if benchmark_x >= 1 & benchmark_x <= 3, lcolor(erose) lpattern(dash) lwidth(thick)) ///
			(rcap h_vpn_purchase_yes l_vpn_purchase_yes treatment_main), ///
			xlabel(1 "Control" 2 "Access" 3 "Access + Encour." 4 "Existing users", noticks) ///
			xtitle("") ///
			ytitle("% purchased any circumvention tool after April 2017") ///
			yscale(r(0 1)) ///
			ylabel(0 "0" 0.2 "20" 0.4 "40" 0.6 "60" 0.8 "80" 1 "100") ///
			xsize(5) ysize(3) ///
			legend(off) ///
			graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
			saving(figure_vpn_purchase_yes_2, replace)
	graph export "figure_vpn_purchase_yes_2.pdf", replace
			
	restore

	
		
*** Table 1
	
	preserve
	keep if panelmerged_wave3 == 1
	gen wave3 = 1
	drop treatment_master
	save panelmerged_wave3_temp, replace
	restore
	
	preserve
	
	// duplicate for attrition test
	append using panelmerged_wave3_temp
	erase panelmerged_wave3_temp.dta
	replace wave3 = 0 if wave3 == .
	
	tempname p
	postfile `p' str50 variable float mu_w1 sd_w1 N_w1 ///
		mu_w3 sd_w3 N_w3 ///
		attrition_pvalue ///
		mu_c sd_c N_c ///
		mu_a sd_a N_a ///
		mu_ce sd_ce N_ce ///
		mu_ae sd_ae N_ae ///
		mu_ex sd_ex N_ex ///
		anova_fstat anova_pvalue ///
		using summstats_main, replace

	foreach X of varlist `var_demog_reg_personal' az_demographics_education az_demographics_english az_demographics_travel az_demographics_household az_preference_risk az_preference_time az_preference_altruism az_preference_reciprocity panelmerged_wave3 treatment_control treatment_vpnonly treatment_nlonly treatment_vpnnl treatment_user {
		qui sum `X' if wave3 == 0
		local N_w1=r(N)
		local mu_w1=r(mean)
		local sd_w1=r(sd)

		qui sum `X' if panelmerged_wave3 == 1 & wave3 == 0
		local N_w3=r(N)
		local mu_w3=r(mean)
		local sd_w3=r(sd)
		
		qui sum `X' if panelmerged_wave3 == 1 & wave3 == 0 & treatment_master == 1
		local N_c=r(N)
		local mu_c=r(mean)
		local sd_c=r(sd)

		qui sum `X' if panelmerged_wave3 == 1 & wave3 == 0 & treatment_master == 2
		local N_a=r(N)
		local mu_a=r(mean)
		local sd_a=r(sd)

		qui sum `X' if panelmerged_wave3 == 1 & wave3 == 0 & treatment_master == 3
		local N_ce=r(N)
		local mu_ce=r(mean)
		local sd_ce=r(sd)

		qui sum `X' if panelmerged_wave3 == 1 & wave3 == 0 & treatment_master == 4
		local N_ae=r(N)
		local mu_ae=r(mean)
		local sd_ae=r(sd)

		qui sum `X' if panelmerged_wave3 == 1 & wave3 == 0 & treatment_master == 5
		local N_ex=r(N)
		local mu_ex=r(mean)
		local sd_ex=r(sd)
		
		qui anova `X' treatment_master if panelmerged_wave3 == 1 & wave3 == 0 & treatment_master != 5
		local anova_fstat=e(F)
		local anova_pvalue=Ftail(e(df_m), e(df_r), e(F))
		
		ttest `X', by(wave3)
		local attrition_pvalue=r(p)
		
		post `p' ("`X'") /// 
		(`mu_w1') (`sd_w1') (`N_w1') ///
		(`mu_w3') (`sd_w3') (`N_w3') ///
		(`attrition_pvalue') ///
		(`mu_c') (`sd_c') (`N_c') ///
		(`mu_ce') (`sd_ce') (`N_ce') /// 
		(`mu_a') (`sd_a') (`N_a') /// 
		(`mu_ae') (`sd_ae') (`N_ae') /// 
		(`mu_ex') (`sd_ex') (`N_ex') ///
		(`anova_fstat') (`anova_pvalue')
	}
		
	postclose `p'
	restore
	
	preserve	
	clear
	use summstats_main
	outsheet using summstats_main.xls, replace
	clear
	erase summstats_main.dta
	restore
	
	
	
*** Table 3

	cap erase outregs_maintreatmenteffects.xls
	cap erase outregs_maintreatmenteffects.txt	
	
	preserve
	keep if panelmerged_wave3 == 1
		
	foreach Y in a b c d e {
		
		// Panel A: reduced form
		qui reg az_overall_`Y'_w3 treatment_vpnonly treatment_nlonly treatment_vpnnl if treatment_user != 1, r
		qui summ az_overall_`Y'_w3 if treatment_user != 1
		local m_all = `r(mean)'
		local sd_all = `r(sd)'
		qui summ az_overall_`Y'_w3 if treatment_control == 1
		local m_control = `r(mean)'
		local sd_control = `r(sd)'
		qui summ az_overall_`Y'_w3 if treatment_user == 1
		local m_user = `r(mean)'
		local sd_user = `r(sd)'		
		local labelY: var label az_overall_`Y'_w3
		local t_treatvpnXnl = _b[treatment_vpnnl]/_se[treatment_vpnnl]
		local p_treatvpnXnl = 2*ttail(e(df_r),abs(`t_treatvpnXnl'))
		outreg2 treatment_vpnonly treatment_nlonly treatment_vpnnl using outregs_maintreatmenteffects.xls, addstat("Mean DV all", `m_all', "Std.Dev.\ DV all", `sd_all', "Mean DV control", `m_control', "Std.Dev.\ DV control", `sd_control', "Mean DV user", `m_user', "Std.Dev.\ DV user", `sd_user', p-value, `p_treatvpnXnl') label dec(3) adec(3) ctitle(`Y', "`labelY'", ) addtext(Specification, Raw diff) br		
				
		}
		
	restore

	// Panel A, addition: attrition bound (Lee 2009)
	preserve
	do "`directory_path'/trimbound.do"	
	gen weights = 1
	
	keep if treatment_control == 1 | treatment_vpnnl == 1
	
	trimbound az_overall_a_w3 treatment_vpnnl weights 
	trimbound az_overall_b_w3 treatment_vpnnl weights 
	trimbound az_overall_c_w3 treatment_vpnnl weights 
	trimbound az_overall_d_w3 treatment_vpnnl weights 
	trimbound az_overall_e_w3 treatment_vpnnl weights
	
	restore
	
	// Panel B: two-stage estimates
	preserve
	keep if panelmerged_wave3 == 1
	keep if treatment_user != 1
	
	ivregress 2sls az_overall_a_w3 (active_user = treatment_vpnonly treatment_nlonly treatment_vpnnl), first
	ivregress 2sls az_overall_b_w3 (active_user = treatment_vpnonly treatment_nlonly treatment_vpnnl), first
	ivregress 2sls az_overall_c_w3 (active_user = treatment_vpnonly treatment_nlonly treatment_vpnnl), first
	ivregress 2sls az_overall_d_w3 (active_user = treatment_vpnonly treatment_nlonly treatment_vpnnl), first
	ivregress 2sls az_overall_e_w3 (active_user = treatment_vpnonly treatment_nlonly treatment_vpnnl), first
	restore
	
	
	
*** Table 4

	cap erase outregs_wave3sociallearning.xls
	cap erase outregs_wave3sociallearning.txt	
	
	preserve
	
	foreach Y in news_c_panamapapers {
		qui reg `Y' soclearning_ownaccess soclearning_rm_new_w2 soclearning_ownXnew_w2 if vpn_roommate_existing == 0 & soclearning_rm_new_w2 < 2
		
		// estimating p(direct learning)
		qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 0)
		local mu_100=r(mean)
		qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 0)
		local mu_000=r(mean)
		local p_all = `mu_100'-`mu_000'
		
		// estimating q based on 1 roommate
		qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 1)
		local mu_101=r(mean)
		qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 1)
		local mu_001=r(mean)
		local q_new_1_0 = (`mu_001'-`mu_000')/`mu_100'
		local q_new_1_1 = (`mu_101'-`mu_100')/`mu_100'

		// predicted and actual based on 2 roommates
		local pred_mu_102 = `mu_100' + (1-(1-`mu_100'*`q_new_1_1')^2)
		local pred_mu_002 = `mu_000' + (1-(1-`mu_100'*`q_new_1_0')^2)
		qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 2)
		local mu_002=r(mean)
		qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 2)
		local mu_102=r(mean)
		
		local labelY: var label `Y'
		outreg2 soclearning_ownaccess soclearning_rm_new_w2 soclearning_ownXnew_w2 using outregs_wave3sociallearning.xls, addstat("p_all", `p_all', "q_new_1_0", `q_new_1_0', "q_new_1_1", `q_new_1_1', "pred_mu_002", `pred_mu_002', "mu_002", `mu_002', "pred_mu_102", `pred_mu_102', "mu_102", `mu_102') label ctitle(`Y', "`labelY'", ) addtext(Specification, Raw Diff) br		
		}

		
	foreach Y in news_c_coalprod news_c_hkceelection news_perccor_cen_all {
		qui reg `Y' soclearning_ownaccess soclearning_rm_new_w3 soclearning_ownXnew_w3 if vpn_roommate_existing == 0 & soclearning_rm_new_w3 < 2
		
		// estimating p(direct learning)
		qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 0)
		local mu_100=r(mean)
		qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 0)
		local mu_000=r(mean)
		local p_all = `mu_100'-`mu_000'
		
		// estimating q based on 1 roommate
		qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 1)
		local mu_101=r(mean)
		local var_
		qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 1)
		local mu_001=r(mean)
		local q_new_1_0 = (`mu_001'-`mu_000')/`mu_100'
		local q_new_1_1 = (`mu_101'-`mu_100')/`mu_100'

		// predicted and actual based on 2 roommates
		local pred_mu_102 = `mu_100' + (1-(1-`mu_100'*`q_new_1_1')^2)
		local pred_mu_002 = `mu_000' + (1-(1-`mu_100'*`q_new_1_0')^2)
		qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 2)
		local mu_002=r(mean)
		qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 2)
		local mu_102=r(mean)
		
		local labelY: var label `Y'
		outreg2 soclearning_ownaccess soclearning_rm_new_w3 soclearning_ownXnew_w3 using outregs_wave3sociallearning.xls, addstat("p_all", `p_all', "q_new_1_0", `q_new_1_0', "q_new_1_1", `q_new_1_1', "pred_mu_002", `pred_mu_002', "mu_002", `mu_002', "pred_mu_102", `pred_mu_102', "mu_102", `mu_102') label ctitle(`Y', "`labelY'", ) addtext(Specification, Raw Diff) br		
		}
		
	restore
	
	// exporting standard error of estimates
	
	preserve
	gen sc_reg_100 = (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 0)
	gen sc_reg_000 = (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 0)
	gen sc_reg_101 = (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 1)
	gen sc_reg_001 = (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 1)

	foreach Y in news_c_panamapapers {
		qui reg `Y' sc_reg_100 sc_reg_000 sc_reg_101 sc_reg_001 if (vpn_roommate_existing == 0 & vpn_roommate_new_w2 < 2), nocons
		nlcom 	(q_new_1_0: (_b[sc_reg_001] - _b[sc_reg_000])/_b[sc_reg_100]) ///
				(q_new_1_1: (_b[sc_reg_101] - _b[sc_reg_100])/_b[sc_reg_100]) ///
				(pred_mu_002: _b[sc_reg_000] + (1-(1-_b[sc_reg_100]*((_b[sc_reg_001] - _b[sc_reg_000])/_b[sc_reg_100]))^2)) ///
				(pred_mu_102: _b[sc_reg_100] + (1-(1-_b[sc_reg_100]*((_b[sc_reg_101] - _b[sc_reg_100])/_b[sc_reg_100]))^2))
		}
		
	restore
	
	preserve
	gen sc_reg_100 = (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 0)
	gen sc_reg_000 = (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 0)
	gen sc_reg_101 = (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 1)
	gen sc_reg_001 = (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 1)

	foreach Y in news_c_coalprod news_c_hkceelection news_perccor_cen_all {
		qui reg `Y' sc_reg_100 sc_reg_000 sc_reg_101 sc_reg_001 if (vpn_roommate_existing == 0 & vpn_roommate_new_w3 < 2), nocons
		nlcom 	(q_new_1_0: (_b[sc_reg_001] - _b[sc_reg_000])/_b[sc_reg_100]) ///
				(q_new_1_1: (_b[sc_reg_101] - _b[sc_reg_100])/_b[sc_reg_100]) ///
				(pred_mu_002: _b[sc_reg_000] + (1-(1-_b[sc_reg_100]*((_b[sc_reg_001] - _b[sc_reg_000])/_b[sc_reg_100]))^2)) ///
				(pred_mu_102: _b[sc_reg_100] + (1-(1-_b[sc_reg_100]*((_b[sc_reg_101] - _b[sc_reg_100])/_b[sc_reg_100]))^2))
		}
		
	restore

	/*
	// calculate bootstrap standard errors for out-of-sample predictions
		
	preserve
	clear
	gen i = .
	gen y = ""
	gen pred_mu_102 = .
	gen pred_mu_002 = .
	save sociallearning_pred_bootstrap, replace
	restore
	
	keep treatment_master vpn_roommate_existing vpn_roommate_new_w2 vpn_roommate_new_w3 `var_knowledge_news_reg_cen_w2' `var_knowledge_news_reg_cen_w3' news_perccor_cen_all
	keep if news_perccor_cen_all != .
	
	forvalues i = 1/1000 {
		
		foreach Y in news_c_panamapapers {
			
			preserve
			
			bsample	
	
			// estimating p(direct learning)
			qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 0)
			local mu_100=r(mean)
			qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 0)
			local mu_000=r(mean)
			local p_all = `mu_100'-`mu_000'
			
			// estimating q based on 1 roommate
			qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 1)
			local mu_101=r(mean)
			local var_
			qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 1)
			local mu_001=r(mean)
			local q_new_1_0 = (`mu_001'-`mu_000')/`mu_100'
			local q_new_1_1 = (`mu_101'-`mu_100')/`mu_100'

			// predicted and actual based on 2 roommates
			local pred_mu_102 = `mu_100' + (1-(1-`mu_100'*`q_new_1_1')^2)
			local pred_mu_002 = `mu_000' + (1-(1-`mu_100'*`q_new_1_0')^2)
			qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 2)
			local mu_002=r(mean)
			qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 2)
			local mu_102=r(mean)
			
			// export estimates
			gen i = `i'
			gen y = "`Y'"
			gen pred_mu_102 = `pred_mu_102'
			gen pred_mu_002 = `pred_mu_002'
			keep i y pred_mu_102 pred_mu_002
			duplicates drop
			append using sociallearning_pred_bootstrap
			save sociallearning_pred_bootstrap, replace
		
			restore
			
			}
			
		
		foreach Y in news_c_coalprod news_c_hkceelection news_perccor_cen_all {
			
			preserve
			
			bsample	
	
			// estimating p(direct learning)
			qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 0)
			local mu_100=r(mean)
			qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 0)
			local mu_000=r(mean)
			local p_all = `mu_100'-`mu_000'
			
			// estimating q based on 1 roommate
			qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 1)
			local mu_101=r(mean)
			local var_
			qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 1)
			local mu_001=r(mean)
			local q_new_1_0 = (`mu_001'-`mu_000')/`mu_100'
			local q_new_1_1 = (`mu_101'-`mu_100')/`mu_100'

			// predicted and actual based on 2 roommates
			local pred_mu_102 = `mu_100' + (1-(1-`mu_100'*`q_new_1_1')^2)
			local pred_mu_002 = `mu_000' + (1-(1-`mu_100'*`q_new_1_0')^2)
			qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 2)
			local mu_002=r(mean)
			qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 2)
			local mu_102=r(mean)
			
			// export estimates
			gen i = `i'
			gen y = "`Y'"
			gen pred_mu_102 = `pred_mu_102'
			gen pred_mu_002 = `pred_mu_002'
			keep i y pred_mu_102 pred_mu_002
			duplicates drop
			append using sociallearning_pred_bootstrap
			save sociallearning_pred_bootstrap, replace
		
			restore
			
			}

		}
		
	use sociallearning_pred_bootstrap, clear
	
	gen pred_mu_102_bse = .
	gen pred_mu_002_bse = .
	
	foreach Y in news_c_coalprod news_c_panamapapers news_c_hkceelection news_perccor_cen_all {
		
		sum pred_mu_102 if y == "`Y'", d
		replace pred_mu_102_bse = `r(sd)' if y == "`Y'"
		sum pred_mu_002 if y == "`Y'", d
		replace pred_mu_002_bse = `r(sd)' if y == "`Y'"
		
		}
		
	keep y pred_mu_102_bse pred_mu_002_bse
	duplicates drop
	export excel using "sociallearning_pred_bootstrap.xlsx", firstrow(variables) replace
	
	*/
	


*** Figure A.8, A.11, A.12, A.13


	preserve
	qui {
	local v = 1
	foreach Y of varlist `var_media_valuation_reg_w3' az_belief_media_value_w3 `var_media_trust_reg_w3' az_belief_media_trust_w3 `var_censor_level_reg_w3' `var_censor_justif_reg_w3' az_belief_media_justif_w3 `var_censor_driver_reg_dom_w3' `var_censor_driver_reg_for_w3' `var_knowledge_news_reg_cen_w3' news_perccor_cen_w3 `var_knowledge_news_reg_unc_w3' news_perccor_unc_w3 `var_knowledge_prot_reg_chi_w3' protest_pcheard_china_w3 `var_knowledge_prot_reg_for_w3' protest_pcheard_foreign_w3 `var_knowledge_prot_reg_fak_w3' `var_knowledge_meta_reg_w3' az_knowledge_meta_w3 `var_econ_guess_reg_cn_perf_w3' az_belief_econ_perf_cn_w3 `var_econ_guess_reg_cn_conf_w3' az_belief_econ_conf_cn_w3 `var_econ_guess_reg_us_perf_w3' az_belief_econ_perf_us_w3 `var_econ_guess_reg_us_conf_w3' az_belief_econ_conf_us_w3 `var_demand_change_reg_w3' az_belief_instchange_w3 `var_trust_inst_reg_govt_w3' az_belief_trust_govt_w3 `var_trust_inst_reg_foreign_w3' az_belief_trust_foreign_w3 `var_eval_govt_reg_w3' az_belief_evalgovt_w3 `var_democracy_reg_w3' `var_willing_fight_reg_w3' az_belief_willing_w3 `var_vpn_purchase_reg_w3' `var_info_ranking_reg_w3' `var_info_freq_reg_w3' `var_socialinteract_reg_w3' `var_polparticipation_reg_w3' `var_stock_invest_reg_w3' `var_planaftergrad_reg_w3' `var_career_sector_reg_w3' `var_career_loc_reg_w3' {
		qui xi: reg `Y' treatment_vpnonly treatment_vpnnl, r
		su `Y' if e(sample)
		gen eb = (`Y' - r(mean))/r(sd)
		qui xi: reg eb treatment_vpnonly treatment_vpnnl if vpn_current_paid_user != 1, r
		gen rcoef_`v' = _b[treatment_vpnnl]
		gen rvpn_`v' = _b[treatment_vpnonly]
		gen rse_`v' = _se[treatment_vpnnl]
		gen rcih_`v' = _b[treatment_vpnnl] + 1.9*_se[treatment_vpnnl]
		gen rcil_`v' = _b[treatment_vpnnl] - 1.9*_se[treatment_vpnnl]
		qui xi: reg eb treatment_user if (treatment_control == 1 | treatment_user == 1), r
		gen ruser_`v' = _b[treatment_user]
		drop eb
		local v = `v' + 1
		}

	keep rcoef_* rvpn_* rse_* rcil_* rcih_* ruser_*
	duplicates drop
	gen i = 1
	reshape long rcoef_ rvpn_ rse_ rcih_ rcil_ ruser_, i(i) j(vnum)
	label drop _all

	expand 2
	gen p = _n > 117
	forval v = 1/117 {
		gen ci`v' = rcil_ if p == 0 & vnum == `v'
		replace ci`v' = rcih_ if p == 1 & vnum == `v'
		replace rcoef_ = . if p == 1 & vnum == `v'
		}

	replace vnum = vnum + 35 if inrange(vnum,110,117) 	// E.7.b
	replace vnum = vnum + 34 if inrange(vnum,102,109) 	// E.7.a
	replace vnum = vnum + 33 if inrange(vnum,97,101) 	// E.6
	replace vnum = vnum + 32 if inrange(vnum,96,96) 	// E.5
	replace vnum = vnum + 31 if inrange(vnum,93,95) 	// E.4
	replace vnum = vnum + 30 if inrange(vnum,91,92) 	// E.3
	replace vnum = vnum + 29 if inrange(vnum,90,90) 	// E.2.b
	replace vnum = vnum + 28 if inrange(vnum,85,89) 	// E.2.a
	replace vnum = vnum + 27 if inrange(vnum,83,84) 	// E.1
	replace vnum = vnum + 25 if inrange(vnum,79,82) 	// D.8
	replace vnum = vnum + 24 if inrange(vnum,78,78) 	// D.6
	replace vnum = vnum + 23 if inrange(vnum,74,77) 	// D.3
	replace vnum = vnum + 22 if inrange(vnum,71,73) 	// D.2.b
	replace vnum = vnum + 21 if inrange(vnum,67,70) 	// D.2.a
	replace vnum = vnum + 20 if inrange(vnum,64,66) 	// D.1
	replace vnum = vnum + 18 if inrange(vnum,61,63) 	// C.4
	replace vnum = vnum + 17 if inrange(vnum,58,60) 	// C.3
	replace vnum = vnum + 16 if inrange(vnum,55,57) 	// C.2
	replace vnum = vnum + 15 if inrange(vnum,52,54) 	// C.1
	replace vnum = vnum + 13 if inrange(vnum,49,51) 	// B.5
	replace vnum = vnum + 12 if inrange(vnum,48,48) 	// B.3.c
	replace vnum = vnum + 11 if inrange(vnum,42,47) 	// B.3.b
	replace vnum = vnum + 10 if inrange(vnum,37,41) 	// B.3.a
	replace vnum = vnum + 9 if inrange(vnum,32,36) 		// B.2.b
	replace vnum = vnum + 8 if inrange(vnum,24,31) 		// B.2.a
	replace vnum = vnum + 6 if inrange(vnum,20,23) 		// A.5.b
	replace vnum = vnum + 5 if inrange(vnum,16,19) 		// A.5.a
	replace vnum = vnum + 4 if inrange(vnum,10,15) 		// A.4
	replace vnum = vnum + 3 if inrange(vnum,9,9) 		// A.3.b
	replace vnum = vnum + 2 if inrange(vnum,8,8) 		// A.3.a
	replace vnum = vnum + 1 if inrange(vnum,4,7) 		// A.2

	label define varname ///
		1 "Willingness to pay for circumvention tool" ///
		2 "Value added of foreign media access" ///
		3 "{it: z-score: valuation of access to foreign media outlets}" ///
		5 "Distrust in domestic state-owned media" ///
		6 "Distrust in domestic privately-owned media" ///
		7 "Trust in foreign media" ///
		8 "{it: z-score: trust in non-domestic media outlets}" ///
		10 "Degree of censorship on domestic news outlets" ///
		12 "Degree of censorship on foreign news outlets" ///
		14 "Unjustified: censoring economic news" ///
		15 "Unjustified: censoring political news" ///
		16 "Unjustified: censoring social news" ///
		17 "Unjustified: censoring foreign news" ///
		18 "Unjustified: censoring pornography" ///
		19 "{it: z-score: censorship unjustified}" ///
		21 "Domestic cens. driven by govt. policies" ///
		22 "Domestic cens. driven by corp. interest" ///
		23 "Domestic cens. driven by media???s ideology" ///
		24 "Domestic cens. driven by readers??? demand" ///
		26 "Foreign cens. driven by govt. policies" ///
		27 "Foreign cens. driven by corp. interest" ///
		28 "Foreign cens. driven by media???s ideology" ///
		29 "Foreign cens. driven by readers??? demand" ///
		32 "Steel production reduction reaches target" ///
		33 "Trump registered trademarks in China" ///
		34 "Jianhua Xiao kidnapped in Hong Kong" ///
		35 "Xinjiang installed GPS on all automobiles" ///
		36 "China and Norway re-normalize ties" ///
		37 "Feminist groups fight women's rights" ///
		38 "Carrie Lam becomes HK Chief Executive" ///
		39 "{it: % quizzes answered correctly: poli. sensitive news}" ///
		41 "China stops importing coal from North Korea" ///
		42 "H7N9 influenza epidemic" ///
		43 "Transnational railway in Ethiopia" ///
		44 "Foreign reserves fall below threshold" ///
		45 "{it: % quizzes answered correctly: nonsensitive news}" ///
		47 "2012 HK Anti-National Curr. Movement" ///
		48 "2014 HK Umbrella Revolution" ///
		49 "2016 HK Mong Kok Revolution" ///
		50 "2014 Taiwan Sunflower Stud. Movement" ///
		51 "{it: % protests in Greater China heard of}" ///
		53 "2014 Ukrainian Euromaidan Revolution" ///
		54 "2010 Arab Spring" ///
		55 "2014 Crimean Status Referendum" ///
		56 "2010 Catalonian Indep. Movement" ///
		57 "2017 Women's March" ///
		58 "{it: % foreign protests heard of}" ///
		60 "2011 Tomorrow Revolution [fake]" ///
		62 "Informedness of issues in China" ///
		63 "Greater informedness than peers" ///
		64 "{it: z-score: self-assessment of knowledge level}" ///
		67 "Guess on GDP growth rate in 2016 China" ///
		68 "Guess on SSCI by end of 2016" ///
		69 "{it: z-score: optimistic belief of Chinese economy}" ///
		71 "Confidence of China GDP guess" ///
		72 "Confidence of SSCI guess" ///
		73 "{it: z-score: confidence of guesses on Chinese economy}" ///
		75 "Guess on GDP growth rate in 2016 US" ///
		76 "Guess on DJI by end of 2016" ///
		77 "{it: z-score: optimistic belief of US economy}" ///
		79 "Confidence of US GDP guess" ///
		80 "Confidence of DJI guess" ///
		81 "{it: z-score: confidence of guesses on US economy}" ///
		84 "Economic system needs changes" ///
		85 "Political system needs changes" ///
		86 "{it: z-score: demand for institutional change}" ///
		88 "Trust in central govt. of China" ///
		89 "Trust in provincial govt. of China" ///
		90 "Trust in local govt. of China" ///
		91 "{it: z-score: trust in Chinese govt.}" ///
		93 "Trust in central govt. of Japan" ///
		94 "Trust in federal govt. of US" ///
		95 "{it: z-score: trust in foreign govt.}" ///
		97 "Satisfaction of economic dev." ///
		98 "Satisfaction of domestic politics" ///
		99 "Satisfaction of diplomatic affairs" ///
		100 "{it: z-score: satisfaction of govt???s performance}" ///
		102 "Living in democracy is not important" ///
		104 "Willing to battle illegal govt. acts" ///
		105 "Willing to report govt. misconduct" ///
		106 "Willing to stand up for the weak" ///
		107 "{it: z-score: willingness to act}" ///
		110 "Purchase discounted tool we offered" ///
		111 "Purchase any tool" ///
		113 "Ranked high: domestic websites" ///
		114 "Ranked high: foreign websites" ///
		115 "Ranked high: domestic social media" ///
		116 "Ranked high: foreign social media" ///
		117 "Ranked high: word of mouth" ///
		119 "Frequency of visiting foreign websites for info." ///
		121 "Frequency of discussing poli. with friends" ///
		122 "Frequency of persuading others" ///
		124 "Protests concerning social issues" ///
		125 "Plan to vote for local PCR" ///
		126 "Complain to school authorities" ///
		128 "Currently invested in Chinese stock mkt." ///
		130 "Plan: grad. school in China" ///
		131 "Plan: master degree abroad" ///
		132 "Plan: PhD degree abroad" ///
		133 "Plan: military in China" ///
		134 "Plan: work right away" ///
		136 "Sector pref.: national civil service" ///
		137 "Sector pref.: local civil service" ///
		138 "Sector pref.: military" ///
		139 "Sector pref.: private firm in China" ///
		140 "Sector pref.: foreign firm in China" ///
		141 "Sector pref.: SOEs" ///
		142 "Sector pref.: inst. organizations" ///
		143 "Sector pref.: entrepreneurship" ///
		145 "Location pref.: Beijing" ///
		146 "Location pref.: Shanghai" ///
		147 "Location pref.: Guangzhou and Shenzhen" ///
		148 "Location pref.: tier 2 cities in central" ///
		149 "Location pref.: other cities in China" ///
		150 "Location pref.: HK and Macau" ///
		151 "Location pref.: Taiwan" ///
		152 "Location pref.: foreign cities" 
		
	label values vnum varname
	scatter vnum ci1, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || ///
	scatter vnum ci2, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || ///
	scatter vnum ci3, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci4, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci5, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci6, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci7, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci8, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci9, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci10, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci11, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci12, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci13, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci14, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci15, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci16, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci17, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci18, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci19, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci20, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci21, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci22, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci23, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci24, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci25, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci26, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci27, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci28, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci29, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci30, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci31, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci32, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci33, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci34, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci35, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci36, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci37, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci38, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci39, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci40, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci41, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci42, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci43, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci44, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci45, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci46, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci47, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci48, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci49, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci50, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci51, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci52, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci53, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci54, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci55, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci56, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci57, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci58, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci59, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci60, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci61, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci62, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci63, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci64, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci65, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci66, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci67, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci68, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci69, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci70, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci71, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci72, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci73, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci74, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci75, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci76, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci77, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci78, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci79, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci80, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci81, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci82, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci83, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci84, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci85, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci86, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci87, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci88, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci89, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci90, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci91, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci92, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci93, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci94, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci95, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci96, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci97, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci98, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci99, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci100, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci101, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci102, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci103, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci104, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci105, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci106, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci107, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci108, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci109, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci110, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci111, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci112, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci113, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci114, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci115, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci116, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci117, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum rcoef_ if vnum==1, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==2, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==3, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==5, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==6, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==7, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==8, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==10, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==12, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==14, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==15, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==16, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==17, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==18, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==19, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==21, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==22, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==23, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==24, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==26, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==27, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==28, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==29, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==32, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==33, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==34, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==35, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==36, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==37, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==38, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==39, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==41, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==42, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==43, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==44, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==45, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==47, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==48, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==49, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==50, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==51, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==53, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==54, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==55, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==56, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==57, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==58, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==60, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==62, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==63, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==64, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==67, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==68, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==69, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==71, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==72, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==73, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==75, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==76, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==77, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==79, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==80, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==81, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==84, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==85, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==86, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==88, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==89, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==90, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==91, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==93, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==94, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==95, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==97, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==98, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==99, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==100, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==102, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==104, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==105, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==106, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==107, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==110, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==111, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==113, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==114, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==115, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==116, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==117, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==119, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==121, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==122, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==124, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==125, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==126, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==128, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==130, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==131, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==132, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==133, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==134, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==136, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==137, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==138, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==139, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==140, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==141, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==142, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==143, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==145, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==146, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==147, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==148, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==149, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==150, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==151, mcolor(cranberry) m(S) || ///
	scatter vnum rcoef_ if vnum==152, mcolor(cranberry) m(S) || ///
	scatter vnum rvpn_ if vnum==1, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==2, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==3, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==5, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==6, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==7, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==8, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==10, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==12, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==14, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==15, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==16, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==17, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==18, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==19, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==21, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==22, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==23, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==24, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==26, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==27, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==28, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==29, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==32, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==33, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==34, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==35, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==36, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==37, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==38, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==39, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==41, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==42, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==43, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==44, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==45, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==47, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==48, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==49, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==50, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==51, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==53, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==54, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==55, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==56, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==57, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==58, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==60, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==62, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==63, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==64, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==67, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==68, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==69, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==71, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==72, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==73, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==75, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==76, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==77, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==79, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==80, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==81, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==84, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==85, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==86, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==88, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==89, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==90, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==91, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==93, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==94, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==95, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==97, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==98, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==99, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==100, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==102, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==104, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==105, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==106, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==107, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==110, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==111, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==113, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==114, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==115, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==116, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==117, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==119, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==121, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==122, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==124, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==125, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==126, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==128, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==130, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==131, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==132, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==133, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==134, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==136, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==137, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==138, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==139, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==140, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==141, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==142, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==143, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==145, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==146, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==147, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==148, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==149, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==150, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==151, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum rvpn_ if vnum==152, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(t) || ///
	scatter vnum ruser_ if vnum==1, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==2, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==3, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==5, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==6, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==7, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==8, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==10, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==12, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==14, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==15, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==16, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==17, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==18, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==19, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==21, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==22, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==23, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==24, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==26, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==27, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==28, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==29, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==32, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==33, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==34, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==35, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==36, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==37, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==38, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==39, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==41, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==42, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==43, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==44, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==45, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==47, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==48, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==49, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==50, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==51, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==53, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==54, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==55, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==56, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==57, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==58, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==60, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==62, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==63, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==64, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==67, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==68, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==69, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==71, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==72, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==73, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==75, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==76, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==77, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==79, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==80, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==81, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==84, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==85, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==86, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==88, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==89, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==90, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==91, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==93, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==94, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==95, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==97, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==98, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==99, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==100, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==102, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==104, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==105, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==106, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==107, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==110, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==111, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==113, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==114, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==115, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==116, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==117, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==119, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==121, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==122, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==124, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==125, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==126, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==128, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==130, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==131, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==132, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==133, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==134, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==136, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==137, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==138, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==139, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==140, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==141, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==142, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==143, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==145, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==146, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==147, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==148, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==149, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==150, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==151, mcolor(emidblue) m(dh) || ///
	scatter vnum ruser_ if vnum==152, mcolor(emidblue) m(dh) ///
		xline(0,   lc(midblue) lp(line)) ///
		yline(4,   lc(gs5)  lp(longdash)) ///
		yline(9,   lc(gs5)  lp(longdash)) ///
		yline(11,  lc(gs10) lp(dash)) ///
		yline(13,  lc(gs5)  lp(longdash)) ///
		yline(20,  lc(gs5)  lp(longdash)) ///
		yline(25,  lc(gs10) lp(dash)) ///
		yline(31,  lc(gs3)  lp(line) lwidth(medthick)) ///
		yline(40,  lc(gs10) lp(dash)) ///
		yline(46,  lc(gs5)  lp(longdash)) ///
		yline(52,  lc(gs10) lp(dash)) ///
		yline(59,  lc(gs10) lp(dash)) ///
		yline(61,  lc(gs5)  lp(longdash)) ///
		yline(66,  lc(gs3)  lp(line) lwidth(medthick)) ///
		yline(70,  lc(gs5)  lp(longdash)) ///
		yline(74,  lc(gs5)  lp(longdash)) ///
		yline(78,  lc(gs5)  lp(longdash)) ///
		yline(83,  lc(gs3)  lp(line) lwidth(medthick)) ///
		yline(87,  lc(gs5)  lp(longdash)) ///
		yline(92,  lc(gs10) lp(dash)) ///
		yline(96,  lc(gs5)  lp(longdash)) ///
		yline(101, lc(gs5)  lp(longdash)) ///
		yline(103, lc(gs5)  lp(longdash)) ///
		yline(109, lc(gs3)  lp(line) lwidth(medthick)) ///
		yline(112, lc(gs5)  lp(longdash)) ///
		yline(118, lc(gs12) lp(dash)) ///
		yline(120, lc(gs5)  lp(longdash)) ///
		yline(123, lc(gs5)  lp(longdash)) ///
		yline(127, lc(gs5)  lp(longdash)) ///
		yline(129, lc(gs5)  lp(longdash)) ///
		yline(135, lc(gs5)  lp(longdash)) ///
		yline(144, lc(gs12) lp(dash)) ///
		ytitle("") yscale(reverse) ///
		xtitle("Standardized means (Control = 0)", size(small)) xlabel(, grid labsize(small)) ///
		ylab(1 2 3 5 6 7 8 10 12 14 15 16 17 18 19 21 22 23 24 26 27 28 29 32 33 34 35 36 37 38 39 41 42 43 44 45 47 48 49 50 51 53 54 55 56 57 58 60 62 63 64 67 68 69 71 72 73 75 76 77 79 80 81 84 85 86 88 89 90 91 93 94 95 97 98 99 100 102 104 105 106 107 110 111 113 114 115 116 117 119 121 122 124 125 126 128 130 131 132 133 134 136 137 138 139 140 141 142 143 145 146 147 148 149 150 151 152, valuelabel angle(0) labsize(vsmall) nogrid) ///
		ysize(20) xsize(4) ///
		legend(order(235 "Access" 118 "Access + Encour." 352 "Existing users") size(vsmall)) ///
		graphregion(color(white))
	}
	
	graph export figure_dotplot_master_w3.pdf, replace
	restore
	
	


*** Figure A.9

	foreach Y in info_freq_website_for plan_grad_foreignmaster cloc_for stock_participation {
	
		preserve
		keep if panelmerged_wave2 == 1 & panelmerged_wave3 == 1
		collapse (mean) `Y'_w3 `Y'_w2 `Y'_w1 (sd) s_`Y'_w3 = `Y'_w3 s_`Y'_w2 = `Y'_w2 s_`Y'_w1 = `Y'_w1 (count) n_`Y'_w3 = `Y'_w3 n_`Y'_w2 = `Y'_w2 n_`Y'_w1 = `Y'_w1, by(treatment_main)
		gen h_`Y'_w3 = `Y'_w3 + invttail(n_`Y'_w3-1,0.05)*(s_`Y'_w3 / sqrt(n_`Y'_w3))
		gen l_`Y'_w3 = `Y'_w3 - invttail(n_`Y'_w3-1,0.05)*(s_`Y'_w3 / sqrt(n_`Y'_w3))	
		gen h_`Y'_w2 = `Y'_w2 + invttail(n_`Y'_w2-1,0.05)*(s_`Y'_w2 / sqrt(n_`Y'_w2))
		gen l_`Y'_w2 = `Y'_w2 - invttail(n_`Y'_w2-1,0.05)*(s_`Y'_w2 / sqrt(n_`Y'_w2))	
		gen h_`Y'_w1 = `Y'_w1 + invttail(n_`Y'_w1-1,0.05)*(s_`Y'_w1 / sqrt(n_`Y'_w1))
		gen l_`Y'_w1 = `Y'_w1 - invttail(n_`Y'_w1-1,0.05)*(s_`Y'_w1 / sqrt(n_`Y'_w1))	
		drop s_* n_*
		drop if treatment_main == .
		
		reshape long `Y'_w h_`Y'_w l_`Y'_w, i(treatment_main) j(wave)
		replace wave = 4 if wave == 3
		
		twoway 	(connected `Y'_w wave if treatment_main == 1, msymbol(o) mcolor(gs9) lpattern(dash) lcolor(gs9) lwidth(medthick)) /// 
				(connected `Y'_w wave if treatment_main == 2, msymbol(d) mcolor(gs3) lpattern(dash) lcolor(gs3) lwidth(medthick)) /// 
				(connected `Y'_w wave if treatment_main == 3, msymbol(o) mcolor(cranberry) lpattern(line) lcolor(cranberry) lwidth(vthick)) /// 
				(connected `Y'_w wave if treatment_main == 4, msymbol(o) mcolor(emidblue) lpattern(line) lcolor(emidblue) lwidth(thick)) /// 
				(rcap h_`Y'_w l_`Y'_w wave), /// 
				xline(1.2) ///
				title("`l_`Y''", size(medium)) ///
				ytitle("") ///
				xtitle("") ///
				xscale(r(0.8 4.2)) ///
				xlabel(1 "Baseline" 2 "Midline" 4 "Endline", noticks) ///
				legend(order(1 "Control" 2 "Access" 3 "Access + Encour." 4 "Existing users") col(4) size(small)) ///
				graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
				saving(figure_panel_`Y', replace)
				
		graph export "figure_panel_`Y'.pdf", replace
		restore
		
	}

	foreach Y in az_belief_media_value az_belief_media_trust az_belief_media_justif_pf az_knowledge_meta az_belief_econ_perf_cn az_belief_econ_conf_cn az_belief_instchange az_belief_trust_govt az_belief_trust_foreign az_belief_evalgovt az_belief_democracy_fp az_belief_willing az_var_socialinteract az_var_polparticipation {
		
		preserve
		keep if panelmerged_wave2 == 1 & panelmerged_wave3 == 1
		collapse (mean) `Y'_w3 `Y'_w2 `Y'_w1 (sd) s_`Y'_w3 = `Y'_w3 s_`Y'_w2 = `Y'_w2 s_`Y'_w1 = `Y'_w1 (count) n_`Y'_w3 = `Y'_w3 n_`Y'_w2 = `Y'_w2 n_`Y'_w1 = `Y'_w1, by(treatment_main)
		gen h_`Y'_w3 = `Y'_w3 + invttail(n_`Y'_w3-1,0.05)*(s_`Y'_w3 / sqrt(n_`Y'_w3))
		gen l_`Y'_w3 = `Y'_w3 - invttail(n_`Y'_w3-1,0.05)*(s_`Y'_w3 / sqrt(n_`Y'_w3))	
		gen h_`Y'_w2 = `Y'_w2 + invttail(n_`Y'_w2-1,0.05)*(s_`Y'_w2 / sqrt(n_`Y'_w2))
		gen l_`Y'_w2 = `Y'_w2 - invttail(n_`Y'_w2-1,0.05)*(s_`Y'_w2 / sqrt(n_`Y'_w2))	
		gen h_`Y'_w1 = `Y'_w1 + invttail(n_`Y'_w1-1,0.05)*(s_`Y'_w1 / sqrt(n_`Y'_w1))
		gen l_`Y'_w1 = `Y'_w1 - invttail(n_`Y'_w1-1,0.05)*(s_`Y'_w1 / sqrt(n_`Y'_w1))	
		drop s_* n_*
		drop if treatment_main == .
		
		reshape long `Y'_w h_`Y'_w l_`Y'_w, i(treatment_main) j(wave)
		replace wave = 4 if wave == 3
		
		twoway 	(connected `Y'_w wave if treatment_main == 1, msymbol(o) mcolor(gs9) lpattern(dash) lcolor(gs9) lwidth(medthick)) /// 
				(connected `Y'_w wave if treatment_main == 2, msymbol(d) mcolor(gs3) lpattern(dash) lcolor(gs3) lwidth(medthick)) /// 
				(connected `Y'_w wave if treatment_main == 3, msymbol(o) mcolor(cranberry) lpattern(line) lcolor(cranberry) lwidth(vthick)) /// 
				(connected `Y'_w wave if treatment_main == 4, msymbol(o) mcolor(emidblue) lpattern(line) lcolor(emidblue) lwidth(thick)) /// 
				(rcap h_`Y'_w l_`Y'_w wave), /// 
				xline(1.2) ///
				title("`l_`Y''", size(medium)) ///
				ytitle("") ///
				xtitle("") ///
				xscale(r(0.8 4.2)) ///
				yscale(r(-1 1)) ///
				ylabel(-1(0.5)1) ///
				xlabel(1 "Baseline" 2 "Midline" 4 "Endline", noticks) ///
				legend(order(1 "Control" 2 "Access" 3 "Access + Encour." 4 "Existing users") col(4) size(small)) ///
				graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
				saving(figure_panel_`Y', replace)
				
		graph export "figure_panel_`Y'.pdf", replace
		restore
		
	}

	foreach Y in bias_domestic bias_foreign {
		
		preserve
		keep if panelmerged_wave2 == 1 & panelmerged_wave3 == 1
		collapse (mean) `Y'_w3 `Y'_w2 `Y'_w1 (sd) s_`Y'_w3 = `Y'_w3 s_`Y'_w2 = `Y'_w2 s_`Y'_w1 = `Y'_w1 (count) n_`Y'_w3 = `Y'_w3 n_`Y'_w2 = `Y'_w2 n_`Y'_w1 = `Y'_w1, by(treatment_main)
		gen h_`Y'_w3 = `Y'_w3 + invttail(n_`Y'_w3-1,0.05)*(s_`Y'_w3 / sqrt(n_`Y'_w3))
		gen l_`Y'_w3 = `Y'_w3 - invttail(n_`Y'_w3-1,0.05)*(s_`Y'_w3 / sqrt(n_`Y'_w3))	
		gen h_`Y'_w2 = `Y'_w2 + invttail(n_`Y'_w2-1,0.05)*(s_`Y'_w2 / sqrt(n_`Y'_w2))
		gen l_`Y'_w2 = `Y'_w2 - invttail(n_`Y'_w2-1,0.05)*(s_`Y'_w2 / sqrt(n_`Y'_w2))	
		gen h_`Y'_w1 = `Y'_w1 + invttail(n_`Y'_w1-1,0.05)*(s_`Y'_w1 / sqrt(n_`Y'_w1))
		gen l_`Y'_w1 = `Y'_w1 - invttail(n_`Y'_w1-1,0.05)*(s_`Y'_w1 / sqrt(n_`Y'_w1))	
		drop s_* n_*
		drop if treatment_main == .
		
		reshape long `Y'_w h_`Y'_w l_`Y'_w, i(treatment_main) j(wave)
		replace wave = 4 if wave == 3
		
		twoway 	(connected `Y'_w wave if treatment_main == 1, msymbol(o) mcolor(gs9) lpattern(dash) lcolor(gs9) lwidth(medthick)) /// 
				(connected `Y'_w wave if treatment_main == 2, msymbol(d) mcolor(gs3) lpattern(dash) lcolor(gs3) lwidth(medthick)) /// 
				(connected `Y'_w wave if treatment_main == 3, msymbol(o) mcolor(cranberry) lpattern(line) lcolor(cranberry) lwidth(vthick)) /// 
				(connected `Y'_w wave if treatment_main == 4, msymbol(o) mcolor(emidblue) lpattern(line) lcolor(emidblue) lwidth(thick)) /// 
				(rcap h_`Y'_w l_`Y'_w wave), /// 
				xline(1.2) ///
				title("`l_`Y''", size(medium)) ///
				ytitle("") ///
				xtitle("") ///
				xscale(r(0.8 4.2)) ///
				yscale(r(4 9)) ///
				ylabel(4(1)9) ///
				xlabel(1 "Baseline" 2 "Midline" 4 "Endline", noticks) ///
				legend(order(1 "Control" 2 "Access" 3 "Access + Encour." 4 "Existing users") col(4) size(small)) ///
				graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
				saving(figure_panel_`Y', replace)
				
		graph export "figure_panel_`Y'.pdf", replace
		restore
		
	}

	foreach Y in news_perccor_cen news_perccor_unc protest_pcheard_china protest_pcheard_foreign {
		
		preserve
		keep if panelmerged_wave2 == 1 & panelmerged_wave3 == 1
		collapse (mean) `Y'_w3 `Y'_w2 `Y'_w1 (sd) s_`Y'_w3 = `Y'_w3 s_`Y'_w2 = `Y'_w2 s_`Y'_w1 = `Y'_w1 (count) n_`Y'_w3 = `Y'_w3 n_`Y'_w2 = `Y'_w2 n_`Y'_w1 = `Y'_w1, by(treatment_main)
		gen h_`Y'_w3 = `Y'_w3 + invttail(n_`Y'_w3-1,0.05)*(s_`Y'_w3 / sqrt(n_`Y'_w3))
		gen l_`Y'_w3 = `Y'_w3 - invttail(n_`Y'_w3-1,0.05)*(s_`Y'_w3 / sqrt(n_`Y'_w3))	
		gen h_`Y'_w2 = `Y'_w2 + invttail(n_`Y'_w2-1,0.05)*(s_`Y'_w2 / sqrt(n_`Y'_w2))
		gen l_`Y'_w2 = `Y'_w2 - invttail(n_`Y'_w2-1,0.05)*(s_`Y'_w2 / sqrt(n_`Y'_w2))	
		gen h_`Y'_w1 = `Y'_w1 + invttail(n_`Y'_w1-1,0.05)*(s_`Y'_w1 / sqrt(n_`Y'_w1))
		gen l_`Y'_w1 = `Y'_w1 - invttail(n_`Y'_w1-1,0.05)*(s_`Y'_w1 / sqrt(n_`Y'_w1))	
		drop s_* n_*
		drop if treatment_main == .
		
		reshape long `Y'_w h_`Y'_w l_`Y'_w, i(treatment_main) j(wave)
		replace wave = 4 if wave == 3
		
		twoway 	(connected `Y'_w wave if treatment_main == 1, msymbol(o) mcolor(gs9) lpattern(dash) lcolor(gs9) lwidth(medthick)) /// 
				(connected `Y'_w wave if treatment_main == 2, msymbol(d) mcolor(gs3) lpattern(dash) lcolor(gs3) lwidth(medthick)) /// 
				(connected `Y'_w wave if treatment_main == 3, msymbol(o) mcolor(cranberry) lpattern(line) lcolor(cranberry) lwidth(vthick)) /// 
				(connected `Y'_w wave if treatment_main == 4, msymbol(o) mcolor(emidblue) lpattern(line) lcolor(emidblue) lwidth(thick)) /// 
				(rcap h_`Y'_w l_`Y'_w wave), /// 
				xline(1.2) ///
				title("`l_`Y''", size(medium)) ///
				ytitle("") ///
				xtitle("") ///
				xscale(r(0.8 4.2)) ///
				yscale(r(0.2 0.8)) ///
				ylabel(0.2(0.2)0.8) ///
				xlabel(1 "Baseline" 2 "Midline" 4 "Endline", noticks) ///
				legend(order(1 "Control" 2 "Access" 3 "Access + Encour." 4 "Existing users") col(4) size(small)) ///
				graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
				saving(figure_panel_`Y', replace)
				
		graph export "figure_panel_`Y'.pdf", replace
		restore
		
	}
	
	
	foreach Y in wtp_vpn {
		
		preserve
		keep if panelmerged_wave2 == 1 & panelmerged_wave3 == 1
		collapse (mean) `Y'_w3 `Y'_w2 `Y'_w1 (sd) s_`Y'_w3 = `Y'_w3 s_`Y'_w2 = `Y'_w2 s_`Y'_w1 = `Y'_w1 (count) n_`Y'_w3 = `Y'_w3 n_`Y'_w2 = `Y'_w2 n_`Y'_w1 = `Y'_w1, by(treatment_main)
		gen h_`Y'_w3 = `Y'_w3 + invttail(n_`Y'_w3-1,0.05)*(s_`Y'_w3 / sqrt(n_`Y'_w3))
		gen l_`Y'_w3 = `Y'_w3 - invttail(n_`Y'_w3-1,0.05)*(s_`Y'_w3 / sqrt(n_`Y'_w3))	
		gen h_`Y'_w2 = `Y'_w2 + invttail(n_`Y'_w2-1,0.05)*(s_`Y'_w2 / sqrt(n_`Y'_w2))
		gen l_`Y'_w2 = `Y'_w2 - invttail(n_`Y'_w2-1,0.05)*(s_`Y'_w2 / sqrt(n_`Y'_w2))	
		gen h_`Y'_w1 = `Y'_w1 + invttail(n_`Y'_w1-1,0.05)*(s_`Y'_w1 / sqrt(n_`Y'_w1))
		gen l_`Y'_w1 = `Y'_w1 - invttail(n_`Y'_w1-1,0.05)*(s_`Y'_w1 / sqrt(n_`Y'_w1))	
		drop s_* n_*
		drop if treatment_main == .
		
		reshape long `Y'_w h_`Y'_w l_`Y'_w, i(treatment_main) j(wave)
		replace wave = 4 if wave == 3
		
		twoway 	(connected `Y'_w wave if treatment_main == 1, msymbol(o) mcolor(gs9) lpattern(dash) lcolor(gs9) lwidth(medthick)) /// 
				(connected `Y'_w wave if treatment_main == 2, msymbol(d) mcolor(gs3) lpattern(dash) lcolor(gs3) lwidth(medthick)) /// 
				(connected `Y'_w wave if treatment_main == 3, msymbol(o) mcolor(cranberry) lpattern(line) lcolor(cranberry) lwidth(vthick)) /// 
				(connected `Y'_w wave if treatment_main == 4, msymbol(o) mcolor(emidblue) lpattern(line) lcolor(emidblue) lwidth(thick)) /// 
				(rcap h_`Y'_w l_`Y'_w wave), /// 
				xline(1.2) ///
				title("`l_`Y''", size(medium)) ///
				ytitle("") ///
				xtitle("") ///
				xscale(r(0.8 4.2)) ///
				yscale(r(15 35)) ///
				ylabel(18 "3" 24 "4" 30 "5" 36 "6") ///				
				xlabel(1 "Baseline" 2 "Midline" 4 "Endline", noticks) ///
				legend(order(1 "Control" 2 "Access" 3 "Access + Encour." 4 "Existing users") col(2) size(small)) ///
				graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
				saving(figure_panel_`Y', replace)
				
		graph export "figure_panel_`Y'.pdf", replace
		restore
		
	}


	foreach Y in protest_2011_tmrw_parade {
		
		preserve
		keep if panelmerged_wave2 == 1 & panelmerged_wave3 == 1
		collapse (mean) `Y'_w3 `Y'_w2 `Y'_w1 (sd) s_`Y'_w3 = `Y'_w3 s_`Y'_w2 = `Y'_w2 s_`Y'_w1 = `Y'_w1 (count) n_`Y'_w3 = `Y'_w3 n_`Y'_w2 = `Y'_w2 n_`Y'_w1 = `Y'_w1, by(treatment_main)
		gen h_`Y'_w3 = `Y'_w3 + invttail(n_`Y'_w3-1,0.05)*(s_`Y'_w3 / sqrt(n_`Y'_w3))
		gen l_`Y'_w3 = `Y'_w3 - invttail(n_`Y'_w3-1,0.05)*(s_`Y'_w3 / sqrt(n_`Y'_w3))	
		gen h_`Y'_w2 = `Y'_w2 + invttail(n_`Y'_w2-1,0.05)*(s_`Y'_w2 / sqrt(n_`Y'_w2))
		gen l_`Y'_w2 = `Y'_w2 - invttail(n_`Y'_w2-1,0.05)*(s_`Y'_w2 / sqrt(n_`Y'_w2))	
		gen h_`Y'_w1 = `Y'_w1 + invttail(n_`Y'_w1-1,0.05)*(s_`Y'_w1 / sqrt(n_`Y'_w1))
		gen l_`Y'_w1 = `Y'_w1 - invttail(n_`Y'_w1-1,0.05)*(s_`Y'_w1 / sqrt(n_`Y'_w1))	
		drop s_* n_*
		drop if treatment_main == .
		
		reshape long `Y'_w h_`Y'_w l_`Y'_w, i(treatment_main) j(wave)
		replace wave = 4 if wave == 3
		
		twoway 	(connected `Y'_w wave if treatment_main == 1, msymbol(o) mcolor(gs9) lpattern(dash) lcolor(gs9) lwidth(medthick)) /// 
				(connected `Y'_w wave if treatment_main == 2, msymbol(d) mcolor(gs3) lpattern(dash) lcolor(gs3) lwidth(medthick)) /// 
				(connected `Y'_w wave if treatment_main == 3, msymbol(o) mcolor(cranberry) lpattern(line) lcolor(cranberry) lwidth(vthick)) /// 
				(connected `Y'_w wave if treatment_main == 4, msymbol(o) mcolor(emidblue) lpattern(line) lcolor(emidblue) lwidth(thick)) /// 
				(rcap h_`Y'_w l_`Y'_w wave), /// 
				xline(1.2) ///
				title("`l_`Y''", size(medium)) ///
				ytitle("") ///
				xtitle("") ///
				xscale(r(0.8 4.2)) ///
				xlabel(1 "Baseline" 2 "Midline" 4 "Endline", noticks) ///
				yscale(r(0 1)) ///
				ylabel(0(0.2)1) ///
				legend(order(1 "Control" 2 "Access" 3 "Access + Encour." 4 "Existing users") col(4) size(small)) ///
				graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
				saving(figure_panel_`Y', replace)
				
		graph export "figure_panel_`Y'.pdf", replace
		restore
		
	}
	
	// combine graphs: wave 1+2+3
	grc1leg figure_panel_info_freq_website_for.gph figure_panel_az_belief_media_value.gph figure_panel_az_belief_media_trust.gph figure_panel_az_belief_media_justif_pf.gph figure_panel_bias_domestic.gph figure_panel_bias_foreign.gph ///
			figure_panel_news_perccor_cen.gph figure_panel_news_perccor_unc.gph figure_panel_protest_pcheard_china.gph figure_panel_protest_pcheard_foreign.gph figure_panel_az_knowledge_meta.gph ///
			figure_panel_az_belief_econ_perf_cn.gph figure_panel_az_belief_econ_conf_cn.gph ///
			figure_panel_az_belief_instchange.gph figure_panel_az_belief_trust_govt.gph figure_panel_az_belief_trust_foreign.gph figure_panel_az_belief_evalgovt.gph figure_panel_az_belief_democracy_fp.gph figure_panel_az_belief_willing.gph /// 
			figure_panel_az_var_socialinteract.gph figure_panel_az_var_polparticipation.gph figure_panel_plan_grad_foreignmaster.gph figure_panel_cloc_for.gph figure_panel_stock_participation.gph, ///
		cols(6) holes(12 15 16 17 18) ///
		scale(0.75) ///
		ysize(17) xsize(20) ///
		legendfrom(figure_panel_az_belief_media_value.gph) ///
		graphregion(fcolor(white) ilcolor(white) lcolor(white))
	graph export "figure_panel_w123.pdf", replace
	




*** Figure A.10, A.14, A.15, A.16

	preserve
	qui {
	local v = 1
	foreach Y of varlist `var_media_valuation_reg_w3' az_belief_media_value_w3 `var_media_trust_reg_w3' az_belief_media_trust_w3 `var_censor_level_reg_w3' `var_censor_justif_reg_w3' az_belief_media_justif_w3 `var_censor_driver_reg_dom_w3' `var_censor_driver_reg_for_w3' `var_knowledge_news_reg_cen_w3' news_perccor_cen_w3 `var_knowledge_news_reg_unc_w3' news_perccor_unc_w3 `var_knowledge_prot_reg_chi_w3' protest_pcheard_china_w3 `var_knowledge_prot_reg_for_w3' protest_pcheard_foreign_w3 `var_knowledge_prot_reg_fak_w3' `var_knowledge_meta_reg_w3' az_knowledge_meta_w3 `var_econ_guess_reg_cn_perf_w3' az_belief_econ_perf_cn_w3 `var_econ_guess_reg_cn_conf_w3' az_belief_econ_conf_cn_w3 `var_econ_guess_reg_us_perf_w3' az_belief_econ_perf_us_w3 `var_econ_guess_reg_us_conf_w3' az_belief_econ_conf_us_w3 `var_demand_change_reg_w3' az_belief_instchange_w3 `var_trust_inst_reg_govt_w3' az_belief_trust_govt_w3 `var_trust_inst_reg_foreign_w3' az_belief_trust_foreign_w3 `var_eval_govt_reg_w3' az_belief_evalgovt_w3  `var_democracy_reg_w3' `var_willing_fight_reg_w3' az_belief_willing_w3 `var_vpn_purchase_reg_w3' `var_info_ranking_reg_w3' `var_info_freq_reg_w3' `var_socialinteract_reg_w3' `var_polparticipation_reg_w3' `var_stock_invest_reg_w3' `var_planaftergrad_reg_w3' `var_career_sector_reg_w3' `var_career_loc_reg_w3' {
		qui xi: reg `Y' treatment_vpnonly treatment_vpnnl, r
		su `Y' if e(sample)
		gen eb = (`Y' - r(mean))/r(sd)
		qui xi: reg eb treatment_nlonly if treatment_control == 1 | treatment_nlonly == 1, r
		gen rnl_`v' = _b[treatment_nlonly]
		gen rse_`v' = _se[treatment_nlonly]
		gen rcih_`v' = _b[treatment_nlonly] + 1.9*_se[treatment_nlonly]
		gen rcil_`v' = _b[treatment_nlonly] - 1.9*_se[treatment_nlonly]
		drop eb
		local v = `v' + 1
		}

	keep rnl_* rse_* rcil_* rcih_*
	duplicates drop
	gen i = 1
	reshape long rnl_ rse_ rcih_ rcil_, i(i) j(vnum)
	label drop _all

	expand 2
	gen p = _n > 117
	forval v = 1/117 {
		gen ci`v' = rcil_ if p == 0 & vnum == `v'
		replace ci`v' = rcih_ if p == 1 & vnum == `v'
		replace rnl_ = . if p == 1 & vnum == `v'
		}

	replace vnum = vnum + 35 if inrange(vnum,110,117) 	// E.7.b
	replace vnum = vnum + 34 if inrange(vnum,102,109) 	// E.7.a
	replace vnum = vnum + 33 if inrange(vnum,97,101) 	// E.6
	replace vnum = vnum + 32 if inrange(vnum,96,96) 	// E.5
	replace vnum = vnum + 31 if inrange(vnum,93,95) 	// E.4
	replace vnum = vnum + 30 if inrange(vnum,91,92) 	// E.3
	replace vnum = vnum + 29 if inrange(vnum,90,90) 	// E.2.b
	replace vnum = vnum + 28 if inrange(vnum,85,89) 	// E.2.a
	replace vnum = vnum + 27 if inrange(vnum,83,84) 	// E.1
	replace vnum = vnum + 25 if inrange(vnum,79,82) 	// D.8
	replace vnum = vnum + 24 if inrange(vnum,78,78) 	// D.6
	replace vnum = vnum + 23 if inrange(vnum,74,77) 	// D.3
	replace vnum = vnum + 22 if inrange(vnum,71,73) 	// D.2.b
	replace vnum = vnum + 21 if inrange(vnum,67,70) 	// D.2.a
	replace vnum = vnum + 20 if inrange(vnum,64,66) 	// D.1
	replace vnum = vnum + 18 if inrange(vnum,61,63) 	// C.4
	replace vnum = vnum + 17 if inrange(vnum,58,60) 	// C.3
	replace vnum = vnum + 16 if inrange(vnum,55,57) 	// C.2
	replace vnum = vnum + 15 if inrange(vnum,52,54) 	// C.1
	replace vnum = vnum + 13 if inrange(vnum,49,51) 	// B.5
	replace vnum = vnum + 12 if inrange(vnum,48,48) 	// B.3.c
	replace vnum = vnum + 11 if inrange(vnum,42,47) 	// B.3.b
	replace vnum = vnum + 10 if inrange(vnum,37,41) 	// B.3.a
	replace vnum = vnum + 9 if inrange(vnum,32,36) 		// B.2.b
	replace vnum = vnum + 8 if inrange(vnum,24,31) 		// B.2.a
	replace vnum = vnum + 6 if inrange(vnum,20,23) 		// A.5.b
	replace vnum = vnum + 5 if inrange(vnum,16,19) 		// A.5.a
	replace vnum = vnum + 4 if inrange(vnum,10,15) 		// A.4
	replace vnum = vnum + 3 if inrange(vnum,9,9) 		// A.3.b
	replace vnum = vnum + 2 if inrange(vnum,8,8) 		// A.3.a
	replace vnum = vnum + 1 if inrange(vnum,4,7) 		// A.2

	label define varname ///
		1 "Willingness to pay for circumvention tool" ///
		2 "Value added of foreign media access" ///
		3 "{it: z-score: valuation of access to foreign media outlets}" ///
		5 "Distrust in domestic state-owned media" ///
		6 "Distrust in domestic privately-owned media" ///
		7 "Trust in foreign media" ///
		8 "{it: z-score: trust in non-domestic media outlets}" ///
		10 "Degree of censorship on domestic news outlets" ///
		12 "Degree of censorship on foreign news outlets" ///
		14 "Unjustified: censoring economic news" ///
		15 "Unjustified: censoring political news" ///
		16 "Unjustified: censoring social news" ///
		17 "Unjustified: censoring foreign news" ///
		18 "Unjustified: censoring pornography" ///
		19 "{it: z-score: censorship unjustified}" ///
		21 "Domestic cens. driven by govt. policies" ///
		22 "Domestic cens. driven by corp. interest" ///
		23 "Domestic cens. driven by media???s ideology" ///
		24 "Domestic cens. driven by readers??? demand" ///
		26 "Foreign cens. driven by govt. policies" ///
		27 "Foreign cens. driven by corp. interest" ///
		28 "Foreign cens. driven by media???s ideology" ///
		29 "Foreign cens. driven by readers??? demand" ///
		32 "Steel production reduction reaches target" ///
		33 "Trump registered trademarks in China" ///
		34 "Jianhua Xiao kidnapped in Hong Kong" ///
		35 "Xinjiang installed GPS on all automobiles" ///
		36 "China and Norway re-normalize ties" ///
		37 "Feminist groups fight women's rights" ///
		38 "Carrie Lam becomes HK Chief Executive" ///
		39 "{it: % quizzes answered correctly: poli. sensitive news}" ///
		41 "China stops importing coal from North Korea" ///
		42 "H7N9 influenza epidemic" ///
		43 "Transnational railway in Ethiopia" ///
		44 "Foreign reserves fall below threshold" ///
		45 "{it: % quizzes answered correctly: nonsensitive news}" ///
		47 "2012 HK Anti-National Curr. Movement" ///
		48 "2014 HK Umbrella Revolution" ///
		49 "2016 HK Mong Kok Revolution" ///
		50 "2014 Taiwan Sunflower Stud. Movement" ///
		51 "{it: % protests in Greater China heard of}" ///
		53 "2014 Ukrainian Euromaidan Revolution" ///
		54 "2010 Arab Spring" ///
		55 "2014 Crimean Status Referendum" ///
		56 "2010 Catalonian Indep. Movement" ///
		57 "2017 Women's March" ///
		58 "{it: % foreign protests heard of}" ///
		60 "2011 Tomorrow Revolution [fake]" ///
		62 "Informedness of issues in China" ///
		63 "Greater informedness than peers" ///
		64 "{it: z-score: self-assessment of knowledge level}" ///
		67 "Guess on GDP growth rate in 2016 China" ///
		68 "Guess on SSCI by end of 2016" ///
		69 "{it: z-score: optimistic belief of Chinese economy}" ///
		71 "Confidence of China GDP guess" ///
		72 "Confidence of SSCI guess" ///
		73 "{it: z-score: confidence of guesses on Chinese economy}" ///
		75 "Guess on GDP growth rate in 2016 US" ///
		76 "Guess on DJI by end of 2016" ///
		77 "{it: z-score: optimistic belief of US economy}" ///
		79 "Confidence of US GDP guess" ///
		80 "Confidence of DJI guess" ///
		81 "{it: z-score: confidence of guesses on US economy}" ///
		84 "Economic system needs changes" ///
		85 "Political system needs changes" ///
		86 "{it: z-score: demand for institutional change}" ///
		88 "Trust in central govt. of China" ///
		89 "Trust in provincial govt. of China" ///
		90 "Trust in local govt. of China" ///
		91 "{it: z-score: trust in Chinese govt.}" ///
		93 "Trust in central govt. of Japan" ///
		94 "Trust in federal govt. of US" ///
		95 "{it: z-score: trust in foreign govt.}" ///
		97 "Satisfaction of economic dev." ///
		98 "Satisfaction of domestic politics" ///
		99 "Satisfaction of diplomatic affairs" ///
		100 "{it: z-score: satisfaction of govt???s performance}" ///
		102 "Living in democracy is not important" ///
		104 "Willing to battle illegal govt. acts" ///
		105 "Willing to report govt. misconduct" ///
		106 "Willing to stand up for the weak" ///
		107 "{it: z-score: willingness to act}" ///
		110 "Purchase discounted tool we offered" ///
		111 "Purchase any tool" ///
		113 "Ranked high: domestic websites" ///
		114 "Ranked high: foreign websites" ///
		115 "Ranked high: domestic social media" ///
		116 "Ranked high: foreign social media" ///
		117 "Ranked high: word of mouth" ///
		119 "Frequency of visiting foreign websites for info." ///
		121 "Frequency of discussing poli. with friends" ///
		122 "Frequency of persuading others" ///
		124 "Protests concerning social issues" ///
		125 "Plan to vote for local PCR" ///
		126 "Complain to school authorities" ///
		128 "Currently invested in Chinese stock mkt." ///
		130 "Plan: grad. school in China" ///
		131 "Plan: master degree abroad" ///
		132 "Plan: PhD degree abroad" ///
		133 "Plan: military in China" ///
		134 "Plan: work right away" ///
		136 "Sector pref.: national civil service" ///
		137 "Sector pref.: local civil service" ///
		138 "Sector pref.: military" ///
		139 "Sector pref.: private firm in China" ///
		140 "Sector pref.: foreign firm in China" ///
		141 "Sector pref.: SOEs" ///
		142 "Sector pref.: inst. organizations" ///
		143 "Sector pref.: entrepreneurship" ///
		145 "Location pref.: Beijing" ///
		146 "Location pref.: Shanghai" ///
		147 "Location pref.: Guangzhou and Shenzhen" ///
		148 "Location pref.: tier 2 cities in central" ///
		149 "Location pref.: other cities in China" ///
		150 "Location pref.: HK and Macau" ///
		151 "Location pref.: Taiwan" ///
		152 "Location pref.: foreign cities" 
		
	label values vnum varname
	scatter vnum ci1, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || ///
	scatter vnum ci2, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || ///
	scatter vnum ci3, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci4, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci5, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci6, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci7, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci8, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci9, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci10, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci11, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci12, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci13, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci14, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci15, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci16, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci17, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci18, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci19, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci20, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci21, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci22, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci23, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci24, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci25, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci26, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci27, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci28, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci29, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci30, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci31, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci32, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci33, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci34, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci35, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci36, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci37, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci38, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci39, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci40, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci41, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci42, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci43, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci44, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci45, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci46, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci47, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci48, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci49, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci50, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci51, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci52, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci53, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci54, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci55, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci56, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci57, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci58, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci59, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci60, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci61, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci62, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci63, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci64, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci65, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci66, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci67, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci68, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci69, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci70, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci71, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci72, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci73, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci74, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci75, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci76, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci77, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci78, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci79, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci80, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci81, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci82, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci83, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci84, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci85, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci86, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci87, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci88, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci89, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci90, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci91, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci92, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci93, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci94, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci95, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci96, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci97, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci98, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci99, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci100, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci101, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci102, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci103, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci104, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci105, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci106, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci107, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci108, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci109, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci110, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci111, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci112, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci113, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci114, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci115, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci116, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum ci117, c(l) mcolor(gs0) lcolor(gs9) lwidth(thin) m(i) || /// 
	scatter vnum rnl_ if vnum==1, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==2, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==3, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==5, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==6, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==7, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==8, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==10, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==12, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==14, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==15, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==16, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==17, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==18, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==19, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==21, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==22, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==23, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==24, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==26, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==27, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==28, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==29, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==32, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==33, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==34, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==35, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==36, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==37, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==38, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==39, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==41, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==42, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==43, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==44, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==45, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==47, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==48, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==49, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==50, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==51, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==53, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==54, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==55, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==56, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==57, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==58, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==60, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==62, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==63, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==64, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==67, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==68, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==69, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==71, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==72, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==73, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==75, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==76, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==77, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==79, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==80, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==81, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==84, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==85, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==86, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==88, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==89, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==90, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==91, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==93, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==94, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==95, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==97, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==98, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==99, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==100, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==102, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==104, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==105, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==106, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==107, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==110, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==111, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==113, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==114, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==115, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==116, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==117, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==119, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==121, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==122, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==124, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==125, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==126, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==128, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==130, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==131, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==132, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==133, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==134, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==136, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==137, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==138, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==139, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==140, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==141, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==142, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==143, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==145, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==146, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==147, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==148, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==149, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==150, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==151, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) || ///
	scatter vnum rnl_ if vnum==152, mcolor(gs11) mlcolor(gs3) mlwidth(medthin) m(S) ///
		xline(0,   lc(midblue) lp(line)) ///
		yline(4,   lc(gs5)  lp(longdash)) ///
		yline(9,   lc(gs5)  lp(longdash)) ///
		yline(11,  lc(gs10) lp(dash)) ///
		yline(20,  lc(gs5)  lp(longdash)) ///
		yline(13,  lc(gs5)  lp(longdash)) ///
		yline(25,  lc(gs10) lp(dash)) ///
		yline(31,  lc(gs3)  lp(line) lwidth(medthick)) ///
		yline(40,  lc(gs10) lp(dash)) ///
		yline(46,  lc(gs5)  lp(longdash)) ///
		yline(52,  lc(gs10) lp(dash)) ///
		yline(59,  lc(gs10) lp(dash)) ///
		yline(61,  lc(gs5)  lp(longdash)) ///
		yline(66,  lc(gs3)  lp(line) lwidth(medthick)) ///
		yline(70,  lc(gs5)  lp(longdash)) ///
		yline(74,  lc(gs5)  lp(longdash)) ///
		yline(78,  lc(gs5)  lp(longdash)) ///
		yline(83,  lc(gs3)  lp(line) lwidth(medthick)) ///
		yline(87,  lc(gs5)  lp(longdash)) ///
		yline(92,  lc(gs10) lp(dash)) ///
		yline(96,  lc(gs5)  lp(longdash)) ///
		yline(101, lc(gs5)  lp(longdash)) ///
		yline(103, lc(gs5)  lp(longdash)) ///
		yline(109, lc(gs3)  lp(line) lwidth(medthick)) ///
		yline(112, lc(gs5)  lp(longdash)) ///
		yline(118, lc(gs12) lp(dash)) ///
		yline(120, lc(gs5)  lp(longdash)) ///
		yline(123, lc(gs5)  lp(longdash)) ///
		yline(127, lc(gs5)  lp(longdash)) ///
		yline(129, lc(gs5)  lp(longdash)) ///
		yline(135, lc(gs5)  lp(longdash)) ///
		yline(144, lc(gs12) lp(dash)) ///
		ytitle("") yscale(reverse) ///
		xtitle("Standardized means (Control = 0)", size(small)) xlabel(, grid labsize(small)) ///
		ylab(1 2 3 5 6 7 8 10 12 14 15 16 17 18 19 21 22 23 24 26 27 28 29 32 33 34 35 36 37 38 39 41 42 43 44 45 47 48 49 50 51 53 54 55 56 57 58 60 62 63 64 67 68 69 71 72 73 75 76 77 79 80 81 84 85 86 88 89 90 91 93 94 95 97 98 99 100 102 104 105 106 107 110 111 113 114 115 116 117 119 121 122 124 125 126 128 130 131 132 133 134 136 137 138 139 140 141 142 143 145 146 147 148 149 150 151 152, valuelabel angle(0) labsize(vsmall) nogrid) ///
		ysize(20) xsize(4) ///
		legend(order(118 "Control + Encour.") size(vsmall)) ///
		graphregion(color(white))
	}
	
	graph export figure_dotplot_master_w3_cce.pdf, replace
	restore
	

	

*** Figure A.17

	keep if panelmerged_wave3 == 1

	foreach Y in az_overall_all_w3 {
		
		preserve
		
		local i = 0
		
		// all subject
		qui reg `Y' treatment_vpnnl if treatment_user != 1, r
		gen b_`i' = _b[treatment_vpnnl]
		gen se_`i' = _se[treatment_vpnnl]
		gen u_`i' = b_`i' + 1.68 * se_`i'
		gen l_`i' = b_`i' - 1.68 * se_`i'		
		local i = `i' + 2
		
		// cut-offs
		foreach h in gender birth_year residence_coastal hukou_urban university_elite hs_track_science department_ssh domestic_english_atleast4 foreign_english_yes travel_hktaiwan travel_foreign_yes father_edu_hsabove work_father_govt father_ccp mother_edu_hsabove work_mother_govt mother_ccp hh_income az_preference_risk az_preference_time az_preference_altruism az_preference_reciprocity az_overall_a_w1 az_overall_b_w1 az_overall_c_w1 az_overall_d_w1 az_overall_e_w1 {
			
			qui reg `Y' treatment_vpnnl if treatment_user != 1 & h_`h' == 0, r
			gen b_`i' = _b[treatment_vpnnl]
			gen se_`i' = _se[treatment_vpnnl]
			gen u_`i' = b_`i' + 1.68 * se_`i'
			gen l_`i' = b_`i' - 1.68 * se_`i'
			local i = `i' + 1
			
			qui reg `Y' treatment_vpnnl if treatment_user != 1 & h_`h' == 1, r
			gen b_`i' = _b[treatment_vpnnl]
			gen se_`i' = _se[treatment_vpnnl]
			gen u_`i' = b_`i' + 1.68 * se_`i'
			gen l_`i' = b_`i' - 1.68 * se_`i'
			local i = `i' + 2
			
			}
			
		keep b_* se_* u_* l_*
		duplicates drop
		
		local b_all = b_0
		
		gen ii = 1
		reshape long b_ se_ u_ l_, i(ii) j(i)
		drop ii
		
		rename b_ b
		rename se_ se
		rename u_ u
		rename l_ l
		
		label def i ///
			0 "All subjects" ///
			2 "Female {it:vs.} male" ///
			5 "Lower class {it:vs.} upper class" ///
			8 "Non-coastal {it:vs.} coastal" ///
			11 "Rural {it:vs.} urban" ///
			14 "2nd tier {it:vs.} elite univ." ///
			17 "Humanities {it:vs.} science track" ///
			20 "Sc/Eng {it:vs.} SocS/Hum major" ///
			23 "Not passed {it:vs.} at least Eng Level 4" ///
			26 "Not taken {it:vs.} taken TOEFL/IELTS" ///
			29 "Not been {it:vs.} been to HK/TW" ///
			32 "Not been {it:vs.} been abroad" ///
			35 "Father below {it:vs.} above hs" ///
			38 "Father not work {it:vs.} work for govt." ///
			41 "Father not {it:vs.} is CCP member" ///
			44 "Mother below {it:vs.} above hs" ///
			47 "Mother not work {it:vs.} work for govt." ///
			50 "Mother not {it:vs.} is CCP member" ///
			53 "HH income < {it:vs.} > median" ///
			56 "Risk pref. < {it:vs.} > median" ///
			59 "Time pref. < {it:vs.} > median" ///
			62 "Altruism < {it:vs.} > median" ///
			65 "Recipro. < {it:vs.} > median" ///
			68 "(A) media-related < {it:vs.} > median" ///
			71 "(B) knowledge < {it:vs.} > median" ///
			74 "(C) economic beliefs < {it:vs.} > median" ///
			77 "(D) pol. attitudes < {it:vs.} > median" ///
			80 "(E) behaviors < {it:vs.} > median"
					
		label value i i
		
		twoway 	(scatter b i if i > 0, msize(large)) ///
				(scatter b i if i == 0, msize(vlarge) msymbol(square) mcolor(cranberry)) ///	
				(rcap u l i), ///
				yline(`b_all', lpattern(dash)) ///
				xline(1, lcolor(gs2) lwidth(medthick)) ///
				xline(13 22 34 55 67, lcolor(gs6) lwidth(medthick)) ///
				xline(4 7 10 16 19 25 28 31 37 40 43 46 49 52 58 61 64 70 73 76 79, lcolor(gs8) lpattern(longdash)) ///
				ytitle("Coefficient on Group-AE effect") ///
				xtitle("") ///
				xlabel(0 2 5 8 11 14 17 20 23 26 29 32 35 38 41 44 47 50 53 56 59 62 65 68 71 74 77 80, notick valuelabel angle(vertical)) ///
				xscale(r(-0.5 81)) ///
				scale(0.8) ///
				xsize(20) ysize(7) ///
				legend(off) ///
				graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
				saving(figure_heterogeneity_`Y', replace)
		graph export "figure_heterogeneity_`Y'.pdf", replace
	
		restore
		
	}	
	


	
*** Figure A.18

	// wave 2 outcomes
	foreach Y in panamapapers tenyearshk stockcrash economistcensor {
		preserve
		rename soclearning_rm_new_w2 	soclearning_rm_new
		keep news_c_`Y' soclearning_ownaccess soclearning_rm_new soclearning_rm_existing
		keep if soclearning_rm_existing == 0 & soclearning_rm_new < 2
		nl (news_c_`Y' = {alpha} + soclearning_ownaccess*{p} + (1-(({alpha}+{p})*(1-(1-soclearning_ownaccess)*{q0}-soclearning_ownaccess*{q1}) + (1-{alpha}-{p}))^soclearning_rm_new)), initial(alpha 0.5 p 0.25 q0 0.2 q1 0.1) nolog
		gen news = "`l_`Y''"
		gen alpha = _b[/alpha]
		gen p = _b[/p]
		gen q0 = _b[/q0]
		gen q1 = _b[/q1]
		keep news alpha p q0 q1
		duplicates drop
		save "slest_`Y'.dta", replace
		restore
		}
		
	// wave 3 outcomes
	foreach Y in coalprod trumpchina xiaojianhua xijiangcar chinanorway womenrights hkceelection {
		preserve
		rename soclearning_rm_new_w3 	soclearning_rm_new
		keep news_c_`Y' soclearning_ownaccess soclearning_rm_new soclearning_rm_existing
		keep if soclearning_rm_existing == 0 & soclearning_rm_new < 2
		nl (news_c_`Y' = {alpha} + soclearning_ownaccess*{p} + (1-(({alpha}+{p})*(1-(1-soclearning_ownaccess)*{q0}-soclearning_ownaccess*{q1}) + (1-{alpha}-{p}))^soclearning_rm_new)), initial(alpha 0.5 p 0.25 q0 0.2 q1 0.1) nolog
		gen news = "`l_`Y''"
		gen alpha = _b[/alpha]
		gen p = _b[/p]
		gen q0 = _b[/q0]
		gen q1 = _b[/q1]
		keep news alpha p q0 q1
		duplicates drop
		save "slest_`Y'.dta", replace
		restore
		}	
	
	// append
	preserve
	clear all
	gen news = ""
	gen alpha = .
	gen p = .
	gen q0 = .
	gen q1 = .
	save "slest_compiled.dta", replace
	restore
	
	preserve
	foreach Y in panamapapers tenyearshk stockcrash economistcensor coalprod trumpchina xiaojianhua xijiangcar chinanorway womenrights hkceelection {
		use "slest_compiled.dta", clear
		append using "slest_`Y'.dta"	
		save "slest_compiled.dta", replace
		erase "slest_`Y'.dta"
		}
	restore
	
	// plot
	preserve
	
	use "slest_compiled.dta", clear
	sort p
	gen i = _n
	
	gen p_l = p-0.003
	gen p_h = p+0.003
	gen p_u = p+0.009
	gen p_v = p+0.012
	
	twoway 	(rcapsym q0 q1 p, m(D) mc(emidblue) lpattern(dash) lcolor(gs8) lwidth(medium)) ///
			(rcapsym q0 q0 p, m(O) mc(cranberry) msize(large)) ///
			(scatter q0 p if i > 6 & i < 9, m(o) mc(none) mlabel(news) mlabposition(12) mlabgap(*25) mlabangle(vertical) mlabcolor(gs6)) ///
			(scatter q1 p if i == 1, m(o) mc(none) mlabel(news) mlabposition(12) mlabgap(*22) mlabangle(vertical) mlabcolor(gs6)) ///
			(scatter q0 p if i == 2, m(o) mc(none) mlabel(news) mlabposition(12) mlabgap(*22) mlabangle(vertical) mlabcolor(gs6)) ///
			(scatter q0 p_u if i == 3, m(o) mc(none) mlabel(news) mlabposition(12) mlabgap(*18) mlabangle(vertical) mlabcolor(gs6)) ///
			(scatter q0 p if i == 4, m(o) mc(none) mlabel(news) mlabposition(12) mlabgap(*28) mlabangle(vertical) mlabcolor(gs6)) ///
			(scatter q0 p_h if i == 5, m(o) mc(none) mlabel(news) mlabposition(12) mlabgap(*25) mlabangle(vertical) mlabcolor(gs6)) ///
			(scatter q0 p if i == 6, m(o) mc(none) mlabel(news) mlabposition(12) mlabgap(*20) mlabangle(vertical) mlabcolor(gs6)) ///
			(scatter q0 p_l if i == 9, m(o) mc(none) mlabel(news) mlabposition(12) mlabgap(*22) mlabangle(vertical) mlabcolor(gs6)) ///
			(scatter q0 p_v if i == 10, m(o) mc(none) mlabel(news) mlabposition(12) mlabgap(*22) mlabangle(vertical) mlabcolor(gs6)) ///
			(scatter q0 p_h if i == 11, m(o) mc(none) mlabel(news) mlabposition(12) mlabgap(*18) mlabangle(vertical) mlabcolor(gs6)) ///
			(scatter q1 p if i == 1, m(D) mc(emidblue)) ///
			(scatter q0 p if i == 1, m(O) mc(cranberry) msize(large)), ///
			yscale(r(-0.1 0.45)) ///
			ylabel(-0.1(0.1)0.4) ///
			ytitle("Social transmission rates: {it:q}") ///
			xscale(r(0.05 0.35)) ///
			xlabel(0.05(0.05)0.35) ///
			xtitle("Direct learning rate: {it:p}") ///
			legend(order(14 "Social transmission rate towards non-user: {it:q(I=0)}" 13 "Social transmission rate towards user: {it:q(I=1)}") col(1) size(small)) ///
			graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
			saving(figure_sociallearning_est, replace)
	graph export "figure_sociallearning_est.pdf", replace	

	restore
	
	

*** Figure A.19

	preserve
	drop if treatment_main == .
	bysort treatment_main: su news_time_lastclick, d
	bysort treatment_main: su news_time_totalclick, d
	bysort treatment_main: su pppr_time_lastclick, d
	bysort treatment_main: su pppr_time_totalclick, d
	restore
	
	// everyone

	preserve
	
	// labels
	label variable news_time_lastclick  	"Time (sec) spent on news quizzes"
	label variable news_time_totalclick 	"# clicks on news quizzes"
	label variable pppr_time_lastclick 		"Time (sec) spent on notable figures & protests"
	label variable pppr_time_totalclick 	"# clicks on notable figures & protests"

	// top and bottom code
	replace news_time_lastclick = 300 	if (news_time_lastclick > 300 & news_time_lastclick != .)
	replace news_time_totalclick = 20 	if (news_time_totalclick > 20 & news_time_totalclick != .)
	replace pppr_time_lastclick = 200 	if (pppr_time_lastclick > 200 & pppr_time_lastclick != .)
	replace pppr_time_lastclick = 20 	if (pppr_time_lastclick < 20 & pppr_time_lastclick != .)
	replace pppr_time_totalclick = 30 	if (pppr_time_totalclick > 30 & pppr_time_totalclick != .)
	
	foreach g in news_time_lastclick news_time_totalclick pppr_time_lastclick pppr_time_totalclick {
		
		local labelX: var label `g'
		
		graph box `g', ///
			over(treatment_main) ///
			ytitle("`labelX'") ///
			graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
			saving(figure_time_`g'_everyone, replace)	
		graph export "figure_time_`g'_everyone.pdf", replace
		
		}	
		
	restore

	// conditional on answering quiz correctly
	
	preserve
	
	// labels
	label variable news_time_lastclick  	"Time (sec) spent on news quizzes"
	label variable news_time_totalclick 	"# clicks on news quizzes"
	label variable pppr_time_lastclick 		"Time (sec) spent on notable figures & protests"
	label variable pppr_time_totalclick 	"# clicks on notable figures & protests"

	// top and bottom code
	replace news_time_lastclick = 300 	if (news_time_lastclick > 300 & news_time_lastclick != .)
	replace news_time_totalclick = 20 	if (news_time_totalclick > 20 & news_time_totalclick != .)
	replace pppr_time_lastclick = 200 	if (pppr_time_lastclick > 200 & pppr_time_lastclick != .)
	replace pppr_time_lastclick = 20 	if (pppr_time_lastclick < 20 & pppr_time_lastclick != .)
	replace pppr_time_totalclick = 30 	if (pppr_time_totalclick > 30 & pppr_time_totalclick != .)
	
	keep if news_totalcorrect_w2 >= 6
	
	foreach g in news_time_lastclick news_time_totalclick {
		
		local labelX: var label `g'
		
		graph box `g', ///
			over(treatment_main) ///
			ytitle("`labelX'") ///
			graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
			saving(figure_time_`g'_correct, replace)	
		graph export "figure_time_`g'_correct.pdf", replace
		
		}	
		
	restore

	
	preserve
	
	// labels
	label variable news_time_lastclick  	"Time (sec) spent on news quizzes"
	label variable news_time_totalclick 	"# clicks on news quizzes"
	label variable pppr_time_lastclick 		"Time (sec) spent on notable figures & protests"
	label variable pppr_time_totalclick 	"# clicks on notable figures & protests"

	// top and bottom code
	replace news_time_lastclick = 300 	if news_time_lastclick > 300
	replace news_time_totalclick = 20 	if news_time_totalclick > 20
	replace pppr_time_lastclick = 200 	if pppr_time_lastclick > 200
	replace pppr_time_lastclick = 20 	if pppr_time_lastclick < 20
	replace pppr_time_totalclick = 30 	if pppr_time_totalclick > 30
	
	keep if people_totalheard_w2 >= 3 & protest_totalheard_w2 >= 3
	
	foreach g in pppr_time_lastclick pppr_time_totalclick {
		
		local labelX: var label `g'
		
		graph box `g', ///
			over(treatment_main) ///
			ytitle("`labelX'") ///
			graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
			saving(figure_time_`g'_correct, replace)	
		graph export "figure_time_`g'_correct.pdf", replace
		
		}	
		
	restore
	
	// combine graphs
	graph combine figure_time_news_time_lastclick_everyone.gph figure_time_news_time_totalclick_everyone.gph, ///
		title("Panel A: among all participants") ///
		col(2) ///
		scale(1.2) ///
		ysize(5) xsize(13) ///
		graphregion(fcolor(white) ilcolor(white) lcolor(white))
	graph export "figure_time_news_everyone.pdf", replace
	
	graph combine figure_time_news_time_lastclick_everyone.gph figure_time_news_time_totalclick_everyone.gph, ///
		title("Panel B: among those who answered > half questions correctly") ///
		col(2) ///
		scale(1.2) ///
		ysize(5) xsize(13) ///
		graphregion(fcolor(white) ilcolor(white) lcolor(white))
	graph export "figure_time_news_correct.pdf", replace
	
	


*** Table A.1

	preserve
	
	tempname p
	postfile `p' str50 variable float mu_0 sd_0 N_0 ///
		mu_1 sd_1 N_1 ///
		mu_2 sd_2 N_2 ///
		mu_3 sd_3 N_3 ///
		mu_4 sd_4 N_4 ///
		mu_5 sd_5 N_5 ///
		f_stat p_val ///
		using summstats_w1, replace

	foreach X of varlist `var_demog_reg_personal' az_demographics_personal `var_demog_reg_education' az_demographics_education `var_demog_reg_english' az_demographics_english `var_demog_reg_travel' az_demographics_travel `var_demog_reg_household' az_demographics_household `var_preference_risk' az_preference_risk `var_preference_time' az_preference_time `var_preference_altruism' az_preference_altruism `var_preference_reciprocity' az_preference_reciprocity treatment_control treatment_vpnonly treatment_nlonly treatment_vpnnl treatment_user {
		qui sum `X'
		local N_0=r(N)
		local mu_0=r(mean)
		local sd_0=r(sd)

		qui sum `X' if treatment_master == 1
		local N_1=r(N)
		local mu_1=r(mean)
		local sd_1=r(sd)

		qui sum `X' if treatment_master == 2
		local N_2=r(N)
		local mu_2=r(mean)
		local sd_2=r(sd)

		qui sum `X' if treatment_master == 3
		local N_3=r(N)
		local mu_3=r(mean)
		local sd_3=r(sd)

		qui sum `X' if treatment_master == 4
		local N_4=r(N)
		local mu_4=r(mean)
		local sd_4=r(sd)

		qui sum `X' if treatment_master == 5
		local N_5=r(N)
		local mu_5=r(mean)
		local sd_5=r(sd)
		
		qui anova `X' treatment_master if treatment_master != 5
		local f_stat=e(F)
		local p_val=Ftail(e(df_m), e(df_r), e(F))
		
		post `p' ("`X'") /// 
		(`mu_0') (`sd_0') (`N_0') ///
		(`mu_1') (`sd_1') (`N_1') ///
		(`mu_2') (`sd_2') (`N_2') /// 
		(`mu_3') (`sd_3') (`N_3') /// 
		(`mu_4') (`sd_4') (`N_4') /// 
		(`mu_5') (`sd_5') (`N_5') ///
		(`f_stat') (`p_val')
	}
		
	postclose `p'
	restore
	
	preserve	
	clear
	use summstats_w1
	outsheet using summstats_w1.xls, replace
	clear
	erase summstats_w1.dta
	restore	
	
	

*** Table A.2

	preserve
	keep if panelmerged_wave2 == 1
	
	tempname p
	postfile `p' str50 variable float mu_0 sd_0 N_0 ///
		mu_1 sd_1 N_1 ///
		mu_2 sd_2 N_2 ///
		mu_3 sd_3 N_3 ///
		mu_4 sd_4 N_4 ///
		mu_5 sd_5 N_5 ///
		f_stat p_val ///
		using summstats_w2, replace

	foreach X of varlist `var_demog_reg_personal' az_demographics_personal `var_demog_reg_education' az_demographics_education `var_demog_reg_english' az_demographics_english `var_demog_reg_travel' az_demographics_travel `var_demog_reg_household' az_demographics_household `var_preference_risk' az_preference_risk `var_preference_time' az_preference_time `var_preference_altruism' az_preference_altruism `var_preference_reciprocity' az_preference_reciprocity treatment_control treatment_vpnonly treatment_nlonly treatment_vpnnl treatment_user {
		qui sum `X'
		local N_0=r(N)
		local mu_0=r(mean)
		local sd_0=r(sd)

		qui sum `X' if treatment_master == 1
		local N_1=r(N)
		local mu_1=r(mean)
		local sd_1=r(sd)

		qui sum `X' if treatment_master == 2
		local N_2=r(N)
		local mu_2=r(mean)
		local sd_2=r(sd)

		qui sum `X' if treatment_master == 3
		local N_3=r(N)
		local mu_3=r(mean)
		local sd_3=r(sd)

		qui sum `X' if treatment_master == 4
		local N_4=r(N)
		local mu_4=r(mean)
		local sd_4=r(sd)

		qui sum `X' if treatment_master == 5
		local N_5=r(N)
		local mu_5=r(mean)
		local sd_5=r(sd)
		
		qui anova `X' treatment_master if treatment_master != 5
		local f_stat=e(F)
		local p_val=Ftail(e(df_m), e(df_r), e(F))
		
		post `p' ("`X'") /// 
		(`mu_0') (`sd_0') (`N_0') ///
		(`mu_1') (`sd_1') (`N_1') ///
		(`mu_2') (`sd_2') (`N_2') /// 
		(`mu_3') (`sd_3') (`N_3') /// 
		(`mu_4') (`sd_4') (`N_4') /// 
		(`mu_5') (`sd_5') (`N_5') ///
		(`f_stat') (`p_val')
	}
		
	postclose `p'
	restore
	
	preserve	
	clear
	use summstats_w2
	outsheet using summstats_w2.xls, replace
	clear
	erase summstats_w2.dta
	restore	

	

*** Table A.3

	preserve
	keep if panelmerged_wave3 == 1
	
	tempname p
	postfile `p' str50 variable float mu_0 sd_0 N_0 ///
		mu_1 sd_1 N_1 ///
		mu_2 sd_2 N_2 ///
		mu_3 sd_3 N_3 ///
		mu_4 sd_4 N_4 ///
		mu_5 sd_5 N_5 ///
		f_stat p_val ///
		using summstats_w3, replace

	foreach X of varlist `var_demog_reg_personal' az_demographics_personal `var_demog_reg_education' az_demographics_education `var_demog_reg_english' az_demographics_english `var_demog_reg_travel' az_demographics_travel `var_demog_reg_household' az_demographics_household `var_preference_risk' az_preference_risk `var_preference_time' az_preference_time `var_preference_altruism' az_preference_altruism `var_preference_reciprocity' az_preference_reciprocity treatment_control treatment_vpnonly treatment_nlonly treatment_vpnnl treatment_user {
		qui sum `X'
		local N_0=r(N)
		local mu_0=r(mean)
		local sd_0=r(sd)

		qui sum `X' if treatment_master == 1
		local N_1=r(N)
		local mu_1=r(mean)
		local sd_1=r(sd)

		qui sum `X' if treatment_master == 2
		local N_2=r(N)
		local mu_2=r(mean)
		local sd_2=r(sd)

		qui sum `X' if treatment_master == 3
		local N_3=r(N)
		local mu_3=r(mean)
		local sd_3=r(sd)

		qui sum `X' if treatment_master == 4
		local N_4=r(N)
		local mu_4=r(mean)
		local sd_4=r(sd)

		qui sum `X' if treatment_master == 5
		local N_5=r(N)
		local mu_5=r(mean)
		local sd_5=r(sd)
		
		qui anova `X' treatment_master if treatment_master != 5
		local f_stat=e(F)
		local p_val=Ftail(e(df_m), e(df_r), e(F))
		
		post `p' ("`X'") /// 
		(`mu_0') (`sd_0') (`N_0') ///
		(`mu_1') (`sd_1') (`N_1') ///
		(`mu_2') (`sd_2') (`N_2') /// 
		(`mu_3') (`sd_3') (`N_3') /// 
		(`mu_4') (`sd_4') (`N_4') /// 
		(`mu_5') (`sd_5') (`N_5') ///
		(`f_stat') (`p_val')
	}
		
	postclose `p'
	restore
	
	preserve	
	clear
	use summstats_w3
	outsheet using summstats_w3.xls, replace
	clear
	erase summstats_w3.dta
	restore	

	

*** Table A.4

	preserve
	keep if panelmerged_wave2 == 1
	gen wave2 = 1
	drop treatment_master
	save panelmerged_wave2_temp, replace
	restore

	preserve
	keep if panelmerged_wave3 == 1
	gen wave3 = 1
	drop treatment_master
	save panelmerged_wave3_temp, replace
	restore
	
	preserve
	
	// duplicate for attrition test
	append using panelmerged_wave2_temp
	erase panelmerged_wave2_temp.dta
	replace wave2 = 0 if wave2 == .
	
	append using panelmerged_wave3_temp
	erase panelmerged_wave3_temp.dta
	replace wave3 = 0 if wave3 == . & wave2 != 1
	
	tempname p
	postfile `p' str50 variable float mu_w1 sd_w1 N_w1 ///
		mu_w2 sd_w2 N_w2 ///
		attrition_w2_pvalue ///
		mu_w3 sd_w3 N_w3 ///
		attrition_w3_pvalue ///
		using attrition_appendix, replace

	foreach X of varlist `var_info_ranking_reg_w1' `var_info_freq_reg_w1' `var_media_valuation_reg_w1' az_belief_media_value_w1 `var_media_trust_reg_w1' az_belief_media_trust_w1 `var_censor_level_reg_w1' `var_censor_justif_reg_w1' az_belief_media_justif_w1 `var_censor_driver_reg_dom_w1' `var_censor_driver_reg_for_w1' `var_knowledge_news_reg_cen_w1' news_perccor_cen_w1 `var_knowledge_news_reg_unc_w1' news_perccor_unc_w1 `var_knowledge_prot_reg_chi_w1' protest_pcheard_china_w1 `var_knowledge_prot_reg_for_w1' protest_pcheard_foreign_w1 `var_knowledge_prot_reg_fak_w1' `var_knowledge_meta_reg_w1' az_knowledge_meta_w1 `var_econ_guess_reg_cn_perf_w1' az_belief_econ_perf_cn_w1 `var_econ_guess_reg_cn_conf_w1' az_belief_econ_conf_cn_w1 `var_demand_change_reg_w1' az_belief_instchange_w1 `var_trust_inst_reg_govt_w1' az_belief_trust_govt_w1 `var_trust_inst_reg_foreign_w1' az_belief_trust_foreign_w1 `var_eval_govt_reg_w1' az_belief_evalgovt_w1 `var_democracy_reg_fp_w1' az_belief_democracy_fp_w1 `var_willing_fight_reg_w1' az_belief_willing_w1 `var_socialinteract_reg_w1' az_var_socialinteract_w1 `var_polparticipation_reg_pf_w1' az_var_polparticipation_w1 `var_stock_invest_reg_w1' `var_planaftergrad_reg_w1' `var_career_sector_reg_w1' `var_career_loc_reg_w1' `var_demog_reg_personal' az_demographics_personal `var_demog_reg_education' az_demographics_education `var_demog_reg_english' az_demographics_english `var_demog_reg_travel' az_demographics_travel `var_demog_reg_household' az_demographics_household `var_preference_risk' az_preference_risk `var_preference_time' az_preference_time `var_preference_altruism' az_preference_altruism `var_preference_reciprocity' az_preference_reciprocity treatment_control treatment_vpnonly treatment_nlonly treatment_vpnnl treatment_user {
		qui sum `X' if wave2 == 0 & wave3 == 0
		local N_w1=r(N)
		local mu_w1=r(mean)
		local sd_w1=r(sd)

		qui sum `X' if panelmerged_wave2 == 1 & wave2 == 0 & wave3 == 0
		local N_w2=r(N)
		local mu_w2=r(mean)
		local sd_w2=r(sd)

		qui sum `X' if panelmerged_wave3 == 1 & wave2 == 0 & wave3 == 0
		local N_w3=r(N)
		local mu_w3=r(mean)
		local sd_w3=r(sd)
		
		ttest `X', by(wave2)
		local attrition_w2_pvalue=r(p)

		ttest `X', by(wave3)
		local attrition_w3_pvalue=r(p)
		
		post `p' ("`X'") /// 
		(`mu_w1') (`sd_w1') (`N_w1') ///
		(`mu_w2') (`sd_w2') (`N_w2') ///
		(`attrition_w2_pvalue') ///
		(`mu_w3') (`sd_w3') (`N_w3') ///
		(`attrition_w3_pvalue')
		}
		
	postclose `p'
	restore
	
	preserve	
	clear
	use attrition_appendix
	outsheet using attrition_appendix.xls, replace
	clear
	erase attrition_appendix.dta
	restore	
	

	

*** Table A.5

	cap erase outregs_attrition_prediction.xls
	cap erase outregs_attrition_prediction.txt	

	preserve

	// generate interactions
	foreach a of varlist az_overall_a_w1 az_overall_b_w1 az_overall_c_w1 az_overall_d_w1 az_overall_e_w1 {
		forvalues i = 1/5 {
			gen `a'_`i' = `a' * (treatment_master == `i')
			}
		}
		
	gen attrited_wave3 = 1 - panelmerged_wave3
	
	qui reg attrited_wave3 az_overall_a_w1 treatment_vpnonly treatment_nlonly treatment_vpnnl treatment_user az_overall_a_w1_2 az_overall_a_w1_3 az_overall_a_w1_4 az_overall_a_w1_5, r
	outreg2 using outregs_attrition_prediction.xls, label dec(3)

	qui reg attrited_wave3 az_overall_b_w1 treatment_vpnonly treatment_nlonly treatment_vpnnl treatment_user az_overall_b_w1_2 az_overall_b_w1_3 az_overall_b_w1_4 az_overall_b_w1_5, r
	outreg2 using outregs_attrition_prediction.xls, label dec(3)

	qui reg attrited_wave3 az_overall_c_w1 treatment_vpnonly treatment_nlonly treatment_vpnnl treatment_user az_overall_c_w1_2 az_overall_c_w1_3 az_overall_c_w1_4 az_overall_c_w1_5, r
	outreg2 using outregs_attrition_prediction.xls, label dec(3)

	qui reg attrited_wave3 az_overall_d_w1 treatment_vpnonly treatment_nlonly treatment_vpnnl treatment_user az_overall_d_w1_2 az_overall_d_w1_3 az_overall_d_w1_4 az_overall_d_w1_5, r
	outreg2 using outregs_attrition_prediction.xls, label dec(3)

	qui reg attrited_wave3 az_overall_e_w1 treatment_vpnonly treatment_nlonly treatment_vpnnl treatment_user az_overall_e_w1_2 az_overall_e_w1_3 az_overall_e_w1_4 az_overall_e_w1_5, r
	outreg2 using outregs_attrition_prediction.xls, label dec(3)

	restore
	



*** Table A.6

	cap erase outregs_heterogeneity_takeup.xls
	cap erase outregs_heterogeneity_takeup.txt	
		
	foreach X in gender birth_year ethnicity_han birthplace_coastal residence_coastal hukou_urban religion_religious ccp_member university_elite hs_track_science department_ssh domestic_english_atleast4 foreign_english_yes travel_hktaiwan travel_foreign_yes siblings_total father_edu_hsabove work_father_govt father_ccp mother_edu_hsabove work_mother_govt mother_ccp hh_income az_preference_risk az_preference_time az_preference_altruism az_preference_reciprocity az_knowledge_news_cens_w1 az_knowledge_news_unce_w1 az_belief_media_value_w1 az_belief_media_trust_w1 {				
		
		// main specification
		qui reg vpn_adopted h_`X', r
		outreg2 h_`X' using outregs_heterogeneity_takeup.xls, label dec(3) addtext(Specification, Main) br		

		// control for AE treatment
		qui reg vpn_adopted h_`X' treatment_vpnnl, r
		outreg2 h_`X' using outregs_heterogeneity_takeup.xls, label dec(3) addtext(Specification, Control for AE) br		

		}

		

*** Table A.7

	cap erase outregs_heterogeneity_activeuser.xls
	cap erase outregs_heterogeneity_activeuser.txt	
	
	preserve
	keep if (treatment_master == 2 | treatment_master == 4)
	
	foreach X in gender birth_year ethnicity_han birthplace_coastal residence_coastal hukou_urban religion_religious ccp_member university_elite hs_track_science department_ssh domestic_english_atleast4 foreign_english_yes travel_hktaiwan travel_foreign_yes siblings_total father_edu_hsabove work_father_govt father_ccp mother_edu_hsabove work_mother_govt mother_ccp hh_income az_preference_risk az_preference_time az_preference_altruism az_preference_reciprocity az_knowledge_news_cens_w1 az_knowledge_news_unce_w1 az_belief_media_value_w1 az_belief_media_trust_w1 {				
		
		// main specification
		qui reg active_user h_`X', r
		outreg2 h_`X' using outregs_heterogeneity_activeuser.xls, label dec(3) addtext(Specification, Main) br		

		// control for AE treatment
		qui reg active_user h_`X' treatment_vpnnl, r
		outreg2 h_`X' using outregs_heterogeneity_activeuser.xls, label dec(3) addtext(Specification, Control for AE) br		

		}
		
	restore



*** Table A.9

	cap erase outregs_wave3_panelregression.xls
	cap erase outregs_wave3_panelregression.txt	

	preserve
	keep if panelmerged_wave3 == 1
	
	foreach Y in `var_info_ranking_fig_w23' `var_info_freq_fig_w123' `var_vpn_purchase_reg_w3' `var_media_valuation_fig_w123' `var_media_trust_fig_w123' `var_censor_level_fig_w123' `var_censor_justif_fig_w123' `var_censor_justif_fig_w23' `var_censor_driver_fig_w123' `var_knowledge_news_reg_cen_w3' `var_knowledge_news_reg_unc_w3' `var_knowledge_news_fig_w123' `var_knowledge_prot_fig_w123' `var_knowledge_prot_fig_w23' `var_knowledge_prot_fig_w3' `var_knowledge_meta_fig_w123' `var_demand_change_fig_w123' `var_trust_inst_fig_w123' `var_eval_govt_fig_w123' `var_democracy_fig_w123' `var_willing_fight_fig_123' `var_socialinteract_fig_w123' `var_polparticipation_fig_w123' `var_stock_invest_fig_w123' `var_planaftergrad_fig_w123' `var_career_sector_fig_w123' `var_career_loc_reg_w123' {
	
	capture{
		
		// raw regression
		qui reg `Y'_w3 treatment_vpnonly treatment_nlonly treatment_vpnnl if treatment_user != 1, r
		qui summ `Y'_w3 if treatment_user != 1
		local m_all = `r(mean)'
		local sd_all = `r(sd)'
		qui summ `Y'_w3 if treatment_control == 1
		local m_control = `r(mean)'
		local sd_control = `r(sd)'
		qui summ `Y'_w3 if treatment_user == 1
		local m_user = `r(mean)'
		local sd_user = `r(sd)'		
		local labelY: var label `Y'_w3
		local t_treatvpnXnl = _b[treatment_vpnnl]/_se[treatment_vpnnl]
		local p_treatvpnXnl = 2*ttail(e(df_r),abs(`t_treatvpnXnl'))
		outreg2 treatment_vpn treatment_newsletter treatvpnXnl using outregs_wave3_panelregression.xls, addstat("Mean DV all", `m_all', "Std.Dev.\ DV all", `sd_all', "Mean DV control", `m_control', "Std.Dev.\ DV control", `sd_control', "Mean DV user", `m_user', "Std.Dev.\ DV user", `sd_user', p-value, `p_treatvpnXnl') label dec(3) adec(3) ctitle(`Y', "`labelY'", ) addtext(Specification, Raw Diff) br		
		
		// control for imbalanced characteristics
		qui reg `Y'_w3 treatment_vpnonly treatment_nlonly treatment_vpnnl `var_demog_reg_imbalance' if treatment_user != 1, r
		qui summ `Y'_w3 if treatment_user != 1
		local m_all = `r(mean)'
		local sd_all = `r(sd)'
		qui summ `Y'_w3 if treatment_control == 1
		local m_control = `r(mean)'
		local sd_control = `r(sd)'
		qui summ `Y'_w3 if treatment_user == 1
		local m_user = `r(mean)'
		local sd_user = `r(sd)'		
		local labelY: var label `Y'_w3
		local t_treatvpnXnl = _b[treatment_vpnnl]/_se[treatment_vpnnl]
		local p_treatvpnXnl = 2*ttail(e(df_r),abs(`t_treatvpnXnl'))
		outreg2 treatment_vpn treatment_newsletter treatvpnXnl using outregs_wave3_panelregression.xls, addstat("Mean DV all", `m_all', "Std.Dev.\ DV all", `sd_all', "Mean DV control", `m_control', "Std.Dev.\ DV control", `sd_control', "Mean DV user", `m_user', "Std.Dev.\ DV user", `sd_user', p-value, `p_treatvpnXnl') label dec(3) adec(3) ctitle(`Y', "`labelY'", ) addtext(Specification, Control) br		
		
		// panel regression
		qui reg `Y'_w3 treatment_vpnonly treatment_nlonly treatment_vpnnl `Y'_w1 if treatment_user != 1, r
		qui summ `Y'_w3 if treatment_user != 1
		local m_all = `r(mean)'
		local sd_all = `r(sd)'
		qui summ `Y'_w3 if treatment_control == 1
		local m_control = `r(mean)'
		local sd_control = `r(sd)'
		qui summ `Y'_w3 if treatment_user == 1
		local m_user = `r(mean)'
		local sd_user = `r(sd)'		
		local labelY: var label `Y'_w3
		local t_treatvpnXnl = _b[treatment_vpnnl]/_se[treatment_vpnnl]
		local p_treatvpnXnl = 2*ttail(e(df_r),abs(`t_treatvpnXnl'))
		outreg2 treatment_vpn treatment_newsletter treatvpnXnl using outregs_wave3_panelregression.xls, addstat("Mean DV all", `m_all', "Std.Dev.\ DV all", `sd_all', "Mean DV control", `m_control', "Std.Dev.\ DV control", `sd_control', "Mean DV user", `m_user', "Std.Dev.\ DV user", `sd_user', p-value, `p_treatvpnXnl') label dec(3) adec(3) ctitle(`Y', "`labelY'", ) addtext(Specification, Panel) br		
			
		}
			
		}
	
	foreach Y in `var_knowledge_news_reg_cen_w3' `var_knowledge_news_reg_unc_w3' {
	
	capture{
		
		// raw regression
		qui reg `Y' treatment_vpnonly treatment_nlonly treatment_vpnnl if treatment_user != 1, r
		qui summ `Y' if treatment_user != 1
		local m_all = `r(mean)'
		local sd_all = `r(sd)'
		qui summ `Y' if treatment_control == 1
		local m_control = `r(mean)'
		local sd_control = `r(sd)'
		qui summ `Y' if treatment_user == 1
		local m_user = `r(mean)'
		local sd_user = `r(sd)'		
		local labelY: var label `Y'
		local t_treatvpnXnl = _b[treatment_vpnnl]/_se[treatment_vpnnl]
		local p_treatvpnXnl = 2*ttail(e(df_r),abs(`t_treatvpnXnl'))
		outreg2 treatment_vpn treatment_newsletter treatvpnXnl using outregs_wave3_panelregression.xls, addstat("Mean DV all", `m_all', "Std.Dev.\ DV all", `sd_all', "Mean DV control", `m_control', "Std.Dev.\ DV control", `sd_control', "Mean DV user", `m_user', "Std.Dev.\ DV user", `sd_user', p-value, `p_treatvpnXnl') label dec(3) adec(3) ctitle(`Y', "`labelY'", ) addtext(Specification, Raw Diff) br		
		
		// control for imbalanced characteristics
		qui reg `Y' treatment_vpnonly treatment_nlonly treatment_vpnnl `var_demog_reg_imbalance' if treatment_user != 1, r
		qui summ `Y' if treatment_user != 1
		local m_all = `r(mean)'
		local sd_all = `r(sd)'
		qui summ `Y' if treatment_control == 1
		local m_control = `r(mean)'
		local sd_control = `r(sd)'
		qui summ `Y' if treatment_user == 1
		local m_user = `r(mean)'
		local sd_user = `r(sd)'		
		local labelY: var label `Y'
		local t_treatvpnXnl = _b[treatment_vpnnl]/_se[treatment_vpnnl]
		local p_treatvpnXnl = 2*ttail(e(df_r),abs(`t_treatvpnXnl'))
		outreg2 treatment_vpn treatment_newsletter treatvpnXnl using outregs_wave3_panelregression.xls, addstat("Mean DV all", `m_all', "Std.Dev.\ DV all", `sd_all', "Mean DV control", `m_control', "Std.Dev.\ DV control", `sd_control', "Mean DV user", `m_user', "Std.Dev.\ DV user", `sd_user', p-value, `p_treatvpnXnl') label dec(3) adec(3) ctitle(`Y', "`labelY'", ) addtext(Specification, Control) br		
				
		}
			
		}

	foreach Y in `var_econ_guess_fig_w123' `var_econ_guess_fig_w23' {
			
	capture{
	
		// raw regression
		qui reg `Y'_w3 treatment_vpnonly treatment_nlonly treatment_vpnnl if treatment_user != 1, r
		qui summ `Y'_w3 if treatment_user != 1
		local m_all = `r(mean)'
		local sd_all = `r(sd)'
		qui summ `Y'_w3 if treatment_control == 1
		local m_control = `r(mean)'
		local sd_control = `r(sd)'
		qui summ `Y'_w3 if treatment_user == 1
		local m_user = `r(mean)'
		local sd_user = `r(sd)'		
		local labelY: var label `Y'_w3
		local t_treatvpnXnl = _b[treatment_vpnnl]/_se[treatment_vpnnl]
		local p_treatvpnXnl = 2*ttail(e(df_r),abs(`t_treatvpnXnl'))
		outreg2 treatment_vpn treatment_newsletter treatvpnXnl using outregs_wave3_panelregression.xls, addstat("Mean DV all", `m_all', "Std.Dev.\ DV all", `sd_all', "Mean DV control", `m_control', "Std.Dev.\ DV control", `sd_control', "Mean DV user", `m_user', "Std.Dev.\ DV user", `sd_user', p-value, `p_treatvpnXnl') label dec(8) adec(8) ctitle(`Y', "`labelY'", ) addtext(Specification, Raw Diff) br		
		
		// control for imbalanced characteristics
		qui reg `Y'_w3 treatment_vpnonly treatment_nlonly treatment_vpnnl `var_demog_reg_imbalance' if treatment_user != 1, r
		qui summ `Y'_w3 if treatment_user != 1
		local m_all = `r(mean)'
		local sd_all = `r(sd)'
		qui summ `Y'_w3 if treatment_control == 1
		local m_control = `r(mean)'
		local sd_control = `r(sd)'
		qui summ `Y'_w3 if treatment_user == 1
		local m_user = `r(mean)'
		local sd_user = `r(sd)'		
		local labelY: var label `Y'_w3
		local t_treatvpnXnl = _b[treatment_vpnnl]/_se[treatment_vpnnl]
		local p_treatvpnXnl = 2*ttail(e(df_r),abs(`t_treatvpnXnl'))
		outreg2 treatment_vpn treatment_newsletter treatvpnXnl using outregs_wave3_panelregression.xls, addstat("Mean DV all", `m_all', "Std.Dev.\ DV all", `sd_all', "Mean DV control", `m_control', "Std.Dev.\ DV control", `sd_control', "Mean DV user", `m_user', "Std.Dev.\ DV user", `sd_user', p-value, `p_treatvpnXnl') label dec(8) adec(8) ctitle(`Y', "`labelY'", ) addtext(Specification, Control) br		

		// panel regression
		qui reg `Y'_w3 treatment_vpnonly treatment_nlonly treatment_vpnnl `Y'_w1 if treatment_user != 1, r
		qui summ `Y'_w3 if treatment_user != 1
		local m_all = `r(mean)'
		local sd_all = `r(sd)'
		qui summ `Y'_w3 if treatment_control == 1
		local m_control = `r(mean)'
		local sd_control = `r(sd)'
		qui summ `Y'_w3 if treatment_user == 1
		local m_user = `r(mean)'
		local sd_user = `r(sd)'		
		local labelY: var label `Y'_w3
		local t_treatvpnXnl = _b[treatment_vpnnl]/_se[treatment_vpnnl]
		local p_treatvpnXnl = 2*ttail(e(df_r),abs(`t_treatvpnXnl'))
		outreg2 treatment_vpn treatment_newsletter treatvpnXnl using outregs_wave3_panelregression.xls, addstat("Mean DV all", `m_all', "Std.Dev.\ DV all", `sd_all', "Mean DV control", `m_control', "Std.Dev.\ DV control", `sd_control', "Mean DV user", `m_user', "Std.Dev.\ DV user", `sd_user', p-value, `p_treatvpnXnl') label dec(8) adec(8) ctitle(`Y', "`labelY'", ) addtext(Specification, Panel) br		

		}
		
		}
		
	restore
	
	

	
*** Table A.11

	cap erase outregs_maintreatmenteffects_robustness.xls
	cap erase outregs_maintreatmenteffects_robustness.txt	
	
	preserve
	keep if panelmerged_wave3 == 1
		
	foreach Y in a b c d e {
		
		// baseline regression (raw)
		qui reg az_overall_`Y'_w3 treatment_vpnonly treatment_nlonly treatment_vpnnl if treatment_user != 1, r
		qui summ az_overall_`Y'_w3 if treatment_user != 1
		local m_all = `r(mean)'
		local sd_all = `r(sd)'
		qui summ az_overall_`Y'_w3 if treatment_control == 1
		local m_control = `r(mean)'
		local sd_control = `r(sd)'
		qui summ az_overall_`Y'_w3 if treatment_user == 1
		local m_user = `r(mean)'
		local sd_user = `r(sd)'		
		local labelY: var label az_overall_`Y'_w3
		local t_treatvpnXnl = _b[treatment_vpnnl]/_se[treatment_vpnnl]
		local p_treatvpnXnl = 2*ttail(e(df_r),abs(`t_treatvpnXnl'))
		outreg2 treatment_vpnonly treatment_nlonly treatment_vpnnl using outregs_maintreatmenteffects_robustness.xls, addstat("Mean DV all", `m_all', "Std.Dev.\ DV all", `sd_all', "Mean DV control", `m_control', "Std.Dev.\ DV control", `sd_control', "Mean DV user", `m_user', "Std.Dev.\ DV user", `sd_user', p-value, `p_treatvpnXnl') label dec(3) adec(3) ctitle(`Y', "`labelY'", ) addtext(Specification, Raw diff) br		

		// control for imbalanced characteristics
		qui reg az_overall_`Y'_w3 treatment_vpnonly treatment_nlonly treatment_vpnnl `var_demog_reg_imbalance' if treatment_user != 1, r
		qui summ az_overall_`Y'_w3 if treatment_user != 1
		local m_all = `r(mean)'
		local sd_all = `r(sd)'
		qui summ az_overall_`Y'_w3 if treatment_control == 1
		local m_control = `r(mean)'
		local sd_control = `r(sd)'
		qui summ az_overall_`Y'_w3 if treatment_user == 1
		local m_user = `r(mean)'
		local sd_user = `r(sd)'		
		local labelY: var label az_overall_`Y'_w3
		local t_treatvpnXnl = _b[treatment_vpnnl]/_se[treatment_vpnnl]
		local p_treatvpnXnl = 2*ttail(e(df_r),abs(`t_treatvpnXnl'))
		outreg2 treatment_vpnonly treatment_nlonly treatment_vpnnl using outregs_maintreatmenteffects_robustness.xls, addstat("Mean DV all", `m_all', "Std.Dev.\ DV all", `sd_all', "Mean DV control", `m_control', "Std.Dev.\ DV control", `sd_control', "Mean DV user", `m_user', "Std.Dev.\ DV user", `sd_user', p-value, `p_treatvpnXnl') label dec(3) adec(3) ctitle(`Y', "`labelY'", ) addtext(Specification, Control for imbalance) br		
		
		// control for baseline level
		qui reg az_overall_`Y'_w3 treatment_vpnonly treatment_nlonly treatment_vpnnl az_overall_`Y'_w1 if treatment_user != 1, r
		qui summ az_overall_`Y'_w3 if treatment_user != 1
		local m_all = `r(mean)'
		local sd_all = `r(sd)'
		qui summ az_overall_`Y'_w3 if treatment_control == 1
		local m_control = `r(mean)'
		local sd_control = `r(sd)'
		qui summ az_overall_`Y'_w3 if treatment_user == 1
		local m_user = `r(mean)'
		local sd_user = `r(sd)'		
		local labelY: var label az_overall_`Y'_w3
		local t_treatvpnXnl = _b[treatment_vpnnl]/_se[treatment_vpnnl]
		local p_treatvpnXnl = 2*ttail(e(df_r),abs(`t_treatvpnXnl'))
		outreg2 treatment_vpnonly treatment_nlonly treatment_vpnnl using outregs_maintreatmenteffects_robustness.xls, addstat("Mean DV all", `m_all', "Std.Dev.\ DV all", `sd_all', "Mean DV control", `m_control', "Std.Dev.\ DV control", `sd_control', "Mean DV user", `m_user', "Std.Dev.\ DV user", `sd_user', p-value, `p_treatvpnXnl') label dec(3) adec(3) ctitle(`Y', "`labelY'", ) addtext(Specification, Control for baseline) br				
		
		// drop those answered 10 on trust of central govt
		qui reg az_overall_`Y'_w3 treatment_vpnonly treatment_nlonly treatment_vpnnl if treatment_user != 1 & trust_central_govt_w3 != 10, r
		qui summ az_overall_`Y'_w3 if treatment_user != 1 & trust_central_govt_w3 != 10
		local m_all = `r(mean)'
		local sd_all = `r(sd)'
		qui summ az_overall_`Y'_w3 if treatment_control == 1 & trust_central_govt_w3 != 10
		local m_control = `r(mean)'
		local sd_control = `r(sd)'
		qui summ az_overall_`Y'_w3 if treatment_user == 1 & trust_central_govt_w3 != 10
		local m_user = `r(mean)'
		local sd_user = `r(sd)'		
		local labelY: var label az_overall_`Y'_w3
		local t_treatvpnXnl = _b[treatment_vpnnl]/_se[treatment_vpnnl]
		local p_treatvpnXnl = 2*ttail(e(df_r),abs(`t_treatvpnXnl'))
		outreg2 treatment_vpnonly treatment_nlonly treatment_vpnnl using outregs_maintreatmenteffects_robustness.xls, addstat("Mean DV all", `m_all', "Std.Dev.\ DV all", `sd_all', "Mean DV control", `m_control', "Std.Dev.\ DV control", `sd_control', "Mean DV user", `m_user', "Std.Dev.\ DV user", `sd_user', p-value, `p_treatvpnXnl') label dec(3) adec(3) ctitle(`Y', "`labelY'", ) addtext(Specification, Drop pol correct) br		
		
		}
		
	restore
	


*** Table A.12

* 	version A: ranked among all participants

	preserve
	
 	// keep wave 3 participants only
	keep if panelmerged_wave3 == 1
	
	tempname p
	postfile `p' str50 variable float p1 p3 change ///
		using outregs_quantilemovement_wave3_all, replace
	
	// quantile movement for panel variables
	foreach Y in info_foreign_website info_freq_website_for az_belief_media_value az_belief_media_trust bias_domestic bias_foreign_r az_belief_media_justif bias_dom_govt_policy_t1 bias_for_govt_policy_t1 news_perccor_cen news_perccor_unc protest_pcheard_china protest_pcheard_foreign protest_2011_tmrw_parade az_knowledge_meta az_belief_econ_perf_cn_r az_belief_econ_conf_cn az_belief_instchange az_belief_trust_govt_r az_belief_trust_foreign az_belief_evalgovt_r importance_live_in_demo_r az_belief_willing frequency_talk_politic frequency_persuade_friends participate_social_protest participate_plan_vote participate_complain_school stock_participation plan_grad_foreignmaster cp_t3_for_firm cloc_for az_overall {
		
		// add random noise to each variable
		gen noise1 = (runiform() /2)
		gen noise3 = (runiform() /2)
		gen `Y'_w1_n = `Y'_w1 + noise1
		gen `Y'_w3_n = `Y'_w3 + noise3
		drop noise1 noise3
		
		xtile `Y'_p1 = `Y'_w1_n, nq(100)
		xtile `Y'_p3 = `Y'_w3_n, nq(100)
		su `Y'_p1 if treatment_vpnnl == 1, d
		local p1 = r(p50)
		su `Y'_p3 if treatment_vpnnl == 1, d
		local p3 = r(p50)
		local change = `p3'-`p1'
		
		post `p' ("`Y'") (`p1') (`p3') (`change')
		
		}
		
	// quantile movement for non-panel variables
	foreach Y in vpn_purchase_wmt_record vpn_purchase_yes az_belief_econ_perf_us_w3 az_belief_econ_conf_us_w3 {
		
		// add random noise to each variable
		gen noise3 = (runiform() /2)
		gen `Y'_n = `Y' + noise3
		drop noise3
		
		xtile `Y'_p = `Y'_n, nq(100)
		su `Y'_p if treatment_control == 1 | treatment_vpnonly == 1 | treatment_nlonly == 1, d
		local p1 = r(p50)
		su `Y'_p if treatment_vpnnl == 1, d
		local p3 = r(p50)
		local change = `p3'-`p1'
		
		post `p' ("`Y'") (`p1') (`p3') (`change')
		
		}
		
	postclose `p'
	restore
	
	preserve	
	clear
	use outregs_quantilemovement_wave3_all
	outsheet using outregs_quantilemovement_wave3_all.xls, replace
	clear
	erase outregs_quantilemovement_wave3_all.dta
	restore	
	
	
* 	version B: ranked among non-existing users

	preserve
	
 	// keep wave 3 participants only
	keep if panelmerged_wave3 == 1 & treatment_user != 1
	
	tempname p
	postfile `p' str50 variable float p1 p3 change ///
		using outregs_quantilemovement_wave3_nonexistinguser, replace
	
	// quantile movement for panel variables
	foreach Y in info_foreign_website info_freq_website_for az_belief_media_value az_belief_media_trust bias_domestic bias_foreign_r az_belief_media_justif bias_dom_govt_policy_t1 bias_for_govt_policy_t1 news_perccor_cen news_perccor_unc protest_pcheard_china protest_pcheard_foreign protest_2011_tmrw_parade az_knowledge_meta az_belief_econ_perf_cn_r az_belief_econ_conf_cn az_belief_instchange az_belief_trust_govt_r az_belief_trust_foreign az_belief_evalgovt_r importance_live_in_demo_r az_belief_willing frequency_talk_politic frequency_persuade_friends participate_social_protest participate_plan_vote participate_complain_school stock_participation plan_grad_foreignmaster cp_t3_for_firm cloc_for az_overall {
		
		// add random noise to each variable
		gen noise1 = (runiform() /500)
		gen noise3 = (runiform() /500)
		gen `Y'_w1_n = `Y'_w1 + noise1
		gen `Y'_w3_n = `Y'_w3 + noise3
		drop noise1 noise3
		
		xtile `Y'_p1 = `Y'_w1_n, nq(100)
		xtile `Y'_p3 = `Y'_w3_n, nq(100)
		su `Y'_p1 if treatment_vpnnl == 1, d
		local p1 = r(p50)
		su `Y'_p3 if treatment_vpnnl == 1, d
		local p3 = r(p50)
		local change = `p3'-`p1'
		
		post `p' ("`Y'") (`p1') (`p3') (`change')
		
		}
		
	// quantile movement for non-panel variables
	foreach Y in vpn_purchase_wmt_record vpn_purchase_yes az_belief_econ_perf_us_w3 az_belief_econ_conf_us_w3 {
		
		// add random noise to each variable
		gen noise3 = (runiform() /500)
		gen `Y'_n = `Y' + noise3
		drop noise3
		
		xtile `Y'_p = `Y'_n, nq(100)
		su `Y'_p if treatment_control == 1 | treatment_vpnonly == 1 | treatment_nlonly == 1, d
		local p1 = r(p50)
		su `Y'_p if treatment_vpnnl == 1, d
		local p3 = r(p50)
		local change = `p3'-`p1'
		
		post `p' ("`Y'") (`p1') (`p3') (`change')
		
		}
		
	postclose `p'
	restore
	
	preserve	
	clear
	use outregs_quantilemovement_wave3_nonexistinguser
	outsheet using outregs_quantilemovement_wave3_nonexistinguser.xls, replace
	clear
	erase outregs_quantilemovement_wave3_nonexistinguser.dta
	restore	
	
	

			
*** Table A.13
	
	cap erase outregs_persuasionrates_wave3.xls
	cap erase outregs_persuasionrates_wave3.txt	
	
	preserve
	
	// keep wave 3 non-existing users only
	keep if panelmerged_wave3 == 1
	keep if treatment_user != 1
	
	// transform into dummy indicators: uncensored belief - above median
	foreach Y in info_foreign_website_w3 info_freq_website_for_w3 az_belief_media_value_w3 az_belief_media_trust_w3 bias_domestic_w3 az_belief_media_justif_w3 news_perccor_cen_w3 news_perccor_unc_w3 protest_pcheard_china_w3 protest_pcheard_foreign_w3 az_knowledge_meta_w3 az_belief_econ_conf_cn_w3 az_belief_econ_perf_us_w3 az_belief_econ_conf_us_w3 az_belief_instchange_w3 az_belief_trust_foreign_w3 az_belief_willing_w3 frequency_talk_politic_w3 frequency_persuade_friends_w3 {
		egen `Y'_m = median(`Y')
		gen `Y'_p = (`Y' > `Y'_m)
		}

	foreach Y in info_foreign_website_w1 info_freq_website_for_w1 az_belief_media_value_w1 az_belief_media_trust_w1 bias_domestic_w1 az_belief_media_justif_w1 news_perccor_cen_w1 news_perccor_unc_w1 protest_pcheard_china_w1 protest_pcheard_foreign_w1 az_knowledge_meta_w1 az_belief_econ_conf_cn_w1 az_belief_instchange_w1 az_belief_trust_foreign_w1 az_belief_willing_w1 frequency_talk_politic_w1 frequency_persuade_friends_w1 {
		egen `Y'_m = median(`Y')
		gen `Y'_p = (`Y' > `Y'_m)
		}

	// transform into dummy indicators: uncensored belief - below median
	foreach Y in bias_foreign_w3 importance_live_in_demo_w3 az_belief_econ_perf_cn_w3 az_belief_trust_govt_w3 az_belief_evalgovt_w3 {
		egen `Y'_m = median(`Y')
		gen `Y'_p = (`Y' < `Y'_m)
		}

	foreach Y in bias_foreign_w1 importance_live_in_demo_w1 az_belief_econ_perf_cn_w1 az_belief_trust_govt_w1 az_belief_evalgovt_w1 {
		egen `Y'_m = median(`Y')
		gen `Y'_p = (`Y' < `Y'_m)
		}

	// already dummy, keep as it is
	rename participate_complain_school_w1 par_complain_school_w1
	rename participate_complain_school_w3 par_complain_school_w3

	foreach Y in vpn_purchase_wmt_record vpn_purchase_yes bias_dom_govt_policy_t1_w3 protest_2011_tmrw_parade_w3 participate_social_protest_w3 participate_plan_vote_w3 par_complain_school_w3 plan_grad_foreignmaster_w3 cp_t3_for_firm_w3 cloc_for_w3 {
		gen `Y'_p = `Y'
		}

	foreach Y in bias_dom_govt_policy_t1_w1 protest_2011_tmrw_parade_w1 participate_social_protest_w1 participate_plan_vote_w1 par_complain_school_w1 plan_grad_foreignmaster_w1 cp_t3_for_firm_w1 cloc_for_w1 {
		gen `Y'_p = `Y'
		}

	// already dummy, flip sign
	foreach Y in bias_for_govt_policy_t1_w3 stock_participation_w3 {
		gen `Y'_p = 1- `Y'
		}

	foreach Y in bias_for_govt_policy_t1_w1 stock_participation_w1 {
		gen `Y'_p = 1- `Y'
		}

	// persuasion rates for panel variables
	foreach Y in info_foreign_website info_freq_website_for az_belief_media_value az_belief_media_trust bias_domestic az_belief_media_justif news_perccor_cen news_perccor_unc protest_pcheard_china protest_pcheard_foreign az_knowledge_meta az_belief_econ_conf_cn az_belief_instchange az_belief_trust_foreign importance_live_in_demo az_belief_willing frequency_talk_politic frequency_persuade_friends bias_foreign az_belief_econ_perf_cn az_belief_trust_govt az_belief_evalgovt bias_dom_govt_policy_t1 protest_2011_tmrw_parade participate_social_protest participate_plan_vote par_complain_school plan_grad_foreignmaster cp_t3_for_firm cloc_for bias_for_govt_policy_t1 stock_participation {
				
		qui reg `Y'_w3_p treatment_vpnonly treatment_nlonly treatment_vpnnl, r
		qui summ `Y'_w3_p
		local m_all = `r(mean)'
		local sd_all = `r(sd)'
		qui summ `Y'_w3_p if treatment_control == 1
		local m_control = `r(mean)'
		local sd_control = `r(sd)'
		local labelY: var label `Y'_w3_p
		qui summ `Y'_w1_p if treatment_vpnnl == 1
		local t_treatvpnXnl = _b[treatment_vpnnl]/_se[treatment_vpnnl]
		local p_treatvpnXnl = 2*ttail(e(df_r),abs(`t_treatvpnXnl'))
		local persuasion = _b[treatment_vpnnl]/(0.646*(1-`r(mean)'))
		outreg2 treatment_vpnnl using outregs_persuasionrates_wave3.xls, addstat("Mean DV all", `m_all', "SD all", `sd_all', "Mean DV control", `m_control', p-value, `p_treatvpnXnl', persuasion, `persuasion') label dec(3) adec(3) ctitle(`Y'_w3_p, "`labelY'", ) addtext(Specification, Raw Diff) br		
		
		}

	// persuasion rates for non-panel variables
	foreach Y in vpn_purchase_wmt_record vpn_purchase_yes az_belief_econ_perf_us_w3 az_belief_econ_conf_us_w3 {
				
		qui reg `Y'_p treatment_vpnonly treatment_nlonly treatment_vpnnl, r
		qui summ `Y'_p
		local m_all = `r(mean)'
		local sd_all = `r(sd)'
		qui summ `Y'_p if treatment_control == 1
		local m_control = `r(mean)'
		local sd_control = `r(sd)'
		local labelY: var label `Y'
		local t_treatvpnXnl = _b[treatment_vpnnl]/_se[treatment_vpnnl]
		local p_treatvpnXnl = 2*ttail(e(df_r),abs(`t_treatvpnXnl'))
		local persuasion = _b[treatment_vpnnl]/(0.646*(1-_b[_cons]))
		outreg2 treatment_vpnnl using outregs_persuasionrates_wave3.xls, addstat("Mean DV all", `m_all', "SD all", `sd_all', "Mean DV control", `m_control', p-value, `p_treatvpnXnl', persuasion, `persuasion') label dec(3) adec(3) ctitle(`Y', "`labelY'", ) addtext(Specification, Raw Diff) br		
		
		}
	
	restore
	


*** Table A.14

	cap erase outregs_heterogeneity_all.xls
	cap erase outregs_heterogeneity_all.txt	
	
	preserve
	
	keep if panelmerged_wave3 == 1
	keep if treatment_user != 1
	
	// generate interaction terms
	foreach h in gender birth_year residence_coastal hukou_urban university_elite hs_track_science department_ssh domestic_english_atleast4 foreign_english_yes travel_hktaiwan travel_foreign_yes father_edu_hsabove work_father_govt father_ccp mother_edu_hsabove work_mother_govt mother_ccp hh_income az_preference_risk az_preference_time az_preference_altruism az_preference_reciprocity az_overall_a_w1 az_overall_b_w1 az_overall_c_w1 az_overall_d_w1 az_overall_e_w1 {
		gen aeX`h' = treatment_vpnnl * h_`h'
		}
		
	foreach Y in az_overall_a_w3 az_overall_b_w3 az_overall_c_w3 az_overall_d_w3 az_overall_e_w3 {
		reg `Y' treatment_vpnnl h_gender h_birth_year h_residence_coastal h_hukou_urban h_university_elite h_hs_track_science h_department_ssh h_domestic_english_atleast4 h_foreign_english_yes h_travel_hktaiwan h_travel_foreign_yes h_father_edu_hsabove h_work_father_govt h_father_ccp h_mother_edu_hsabove h_work_mother_govt h_mother_ccp h_hh_income h_az_preference_risk h_az_preference_time h_az_preference_altruism h_az_preference_reciprocity h_az_overall_a_w1 h_az_overall_b_w1 h_az_overall_c_w1 h_az_overall_d_w1 h_az_overall_e_w1 aeXgender aeXbirth_year aeXresidence_coastal aeXhukou_urban aeXuniversity_elite aeXhs_track_science aeXdepartment_ssh aeXdomestic_english_atleast4 aeXforeign_english_yes aeXtravel_hktaiwan aeXtravel_foreign_yes aeXfather_edu_hsabove aeXwork_father_govt aeXfather_ccp aeXmother_edu_hsabove aeXwork_mother_govt aeXmother_ccp aeXhh_income aeXaz_preference_risk aeXaz_preference_time aeXaz_preference_altruism aeXaz_preference_reciprocity aeXaz_overall_a_w1 aeXaz_overall_b_w1 aeXaz_overall_c_w1 aeXaz_overall_d_w1 aeXaz_overall_e_w1
		outreg2 using outregs_heterogeneity_all.xls, dec(3)
	
		}
	
	restore
	
	



*** Table A.15

	cap erase outregs_wave3sociallearning_full.xls
	cap erase outregs_wave3sociallearning_full.txt	
	
	preserve
	
	foreach Y in `var_knowledge_news_reg_cen_w2' {
		qui reg `Y' soclearning_ownaccess soclearning_rm_new_w2 soclearning_ownXnew_w2 if vpn_roommate_existing == 0 & soclearning_rm_new_w2 < 2
		
		// estimating p(direct learning)
		qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 0)
		local mu_100=r(mean)
		qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 0)
		local mu_000=r(mean)
		local p_all = `mu_100'-`mu_000'
		
		// estimating q based on 1 roommate
		qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 1)
		local mu_101=r(mean)
		qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 1)
		local mu_001=r(mean)
		local q_new_1_0 = (`mu_001'-`mu_000')/`mu_100'
		local q_new_1_1 = (`mu_101'-`mu_100')/`mu_100'

		// predicted and actual based on 2 roommates
		local pred_mu_102 = `mu_100' + (1-(1-`mu_100'*`q_new_1_1')^2)
		local pred_mu_002 = `mu_000' + (1-(1-`mu_100'*`q_new_1_0')^2)
		qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 2)
		local mu_002=r(mean)
		qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 2)
		local mu_102=r(mean)
		
		local labelY: var label `Y'
		outreg2 soclearning_ownaccess soclearning_rm_new_w2 soclearning_ownXnew_w2 using outregs_wave3sociallearning_full.xls, addstat("p_all", `p_all', "q_new_1_0", `q_new_1_0', "q_new_1_1", `q_new_1_1', "pred_mu_002", `pred_mu_002', "mu_002", `mu_002', "pred_mu_102", `pred_mu_102', "mu_102", `mu_102') label ctitle(`Y', "`labelY'", ) addtext(Specification, Raw Diff) br		
		}

		
	foreach Y in `var_knowledge_news_reg_cen_w3' news_perccor_cen_all {
		qui reg `Y' soclearning_ownaccess soclearning_rm_new_w3 soclearning_ownXnew_w3 if vpn_roommate_existing == 0 & soclearning_rm_new_w3 < 2
		
		// estimating p(direct learning)
		qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 0)
		local mu_100=r(mean)
		qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 0)
		local mu_000=r(mean)
		local p_all = `mu_100'-`mu_000'
		
		// estimating q based on 1 roommate
		qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 1)
		local mu_101=r(mean)
		local var_
		qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 1)
		local mu_001=r(mean)
		local q_new_1_0 = (`mu_001'-`mu_000')/`mu_100'
		local q_new_1_1 = (`mu_101'-`mu_100')/`mu_100'

		// predicted and actual based on 2 roommates
		local pred_mu_102 = `mu_100' + (1-(1-`mu_100'*`q_new_1_1')^2)
		local pred_mu_002 = `mu_000' + (1-(1-`mu_100'*`q_new_1_0')^2)
		qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 2)
		local mu_002=r(mean)
		qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 2)
		local mu_102=r(mean)
		
		local labelY: var label `Y'
		outreg2 soclearning_ownaccess soclearning_rm_new_w3 soclearning_ownXnew_w3 using outregs_wave3sociallearning_full.xls, addstat("p_all", `p_all', "q_new_1_0", `q_new_1_0', "q_new_1_1", `q_new_1_1', "pred_mu_002", `pred_mu_002', "mu_002", `mu_002', "pred_mu_102", `pred_mu_102', "mu_102", `mu_102') label ctitle(`Y', "`labelY'", ) addtext(Specification, Raw Diff) br		
		}
		
	restore
	
	// exporting standard error of estimates
	
	preserve
	gen sc_reg_100 = (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 0)
	gen sc_reg_000 = (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 0)
	gen sc_reg_101 = (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 1)
	gen sc_reg_001 = (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 1)

	foreach Y in `var_knowledge_news_reg_cen_w2' {
		qui reg `Y' sc_reg_100 sc_reg_000 sc_reg_101 sc_reg_001 if (vpn_roommate_existing == 0 & vpn_roommate_new_w2 < 2), nocons
		nlcom 	(q_new_1_0: (_b[sc_reg_001] - _b[sc_reg_000])/_b[sc_reg_100]) ///
				(q_new_1_1: (_b[sc_reg_101] - _b[sc_reg_100])/_b[sc_reg_100]) ///
				(pred_mu_002: _b[sc_reg_000] + (1-(1-_b[sc_reg_100]*((_b[sc_reg_001] - _b[sc_reg_000])/_b[sc_reg_100]))^2)) ///
				(pred_mu_102: _b[sc_reg_100] + (1-(1-_b[sc_reg_100]*((_b[sc_reg_101] - _b[sc_reg_100])/_b[sc_reg_100]))^2))
		}
		
	restore
	
	preserve
	gen sc_reg_100 = (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 0)
	gen sc_reg_000 = (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 0)
	gen sc_reg_101 = (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 1)
	gen sc_reg_001 = (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 1)

	foreach Y in `var_knowledge_news_reg_cen_w3' news_perccor_cen_all {
		qui reg `Y' sc_reg_100 sc_reg_000 sc_reg_101 sc_reg_001 if (vpn_roommate_existing == 0 & vpn_roommate_new_w3 < 2), nocons
		nlcom 	(q_new_1_0: (_b[sc_reg_001] - _b[sc_reg_000])/_b[sc_reg_100]) ///
				(q_new_1_1: (_b[sc_reg_101] - _b[sc_reg_100])/_b[sc_reg_100]) ///
				(pred_mu_002: _b[sc_reg_000] + (1-(1-_b[sc_reg_100]*((_b[sc_reg_001] - _b[sc_reg_000])/_b[sc_reg_100]))^2)) ///
				(pred_mu_102: _b[sc_reg_100] + (1-(1-_b[sc_reg_100]*((_b[sc_reg_101] - _b[sc_reg_100])/_b[sc_reg_100]))^2))
		}
		
	restore

	
	/*
	// calculate bootstrap standard errors for out-of-sample predictions
		
	preserve
	clear
	gen i = .
	gen y = ""
	gen pred_mu_102 = .
	gen pred_mu_002 = .
	save "`directory_path'/Outregs/sociallearning_pred_bootstrap_full.dta", replace
	restore
	
	keep treatment_master vpn_roommate_existing vpn_roommate_new_w2 vpn_roommate_new_w3 `var_knowledge_news_reg_cen_w2' `var_knowledge_news_reg_cen_w3' news_perccor_cen_all
	keep if news_perccor_cen_all != .
	
	forvalues i = 1/1000 {
		
		foreach Y in `var_knowledge_news_reg_cen_w2' {
			
			preserve
			
			bsample	
	
			// estimating p(direct learning)
			qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 0)
			local mu_100=r(mean)
			qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 0)
			local mu_000=r(mean)
			local p_all = `mu_100'-`mu_000'
			
			// estimating q based on 1 roommate
			qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 1)
			local mu_101=r(mean)
			local var_
			qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 1)
			local mu_001=r(mean)
			local q_new_1_0 = (`mu_001'-`mu_000')/`mu_100'
			local q_new_1_1 = (`mu_101'-`mu_100')/`mu_100'

			// predicted and actual based on 2 roommates
			local pred_mu_102 = `mu_100' + (1-(1-`mu_100'*`q_new_1_1')^2)
			local pred_mu_002 = `mu_000' + (1-(1-`mu_100'*`q_new_1_0')^2)
			qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 2)
			local mu_002=r(mean)
			qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w2 == 2)
			local mu_102=r(mean)
			
			// export estimates
			gen i = `i'
			gen y = "`Y'"
			gen pred_mu_102 = `pred_mu_102'
			gen pred_mu_002 = `pred_mu_002'
			keep i y pred_mu_102 pred_mu_002
			duplicates drop
			append using "`directory_path'/Outregs/sociallearning_pred_bootstrap_full.dta"
			save "`directory_path'/Outregs/sociallearning_pred_bootstrap_full.dta", replace
		
			restore
			
			}
			
		
		foreach Y in `var_knowledge_news_reg_cen_w3' news_perccor_cen_all {
			
			preserve
			
			bsample	
	
			// estimating p(direct learning)
			qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 0)
			local mu_100=r(mean)
			qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 0)
			local mu_000=r(mean)
			local p_all = `mu_100'-`mu_000'
			
			// estimating q based on 1 roommate
			qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 1)
			local mu_101=r(mean)
			local var_
			qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 1)
			local mu_001=r(mean)
			local q_new_1_0 = (`mu_001'-`mu_000')/`mu_100'
			local q_new_1_1 = (`mu_101'-`mu_100')/`mu_100'

			// predicted and actual based on 2 roommates
			local pred_mu_102 = `mu_100' + (1-(1-`mu_100'*`q_new_1_1')^2)
			local pred_mu_002 = `mu_000' + (1-(1-`mu_100'*`q_new_1_0')^2)
			qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 2)
			local mu_002=r(mean)
			qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 2)
			local mu_102=r(mean)
			
			// export estimates
			gen i = `i'
			gen y = "`Y'"
			gen pred_mu_102 = `pred_mu_102'
			gen pred_mu_002 = `pred_mu_002'
			keep i y pred_mu_102 pred_mu_002
			duplicates drop
			append using "`directory_path'/Outregs/sociallearning_pred_bootstrap_full.dta"
			save "`directory_path'/Outregs/sociallearning_pred_bootstrap_full.dta", replace
		
			restore
			
			}

		}
		
	use "`directory_path'/Outregs/sociallearning_pred_bootstrap_full.dta", clear
	
	gen pred_mu_102_bse = .
	gen pred_mu_002_bse = .
	
	foreach Y in `var_knowledge_news_reg_cen_w2' `var_knowledge_news_reg_cen_w3' news_perccor_cen_all {
		
		sum pred_mu_102 if y == "`Y'", d
		replace pred_mu_102_bse = `r(sd)' if y == "`Y'"
		sum pred_mu_002 if y == "`Y'", d
		replace pred_mu_002_bse = `r(sd)' if y == "`Y'"
		
		}
		
	keep y pred_mu_102_bse pred_mu_002_bse
	duplicates drop
	export excel using "sociallearning_pred_bootstrap_full.xlsx", firstrow(variables) replace

	*/
	
	

	
*** Table A.16

	cap erase outregs_wave3sociallearning_additional.xls
	cap erase outregs_wave3sociallearning_additional.txt	
	
	preserve
	keep if panelmerged_wave3 == 1
	
	foreach Y in az_overall_a_w3 az_overall_b_w3 az_overall_c_w3 az_overall_d_w3 az_overall_e_w3 {
		qui reg `Y' soclearning_ownaccess soclearning_rm_new_w3 soclearning_ownXnew_w3 if vpn_roommate_existing == 0 & soclearning_rm_new_w3 < 2
		
		// estimating p(direct learning)
		qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 0)
		local mu_100=r(mean)
		qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 0)
		local mu_000=r(mean)
		local p_all = `mu_100'-`mu_000'
		
		// estimating q based on 1 roommate
		qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 1)
		local mu_101=r(mean)
		local var_
		qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 1)
		local mu_001=r(mean)
		local q_new_1_0 = (`mu_001'-`mu_000')/`mu_100'
		local q_new_1_1 = (`mu_101'-`mu_100')/`mu_100'

		// predicted and actual basedon 2 roommates
		local pred_mu_102 = `mu_100' + (1-(1-`mu_100'*`q_new_1_1')^2)
		local pred_mu_002 = `mu_000' + (1-(1-`mu_100'*`q_new_1_0')^2)
		qui sum `Y' if (treatment_master < 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 2)
		local mu_002=r(mean)
		qui sum `Y' if (treatment_master >= 4 & vpn_roommate_existing == 0 & vpn_roommate_new_w3 == 2)
		local mu_102=r(mean)
		
		local labelY: var label `Y'
		outreg2 soclearning_ownaccess soclearning_rm_new_w3 soclearning_ownXnew_w3 using outregs_wave3sociallearning_additional.xls, addstat("p_all", `p_all', "q_new_1_0", `q_new_1_0', "q_new_1_1", `q_new_1_1', "pred_mu_002", `pred_mu_002', "mu_002", `mu_002', "pred_mu_102", `pred_mu_102', "mu_102", `mu_102') label ctitle(`Y', "`labelY'", ) addtext(Specification, Raw Diff) br		
		}
		
	restore
	
	
	

*** Table A.17

	cap erase outregs_wave2_maineffects.xls
	cap erase outregs_wave2_maineffects.txt	
	
	preserve
	keep if panelmerged_wave2 == 1
	
	foreach Y in `var_info_ranking_reg_w2' `var_info_freq_reg_w2' `var_media_valuation_reg_w2' az_belief_media_value_w2 `var_media_trust_reg_w2' az_belief_media_trust_w2 `var_censor_level_reg_w2' `var_censor_justif_reg_w2' az_belief_media_justif_w2  `var_censor_driver_reg_dom_w2' `var_censor_driver_reg_for_w2' `var_percmediabias_reg_ce_cn_w2' az_belief_media_cens_cn `var_percmediabias_reg_ce_us_w2' az_belief_media_cens_us `var_percmediabias_reg_di_cn_w2' az_belief_media_bias_cn `var_percmediabias_reg_di_us_w2' az_belief_media_bias_us `var_knowledge_news_reg_quiz' news_perccor_qui `var_knowledge_news_reg_cen_w2' news_perccor_cen_w2 `var_knowledge_news_reg_unc_w2' news_perccor_unc_w2 `var_knowledge_prot_reg_chi_w2' protest_pcheard_china_w2 `var_knowledge_prot_reg_for_w2' protest_pcheard_foreign_w2 `var_knowledge_prot_reg_fak_w2' `var_knowledge_people_reg_tocens' people_perchd_unctocen `var_knowledge_people_reg_censor' people_perchd_censored `var_knowledge_people_reg_uncens' people_perchd_uncensor `var_knowledge_people_reg_fake' `var_knowledge_meta_reg_w2' az_knowledge_meta_w2 `var_econ_guess_reg_cn_perf_w2' az_belief_econ_perf_cn_w2 `var_econ_guess_reg_cn_conf_w2' az_belief_econ_conf_cn_w2 `var_econ_guess_reg_us_perf_w2' az_belief_econ_perf_us_w2 `var_econ_guess_reg_us_conf_w2' az_belief_econ_conf_us_w2 `var_demand_change_reg_w2' az_belief_instchange_w2 `var_trust_inst_reg_govt_w2' az_belief_trust_govt_w2 `var_trust_inst_reg_copo_w2' az_belief_trust_copo_w2 trust_financial_domestic_w2 `var_trust_inst_reg_foreign_w2' az_belief_trust_foreign_w2 trust_financial_foreign_w2 trust_ngo_w2 `var_eval_govt_reg_w2' az_belief_evalgovt_w2 `var_eval_criteria_reg_w2' `var_severity_reg_w2' az_belief_severity_w2 `var_democracy_reg_w2' az_belief_democracy_w2 `var_contro_justi_reg_policy_w2' az_belief_justify_policy_w2 `var_contro_justi_reg_liberal_w2' az_belief_justify_liberal_w2 `var_willing_fight_reg_w2' az_belief_willing_w2 `var_interest_reg_w2' az_belief_interest_w2 `var_patriotism_reg_w2' `var_fear_critgovt_reg_w2' `var_socialinteract_reg_w2' `var_polparticipation_reg_w2' `var_stock_invest_reg_w2' `var_planaftergrad_reg_w2' `var_career_sector_reg_w2' `var_career_loc_reg_w2' {
		
		capture{
		
		// raw regression
		qui reg `Y' treatment_vpnonly treatment_nlonly treatment_vpnnl if treatment_user != 1, r
		qui summ `Y' if treatment_user != 1
		local m_all = `r(mean)'
		local sd_all = `r(sd)'
		qui summ `Y' if treatment_control == 1
		local m_control = `r(mean)'
		local sd_control = `r(sd)'
		qui summ `Y' if treatment_user == 1
		local m_user = `r(mean)'
		local sd_user = `r(sd)'		
		local labelY: var label `Y'
		local t_treatvpnXnl = _b[treatment_vpnnl]/_se[treatment_vpnnl]
		local p_treatvpnXnl = 2*ttail(e(df_r),abs(`t_treatvpnXnl'))
		outreg2 treatment_vpn treatment_newsletter treatvpnXnl using outregs_wave2_maineffects.xls, addstat("Mean DV all", `m_all', "Std.Dev.\ DV all", `sd_all', "Mean DV control", `m_control', "Std.Dev.\ DV control", `sd_control', "Mean DV user", `m_user', "Std.Dev.\ DV user", `sd_user', p-value, `p_treatvpnXnl') label dec(3) adec(3) ctitle(`Y', "`labelY'", ) addtext(Specification, Raw Diff) br		
			}
			
		}

		
	foreach Y in `var_econ_guess_reg_cn_perf_w2' `var_econ_guess_reg_us_perf_w2' {
			
		capture{
		
		// raw regression
		qui reg `Y' treatment_vpnonly treatment_nlonly treatment_vpnnl if treatment_user != 1, r
		qui summ `Y' if treatment_user != 1
		local m_all = `r(mean)'
		local sd_all = `r(sd)'
		qui summ `Y' if treatment_control == 1
		local m_control = `r(mean)'
		local sd_control = `r(sd)'
		qui summ `Y' if treatment_user == 1
		local m_user = `r(mean)'
		local sd_user = `r(sd)'		
		local labelY: var label `Y'
		local t_treatvpnXnl = _b[treatment_vpnnl]/_se[treatment_vpnnl]
		local p_treatvpnXnl = 2*ttail(e(df_r),abs(`t_treatvpnXnl'))
		outreg2 treatment_vpn treatment_newsletter treatvpnXnl using outregs_wave2_maineffects.xls, addstat("Mean DV all", `m_all', "Std.Dev.\ DV all", `sd_all', "Mean DV control", `m_control', "Std.Dev.\ DV control", `sd_control', "Mean DV user", `m_user', "Std.Dev.\ DV user", `sd_user', p-value, `p_treatvpnXnl') label dec(8) adec(8) ctitle(`Y', "`labelY'", ) addtext(Specification, Raw Diff) br		
			}
			
		}

	restore
