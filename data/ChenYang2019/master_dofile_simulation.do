*** Do-file: replication of Chen and Yang (2018)
*** Part 3: simulation exercises
*** December 2018


clear all
set more off
set maxvar 120000


*** Define root name for data folder
	local directory_path 	`""'
	cd "`directory_path'/Outregs"



	
*** Simulation of knowledge on Panama Papers: clustering before 21 but not after
	
	clear all
	
	// primitives
	set obs 10000
	gen id = _n
	
	// dorm ids (4 people per room)
	gen dorm = _n
	replace dorm = dorm - 7500 	if dorm > 7500
	replace dorm = dorm - 5000 	if dorm > 5000
	replace dorm = dorm - 2500 	if dorm > 2500
	
	// learning rates
	gen alpha = 0.560
	gen direct = 0.332
	gen indirect_0 = 0.246
	gen indirect_1 = 0.121
	
	// assigning access: completely randomly across population
	gen access_random_1 = rnormal()
	egen access_random = rank(access_random_1)
	drop access_random_1
	
	// social network: prior to 21% with clustering, after 21% without clustering
	
	preserve
	clear all
	set obs 10000
	gen id = _n
	gen access_newdorm_1 = runiform() * 12
	gen access_newdorm = (access_newdorm_1 > 11)
	drop access_newdorm_1
	rename id newdorm_id
	save random_access_newdorm, replace
	restore
	
	gen access_cluster = 1 	if id == 1
	
	forvalue i = 2/2100 {
		bysort access_cluster: egen dorm_yes1 = sum(access_cluster)
		bysort dorm: egen dorm_yes2 = sum(dorm_yes1)
		gen dorm_yes = (dorm_yes2 > 0)
		gen newdorm_id = `i'
		
		preserve
		use random_access_newdorm, clear
		keep if newdorm_id == `i'
		save random_access_newdorm_temp, replace
		restore
		merge m:1 newdorm_id using random_access_newdorm_temp
		keep if _merge == 3
		drop _merge
		
		bysort dorm: egen total_dorm = count(access_cluster)
		gen total_dorm_existing1 = (total_dorm < 4) 	if dorm_yes == 1
		egen total_dorm_existing = sum(total_dorm_existing)
		
		gen access_cluster_10 = rnormal() if (dorm_yes == 0)
		gen access_cluster_11 = rnormal() if (dorm_yes == 1)
		gen access_cluster_1 = access_cluster_10 		if (access_newdorm == 1 | total_dorm_existing == 0)
		replace access_cluster_1 = access_cluster_11 	if (access_newdorm == 0 & total_dorm_existing > 0)
		replace access_cluster_1 = . if access_cluster != .
		egen access_cluster_max = max(access_cluster_1)
		replace access_cluster = `i' if access_cluster_1 == access_cluster_max
		drop dorm_yes1 dorm_yes2 dorm_yes newdorm_id access_cluster_10 access_cluster_11 access_cluster_1 access_cluster_max total_dorm total_dorm_existing1 total_dorm_existing
		erase random_access_newdorm_temp.dta
		}

	gen access_cluster_r = rnormal()  						if access_cluster == .
	egen access_cluster_rr = rank(access_cluster_r) 		if access_cluster == .
	replace access_cluster = access_cluster_rr + 2100 		if access_cluster_rr != .
	drop access_cluster_r access_cluster_rr
	
	// assign access and roommate status
	forvalue i = 1/10000 {
		gen soclearning_ownaccess_`i' = (access_cluster <= `i')
		bysort dorm: gen soclearning_rm_new_`i' = sum(soclearning_ownaccess_`i')
		replace soclearning_rm_new_`i' = soclearning_rm_new_`i' - 1 	if soclearning_rm_new_`i' > 1
		}
		
	// learning probability, step 1: pure direct learning
	forvalue i = 1/10000 {
		gen quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i'
		gen quiz_random_`i' = runiform()
		gen quiz_yes_`i' = (quiz_random_`i' < quiz_rate_`i')
		drop quiz_rate_`i' quiz_random_`i'
		}
		
	// total number of people per dorm who know the info from direct learning (with access)
	forvalue i = 1/10000 {
		gen quiz_yesdrac_`i' = quiz_yes_`i' 	if soclearning_ownaccess_`i' == 1
		bysort dorm: egen quiz_yesdracd_`i' = sum(quiz_yesdrac_`i')
		replace quiz_yesdracd_`i' = quiz_yesdracd_`i' - quiz_yes_`i'  if soclearning_ownaccess_`i' == 1
		drop quiz_yesdrac_`i'
		}
	
	// learning probability, step 2: first degree social learning
	forvalue i = 1/10000 {
		gen quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i' + (1-((1-indirect_0)*(1-soclearning_ownaccess_`i') + (1-indirect_1)*(soclearning_ownaccess_`i'))^quiz_yesdracd_`i')	if quiz_yesdracd_`i' > 0
		replace quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i' 	if quiz_yesdracd_`i' == 0
		gen quiz_random_`i' = runiform()
		gen quiz_yessl_`i' = (quiz_random_`i' < quiz_rate_`i')
		replace quiz_yes_`i' = quiz_yessl_`i' if (quiz_yesdracd_`i' > 0 & quiz_yes_`i' == 0)
		drop quiz_rate_`i' quiz_random_`i' soclearning_ownaccess_`i' soclearning_rm_new_`i' quiz_yessl_`i'
		}
		
	// actual learning rates
	forvalue i = 1/10000 {
		egen quiz_total_`i' = sum(quiz_yes_`i')
		drop quiz_yes_`i'
		}
	keep quiz_total_*	
	duplicates drop	
	
	// reshape (split into 4 parts)
	gen i = 1

	preserve
	forvalue i = 1/10000 {
		if `i' > 2500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_1, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 2500 | `i' > 5000 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_2, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 5000 | `i' > 7500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_3, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 7500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_4, replace
	restore

	// append and save
	use simulation_temp_1, clear
	append using simulation_temp_2
	append using simulation_temp_3
	append using simulation_temp_4
	erase simulation_temp_1.dta
	erase simulation_temp_2.dta
	erase simulation_temp_3.dta
	erase simulation_temp_4.dta
	
	drop i
	rename quiz_total_ quiz_total
	gen quiz_rate = quiz_total / 10000
	save simulation_access_cluster_withsociallearning, replace
	


	
*** Simulation of knowledge on Panama Papers: clustering before 21 but not after, no social learning

	clear all
	
	// primitives
	set obs 10000
	gen id = _n
	
	// dorm ids (4 people per room)
	gen dorm = _n
	replace dorm = dorm - 7500 	if dorm > 7500
	replace dorm = dorm - 5000 	if dorm > 5000
	replace dorm = dorm - 2500 	if dorm > 2500
	
	// learning rates
	gen alpha = 0.560
	gen direct = 0.332
	gen indirect_0 = 0.246
	gen indirect_1 = 0.121
	
	// assigning access: completely randomly across population
	gen access_random_1 = rnormal()
	egen access_random = rank(access_random_1)
	drop access_random_1
	
	// social network: prior to 21% with clustering, after 21% without clustering
	
	preserve
	clear all
	set obs 10000
	gen id = _n
	gen access_newdorm_1 = runiform() * 12
	gen access_newdorm = (access_newdorm_1 > 11)
	drop access_newdorm_1
	rename id newdorm_id
	save random_access_newdorm, replace
	restore
	
	gen access_cluster = 1 	if id == 1
	
	forvalue i = 2/2100 {
		bysort access_cluster: egen dorm_yes1 = sum(access_cluster)
		bysort dorm: egen dorm_yes2 = sum(dorm_yes1)
		gen dorm_yes = (dorm_yes2 > 0)
		gen newdorm_id = `i'
		
		preserve
		use random_access_newdorm, clear
		keep if newdorm_id == `i'
		save random_access_newdorm_temp, replace
		restore
		merge m:1 newdorm_id using random_access_newdorm_temp
		keep if _merge == 3
		drop _merge
		
		bysort dorm: egen total_dorm = count(access_cluster)
		gen total_dorm_existing1 = (total_dorm < 4) 	if dorm_yes == 1
		egen total_dorm_existing = sum(total_dorm_existing)
		
		gen access_cluster_10 = rnormal() if (dorm_yes == 0)
		gen access_cluster_11 = rnormal() if (dorm_yes == 1)
		gen access_cluster_1 = access_cluster_10 		if (access_newdorm == 1 | total_dorm_existing == 0)
		replace access_cluster_1 = access_cluster_11 	if (access_newdorm == 0 & total_dorm_existing > 0)
		replace access_cluster_1 = . if access_cluster != .
		egen access_cluster_max = max(access_cluster_1)
		replace access_cluster = `i' if access_cluster_1 == access_cluster_max
		drop dorm_yes1 dorm_yes2 dorm_yes newdorm_id access_cluster_10 access_cluster_11 access_cluster_1 access_cluster_max total_dorm total_dorm_existing1 total_dorm_existing
		erase random_access_newdorm_temp.dta
		}

	gen access_cluster_r = rnormal()  						if access_cluster == .
	egen access_cluster_rr = rank(access_cluster_r) 		if access_cluster == .
	replace access_cluster = access_cluster_rr + 2100 		if access_cluster_rr != .
	drop access_cluster_r access_cluster_rr
	
	// assign access and roommate status
	forvalue i = 1/10000 {
		gen soclearning_ownaccess_`i' = (access_cluster <= `i')
		bysort dorm: gen soclearning_rm_new_`i' = sum(soclearning_ownaccess_`i')
		replace soclearning_rm_new_`i' = soclearning_rm_new_`i' - 1 	if soclearning_rm_new_`i' > 1
		}
		
	// learning probability
	forvalue i = 1/10000 {
		gen quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i'
		gen quiz_random_`i' = runiform()
		gen quiz_yes_`i' = (quiz_random_`i' < quiz_rate_`i')
		drop quiz_rate_`i' quiz_random_`i' soclearning_ownaccess_`i' soclearning_rm_new_`i'
		}
		
	// actual learning rates
	forvalue i = 1/10000 {
		egen quiz_total_`i' = sum(quiz_yes_`i')
		drop quiz_yes_`i'
		}
	keep quiz_total_*	
	duplicates drop
	
	// reshape (split into 4 parts)
	gen i = 1

	preserve
	forvalue i = 1/10000 {
		if `i' > 2500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_1, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 2500 | `i' > 5000 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_2, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 5000 | `i' > 7500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_3, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 7500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_4, replace
	restore

	// append and save
	use simulation_temp_1, clear
	append using simulation_temp_2
	append using simulation_temp_3
	append using simulation_temp_4
	erase simulation_temp_1.dta
	erase simulation_temp_2.dta
	erase simulation_temp_3.dta
	erase simulation_temp_4.dta
	
	drop i
	rename quiz_total_ quiz_total
	gen quiz_rate = quiz_total / 10000
	save simulation_access_cluster_nosociallearning, replace
	
		
	
	
*** Simulation of knowledge on Panama Papers: clustering throughout

	clear all
	
	// primitives
	set obs 10000
	gen id = _n
	
	// dorm ids (4 people per room)
	gen dorm = _n
	replace dorm = dorm - 7500 	if dorm > 7500
	replace dorm = dorm - 5000 	if dorm > 5000
	replace dorm = dorm - 2500 	if dorm > 2500
	
	// learning rates
	gen alpha = 0.560
	gen direct = 0.332
	gen indirect_0 = 0.246
	gen indirect_1 = 0.121
	
	// assigning access: completely randomly across population
	gen access_random_1 = rnormal()
	egen access_random = rank(access_random_1)
	drop access_random_1
	
	// social network: prior to 21% with clustering, after 21% without clustering
	
	preserve
	clear all
	set obs 10000
	gen id = _n
	gen access_newdorm_1 = runiform() * 12
	gen access_newdorm = (access_newdorm_1 > 11)
	drop access_newdorm_1
	rename id newdorm_id
	save random_access_newdorm, replace
	restore
	
	gen access_cluster = 1 	if id == 1
	
	forvalue i = 2/10000 {
		bysort access_cluster: egen dorm_yes1 = sum(access_cluster)
		bysort dorm: egen dorm_yes2 = sum(dorm_yes1)
		gen dorm_yes = (dorm_yes2 > 0)
		gen newdorm_id = `i'
		
		preserve
		use random_access_newdorm, clear
		keep if newdorm_id == `i'
		save random_access_newdorm_temp, replace
		restore
		merge m:1 newdorm_id using random_access_newdorm_temp
		keep if _merge == 3
		drop _merge
		
		bysort dorm: egen total_dorm = count(access_cluster)
		gen total_dorm_existing1 = (total_dorm < 4) 	if dorm_yes == 1
		egen total_dorm_existing = sum(total_dorm_existing)
		
		gen access_cluster_10 = rnormal() if (dorm_yes == 0)
		gen access_cluster_11 = rnormal() if (dorm_yes == 1)
		gen access_cluster_1 = access_cluster_10 		if (access_newdorm == 1 | total_dorm_existing == 0)
		replace access_cluster_1 = access_cluster_11 	if (access_newdorm == 0 & total_dorm_existing > 0)
		replace access_cluster_1 = . if access_cluster != .
		egen access_cluster_max = max(access_cluster_1)
		replace access_cluster = `i' if access_cluster_1 == access_cluster_max
		drop dorm_yes1 dorm_yes2 dorm_yes newdorm_id access_cluster_10 access_cluster_11 access_cluster_1 access_cluster_max total_dorm total_dorm_existing1 total_dorm_existing
		erase random_access_newdorm_temp.dta
		}

	// assign access and roommate status
	forvalue i = 1/10000 {
		gen soclearning_ownaccess_`i' = (access_cluster <= `i')
		bysort dorm: gen soclearning_rm_new_`i' = sum(soclearning_ownaccess_`i')
		replace soclearning_rm_new_`i' = soclearning_rm_new_`i' - 1 	if soclearning_rm_new_`i' > 1
		}
		
	// learning probability, step 1: pure direct learning
	forvalue i = 1/10000 {
		gen quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i'
		gen quiz_random_`i' = runiform()
		gen quiz_yes_`i' = (quiz_random_`i' < quiz_rate_`i')
		drop quiz_rate_`i' quiz_random_`i'
		}
		
	// total number of people per dorm who know the info from direct learning (with access)
	forvalue i = 1/10000 {
		gen quiz_yesdrac_`i' = quiz_yes_`i' 	if soclearning_ownaccess_`i' == 1
		bysort dorm: egen quiz_yesdracd_`i' = sum(quiz_yesdrac_`i')
		replace quiz_yesdracd_`i' = quiz_yesdracd_`i' - quiz_yes_`i'  if soclearning_ownaccess_`i' == 1
		drop quiz_yesdrac_`i'
		}
	
	// learning probability, step 2: first degree social learning
	forvalue i = 1/10000 {
		gen quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i' + (1-((1-indirect_0)*(1-soclearning_ownaccess_`i') + (1-indirect_1)*(soclearning_ownaccess_`i'))^quiz_yesdracd_`i')	if quiz_yesdracd_`i' > 0
		replace quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i' 	if quiz_yesdracd_`i' == 0
		gen quiz_random_`i' = runiform()
		gen quiz_yessl_`i' = (quiz_random_`i' < quiz_rate_`i')
		replace quiz_yes_`i' = quiz_yessl_`i' if (quiz_yesdracd_`i' > 0 & quiz_yes_`i' == 0)
		drop quiz_rate_`i' quiz_random_`i' soclearning_ownaccess_`i' soclearning_rm_new_`i' quiz_yessl_`i'
		}
		
	// actual learning rates
	forvalue i = 1/10000 {
		egen quiz_total_`i' = sum(quiz_yes_`i')
		drop quiz_yes_`i'
		}
	keep quiz_total_*	
	duplicates drop
	
	// reshape (split into 4 parts)
	gen i = 1

	preserve
	forvalue i = 1/10000 {
		if `i' > 2500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_1, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 2500 | `i' > 5000 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_2, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 5000 | `i' > 7500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_3, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 7500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_4, replace
	restore

	// append and save
	use simulation_temp_1, clear
	append using simulation_temp_2
	append using simulation_temp_3
	append using simulation_temp_4
	erase simulation_temp_1.dta
	erase simulation_temp_2.dta
	erase simulation_temp_3.dta
	erase simulation_temp_4.dta
	
	drop i
	rename quiz_total_ quiz_total
	gen quiz_rate = quiz_total / 10000
	save simulation_access_allcluster_withsociallearning, replace
	
	

	
*** Simulation of knowledge on Panama Papers: 8 people dorm

	clear all
	
	// primitives
	set obs 10000
	gen id = _n
	
	// dorm ids (8 people per room)
	gen dorm = _n
	replace dorm = dorm - 8750 	if dorm > 8750
	replace dorm = dorm - 7500 	if dorm > 7500
	replace dorm = dorm - 6250 	if dorm > 6250
	replace dorm = dorm - 5000 	if dorm > 5000
	replace dorm = dorm - 3750 	if dorm > 3750
	replace dorm = dorm - 2500 	if dorm > 2500
	replace dorm = dorm - 1250 	if dorm > 1250

	// learning rates
	gen alpha = 0.560
	gen direct = 0.332
	gen indirect_0 = 0.246
	gen indirect_1 = 0.121
	
	// assigning access: completely randomly across population
	gen access_random_1 = rnormal()
	egen access_random = rank(access_random_1)
	drop access_random_1
	
	// social network: prior to 21% with clustering, after 21% without clustering
	
	preserve
	clear all
	set obs 10000
	gen id = _n
	gen access_newdorm_1 = runiform() * 12
	gen access_newdorm = (access_newdorm_1 > 11)
	drop access_newdorm_1
	rename id newdorm_id
	save random_access_newdorm, replace
	restore
	
	gen access_cluster = 1 	if id == 1
	
	forvalue i = 2/2100 {
		bysort access_cluster: egen dorm_yes1 = sum(access_cluster)
		bysort dorm: egen dorm_yes2 = sum(dorm_yes1)
		gen dorm_yes = (dorm_yes2 > 0)
		gen newdorm_id = `i'
		
		preserve
		use random_access_newdorm, clear
		keep if newdorm_id == `i'
		save random_access_newdorm_temp, replace
		restore
		merge m:1 newdorm_id using random_access_newdorm_temp
		keep if _merge == 3
		drop _merge
		
		bysort dorm: egen total_dorm = count(access_cluster)
		gen total_dorm_existing1 = (total_dorm < 4) 	if dorm_yes == 1
		egen total_dorm_existing = sum(total_dorm_existing)
		
		gen access_cluster_10 = rnormal() if (dorm_yes == 0)
		gen access_cluster_11 = rnormal() if (dorm_yes == 1)
		gen access_cluster_1 = access_cluster_10 		if (access_newdorm == 1 | total_dorm_existing == 0)
		replace access_cluster_1 = access_cluster_11 	if (access_newdorm == 0 & total_dorm_existing > 0)
		replace access_cluster_1 = . if access_cluster != .
		egen access_cluster_max = max(access_cluster_1)
		replace access_cluster = `i' if access_cluster_1 == access_cluster_max
		drop dorm_yes1 dorm_yes2 dorm_yes newdorm_id access_cluster_10 access_cluster_11 access_cluster_1 access_cluster_max total_dorm total_dorm_existing1 total_dorm_existing
		erase random_access_newdorm_temp.dta
		}

	gen access_cluster_r = rnormal()  						if access_cluster == .
	egen access_cluster_rr = rank(access_cluster_r) 		if access_cluster == .
	replace access_cluster = access_cluster_rr + 2100 		if access_cluster_rr != .
	drop access_cluster_r access_cluster_rr
	
	// assign access and roommate status
	forvalue i = 1/10000 {
		gen soclearning_ownaccess_`i' = (access_cluster <= `i')
		bysort dorm: gen soclearning_rm_new_`i' = sum(soclearning_ownaccess_`i')
		replace soclearning_rm_new_`i' = soclearning_rm_new_`i' - 1 	if soclearning_rm_new_`i' > 1
		}
		
	// learning probability, step 1: pure direct learning
	forvalue i = 1/10000 {
		gen quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i'
		gen quiz_random_`i' = runiform()
		gen quiz_yes_`i' = (quiz_random_`i' < quiz_rate_`i')
		drop quiz_rate_`i' quiz_random_`i'
		}
		
	// total number of people per dorm who know the info from direct learning (with access)
	forvalue i = 1/10000 {
		gen quiz_yesdrac_`i' = quiz_yes_`i' 	if soclearning_ownaccess_`i' == 1
		bysort dorm: egen quiz_yesdracd_`i' = sum(quiz_yesdrac_`i')
		replace quiz_yesdracd_`i' = quiz_yesdracd_`i' - quiz_yes_`i'  if soclearning_ownaccess_`i' == 1
		drop quiz_yesdrac_`i'
		}
	
	// learning probability, step 2: first degree social learning
	forvalue i = 1/10000 {
		gen quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i' + (1-((1-indirect_0)*(1-soclearning_ownaccess_`i') + (1-indirect_1)*(soclearning_ownaccess_`i'))^quiz_yesdracd_`i')	if quiz_yesdracd_`i' > 0
		replace quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i' 	if quiz_yesdracd_`i' == 0
		gen quiz_random_`i' = runiform()
		gen quiz_yessl_`i' = (quiz_random_`i' < quiz_rate_`i')
		replace quiz_yes_`i' = quiz_yessl_`i' if (quiz_yesdracd_`i' > 0 & quiz_yes_`i' == 0)
		drop quiz_rate_`i' quiz_random_`i' soclearning_ownaccess_`i' soclearning_rm_new_`i' quiz_yessl_`i'
		}
		
	// actual learning rates
	forvalue i = 1/10000 {
		egen quiz_total_`i' = sum(quiz_yes_`i')
		drop quiz_yes_`i'
		}
	keep quiz_total_*	
	duplicates drop
	
	// reshape (split into 4 parts)
	gen i = 1

	preserve
	forvalue i = 1/10000 {
		if `i' > 2500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_1, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 2500 | `i' > 5000 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_2, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 5000 | `i' > 7500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_3, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 7500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_4, replace
	restore

	// append and save
	use simulation_temp_1, clear
	append using simulation_temp_2
	append using simulation_temp_3
	append using simulation_temp_4
	erase simulation_temp_1.dta
	erase simulation_temp_2.dta
	erase simulation_temp_3.dta
	erase simulation_temp_4.dta
	
	drop i
	rename quiz_total_ quiz_total
	gen quiz_rate = quiz_total / 10000
	save simulation_access_cluster_withsociallearning_doubledormsize, replace
	
	
	
	
*** Simulation of knowledge on Panama Papers: baseline, double social learning

	clear all
	
	// primitives
	set obs 10000
	gen id = _n
	
	// dorm ids (4 people per room)
	gen dorm = _n
	replace dorm = dorm - 7500 	if dorm > 7500
	replace dorm = dorm - 5000 	if dorm > 5000
	replace dorm = dorm - 2500 	if dorm > 2500
	
	// learning rates
	gen alpha = 0.560
	gen direct = 0.332
	gen indirect_0 = 0.246 * 2
	gen indirect_1 = 0.121 * 2
	
	// assigning access: completely randomly across population
	gen access_random_1 = rnormal()
	egen access_random = rank(access_random_1)
	drop access_random_1
	
	// social network: prior to 21% with clustering, after 21% without clustering
	
	preserve
	clear all
	set obs 10000
	gen id = _n
	gen access_newdorm_1 = runiform() * 12
	gen access_newdorm = (access_newdorm_1 > 11)
	drop access_newdorm_1
	rename id newdorm_id
	save random_access_newdorm, replace
	restore
	
	gen access_cluster = 1 	if id == 1
	
	forvalue i = 2/2100 {
		bysort access_cluster: egen dorm_yes1 = sum(access_cluster)
		bysort dorm: egen dorm_yes2 = sum(dorm_yes1)
		gen dorm_yes = (dorm_yes2 > 0)
		gen newdorm_id = `i'
		
		preserve
		use random_access_newdorm, clear
		keep if newdorm_id == `i'
		save random_access_newdorm_temp, replace
		restore
		merge m:1 newdorm_id using random_access_newdorm_temp
		keep if _merge == 3
		drop _merge
		
		bysort dorm: egen total_dorm = count(access_cluster)
		gen total_dorm_existing1 = (total_dorm < 4) 	if dorm_yes == 1
		egen total_dorm_existing = sum(total_dorm_existing)
		
		gen access_cluster_10 = rnormal() if (dorm_yes == 0)
		gen access_cluster_11 = rnormal() if (dorm_yes == 1)
		gen access_cluster_1 = access_cluster_10 		if (access_newdorm == 1 | total_dorm_existing == 0)
		replace access_cluster_1 = access_cluster_11 	if (access_newdorm == 0 & total_dorm_existing > 0)
		replace access_cluster_1 = . if access_cluster != .
		egen access_cluster_max = max(access_cluster_1)
		replace access_cluster = `i' if access_cluster_1 == access_cluster_max
		drop dorm_yes1 dorm_yes2 dorm_yes newdorm_id access_cluster_10 access_cluster_11 access_cluster_1 access_cluster_max total_dorm total_dorm_existing1 total_dorm_existing
		erase random_access_newdorm_temp.dta
		}

	gen access_cluster_r = rnormal()  						if access_cluster == .
	egen access_cluster_rr = rank(access_cluster_r) 		if access_cluster == .
	replace access_cluster = access_cluster_rr + 2100 		if access_cluster_rr != .
	drop access_cluster_r access_cluster_rr
	
	// assign access and roommate status
	forvalue i = 1/10000 {
		gen soclearning_ownaccess_`i' = (access_cluster <= `i')
		bysort dorm: gen soclearning_rm_new_`i' = sum(soclearning_ownaccess_`i')
		replace soclearning_rm_new_`i' = soclearning_rm_new_`i' - 1 	if soclearning_rm_new_`i' > 1
		}
		
	// learning probability, step 1: pure direct learning
	forvalue i = 1/10000 {
		gen quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i'
		gen quiz_random_`i' = runiform()
		gen quiz_yes_`i' = (quiz_random_`i' < quiz_rate_`i')
		drop quiz_rate_`i' quiz_random_`i'
		}
		
	// total number of people per dorm who know the info from direct learning (with access)
	forvalue i = 1/10000 {
		gen quiz_yesdrac_`i' = quiz_yes_`i' 	if soclearning_ownaccess_`i' == 1
		bysort dorm: egen quiz_yesdracd_`i' = sum(quiz_yesdrac_`i')
		replace quiz_yesdracd_`i' = quiz_yesdracd_`i' - quiz_yes_`i'  if soclearning_ownaccess_`i' == 1
		drop quiz_yesdrac_`i'
		}
	
	// learning probability, step 2: first degree social learning
	forvalue i = 1/10000 {
		gen quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i' + (1-((1-indirect_0)*(1-soclearning_ownaccess_`i') + (1-indirect_1)*(soclearning_ownaccess_`i'))^quiz_yesdracd_`i')	if quiz_yesdracd_`i' > 0
		replace quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i' 	if quiz_yesdracd_`i' == 0
		gen quiz_random_`i' = runiform()
		gen quiz_yessl_`i' = (quiz_random_`i' < quiz_rate_`i')
		replace quiz_yes_`i' = quiz_yessl_`i' if (quiz_yesdracd_`i' > 0 & quiz_yes_`i' == 0)
		drop quiz_rate_`i' quiz_random_`i' soclearning_ownaccess_`i' soclearning_rm_new_`i' quiz_yessl_`i'
		}
		
	// actual learning rates
	forvalue i = 1/10000 {
		egen quiz_total_`i' = sum(quiz_yes_`i')
		drop quiz_yes_`i'
		}
	keep quiz_total_*	
	duplicates drop
	
	// reshape (split into 4 parts)
	gen i = 1

	preserve
	forvalue i = 1/10000 {
		if `i' > 2500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_1, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 2500 | `i' > 5000 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_2, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 5000 | `i' > 7500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_3, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 7500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_4, replace
	restore

	// append and save
	use simulation_temp_1, clear
	append using simulation_temp_2
	append using simulation_temp_3
	append using simulation_temp_4
	erase simulation_temp_1.dta
	erase simulation_temp_2.dta
	erase simulation_temp_3.dta
	erase simulation_temp_4.dta
	
	drop i
	rename quiz_total_ quiz_total
	gen quiz_rate = quiz_total / 10000
	save simulation_access_cluster_withsociallearning_doublesociallearning, replace
	

	
*** Simulation of knowledge on Panama Papers: baseline, double direct learning

	clear all
	
	// primitives
	set obs 10000
	gen id = _n
	
	// dorm ids (4 people per room)
	gen dorm = _n
	replace dorm = dorm - 7500 	if dorm > 7500
	replace dorm = dorm - 5000 	if dorm > 5000
	replace dorm = dorm - 2500 	if dorm > 2500
	
	// learning rates
	gen alpha = 0.560
	gen direct = 0.332 * 2
	gen indirect_0 = 0.246
	gen indirect_1 = 0.121
	
	// assigning access: completely randomly across population
	gen access_random_1 = rnormal()
	egen access_random = rank(access_random_1)
	drop access_random_1
	
	// social network: prior to 21% with clustering, after 21% without clustering
	
	preserve
	clear all
	set obs 10000
	gen id = _n
	gen access_newdorm_1 = runiform() * 12
	gen access_newdorm = (access_newdorm_1 > 11)
	drop access_newdorm_1
	rename id newdorm_id
	save random_access_newdorm, replace
	restore
	
	gen access_cluster = 1 	if id == 1
	
	forvalue i = 2/2100 {
		bysort access_cluster: egen dorm_yes1 = sum(access_cluster)
		bysort dorm: egen dorm_yes2 = sum(dorm_yes1)
		gen dorm_yes = (dorm_yes2 > 0)
		gen newdorm_id = `i'
		
		preserve
		use random_access_newdorm, clear
		keep if newdorm_id == `i'
		save random_access_newdorm_temp, replace
		restore
		merge m:1 newdorm_id using random_access_newdorm_temp
		keep if _merge == 3
		drop _merge
		
		bysort dorm: egen total_dorm = count(access_cluster)
		gen total_dorm_existing1 = (total_dorm < 4) 	if dorm_yes == 1
		egen total_dorm_existing = sum(total_dorm_existing)
		
		gen access_cluster_10 = rnormal() if (dorm_yes == 0)
		gen access_cluster_11 = rnormal() if (dorm_yes == 1)
		gen access_cluster_1 = access_cluster_10 		if (access_newdorm == 1 | total_dorm_existing == 0)
		replace access_cluster_1 = access_cluster_11 	if (access_newdorm == 0 & total_dorm_existing > 0)
		replace access_cluster_1 = . if access_cluster != .
		egen access_cluster_max = max(access_cluster_1)
		replace access_cluster = `i' if access_cluster_1 == access_cluster_max
		drop dorm_yes1 dorm_yes2 dorm_yes newdorm_id access_cluster_10 access_cluster_11 access_cluster_1 access_cluster_max total_dorm total_dorm_existing1 total_dorm_existing
		erase random_access_newdorm_temp.dta
		}

	gen access_cluster_r = rnormal()  						if access_cluster == .
	egen access_cluster_rr = rank(access_cluster_r) 		if access_cluster == .
	replace access_cluster = access_cluster_rr + 2100 		if access_cluster_rr != .
	drop access_cluster_r access_cluster_rr
	
	// assign access and roommate status
	forvalue i = 1/10000 {
		gen soclearning_ownaccess_`i' = (access_cluster <= `i')
		bysort dorm: gen soclearning_rm_new_`i' = sum(soclearning_ownaccess_`i')
		replace soclearning_rm_new_`i' = soclearning_rm_new_`i' - 1 	if soclearning_rm_new_`i' > 1
		}
		
	// learning probability, step 1: pure direct learning
	forvalue i = 1/10000 {
		gen quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i'
		gen quiz_random_`i' = runiform()
		gen quiz_yes_`i' = (quiz_random_`i' < quiz_rate_`i')
		drop quiz_rate_`i' quiz_random_`i'
		}
		
	// total number of people per dorm who know the info from direct learning (with access)
	forvalue i = 1/10000 {
		gen quiz_yesdrac_`i' = quiz_yes_`i' 	if soclearning_ownaccess_`i' == 1
		bysort dorm: egen quiz_yesdracd_`i' = sum(quiz_yesdrac_`i')
		replace quiz_yesdracd_`i' = quiz_yesdracd_`i' - quiz_yes_`i'  if soclearning_ownaccess_`i' == 1
		drop quiz_yesdrac_`i'
		}
	
	// learning probability, step 2: first degree social learning
	forvalue i = 1/10000 {
		gen quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i' + (1-((1-indirect_0)*(1-soclearning_ownaccess_`i') + (1-indirect_1)*(soclearning_ownaccess_`i'))^quiz_yesdracd_`i')	if quiz_yesdracd_`i' > 0
		replace quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i' 	if quiz_yesdracd_`i' == 0
		gen quiz_random_`i' = runiform()
		gen quiz_yessl_`i' = (quiz_random_`i' < quiz_rate_`i')
		replace quiz_yes_`i' = quiz_yessl_`i' if (quiz_yesdracd_`i' > 0 & quiz_yes_`i' == 0)
		drop quiz_rate_`i' quiz_random_`i' soclearning_ownaccess_`i' soclearning_rm_new_`i' quiz_yessl_`i'
		}
		
	// actual learning rates
	forvalue i = 1/10000 {
		egen quiz_total_`i' = sum(quiz_yes_`i')
		drop quiz_yes_`i'
		}
	keep quiz_total_*	
	duplicates drop
	
	// reshape (split into 4 parts)
	gen i = 1

	preserve
	forvalue i = 1/10000 {
		if `i' > 2500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_1, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 2500 | `i' > 5000 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_2, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 5000 | `i' > 7500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_3, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 7500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_4, replace
	restore

	// append and save
	use simulation_temp_1, clear
	append using simulation_temp_2
	append using simulation_temp_3
	append using simulation_temp_4
	erase simulation_temp_1.dta
	erase simulation_temp_2.dta
	erase simulation_temp_3.dta
	erase simulation_temp_4.dta
	
	drop i
	rename quiz_total_ quiz_total
	gen quiz_rate = quiz_total / 10000
	save simulation_access_cluster_withsociallearning_doubledirectlearning, replace
	

	
	
*** Simulation of knowledge on Panama Papers: clustering before 21, targeted diffusion after
	
	clear all
	
	// primitives
	set obs 10000
	gen id = _n
	
	// dorm ids (4 people per room)
	gen dorm = _n
	replace dorm = dorm - 7500 	if dorm > 7500
	replace dorm = dorm - 5000 	if dorm > 5000
	replace dorm = dorm - 2500 	if dorm > 2500
	
	// learning rates
	gen alpha = 0.560
	gen direct = 0.332
	gen indirect_0 = 0.246
	gen indirect_1 = 0.121
	
	// assigning access: completely randomly across population
	gen access_random_1 = rnormal()
	egen access_random = rank(access_random_1)
	drop access_random_1
	
	// social network: prior to 21% with clustering, after 21% without clustering
	
	preserve
	clear all
	set obs 10000
	gen id = _n
	gen access_newdorm_1 = runiform() * 12
	gen access_newdorm = (access_newdorm_1 > 11)
	drop access_newdorm_1
	rename id newdorm_id
	save random_access_newdorm, replace
	restore
	
	gen access_cluster = 1 	if id == 1
	
	forvalue i = 2/210 {
		bysort access_cluster: egen dorm_yes1 = sum(access_cluster)
		bysort dorm: egen dorm_yes2 = sum(dorm_yes1)
		gen dorm_yes = (dorm_yes2 > 0)
		gen newdorm_id = `i'
		
		preserve
		use random_access_newdorm, clear
		keep if newdorm_id == `i'
		save random_access_newdorm_temp, replace
		restore
		merge m:1 newdorm_id using random_access_newdorm_temp
		keep if _merge == 3
		drop _merge
		
		bysort dorm: egen total_dorm = count(access_cluster)
		gen total_dorm_existing1 = (total_dorm < 4) 	if dorm_yes == 1
		egen total_dorm_existing = sum(total_dorm_existing)
		
		gen access_cluster_10 = rnormal() if (dorm_yes == 0)
		gen access_cluster_11 = rnormal() if (dorm_yes == 1)
		gen access_cluster_1 = access_cluster_10 		if (access_newdorm == 1 | total_dorm_existing == 0)
		replace access_cluster_1 = access_cluster_11 	if (access_newdorm == 0 & total_dorm_existing > 0)
		replace access_cluster_1 = . if access_cluster != .
		egen access_cluster_max = max(access_cluster_1)
		replace access_cluster = `i' if access_cluster_1 == access_cluster_max
		drop dorm_yes1 dorm_yes2 dorm_yes newdorm_id access_cluster_10 access_cluster_11 access_cluster_1 access_cluster_max total_dorm total_dorm_existing1 total_dorm_existing
		erase random_access_newdorm_temp.dta
		}

	gen access_cluster_r = rnormal()  						if access_cluster == .
	egen access_cluster_rr = rank(access_cluster_r) 		if access_cluster == .
	replace access_cluster = access_cluster_rr + 2100 		if access_cluster_rr != .
	drop access_cluster_r access_cluster_rr

	// generate targetted dorm structure
	gen access_target_r = rnormal() 							if access_cluster > 2100
	bysort dorm: egen access_target_rr = rank(access_target_r) 	if access_cluster > 2100
	sort access_target_rr access_cluster
	gen access_target_rrr = _n + 2100 							if access_cluster > 2100
	gen access_target = access_cluster 							if access_cluster <= 2100
	replace access_target = access_target_rrr  					if access_cluster > 2100
	drop access_target_r access_target_rr access_target_rrr
	
	// assign access and roommate status
	forvalue i = 1/10000 {
		gen soclearning_ownaccess_`i' = (access_target <= `i')
		bysort dorm: gen soclearning_rm_new_`i' = sum(soclearning_ownaccess_`i')
		replace soclearning_rm_new_`i' = soclearning_rm_new_`i' - 1 	if soclearning_rm_new_`i' > 1
		}
		
	// learning probability, step 1: pure direct learning
	forvalue i = 1/10000 {
		gen quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i'
		gen quiz_random_`i' = runiform()
		gen quiz_yes_`i' = (quiz_random_`i' < quiz_rate_`i')
		drop quiz_rate_`i' quiz_random_`i'
		}
		
	// total number of people per dorm who know the info from direct learning (with access)
	forvalue i = 1/10000 {
		gen quiz_yesdrac_`i' = quiz_yes_`i' 	if soclearning_ownaccess_`i' == 1
		bysort dorm: egen quiz_yesdracd_`i' = sum(quiz_yesdrac_`i')
		replace quiz_yesdracd_`i' = quiz_yesdracd_`i' - quiz_yes_`i'  if soclearning_ownaccess_`i' == 1
		drop quiz_yesdrac_`i'
		}
	
	// learning probability, step 2: first degree social learning
	forvalue i = 1/10000 {
		gen quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i' + (1-((1-indirect_0)*(1-soclearning_ownaccess_`i') + (1-indirect_1)*(soclearning_ownaccess_`i'))^quiz_yesdracd_`i')	if quiz_yesdracd_`i' > 0
		replace quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i' 	if quiz_yesdracd_`i' == 0
		gen quiz_random_`i' = runiform()
		gen quiz_yessl_`i' = (quiz_random_`i' < quiz_rate_`i')
		replace quiz_yes_`i' = quiz_yessl_`i' if (quiz_yesdracd_`i' > 0 & quiz_yes_`i' == 0)
		drop quiz_rate_`i' quiz_random_`i' soclearning_ownaccess_`i' soclearning_rm_new_`i' quiz_yessl_`i'
		}
		
	// actual learning rates
	forvalue i = 1/10000 {
		egen quiz_total_`i' = sum(quiz_yes_`i')
		drop quiz_yes_`i'
		}
	keep quiz_total_*	
	duplicates drop
	
	// reshape (split into 4 parts)
	gen i = 1

	preserve
	forvalue i = 1/10000 {
		if `i' > 2500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_1, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 2500 | `i' > 5000 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_2, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 5000 | `i' > 7500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_3, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 7500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_4, replace
	restore

	// append and save
	use simulation_temp_1, clear
	append using simulation_temp_2
	append using simulation_temp_3
	append using simulation_temp_4
	erase simulation_temp_1.dta
	erase simulation_temp_2.dta
	erase simulation_temp_3.dta
	erase simulation_temp_4.dta
	
	drop i
	rename quiz_total_ quiz_total
	gen quiz_rate = quiz_total / 10000
	save simulation_access_cluster_targeteddiffusion, replace

	

	
*** Simulation of knowledge on Panama Papers: clustering before 21 but not after; allow 2 degrees of transmission
	
	clear all
	
	// primitives
	set obs 10000
	gen id = _n
	
	// dorm ids (4 people per room)
	gen dorm = _n
	replace dorm = dorm - 7500 	if dorm > 7500
	replace dorm = dorm - 5000 	if dorm > 5000
	replace dorm = dorm - 2500 	if dorm > 2500
	
	// learning rates
	gen alpha = 0.560
	gen direct = 0.332
	gen indirect_0 = 0.246
	gen indirect_1 = 0.121
	
	// assigning access: completely randomly across population
	gen access_random_1 = rnormal()
	egen access_random = rank(access_random_1)
	drop access_random_1
	
	// social network: prior to 21% with clustering, after 21% without clustering
	
	preserve
	clear all
	set obs 10000
	gen id = _n
	gen access_newdorm_1 = runiform() * 12
	gen access_newdorm = (access_newdorm_1 > 11)
	drop access_newdorm_1
	rename id newdorm_id
	save random_access_newdorm, replace
	restore
	
	gen access_cluster = 1 	if id == 1
	
	forvalue i = 2/2100 {
		bysort access_cluster: egen dorm_yes1 = sum(access_cluster)
		bysort dorm: egen dorm_yes2 = sum(dorm_yes1)
		gen dorm_yes = (dorm_yes2 > 0)
		gen newdorm_id = `i'
		
		preserve
		use random_access_newdorm, clear
		keep if newdorm_id == `i'
		save random_access_newdorm_temp, replace
		restore
		merge m:1 newdorm_id using random_access_newdorm_temp
		keep if _merge == 3
		drop _merge
		
		bysort dorm: egen total_dorm = count(access_cluster)
		gen total_dorm_existing1 = (total_dorm < 4) 	if dorm_yes == 1
		egen total_dorm_existing = sum(total_dorm_existing)
		
		gen access_cluster_10 = rnormal() if (dorm_yes == 0)
		gen access_cluster_11 = rnormal() if (dorm_yes == 1)
		gen access_cluster_1 = access_cluster_10 		if (access_newdorm == 1 | total_dorm_existing == 0)
		replace access_cluster_1 = access_cluster_11 	if (access_newdorm == 0 & total_dorm_existing > 0)
		replace access_cluster_1 = . if access_cluster != .
		egen access_cluster_max = max(access_cluster_1)
		replace access_cluster = `i' if access_cluster_1 == access_cluster_max
		drop dorm_yes1 dorm_yes2 dorm_yes newdorm_id access_cluster_10 access_cluster_11 access_cluster_1 access_cluster_max total_dorm total_dorm_existing1 total_dorm_existing
		erase random_access_newdorm_temp.dta
		}

	gen access_cluster_r = rnormal()  						if access_cluster == .
	egen access_cluster_rr = rank(access_cluster_r) 		if access_cluster == .
	replace access_cluster = access_cluster_rr + 2100 		if access_cluster_rr != .
	drop access_cluster_r access_cluster_rr
	
	// assign access and roommate status
	forvalue i = 1/10000 {
		gen soclearning_ownaccess_`i' = (access_cluster <= `i')
		bysort dorm: gen soclearning_rm_new_`i' = sum(soclearning_ownaccess_`i')
		replace soclearning_rm_new_`i' = soclearning_rm_new_`i' - 1 	if soclearning_rm_new_`i' > 1
		}
		
	// learning probability, step 1: pure direct learning
	forvalue i = 1/10000 {
		gen quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i'
		gen quiz_random_`i' = runiform()
		gen quiz_yes_`i' = (quiz_random_`i' < quiz_rate_`i')
		drop quiz_rate_`i' quiz_random_`i'
		}
		
	// total number of people per dorm who know the info from direct learning (with access)
	forvalue i = 1/10000 {
		gen quiz_yesdrac_`i' = quiz_yes_`i' 	if soclearning_ownaccess_`i' == 1
		bysort dorm: egen quiz_yesdracd_`i' = sum(quiz_yesdrac_`i')
		replace quiz_yesdracd_`i' = quiz_yesdracd_`i' - quiz_yes_`i'  if soclearning_ownaccess_`i' == 1
		drop quiz_yesdrac_`i'
		}
	
	// learning probability, step 2: first degree social learning
	forvalue i = 1/10000 {
		gen quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i' + (1-((1-indirect_0)*(1-soclearning_ownaccess_`i') + (1-indirect_1)*(soclearning_ownaccess_`i'))^quiz_yesdracd_`i')	if quiz_yesdracd_`i' > 0
		replace quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i' 	if quiz_yesdracd_`i' == 0
		gen quiz_random_`i' = runiform()
		gen quiz_yessl_`i' = (quiz_random_`i' < quiz_rate_`i')
		replace quiz_yes_`i' = quiz_yessl_`i' if (quiz_yesdracd_`i' > 0 & quiz_yes_`i' == 0)
		drop quiz_rate_`i' quiz_random_`i' quiz_yessl_`i'
		}
		
	// total number of people per dorm who know the info, after 1st degree of transmission
	forvalue i = 1/10000 {
		gen quiz_yesdrac2_`i' = quiz_yes_`i' 	if soclearning_ownaccess_`i' == 1
		bysort dorm: egen quiz_yesdracd2_`i' = sum(quiz_yesdrac2_`i')
		replace quiz_yesdracd2_`i' = quiz_yesdracd2_`i' - quiz_yes_`i'  if soclearning_ownaccess_`i' == 1
		drop quiz_yesdrac2_`i'
		}
	
	// learning probability, step 3: second degree social learning
	forvalue i = 1/10000 {
		gen quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i' + (1-((1-indirect_0)*(1-soclearning_ownaccess_`i') + (1-indirect_1)*(soclearning_ownaccess_`i'))^quiz_yesdracd2_`i')	if quiz_yesdracd2_`i' > 0
		replace quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i' 	if quiz_yesdracd2_`i' == 0
		gen quiz_random_`i' = runiform()
		gen quiz_yessl_`i' = (quiz_random_`i' < quiz_rate_`i')
		replace quiz_yes_`i' = quiz_yessl_`i' if ((quiz_yesdracd2_`i' > quiz_yesdracd_`i') & quiz_yes_`i' == 0)
		drop quiz_rate_`i' quiz_random_`i' quiz_yessl_`i'
		}
		
	// actual learning rates
	forvalue i = 1/10000 {
		egen quiz_total_`i' = sum(quiz_yes_`i')
		drop quiz_yes_`i'
		}
	keep quiz_total_*	
	duplicates drop
	
	// reshape (split into 4 parts)
	gen i = 1

	preserve
	forvalue i = 1/10000 {
		if `i' > 2500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_1, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 2500 | `i' > 5000 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_2, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 5000 | `i' > 7500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_3, replace
	restore

	preserve
	forvalue i = 1/10000 {
		if `i' <= 7500 {
		drop quiz_total_`i'
		}
		}
	reshape long quiz_total_, i(i) j(count)
	save simulation_temp_4, replace
	restore

	// append and save
	use simulation_temp_1, clear
	append using simulation_temp_2
	append using simulation_temp_3
	append using simulation_temp_4
	erase simulation_temp_1.dta
	erase simulation_temp_2.dta
	erase simulation_temp_3.dta
	erase simulation_temp_4.dta
	
	drop i
	rename quiz_total_ quiz_total
	gen quiz_rate = quiz_total / 10000
	save simulation_access_cluster_2degreesociallearning, replace




*** Simulation with network structure from Harvard undergraduate, sequential access rollout before 210, various links

	// average conversation links = 3.19
	// average size of conversation networks = 12.60
	
	forvalue k = 1/13 {
	
	clear all
	
	// primitives
	set obs 1000
	gen id = _n
	
	// assigning conversation networks id
	gen connetwork1 = _n
	egen connetwork = cut(connetwork1), at(1(13)1000)
	drop connetwork1
	replace connetwork = connetwork[_n-1] if connetwork == .
	
	// assigning social connection nodes
	forvalue j = 1/1000 {
		gen connected_`j' = 0
		}
		
	forvalue i = 1/1000 {
		local connet = connetwork[`i']
		local connet1 = `connet' + 12
		forvalue j = `connet'/`connet1' {
			gen connected_`j'_r = runiform()
			replace connected_`j' = 1 in `i'	if connected_`j'_r <= (1/13)*`k'
			drop connected_`j'_r
			}
		}
		
	gen totalconnect = 0
	forvalue i = 1/1000 {
		replace totalconnect = totalconnect + connected_`i'
		}

	// learning rates
	gen alpha = 0.560
	gen direct = 0.332
	gen indirect_0 = 0.246
	gen indirect_1 = 0.121
	
	// assigning access: with complementarity (sequential rollout prior to 210)
	gen access_random = id
	gen access_random_1 = rnormal() if id > 210
	egen access_random_11 = rank(access_random_1) if id > 210
	replace access_random = access_random_11 + 210 if id > 210 
	drop access_random_1 access_random_11
	
	forvalue i = 1/1000 {
		gen soclearning_ownaccess_`i' = (access_random <= `i')
		}

	// learning probability, step 1: pure direct learning
	forvalue i = 1/1000 {
		gen quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i'
		gen quiz_random_`i' = runiform()
		gen quiz_yes_`i' = (quiz_random_`i' < quiz_rate_`i')
		drop quiz_rate_`i' quiz_random_`i'
		}
	
	// total number of people in conversation network who know the info from direct learning (with access)
	forvalue n = 1/1000 {
		gen quiz_yesdr_`n' = 0
		replace quiz_yesdr_`n' = quiz_yes_`n' if soclearning_ownaccess_`n' == 1
		
		gen quiz_yesdracd_`n' = 0
		
		forvalue i = 1/1000 {
		
			di "`n'" "-" "`i'"
			local connet = connetwork[`i']
			local connet1 = `connet' + 12
			
			forvalue j = `connet'/`connet1' {
				qui: replace quiz_yesdracd_`n' = quiz_yesdracd_`n' + quiz_yesdr_`n'[`j'] in `i'  if connected_`j' == 1
				}
			}
		}
		
	// learning probability, step 2: first degree social learning
	forvalue i = 1/1000 {
		gen quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i' + (1-((1-indirect_0)*(1-soclearning_ownaccess_`i') + (1-indirect_1)*(soclearning_ownaccess_`i'))^quiz_yesdracd_`i')	if quiz_yesdracd_`i' > 0
		replace quiz_rate_`i' = alpha + direct * soclearning_ownaccess_`i' 	if quiz_yesdracd_`i' == 0
		gen quiz_random_`i' = runiform()
		gen quiz_yessl_`i' = (quiz_random_`i' < quiz_rate_`i')
		replace quiz_yes_`i' = quiz_yessl_`i' if (quiz_yesdracd_`i' > 0 & quiz_yes_`i' == 0)
		drop quiz_rate_`i' quiz_random_`i' soclearning_ownaccess_`i' quiz_yessl_`i'
		}
	
	// actual learning rates
	forvalue i = 1/1000 {
		egen quiz_total_`i' = sum(quiz_yes_`i')
		drop quiz_yes_`i'
		}
	keep quiz_total_*	
	duplicates drop	
	
	// reshape
	gen i = 1
	reshape long quiz_total_, i(i) j(count)	
	drop i
	rename quiz_total_ quiz_total
	gen quiz_rate = quiz_total / 1000
	save simulation_access_cluster_withsociallearning_harvardnetwork_link_`k', replace

	}




*** Merge simulation results

	clear all
	use simulation_access_cluster_withsociallearning, clear
	rename quiz_rate 	quiz_rate_base
	save simulation_compiled, replace

	use simulation_access_cluster_nosociallearning, clear
	rename quiz_rate 	quiz_rate_cl_nosl
	merge 1:1 count using simulation_compiled
	drop _merge
	save simulation_compiled, replace

	use simulation_access_allcluster_withsociallearning, clear
	rename quiz_rate 	quiz_rate_allcl_sl
	merge 1:1 count using simulation_compiled
	drop _merge
	save simulation_compiled, replace
	
	use simulation_access_cluster_withsociallearning_doubledormsize, clear
	rename quiz_rate 	quiz_rate_base_2dorm
	merge 1:1 count using simulation_compiled
	drop _merge
	save simulation_compiled, replace
	
	use simulation_access_cluster_withsociallearning_doublesociallearning, clear
	rename quiz_rate 	quiz_rate_base_2sl
	merge 1:1 count using simulation_compiled
	drop _merge
	save simulation_compiled, replace
	
	use simulation_access_cluster_withsociallearning_doubledirectlearning, clear
	rename quiz_rate 	quiz_rate_base_2dl
	merge 1:1 count using simulation_compiled
	drop _merge
	save simulation_compiled, replace
	
	use simulation_access_cluster_targeteddiffusion, clear
	rename quiz_rate 	quiz_rate_targetted
	merge 1:1 count using simulation_compiled
	drop _merge
	save simulation_compiled, replace

	use simulation_access_cluster_2degreesociallearning, clear
	rename quiz_rate 	quiz_rate_2degree
	merge 1:1 count using simulation_compiled
	drop _merge
	save simulation_compiled, replace

	drop quiz_total
	
	// generate bins
	egen count_bin = cut(count), at(0(80)10000)

	foreach v in quiz_rate_base quiz_rate_cl_nosl quiz_rate_allcl_sl quiz_rate_base_2dorm quiz_rate_base_2sl quiz_rate_base_2dl quiz_rate_targetted quiz_rate_2degree {
		bysort count_bin: egen `v'_b = mean(`v')
		}

	save simulation_compiled, replace

	

*** Figure 4

	twoway 	(line quiz_rate_base_b count_bin, mcolor(navy) lpattern(line) lcolor(navy) lwidth(vthick)) /// 
			(line quiz_rate_cl_nosl_b count_bin, mcolor(gs7) lpattern(dash) lcolor(gs7) lwidth(thick)), /// 
			xline(2100, lpattern(solid) lcolor(black)) ///
			xline(2350, lpattern(dash) lcolor(gs9)) ///
			xline(5770, lpattern(dash) lcolor(gs9)) ///
			xline(7170, lpattern(dash) lcolor(gs9)) ///
			text(1.03 1250 "Pre-treatment" "proportion", size(small)) ///
			text(1.03 3200 "If cost = 0" "(current demand)", size(small)) ///
			text(1.03 5500 "If cost > 0" "(encouragement)", size(small)) ///
			text(1.03 7400 "If cost = 0" "(encouragement)", size(small)) ///
			yscale(r(0.5 1.03)) ylabel(0.5 "50" 0.6 "60" 0.7 "70" 0.8 "80" 0.9 "90" 1 "100", grid) ///
			xlabel(0 "0" 2000 "20" 4000 "40" 6000 "60" 8000 "80" 10000 "100") ///
			ytitle("% students able to answer quiz correctly") ///
			xtitle("Proportion of students actively visit foreign news websites") ///
			legend(order(1 "Baseline" 2 "If no social learning")) ///
			graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
			saving(figure_simulation_social_5, replace)		
	graph export "figure_simulation_social_5.pdf", replace
	
	
	
*** Figure A.21

	twoway 	(line quiz_rate_base_b count_bin, mcolor(navy) lpattern(line) lcolor(navy) lwidth(vthick)) /// 
			(line quiz_rate_base_2dl_b count_bin, mcolor(gs5) lpattern(shortdash dot) lcolor(gs5) lwidth(thick)) /// 
			(line quiz_rate_cl_nosl_b count_bin, mcolor(gs7) lpattern(dash) lcolor(gs7) lwidth(thick)) /// 
			(line quiz_rate_base_2sl_b count_bin, mcolor(gs3) lpattern(longdash dot) lcolor(gs3) lwidth(thick)), /// 
			xline(2100, lpattern(solid) lcolor(black)) ///
			xline(2350, lpattern(dash) lcolor(gs9)) ///
			xline(5770, lpattern(dash) lcolor(gs9)) ///
			xline(7170, lpattern(dash) lcolor(gs9)) ///
			text(1.03 1250 "Pre-treatment" "proportion", size(small)) ///
			text(1.03 3200 "If cost = 0" "(current demand)", size(small)) ///
			text(1.03 5500 "If cost > 0" "(encouragement)", size(small)) ///
			text(1.03 7400 "If cost = 0" "(encouragement)", size(small)) ///
			yscale(r(0.5 1.03)) ylabel(0.5 "50" 0.6 "60" 0.7 "70" 0.8 "80" 0.9 "90" 1 "100", grid) ///
			xlabel(0 "0" 2000 "20" 4000 "40" 6000 "60" 8000 "80" 10000 "100") ///
			ytitle("% students able to answer quiz correctly") ///
			xtitle("Proportion of students actively visit foreign news websites") ///
			legend(order(1 "Baseline" 2 "If direct learning rate doubles" 3 "If no social learning" 4 "If social transmission rate doubles")) ///
			graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
			saving(figure_simulation_learningrates, replace)		
	graph export "figure_simulation_learningrates.pdf", replace
	

	
*** Figure A.22

	twoway 	(line quiz_rate_base_b count_bin, mcolor(navy) lpattern(line) lcolor(navy) lwidth(vthick)) /// 
			(line quiz_rate_cl_nosl_b count_bin, mcolor(gs7) lpattern(dash) lcolor(gs7) lwidth(thick)) /// 
			(line quiz_rate_allcl_sl_b count_bin, mcolor(gs3) lpattern(longdash dot) lcolor(gs3) lwidth(thick)) /// 
			(line quiz_rate_base_2dorm_b count_bin, mcolor(gs5) lpattern(shortdash dot) lcolor(gs5) lwidth(thick)) /// 
			(line quiz_rate_targetted_b count_bin, mcolor(gs10) lpattern(longdash) lcolor(gs10) lwidth(thick)) /// 
			(line quiz_rate_2degree_b count_bin, mcolor(gs10) lpattern(shortdash) lcolor(gs13) lwidth(thick)), /// 
			xline(2100, lpattern(solid) lcolor(black)) ///
			xline(2350, lpattern(dash) lcolor(gs9)) ///
			xline(5770, lpattern(dash) lcolor(gs9)) ///
			xline(7170, lpattern(dash) lcolor(gs9)) ///
			text(1.03 1250 "Pre-treatment" "proportion", size(small)) ///
			text(1.03 3200 "If cost = 0" "(current demand)", size(small)) ///
			text(1.03 5500 "If cost > 0" "(encouragement)", size(small)) ///
			text(1.03 7400 "If cost = 0" "(encouragement)", size(small)) ///
			yscale(r(0.5 1.03)) ylabel(0.5 "50" 0.6 "60" 0.7 "70" 0.8 "80" 0.9 "90" 1 "100", grid) ///
			xlabel(0 "0" 2000 "20" 4000 "40" 6000 "60" 8000 "80" 10000 "100") ///
			ytitle("% students able to answer quiz correctly") ///
			xtitle("Proportion of students actively visit foreign news websites") ///
			legend(order(1 "Baseline" 2 "If no social learning" 3 "If diffuse in clusters" 4 "If dorm size doubles" 5 "If targeted diffusion" 6 "2nd degree transmission")) ///
			graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
			saving(figure_simulation_social, replace)		
	graph export "figure_simulation_social.pdf", replace
		



*** Figure A.23		

	clear all
	use simulation_access_cluster_withsociallearning_harvardnetwork_link_1, clear
	rename quiz_rate 	quiz_rate_harvard_link_1
	save simulation_compiled_harvard_links, replace
	
	forvalue k = 2/13 {
		use simulation_access_cluster_withsociallearning_harvardnetwork_link_`k', clear
		rename quiz_rate 	quiz_rate_harvard_link_`k'
		merge 1:1 count using simulation_compiled_harvard_links
		drop _merge
		save simulation_compiled_harvard_links, replace
		}

	use simulation_access_cluster_withsociallearning, clear
	rename quiz_rate 	quiz_rate_base
	replace count = count / 10
	merge 1:1 count using simulation_compiled_harvard_links
	keep if _merge == 3
	drop _merge	
	save simulation_compiled_harvard_links, replace
		
	// generate bins
	egen count_bin = cut(count), at(0(80)1000)

	forvalue k = 1/13 {
		foreach v in quiz_rate_harvard_link_`k' {
			bysort count_bin: egen `v'_b = mean(`v')
			}
		}
		
	bysort count_bin: egen quiz_rate_base_b = mean(quiz_rate_base)

	twoway 	(line quiz_rate_harvard_link_1_b count_bin, mcolor(gs1) lpattern(dash) lcolor(gs1) lwidth(median)) /// 
			(line quiz_rate_harvard_link_2_b count_bin, mcolor(gs2) lpattern(dash) lcolor(gs2) lwidth(median)) ///
			(line quiz_rate_harvard_link_3_b count_bin, mcolor(orange) lpattern(longdash) lcolor(orange) lwidth(thick)) /// 
			(line quiz_rate_harvard_link_4_b count_bin, mcolor(gs4) lpattern(dash) lcolor(gs4) lwidth(median)) /// 
			(line quiz_rate_harvard_link_5_b count_bin, mcolor(gs5) lpattern(dash) lcolor(gs5) lwidth(median)) /// 
			(line quiz_rate_harvard_link_6_b count_bin, mcolor(gs6) lpattern(dash) lcolor(gs6) lwidth(median)) /// 
			(line quiz_rate_harvard_link_7_b count_bin, mcolor(gs7) lpattern(dash) lcolor(gs7) lwidth(median)) /// 
			(line quiz_rate_harvard_link_8_b count_bin, mcolor(gs8) lpattern(dash) lcolor(gs8) lwidth(median)) /// 
			(line quiz_rate_harvard_link_9_b count_bin, mcolor(gs9) lpattern(dash) lcolor(gs9) lwidth(median)) /// 
			(line quiz_rate_harvard_link_10_b count_bin, mcolor(gs10) lpattern(dash) lcolor(gs10) lwidth(median)) /// 
			(line quiz_rate_harvard_link_11_b count_bin, mcolor(gs11) lpattern(dash) lcolor(gs11) lwidth(median)) /// 
			(line quiz_rate_harvard_link_12_b count_bin, mcolor(gs12) lpattern(dash) lcolor(gs12) lwidth(median)) /// 
			(line quiz_rate_harvard_link_13_b count_bin, mcolor(gs13) lpattern(dash) lcolor(gs13) lwidth(median)) /// 
			(line quiz_rate_base_b count_bin, mcolor(navy) lpattern(line) lcolor(navy) lwidth(thick)), /// 
			xline(210, lpattern(solid) lcolor(black)) ///
			xline(235, lpattern(dash) lcolor(gs9)) ///
			xline(577, lpattern(dash) lcolor(gs9)) ///
			xline(717, lpattern(dash) lcolor(gs9)) ///
			text(1.03 125 "Pre-treatment" "proportion", size(small)) ///
			text(1.03 320 "If cost = 0" "(current demand)", size(small)) ///
			text(1.03 550 "If cost > 0" "(encouragement)", size(small)) ///
			text(1.03 740 "If cost = 0" "(encouragement)", size(small)) ///
			yscale(r(0.5 1.03)) ylabel(0.5 "50" 0.6 "60" 0.7 "70" 0.8 "80" 0.9 "90" 1 "100", grid) ///
			xlabel(0 "0" 200 "20" 400 "40" 600 "60" 800 "80" 1000 "100") ///
			ytitle("% students able to answer quiz correctly") ///
			xtitle("Proportion of students actively visit foreign news websites") ///
			legend(order(14 "Baseline" 1 "1 link" 2 "2 links" 3 "3 links" 4 "4 links" 5 "5 links" 6 "6 links" 7 "7 links" 8 "8 links" 9 "9 links" 10 "10 links" 11 "11 links" 12 "12 links" 13 "13 links") col(4)) ///
			graphregion(fcolor(white) ilcolor(white) lcolor(white))	///
			saving(figure_simulation_harvard_links, replace)		
	graph export "figure_simulation_harvard_links.pdf", replace
		
