**********************
* OVERVIEW
*   This script generates selected tables for the paper:
*       "IDENTIFYING THE EFFECT OF PERSUASION" (Sung Jae Jun and Sokbae Lee)
*   All data are stored in /data
*   All results are outputted to /results
* 
* SOFTWARE REQUIREMENTS
*   Stata version 16 or newer
*
* TO PERFORM A CLEAN RUN, DELETE THE FOLLOWING FOLDER:
*   $Persuasion/results
*
* ACKNOWLEDGMENT
*	The structure of the replication files is forked from 
*   Illinois Workplace Wellness Study: Public Use Data Repository
*   at https://github.com/reifjulian/illinois-wellness-data
**********************

**********************
* Parameters defined by user
**********************

* "Persuasion" points to the root directory, which contains the subfolders "data" and "scripts"
global Persuasion "TO BE ADDED"

* Number of bootstraps for Table 2
* The paper used nbt = 10000. Runtime is approximately 3 hours on an iMac (M1, 2021).
* Use nbt = 10 to check whether replication code runs without an error. 
global nbt = 10

* set seed for replicability 
set seed 987975	

**********************
**********************

* Confirm that the global for the root directory has been defined
assert !missing("$Persuasion")

* Initialize log and record system parameters
clear 
set more off
cap mkdir "$Persuasion/scripts/logs"
cap log close
local datetime : di %tcCCYY.NN.DD!-HH.MM.SS `=clock("$S_DATE $S_TIME", "DMYhms")'
local logfile "$Persuasion/scripts/logs/`datetime'.log.txt"
log using "`logfile'", text

di "Begin date and time: $S_DATE $S_TIME"
di "Stata version: `c(stata_version)'"
di "Updated as of: `c(born_date)'"
di "Variant:       `=cond( c(MP),"MP",cond(c(SE),"SE",c(flavor)) )'"
di "Processors:    `c(processors)'"
di "OS:            `c(os)' `c(osdtl)'"
di "Machine type:  `c(machine_type)'"

* All required Stata packages are available in the /libraries folder
adopath ++ "$Persuasion/scripts/libraries"
mata: mata mlib index

* Stata version control
version 16

* Create directories for output files
cap mkdir "$Persuasion/results"

**********************
* Run analysis
**********************
* TABLE 1. Persuasion Rates: Papers on Voter Turnout
do "$Persuasion/scripts/table1.do"

* Effects of Uncensored Media: Revisiting Chen and Yang (2019)
*   To generate the dataset "ChenYang2019.dta", it is necessary to run CY19_data.do 
*   after downloading the original dataset from 
*   the AER webpage at https://www.aeaweb.org/articles?id=10.1257/aer.20171765
*   and storing them at "$Persuasion/data/ChenYang2019".
*   As the dataset is already created and stored at "$Persuasion/data/ChenYang2019", 
*   the do file is commented below.
*do "$Persuasion/scripts/CY19_data.do"
*   TABLE 2. Persuasion Rates of Exposure to Uncensored Internet
do "$Persuasion/scripts/table2.do"

* Effects of Political News: Revisting Gerber, Karlan, and Bergan (2009) 
*   TABLE 3. Summary Statistics of the GKB Data
do "$Persuasion/scripts/table3.do"
*   TABLE 4. Estimates of the Key Parameters
do "$Persuasion/scripts/table4.do"

*	TABLE D1. Persuasion Rates: Fox News Effects
do "$Persuasion/scripts/tableD1.do"

*	FIGURE D1. Estimates of Marginal and Average Persuasion Rates
* 	The following do file generates an excel file that will be used to 
* 	create Figure D1 by separate R code.
do "$Persuasion/scripts/figureD1input.do"

*	TABLE E1. Persuasive Effect by Treatment in Landry et al. (2006)
do "$Persuasion/scripts/tableE1.do"
*	TABLE E2. Persuasive Effect by Treatment in DLM
do "$Persuasion/scripts/tableE2.do"

*	TABLE H1. Persuasion Rates: NTV Effects Using a Binary Instrument
do "$Persuasion/scripts/tableH1.do"

* End log
di "End date and time: $S_DATE $S_TIME"
log close

** EOF
