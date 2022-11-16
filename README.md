# Replication: Jun and Lee (2022)

## Overview

This folder contains files to replicate results from the following paper:

Sung Jae Jun and Sokbae Lee. 2022. Identifying the Effect of Persuasion. https://arxiv.org/abs/1812.02276.

## Replication folder structure 

The following diagram summarizes the organization of the replication files.

```
replications-JunLee-JPE         # Replication files for Jun and Lee (2022)
    ├── data                    # 	Datasets
    |   ├── ChenYang2019                        #	Dataset from Chen and Yang (2019)
    |   ├── DellaVignaKaplan2007                #	Dataset from DellaVigna and Kaplan (2007)
    |   ├── DellaVignaListMalmendier2012        #	Dataset from DellaVigna, List, and Malmendier (2012)
    |   ├── EnikolopovPetrovaZhuravskaya2011	#	Dataset from Enikolopov, Petrova, and Zhuravskaya (2011)
    |   ├── GerberKarlanBergan2009              #	Dataset from Gerber, Karlan, and Bergan (2009)
    ├── scripts                 # 	Replication scripts
    |   ├── libraries           #       	Stata module `persuasio`
    |   ├── logs                #       	Stata log files are saved here
    |   ├── table1.do           #       	Make Table 1
    |   ├── table2.do           #       	Make Table 2
    |   ├── table3.do           #       	Make Table 3
    |   ├── table4.do           #       	Make Table 4
    |   ├── tableD1.do          #       	Make Table D1
    |   |── figureD1input.do    #       	Prepare input for Figure D1
    |   ├── FigureD1.R          #       	Make Figure D1
    |   ├── tableE1.do          #       	Make Table E1
    |   ├── tableE2.do          #       	Make Table E2
    |   ├── tableH1.do          #       	Make Table H1
    |   ├── CY19_data.do        #		Prepare input for Table 2
    ├── results                 #	Replication Figure D1 and all the tables
    └── run.do                  #	Master script that calls the replication scripts (except FigureD1.R)
```

## Instructions

First, [download this repository](https://github.com/persuasio/replication-JunLee-JPE/archive/main.zip). Then, run the following Stata script: `run.do`. This master script runs all the scripts required to replicate the main results except Figure D1, which can be generated by running `FigureD1.R`. Other scripts are contained in the `scripts` directory and are organized as follows:

- `scripts/libraries` contains user-built Stata commands used in the analysis.
- To generate the dataset `ChenYang2019.dta`, it is necessary to run `CY19_data.do` after downloading the original dataset from Chen and Yang (2019) at https://www.aeaweb.org/articles?id=10.1257/aer.20171765 and storing them at `data/ChenYang2019`. This dataset is already created and stored at `data/ChenYang2019`. 

Before running `run.do`, the user must define the global macro `Persuasion` at the top of that script. Also, it is necessary to install the Stata module `outreg` by `ssc install outreg, replace` before running the do file.

## Software, memory, and runtime requirements

Running this analysis requires Stata version 16 or higher. Add-on packages are included in `scripts/libraries` and do not need to be installed by the user. In particular, the package `persuasio` is a Stata module that estimates the effect of persuasion and conducts inference. 

Memory requirements are minimal. Runtime is approximately 3 hours on an iMac (M1, 2021) when running the default specification of `nbt = 10000` bootstraps. 

## Acknowledgment

The structure of the replication files is forked from Illinois Workplace Wellness Study: Public Use Data Repository at https://github.com/reifjulian/illinois-wellness-data.
