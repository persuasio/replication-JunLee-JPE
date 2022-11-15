# Replication: Jun and Lee (2022)

## Overview

This folder contains files to replicate results from the following paper:

Sung Jae Jun and Sokbae Lee. 2022. Identifying the Effect of Persuasion. https://arxiv.org/abs/1812.02276.

## Replication folder structure 

The following diagram summarizes the organization of the replication files.

```
replications-JunLee-JPE         #	Replication files for Jun and Lee (2022)
    ├── data                    #   	Datasets
    |   ├── ChenYang2019                        #	Dataset from Chen and Yang (2019)
    |   ├── DellaVignaKaplan2007                #	Dataset from DellaVigna and Kaplan (2007)
    |   ├── DellaVignaListMalmendier2012        #	Dataset from DellaVigna, List, and Malmendier (2012)
    |   ├── EnikolopovPetrovaZhuravskaya2011	#	Dataset from Enikolopov, Petrova, and Zhuravskaya (2011)
    |   ├── GerberKarlanBergan2009              #	Dataset from Gerber, Karlan, and Bergan (2009)
    ├── scripts                 #		Replication scripts
    |   ├── libraries           #       	Add-on Stata packages
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
    |   ├── CY19_data.do        #           Prepare input for Table 2
    ├── results                 #	Replication Figure D1 and all the tables
    └── run.do                  #	Master script that calls the replication scripts (except FigureD1.R)
```



## Acknowledgment



