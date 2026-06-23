# Replication: Jun and Lee (2022)

## Overview

The files in this folder replicate results from the following paper:

Sung Jae Jun and Sokbae Lee. 2023. Identifying the Effect of Persuasion. <https://doi.org/10.1086/724114>. Available in the _Journal of Political Economy_.

Specifically, this folder provides R codes that match with the paper's accompanying Stata replication package. 

## Replication folder structure 

The replication files have the following structure:

```
replications-JunLee-JPE-R                # R replication of Jun and Lee (2022)
    ├── data                            # Datasets
    │   ├── ChenYang2019                #   Dataset from Chen and Yang (2019)
    │   ├── DellaVignaKaplan2007         #   Dataset from DellaVigna and Kaplan (2007)
    │   ├── DellaVignaListMalmendier2012 #   Dataset from DellaVigna, List, and Malmendier (2012)
    │   ├── EnikolopovPetrovaZhuravskaya2011 # Dataset from Enikolopov, Petrova, and Zhuravskaya (2011)
    │   ├── GerberKarlanBergan2009       #   Dataset from Gerber, Karlan, and Bergan (2009)
    ├── R                               # Replication scripts
    │   ├── renv.lock                   #   Locked package versions (replaces "libraries")
    │   ├── logs                        #   R session/log output saved here
    │   ├── table1.R                    #   Make Table 1
    │   ├── table2.R                    #   Make Table 2
    │   ├── table3.R                    #   Make Table 3
    │   ├── table4.R                    #   Make Table 4
    │   ├── tableD1.R                   #   Make Table D1
    │   ├── figureD1input.R             #   Prepare input for Figure D1
    │   ├── FigureD1.R                  #   Make Figure D1
    │   ├── tableE1.R                   #   Make Table E1
    │   ├── tableE2.R                   #   Make Table E2
    │   ├── tableH1.R                   #   Make Table H1
    │   ├── CY19_data.R                 #   Prepare input for Table 2
    ├── results                         # Replication Figure D1 and all the tables
    └── run.R                           # Master script that calls all replication scripts
```

## Instructions

First, [download this repository](https://github.com/persuasio/replication-JunLee-JPE/archive/main.zip). Then, run the following R script: `run.R`. This master script runs all the scripts required to replicate the main results .

- `R/renv.lock` records the package versions used in the analysis; running `renv::restore()` at the top of run.R will install them.
- To generate the dataset `ChenYang2019.rds`, it is necessary to run `CY19_data.R` after downloading the original dataset from Chen and Yang (2019) at https://www.aeaweb.org/articles?id=10.1257/aer.20171765 and storing it at `data/ChenYang2019`. This dataset is already created and stored at `data/ChenYang2019`. 

Before running `run.R`, the user must define the variable `persuasion_dir` at the top of that script, pointing to the repository's root directory. Also, it is necessary to install the package `persuasio` from its source (included in `R/persuasio.R`) and the stargazer package (`install.packages("stargazer")`) before running the script, as it is used in place of Stata's `outreg` for formatted regression tables.


## Software, memory, and runtime requirements

Running this analysis requires [placeholder] or higher. Add-on packages are included in `[placeholder] and do not need to be installed by the user. In particular, the package `persuasio` is a R module that estimates the effect of persuasion and conducts inference. 

Memory requirements are minimal. Runtime is approximately [placeholder] on an [placeholder] when running the default specification of `nbt = 10000` bootstraps. 

## Acknowledgment

The structure of the replication files is forked from Illinois Workplace Wellness Study: Public Use Data Repository at <https://github.com/reifjulian/illinois-wellness-data>.
