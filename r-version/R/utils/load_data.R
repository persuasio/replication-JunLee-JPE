# This script loads data files from the replication package and makes them workable in R.

library(haven)

load_data <- function() {
  list(
    cy  = zap_labels(read_dta("data/ChenYang2019/ChenYang2019.dta")),
    dk  = zap_labels(read_dta("data/DellaVignaKaplan2007/DellaVignaKaplan2007.dta")),
    dlm = zap_labels(read_dta("data/DellaVignaListMalmendier2012/DellaVignaListMalmendier2012.dta")),
    epz = zap_labels(read_dta("data/EnikolopovPetrovaZhuravskaya2011/EnikolopovPetrovaZhuravskaya2011.dta")),
    gkb = zap_labels(read_dta("data/GerberKarlanBergan2009/GerberKarlanBergan2009.dta"))
  )
}
