# ******************************************************************************
# SCRIPT:   CY19_data.R PURPOSE:  Creates a dataset for Table 2
#
# ACKNOWLEDGMENT: This script is a slightly modified version of the do file
# "master_dofile_panelsurvey.do" included in the replication files for Chen and
# Yang (2019, AER). Before running this code, it is necessary to download the
# original dataset from the AER webpage at
# https://www.aeaweb.org/articles?id=10.1257/aer.20171765 and store them at
# "$Persuasion/data/ChenYang2019"

# R file: replication of Chen and Yang (2019)
# Part 1: based on student survey data

library(haven)
library(dplyr)
library(tidyr)

data_path <- here::here("data/ChenYang2019")

#------------------------------------------------------------------------------
# Set directory for output
#------------------------------------------------------------------------------

output_dir <- here::here("output")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# HELPER: Anderson (2008) GLS z-score -------------------------------------


andersonz <- function(data, varlist, subset = NULL, generate) {

  n <- nrow(data)

  # `if_mask` (used where Stata used `if`) and
  # `full_mask` (used where Stata omitted it, i.e. always all rows).
  if (is.null(subset)) {
    if_mask <- rep(TRUE, n)
  } else {
    if_mask <- subset
  }
  full_mask <- rep(TRUE, n)

  # Convert to effect sizes

  varliste <- character(length(varlist))

  for (i in seq_along(varlist)) {
    X <- varlist[i]
    Xe <- paste0(X, "e")

    x_vals <- data[[X]]

    # qui sum `X' `if'  -> sd computed on the `if' subset
    csd <- sd(x_vals[if_mask], na.rm = TRUE)

    Xe_vals <- rep(NA_real_, n)
    Xe_vals[if_mask] <- x_vals[if_mask] / csd

    # qui sum `X'   -> mean is over FULL sample
    full_mean <- mean(x_vals[full_mask], na.rm = TRUE)

    # replace `X'e = `X'e - r(mean)/`csd' `if'
    Xe_vals[if_mask] <- Xe_vals[if_mask] - full_mean / csd

    data[[Xe]] <- Xe_vals
    varliste[i] <- Xe
  }

  # Generate GLS weighting matrix
  Xe_mat_full <- as.matrix(data[, varliste, drop = FALSE])
  complete_rows <- stats::complete.cases(Xe_mat_full)
  Xe_complete <- Xe_mat_full[complete_rows, , drop = FALSE]
  Ncases <- nrow(Xe_complete)

  # deviations from column means, cross product, divided by N
  centered <- scale(Xe_complete, center = TRUE, scale = FALSE)
  R <- (t(centered) %*% centered) / Ncases

  # matrix R = syminv(R)  -> symmetric (generalized) inverse using solve()
  R_inv <- solve(R)  # R is symmetric PD in the standard Anderson setup;

  # Calculate weights
  k <- ncol(R_inv)
  J <- matrix(1, nrow = k, ncol = 1)
  weights <- as.numeric(R_inv %*% J)  # weight[i] = sum_j R_inv[i, j]
  names(weights) <- varlist

  # Generate total number of variables per obs
  # Consistent with the Stata code, as opposed to Anderson's code, do not replace missing with zeros

  sample_wt <- rep(0, n)
  for (i in seq_along(varlist)) {
    X <- varlist[i]
    nonmissing_X <- !is.na(data[[X]])
    sample_wt[nonmissing_X] <- sample_wt[nonmissing_X] + weights[i]
  }
  # Apply weights to outcomes

  out <- rep(NA_real_, n)
  out[if_mask] <- 0

  for (i in seq_along(varlist)) {
    Xe <- varliste[i]
    out[if_mask] <- data[[Xe]][if_mask] * weights[i] + out[if_mask]
  }

  out[if_mask] <- out[if_mask] / sample_wt[if_mask]

  # Drop the helper columns
  data[varliste] <- NULL

  # Re-standardize aggregate z-score
  z_out <- as.numeric(scale(out))  # scale() uses sd() with N-1 denominator,
  # NA-handling: by default scale() does not drop NAs internally, so compute explicitly below instead.

  out_mean <- mean(out, na.rm = TRUE)
  out_sd   <- sd(out, na.rm = TRUE)
  z_out <- (out - out_mean) / out_sd

  data[[generate]] <- z_out

  data
}

# HELPER: Lee (2009) attrition bounds --------------------------------------

trimbound <- function(data, y, d, wgt, sv = NULL) {

  # Stata: preserve
  # R: work on a local copy.
  dat <- data %>%
    filter(!is.na(.data[[d]])) %>%
    mutate(
      s = as.integer(!is.na(.data[[y]]))
    )

  if (is.null(sv)) {
    # Stata allows a sortvar to break ties. If none is supplied,
    # use row order as the tie-breaker.
    dat <- dat %>% mutate(.sortvar = row_number())
    sv <- ".sortvar"
  }

  wmean <- function(x, w) {
    weighted.mean(x, w, na.rm = TRUE)
  }

  wvar <- function(x, w) {
    keep <- !is.na(x) & !is.na(w)
    x <- x[keep]
    w <- w[keep]
    mu <- weighted.mean(x, w)
    sum(w * (x - mu)^2) / sum(w)
  }

  # Stata: reg s d [aw=wgt]
  sel_fit <- lm(
    reformulate(d, response = "s"),
    data = dat,
    weights = dat[[wgt]]
  )

  b_d <- coef(sel_fit)[[d]]

  # If treatment has lower selection, flip treatment status exactly
  # as in the Stata code.
  trim_control <- b_d < 0

  if (trim_control) {
    dat[[d]] <- 1 - dat[[d]]

    sel_fit <- lm(
      reformulate(d, response = "s"),
      data = dat,
      weights = dat[[wgt]]
    )
  }

  alpha <- coef(sel_fit)[["(Intercept)"]]
  b_d <- coef(sel_fit)[[d]]
  p <- b_d / (b_d + alpha)

  sumw <- dat %>%
    filter(.data[[d]] == 1, !is.na(.data[[y]])) %>%
    summarise(x = sum(.data[[wgt]], na.rm = TRUE)) %>%
    pull(x)

  dat <- dat %>%
    arrange(.data[[d]], .data[[y]], .data[[sv]]) %>%
    group_by(.data[[d]]) %>%
    mutate(
      swgt = cumsum(.data[[wgt]]),
      upper = as.integer(swgt > p * sumw),
      lower = as.integer(swgt < (1 - p) * sumw)
    ) %>%
    ungroup()

  controlN <- dat %>% filter(.data[[d]] == 0) %>% nrow()
  treatN   <- dat %>% filter(.data[[d]] == 1) %>% nrow()

  control <- dat %>% filter(.data[[d]] == 0)
  treat   <- dat %>% filter(.data[[d]] == 1)

  yc <- wmean(control[[y]], control[[wgt]])
  yt <- wmean(treat[[y]], treat[[wgt]])

  ycvar <- wvar(control[[y]], control[[wgt]]) / sum(!is.na(control[[y]]))
  ytvar <- wvar(treat[[y]], treat[[wgt]]) / sum(!is.na(treat[[y]]))

  UtreatN <- sum(!is.na(treat[[y]]))

  upper_treat <- dat %>% filter(.data[[d]] == 1, upper == 1)
  lower_treat <- dat %>% filter(.data[[d]] == 1, lower == 1)

  uyt <- wmean(upper_treat[[y]], upper_treat[[wgt]])
  lyt <- wmean(lower_treat[[y]], lower_treat[[wgt]])

  alphat <- alpha / (1 - p)

  pvar <- ((1 - p)^2) *
    (((1 - alphat) / (treatN * alphat)) +
       ((1 - alpha) / (controlN * alpha)))

  ucomp1 <- wvar(upper_treat[[y]], upper_treat[[wgt]]) / sum(!is.na(upper_treat[[y]]))
  ucomp2 <- ((min(upper_treat[[y]], na.rm = TRUE) - uyt)^2) *
    p / ((1 - p) * UtreatN)
  ucomp3 <- (((min(upper_treat[[y]], na.rm = TRUE) - uyt)^2) /
               ((1 - p)^2)) * pvar
  ubse <- sqrt(ucomp1 + ucomp2 + ucomp3)

  lcomp1 <- wvar(lower_treat[[y]], lower_treat[[wgt]]) / sum(!is.na(lower_treat[[y]]))
  lcomp2 <- ((max(lower_treat[[y]], na.rm = TRUE) - lyt)^2) *
    p / ((1 - p) * UtreatN)
  lcomp3 <- (((max(lower_treat[[y]], na.rm = TRUE) - lyt)^2) /
               ((1 - p)^2)) * pvar
  lbse <- sqrt(lcomp1 + lcomp2 + lcomp3)

  if (!trim_control) {

    tibble(
      outcome = y,
      trimmed_group = "treatment",
      p = p,
      p_se = sqrt(pvar),
      lower_cutoff = min(upper_treat[[y]], na.rm = TRUE),
      upper_cutoff = max(lower_treat[[y]], na.rm = TRUE),

      control_mean = yc,
      control_se = sqrt(ycvar),
      treatment_mean = yt,
      treatment_se = sqrt(ytvar),
      untrimmed_effect = yt - yc,
      untrimmed_effect_se = sqrt(ytvar + ycvar),

      upper_bound_mean = uyt,
      upper_bound_se = ubse,
      lower_bound_mean = lyt,
      lower_bound_se = lbse,

      upper_bound_effect = uyt - yc,
      upper_bound_effect_se = sqrt(ycvar + ubse^2),
      lower_bound_effect = lyt - yc,
      lower_bound_effect_se = sqrt(ycvar + lbse^2)
    )

  } else {

    tibble(
      outcome = y,
      trimmed_group = "control",
      p = p,
      p_se = sqrt(pvar),
      lower_cutoff = min(upper_treat[[y]], na.rm = TRUE),
      upper_cutoff = max(lower_treat[[y]], na.rm = TRUE),

      control_mean = yt,
      control_se = sqrt(ytvar),
      treatment_mean = yc,
      treatment_se = sqrt(ycvar),
      untrimmed_effect = yc - yt,
      untrimmed_effect_se = sqrt(ytvar + ycvar),

      upper_bound_mean = uyt,
      upper_bound_se = ubse,
      lower_bound_mean = lyt,
      lower_bound_se = lbse,

      upper_bound_effect = yc - lyt,
      upper_bound_effect_se = sqrt(ycvar + lbse^2),
      lower_bound_effect = yc - uyt,
      lower_bound_effect_se = sqrt(ycvar + ubse^2)
    )
  }
}

# Load Data ---------------------------------------------------------------

library(haven)
library(dplyr)
library(tidyr)
library(stringr)

cy <- zap_labels(read_dta(file.path(data_path, "panelsurvey_raw.dta")))

# Define Variables --------------------------------------------------------
# The code keeps variable names identical to the original Stata code.

var_identity_w1                <- c("responseID_wave1")
var_university                 <- c("university")
var_treatment                  <- c("treatment_vpn", "treatment_newsletter", "vpn_current_paid_user", "treatment_vpnexpiration")


# 	A: Media related beliefs, attitudes, and behaviors
var_info_ranking_reg_w1         <- c("info_domestic_website_w1", "info_foreign_website_w1", "info_social_media_dom_w1", "info_word_of_mouth_w1")
var_info_ranking_reg_w2         <- c("info_domestic_website_w2", "info_foreign_website_w2", "info_social_media_dom_w2", "info_social_media_for_w2", "info_word_of_mouth_w2")
var_info_ranking_reg_w3         <- c("info_domestic_website_w3", "info_foreign_website_w3", "info_social_media_dom_w3", "info_social_media_for_w3", "info_word_of_mouth_w3")
var_info_ranking_fig_w23        <- c("info_domestic_website", "info_foreign_website", "info_social_media_dom", "info_word_of_mouth", "info_social_media_for")

var_info_freq_reg_w1            <- c("info_freq_website_for_w1")
var_info_freq_reg_w2            <- c("info_freq_website_for_w2")
var_info_freq_reg_w3            <- c("info_freq_website_for_w3")
var_info_freq_fig_w123          <- c("info_freq_website_for")

var_vpnusage_own                <- c("vpn_usage", "vpn_neveruser", "vpn_payment_now", "vpn_stop_when", "vpn_stop_why", "vpn_stop_why_text")

var_vpnusage_rm_all             <- c("vpn_roommate_usage")

var_vpn_purchase_w3             <- c("vpn_purchase", "vpn_purchase_yes", "vpn_purchase_wmt", "vpn_purchase_premium")
var_vpn_purchase_reg_w3         <- c("vpn_purchase_wmt_record", "vpn_purchase_yes")

var_media_valuation_reg_w1      <- c("wtp_vpn_w1", "added_value_foreign_media_w1")
var_media_valuation_reg_w2      <- c("wtp_vpn_w2", "added_value_foreign_media_w2")
var_media_valuation_reg_w3      <- c("wtp_vpn_w3", "added_value_foreign_media_w3")
var_media_valuation_fig_w123    <- c("wtp_vpn", "added_value_foreign_media", "az_belief_media_value")

var_media_trust_reg_w1          <- c("trust_media_dom_state_w1", "trust_media_dom_private_w1", "trust_media_foreign_w1")
var_media_trust_reg_w2          <- c("trust_media_dom_state_w2", "trust_media_dom_private_w2", "trust_media_foreign_w2")
var_media_trust_reg_w3          <- c("trust_media_dom_state_w3", "trust_media_dom_private_w3", "trust_media_foreign_w3")
var_media_trust_fig_w123        <- c("trust_media_dom_state", "trust_media_dom_private", "trust_media_foreign", "az_belief_media_trust")

var_percmediabias_reg_di_cn_w2  <- c("distneutral_cn_neg_cn_w2", "distneutral_cn_pos_cn_w2", "distneutral_us_neg_cn_w2", "distneutral_us_pos_cn_w2")
var_percmediabias_reg_di_us_w2  <- c("distneutral_cn_neg_us_w2", "distneutral_cn_pos_us_w2", "distneutral_us_neg_us_w2", "distneutral_us_pos_us_w2")
var_percmediabias_reg_ce_cn_w2  <- c("bias_cn_neg_cn_cens_w2", "bias_cn_pos_cn_cens_w2", "bias_us_neg_cn_cens_w2", "bias_us_pos_cn_cens_w2")
var_percmediabias_reg_ce_us_w2  <- c("bias_cn_neg_us_cens_w2", "bias_cn_pos_us_cens_w2", "bias_us_neg_us_cens_w2", "bias_us_pos_us_cens_w2")
var_percmediabias_fig_w12       <- c("distneutral_cn_neg_cn", "distneutral_cn_pos_cn", "distneutral_us_neg_cn", "distneutral_us_pos_cn",
                                     "distneutral_cn_neg_us", "distneutral_cn_pos_us", "distneutral_us_neg_us", "distneutral_us_pos_us",
                                     "bias_cn_neg_cn_cens", "bias_cn_pos_cn_cens", "bias_us_neg_cn_cens", "bias_us_pos_cn_cens",
                                     "bias_cn_neg_us_cens", "bias_cn_pos_us_cens", "bias_us_neg_us_cens", "bias_us_pos_us_cens")

var_censor_justif_reg_w1        <- c("censor_just_dom_economic_w1", "censor_just_dom_political_w1", "censor_just_dom_social_w1", "censor_just_for_w1")
var_censor_justif_reg_w2        <- c("censor_just_dom_economic_w2", "censor_just_dom_political_w2", "censor_just_dom_social_w2", "censor_just_for_w2", "censor_just_porn_w2")
var_censor_justif_reg_w3        <- c("censor_just_dom_economic_w3", "censor_just_dom_political_w3", "censor_just_dom_social_w3", "censor_just_for_w3", "censor_just_porn_w3")
var_censor_justif_reg_fp_w1     <- c("censor_just_dom_economic_w1", "censor_just_dom_political_w1", "censor_just_dom_social_w1", "censor_just_for_w1")
var_censor_justif_reg_fp_w2     <- c("censor_just_dom_economic_w2", "censor_just_dom_political_w2", "censor_just_dom_social_w2", "censor_just_for_w2")
var_censor_justif_reg_fp_w3     <- c("censor_just_dom_economic_w3", "censor_just_dom_political_w3", "censor_just_dom_social_w3", "censor_just_for_w3")
var_censor_justif_fig_w123      <- c("censor_just_dom_economic", "censor_just_dom_political", "censor_just_dom_social", "censor_just_for", "az_belief_media_justif_pf")
var_censor_justif_fig_w23       <- c("censor_just_porn")

var_censor_level_reg_w1         <- c("bias_domestic_w1", "bias_foreign_w1")
var_censor_level_reg_w2         <- c("bias_domestic_w2", "bias_foreign_w2")
var_censor_level_reg_w3         <- c("bias_domestic_w3", "bias_foreign_w3")
var_censor_driver_reg_dom_w1    <- c("bias_dom_govt_policy_t1_w1", "bias_dom_firm_interest_t1_w1", "bias_dom_media_pref_t1_w1", "bias_dom_reader_demand_t1_w1")
var_censor_driver_reg_dom_w2    <- c("bias_dom_govt_policy_t1_w2", "bias_dom_firm_interest_t1_w2", "bias_dom_media_pref_t1_w2", "bias_dom_reader_demand_t1_w2")
var_censor_driver_reg_dom_w3    <- c("bias_dom_govt_policy_t1_w3", "bias_dom_firm_interest_t1_w3", "bias_dom_media_pref_t1_w3", "bias_dom_reader_demand_t1_w3")
var_censor_driver_reg_for_w1    <- c("bias_for_govt_policy_t1_w1", "bias_for_firm_interest_t1_w1", "bias_for_media_pref_t1_w1", "bias_for_reader_demand_t1_w1")
var_censor_driver_reg_for_w2    <- c("bias_for_govt_policy_t1_w2", "bias_for_firm_interest_t1_w2", "bias_for_media_pref_t1_w2", "bias_for_reader_demand_t1_w2")
var_censor_driver_reg_for_w3    <- c("bias_for_govt_policy_t1_w3", "bias_for_firm_interest_t1_w3", "bias_for_media_pref_t1_w3", "bias_for_reader_demand_t1_w3")
var_censor_level_fig_w123       <- c("bias_domestic", "bias_foreign")
var_censor_driver_fig_w123      <- c("bias_dom_govt_policy_t1", "bias_dom_firm_interest_t1", "bias_dom_media_pref_t1", "bias_dom_reader_demand_t1",
                                     "bias_for_govt_policy_t1", "bias_for_firm_interest_t1", "bias_for_media_pref_t1", "bias_for_reader_demand_t1")


# 	B: Knowledge
var_knowledge_news_reg_w1       <- c("news_c_stock_rise", "news_c_exchange_depreciate", "news_c_train_brazil_peru", "news_c_army_reduce",
                                     "news_c_taiwan_election", "news_c_china_us_network_coop", "news_c_lijiacheng_china", "news_c_nanjing_memory_list")
var_knowledge_news_reg_w2       <- c("news_c_panamapapers", "news_c_stockcrash", "news_c_topinequality", "news_c_cpicensorship",
                                     "news_c_laborunrest", "news_c_tenyearshk", "news_c_economistcensor", "news_c_waterpollution",
                                     "news_c_applefbi", "news_c_taiwanelection", "news_c_yihehotel", "news_perccorrect_w2",
                                     "news_perccor_cen_w2", "news_perccor_unc_w2", "news_perccor_qui")
var_knowledge_news_reg_quiz      <- c("news_c_topinequality", "news_c_cpicensorship", "news_c_laborunrest", "news_c_waterpollution")
var_knowledge_news_reg_cen_w1    <- c("news_c_stock_rise", "news_c_train_brazil_peru", "news_c_taiwan_election", "news_c_china_us_network_coop")
var_knowledge_news_reg_cen_w2    <- c("news_c_panamapapers", "news_c_tenyearshk", "news_c_stockcrash", "news_c_economistcensor")
var_knowledge_news_reg_cen_w3    <- c("news_c_coalprod", "news_c_trumpchina", "news_c_xiaojianhua", "news_c_xijiangcar",
                                      "news_c_chinanorway", "news_c_womenrights", "news_c_hkceelection")
var_knowledge_news_reg_unc_w1    <- c("news_c_exchange_depreciate", "news_c_army_reduce", "news_c_lijiacheng_china", "news_c_nanjing_memory_list")
var_knowledge_news_reg_unc_w2    <- c("news_c_applefbi", "news_c_taiwanelection", "news_c_yihehotel")
var_knowledge_news_reg_unc_w3    <- c("news_c_northkoreacoal", "news_c_birdflu", "news_c_ethiopiatrain", "news_c_foreignreserve")
var_knowledge_news_fig_w123      <- c("news_perccor_cen", "news_perccor_unc")

var_knowledge_people_reg_tocens  <- c("people_puzhiqiang_w2", "people_renzhiqiang_w2", "people_huangzhifeng_w2")
var_knowledge_people_reg_censor  <- c("people_lizehou_w2", "people_chenguangcheng_w2", "people_lixiaolin_w2")
var_knowledge_people_reg_uncens  <- c("people_maoyushi_w2", "people_honghuang_w2", "people_liuqiangdong_w2")
var_knowledge_people_reg_fake    <- c("people_jialequn_w2")
var_knowledge_people_censor_w1   <- c("people_lizehou_w1", "people_huangzhifeng_w1", "people_chenguangcheng_w1", "people_lixiaolin_w1")
var_knowledge_people_uncens_w1   <- c("people_puzhiqiang_w1", "people_renzhiqiang_w1", "people_maoyushi_w1", "people_honghuang_w1", "people_liuqiangdong_w1")
var_knowledge_people_fig_w12     <- c("people_puzhiqiang", "people_lizehou", "people_huangzhifeng", "people_chenguangcheng", "people_lixiaolin",
                                      "people_renzhiqiang", "people_maoyushi", "people_honghuang", "people_liuqiangdong", "people_jialequn")

var_knowledge_prot_reg_w2        <- c("protest_2012_hk_curriculum_w2", "protest_2014_umbrella_w2", "protest_2016_mongkok_riot_w2",
                                      "protest_2014_sun_flower_w2", "protest_2014_europe_square_w2", "protest_2010_arabic_spring_w2",
                                      "protest_2014_crimea_vote_w2", "protest_2010_catal_indep_w2", "protest_2011_tmrw_parade_w2",
                                      "protest_pcheard_total_w2", "protest_pcheard_china_w2", "protest_pcheard_foreign_w2")
var_knowledge_prot_reg_w3        <- c("protest_2012_hk_curriculum_w3", "protest_2014_umbrella_w3", "protest_2016_mongkok_riot_w3",
                                      "protest_2014_sun_flower_w3", "protest_2014_europe_square_w3", "protest_2010_arabic_spring_w3",
                                      "protest_2014_crimea_vote_w3", "protest_2010_catal_indep_w3", "protest_2011_tmrw_parade_w3",
                                      "protest_2017_women_march_w3", "protest_pcheard_total_w3", "protest_pcheard_china_w3", "protest_pcheard_foreign_w3")
var_knowledge_prot_reg_chi_w1    <- c("protest_2012_hk_curriculum_w1", "protest_2014_umbrella_w1", "protest_2014_sun_flower_w1")
var_knowledge_prot_reg_chi_w2    <- c("protest_2012_hk_curriculum_w2", "protest_2014_umbrella_w2", "protest_2016_mongkok_riot_w2", "protest_2014_sun_flower_w2")
var_knowledge_prot_reg_chi_w3    <- c("protest_2012_hk_curriculum_w3", "protest_2014_umbrella_w3", "protest_2016_mongkok_riot_w3", "protest_2014_sun_flower_w3")
var_knowledge_prot_reg_for_w1    <- c("protest_2014_europe_square_w1", "protest_2010_arabic_spring_w1", "protest_2014_crimea_vote_w1", "protest_2010_catal_indep_w1")
var_knowledge_prot_reg_for_w2    <- c("protest_2014_europe_square_w2", "protest_2010_arabic_spring_w2", "protest_2014_crimea_vote_w2", "protest_2010_catal_indep_w2")
var_knowledge_prot_reg_for_w3    <- c("protest_2014_europe_square_w3", "protest_2010_arabic_spring_w3", "protest_2014_crimea_vote_w3",
                                      "protest_2010_catal_indep_w3", "protest_2017_women_march_w3")
var_knowledge_prot_reg_fak_w1    <- c("protest_2011_tmrw_parade_w1")
var_knowledge_prot_reg_fak_w2    <- c("protest_2011_tmrw_parade_w2")
var_knowledge_prot_reg_fak_w3    <- c("protest_2011_tmrw_parade_w3")
var_knowledge_prot_fig_w123      <- c("protest_2014_europe_square", "protest_2014_sun_flower", "protest_2010_arabic_spring", "protest_2014_crimea_vote",
                                      "protest_2012_hk_curriculum", "protest_2010_catal_indep", "protest_2014_umbrella", "protest_2011_tmrw_parade",
                                      "protest_pcheard_china", "protest_pcheard_foreign")
var_knowledge_prot_fig_w23       <- c("protest_2016_mongkok_riot")
var_knowledge_prot_fig_w3        <- c("protest_2017_women_march")

var_knowledge_time_w2            <- c("news_time_firstclick", "news_time_lastclick", "news_time_submit", "news_time_totalclick",
                                      "pppr_time_firstclick", "pppr_time_lastclick", "pppr_time_submit", "pppr_time_totalclick")

var_knowledge_past_w1            <- c("event_1994_beijing_jianguomen", "event_1994_kelamayi_fire", "event_2003_sunzhigang",
                                      "event_2008_guizhou_wengan", "event_2008_sanlu_milk_powder", "event_2010_hefei_wetland",
                                      "event_2011_guangdong_wukan", "event_2013_fudan_poison", "event_2014_caixin_xuxiao",
                                      "event_2015_sichuan_linshui_road")

var_knowledge_meta_reg_w1        <- c("familiar_china_issues_self_w1", "familiar_china_others_w1")
var_knowledge_meta_reg_w2        <- c("familiar_china_issues_self_w2", "familiar_china_others_w2")
var_knowledge_meta_reg_w3        <- c("familiar_china_issues_self_w3", "familiar_china_others_w3")
var_knowledge_meta_fig_w123      <- c("familiar_china_issues_self", "familiar_china_others", "az_knowledge_meta")


# 	C: Economic beliefs
var_econ_guess_reg_cn_perf_w1   <- c("guess_gdp_growth_china_w1", "guess_stock_index_sh_w1")
var_econ_guess_reg_cn_perf_w2   <- c("guess_gdp_growth_china_w2", "guess_stock_index_sh_w2")
var_econ_guess_reg_cn_perf_w3   <- c("guess_gdp_growth_china_w3", "guess_stock_index_sh_w3")
var_econ_guess_reg_cn_conf_w1   <- c("guess_gdp_growth_china_con_w1", "guess_stock_index_sh_con_w1")
var_econ_guess_reg_cn_conf_w2   <- c("guess_gdp_growth_china_con_w2", "guess_stock_index_sh_con_w2")
var_econ_guess_reg_cn_conf_w3   <- c("guess_gdp_growth_china_con_w3", "guess_stock_index_sh_con_w3")
var_econ_guess_reg_us_perf_w2   <- c("guess_gdp_growth_us_w2", "guess_stock_index_dj_w2")
var_econ_guess_reg_us_perf_w3   <- c("guess_gdp_growth_us_w3", "guess_stock_index_dj_w3")
var_econ_guess_reg_us_conf_w2   <- c("guess_gdp_growth_us_con_w2", "guess_stock_index_dj_con_w2")
var_econ_guess_reg_us_conf_w3   <- c("guess_gdp_growth_us_con_w3", "guess_stock_index_dj_con_w3")
var_econ_guess_fig_w123         <- c("guess_gdp_growth_china", "guess_gdp_growth_china_con", "guess_stock_index_sh", "guess_stock_index_sh_con",
                                     "az_belief_econ_perf_cn", "az_belief_econ_conf_cn")
var_econ_guess_fig_w23          <- c("guess_gdp_growth_us", "guess_gdp_growth_us_con", "guess_stock_index_dj", "guess_stock_index_dj_con",
                                     "az_belief_econ_perf_us", "az_belief_econ_conf_us")


# 	D: Political attitudes
var_demand_change_reg_w1        <- c("inst_change_econ_w1", "inst_change_poli_w1")
var_demand_change_reg_w2        <- c("inst_change_econ_w2", "inst_change_poli_w2")
var_demand_change_reg_w3        <- c("inst_change_econ_w3", "inst_change_poli_w3")
var_demand_change_fig_w123      <- c("inst_change_econ", "inst_change_poli", "az_belief_instchange")

var_trust_inst_reg_govt_w1      <- c("trust_central_govt_w1", "trust_provincial_govt_w1", "trust_local_govt_w1")
var_trust_inst_reg_govt_w2      <- c("trust_central_govt_w2", "trust_provincial_govt_w2", "trust_local_govt_w2")
var_trust_inst_reg_govt_w3      <- c("trust_central_govt_w3", "trust_provincial_govt_w3", "trust_local_govt_w3")
var_trust_inst_reg_foreign_w1   <- c("trust_japan_govt_w1", "trust_us_govt_w1")
var_trust_inst_reg_foreign_w2   <- c("trust_japan_govt_w2", "trust_us_govt_w2")
var_trust_inst_reg_foreign_w3   <- c("trust_japan_govt_w3", "trust_us_govt_w3")
var_trust_inst_reg_finance_w2   <- c("trust_financial_domestic_w2", "trust_financial_foreign_w2")
var_trust_inst_reg_ngo_w2       <- c("trust_ngo_w2")
var_trust_inst_reg_copo_w2      <- c("trust_court_w2", "trust_police_w2")
var_trust_inst_fig_w123         <- c("trust_central_govt", "trust_provincial_govt", "trust_local_govt", "trust_japan_govt", "trust_us_govt",
                                     "az_belief_trust_govt", "az_belief_trust_foreign")
var_trust_inst_fig_w12          <- c("trust_court", "trust_police", "trust_ngo", "trust_financial_domestic", "trust_financial_foreign")

var_eval_govt_reg_w1            <- c("eval_govt_economic_w1", "eval_govt_dom_politics_w1", "eval_govt_for_relations_w1")
var_eval_govt_reg_w2            <- c("eval_govt_economic_w2", "eval_govt_dom_politics_w2", "eval_govt_for_relations_w2")
var_eval_govt_reg_w3            <- c("eval_govt_economic_w3", "eval_govt_dom_politics_w3", "eval_govt_for_relations_w3")
var_eval_govt_fig_w123          <- c("eval_govt_economic", "eval_govt_dom_politics", "eval_govt_for_relations", "az_belief_evalgovt")

var_eval_criteria_reg_w2        <- c("revalgovt_election_w2", "revalgovt_economy_w2", "revalgovt_equality_w2", "revalgovt_ruleoflaw_w2",
                                     "revalgovt_human_rights_w2", "revalgovt_freedom_speech_w2", "revalgovt_global_power_w2", "revalgovt_fair_history_w2")
var_eval_criteria_fig_w12       <- c("revalgovt_election", "revalgovt_economy", "revalgovt_equality", "revalgovt_ruleoflaw",
                                     "revalgovt_human_rights", "revalgovt_freedom_speech", "revalgovt_global_power", "revalgovt_fair_history")

var_severity_reg_w2             <- c("severity_welfare_w2", "severity_employment_w2", "severity_pollution_w2",
                                     "severity_inequality_w2", "severity_corruption_w2", "severity_dscrm_minority_w2")
var_severity_fig_w12            <- c("severity_welfare", "severity_employment", "severity_pollution",
                                     "severity_inequality", "severity_corruption", "severity_dscrm_minority")

var_democracy_reg_w2            <- c("china_interest_group_w2", "china_rate_democracy_w2", "china_rate_humanrights_w2", "importance_live_in_demo_w2")
var_democracy_reg_w3            <- c("importance_live_in_demo_w3")
var_democracy_reg_fp_w1         <- c("importance_live_in_demo_w1")
var_democracy_reg_fp_w2         <- c("importance_live_in_demo_w2")
var_democracy_reg_fp_w3         <- c("importance_live_in_demo_w3")
var_democracy_fig_w123          <- c("importance_live_in_demo", "az_belief_democracy_fp")
var_democracy_fig_w12           <- c("china_interest_group", "china_rate_democracy", "china_rate_humanrights")

var_contro_justi_reg_policy_w2  <- c("justify_minority_policy_w2", "justify_hukou_w2", "justify_one_child_w2", "justify_hongkong_policy_w2",
                                     "justify_taiwan_policy_w2", "justify_violence_stability_w2", "justify_receive_refugee_w2",
                                     "justify_reduce_pollution_w2", "justify_gaokao_w2", "justify_soe_privatize_w2")
var_contro_justi_reg_liberal_w2 <- c("justify_homo_marriage_w2", "justify_legal_prostitute_w2", "justify_abortion_w2",
                                     "justify_exmaritalsex", "justify_transgene_w2", "justify_soft_drug_w2")
var_contro_justi_fig_w12        <- c("justify_minority_policy", "justify_reduce_pollution", "justify_hukou", "justify_one_child",
                                     "justify_gaokao", "justify_hongkong_policy", "justify_taiwan_policy", "justify_transgene",
                                     "justify_receive_refugee", "justify_soe_privatize", "justify_homo_marriage", "justify_legal_prostitute",
                                     "justify_abortion", "justify_soft_drug", "justify_violence_stability")

var_willing_fight_reg_w1        <- c("willing_against_illi_govt_w1", "willing_report_mis_w1", "willing_protect_weak_w1")
var_willing_fight_reg_w2        <- c("willing_against_illi_govt_w2", "willing_report_mis_w2", "willing_protect_weak_w2")
var_willing_fight_reg_w3        <- c("willing_against_illi_govt_w3", "willing_report_mis_w3", "willing_protect_weak_w3")
var_willing_fight_fig_123       <- c("willing_against_illi_govt", "willing_report_mis", "willing_protect_weak", "az_belief_willing")

var_interest_reg_w2             <- c("interest_economic_w2", "interest_politics_w2")
var_interest_fig_w12            <- c("interest_economic", "interest_politics")

var_patriotism_reg_w2           <- c("proud_being_chinese_w2")
var_patriotism_fig_w12          <- c("proud_being_chinese")

var_fear_critgovt_reg_w2        <- c("fear_critic_govt_self_w2")
var_fear_critgovt_fig_w12       <- c("fear_critic_govt_self")


# 	E: Behaviors
var_socialinteract_reg_w1       <- c("frequency_talk_politic_w1", "frequency_persuade_friends_w1")
var_socialinteract_reg_w2       <- c("frequency_talk_politic_w2", "frequency_persuade_friends_w2")
var_socialinteract_reg_w3       <- c("frequency_talk_politic_w3", "frequency_persuade_friends_w3")
var_socialinteract_fig_w123     <- c("frequency_talk_politic", "frequency_persuade_friends", "az_var_socialinteract")

var_socialknowledge_w1          <- c("know_attitudes_relatives", "know_attitudes_schoolmates", "know_attitudes_outschool_friends")

var_polparticipation_reg_w2     <- c("participate_ngo_w2", "participate_social_protest_w2", "participate_plan_vote_w2", "participate_complain_school_w2")
var_polparticipation_reg_w3     <- c("participate_social_protest_w3", "participate_plan_vote_w3", "participate_complain_school_w3")
var_polparticipation_reg_pf_w1  <- c("participate_social_protest_w1", "participate_plan_vote_w1", "participate_complain_school_w1")
var_polparticipation_reg_pf_w2  <- c("participate_social_protest_w2", "participate_plan_vote_w2", "participate_complain_school_w2")
var_polparticipation_reg_pf_w3  <- c("participate_social_protest_w3", "participate_plan_vote_w3", "participate_complain_school_w3")
var_polparticipation_fig_w123   <- c("participate_social_protest", "participate_plan_vote", "participate_complain_school", "az_var_polparticipation")
var_polparticipation_fig_w12    <- c("participate_ngo")

var_polengagement_w1            <- c("ccp_member", "ccp_join_year", "ccp_prep_member", "participate_stuunion",
                                     "participate_stuunion_ever", "participate_tuanwei", "participate_tuanwei_ever")

var_planaftergrad_reg_w1        <- c("plan_grad_gradschool_dom_w1", "plan_grad_foreignmaster_w1", "plan_grad_foreignphd_w1",
                                     "plan_grad_military_w1", "plan_grad_work_w1")
var_planaftergrad_reg_w2        <- c("plan_grad_gradschool_dom_w2", "plan_grad_foreignmaster_w2", "plan_grad_foreignphd_w2",
                                     "plan_grad_military_w2", "plan_grad_work_w2")
var_planaftergrad_reg_w3        <- c("plan_grad_gradschool_dom_w3", "plan_grad_foreignmaster_w3", "plan_grad_foreignphd_w3",
                                     "plan_grad_military_w3", "plan_grad_work_w3")
var_planaftergrad_fig_w123      <- c("plan_grad_gradschool_dom", "plan_grad_foreignmaster", "plan_grad_foreignphd",
                                     "plan_grad_military", "plan_grad_work")

var_career_sector_reg_w1        <- c("cp_t3_national_civil_w1", "cp_t3_local_civil_w1", "cp_t3_military_w1", "cp_t3_chinese_private_w1",
                                     "cp_t3_for_firm_w1", "cp_t3_soe_w1", "cp_t3_institutional_w1", "cp_t3_entrepreneur_w1")
var_career_sector_reg_w2        <- c("cp_t3_national_civil_w2", "cp_t3_local_civil_w2", "cp_t3_military_w2", "cp_t3_chinese_private_w2",
                                     "cp_t3_for_firm_w2", "cp_t3_soe_w2", "cp_t3_institutional_w2", "cp_t3_entrepreneur_w2")
var_career_sector_reg_w3        <- c("cp_t3_national_civil_w3", "cp_t3_local_civil_w3", "cp_t3_military_w3", "cp_t3_chinese_private_w3",
                                     "cp_t3_for_firm_w3", "cp_t3_soe_w3", "cp_t3_institutional_w3", "cp_t3_entrepreneur_w3")
var_career_sector_fig_w123      <- c("cp_t1_national_civil", "cp_t1_local_civil", "cp_t1_military", "cp_t1_chinese_private",
                                     "cp_t1_for_firm", "cp_t1_soe", "cp_t1_institutional", "cp_t1_entrepreneur")

var_career_loc_reg_w1           <- c("cloc_beijing_w1", "cloc_shanghai_w1", "cloc_gzsz_w1", "cloc_tjcq_w1",
                                     "cloc_hkmc_w1", "cloc_taiwan_w1", "cloc_dom_w1", "cloc_for_w1")
var_career_loc_reg_w2           <- c("cloc_beijing_w2", "cloc_shanghai_w2", "cloc_gzsz_w2", "cloc_tjcq_w2",
                                     "cloc_hkmc_w2", "cloc_taiwan_w2", "cloc_dom_w2", "cloc_for_w2")
# NOTE: kept as-is — original Stata line has "cloc_dom_w2" (not "_w3") inside
# the wave-3 list; this looks like a typo in the source code, but per the
# instruction to translate faithfully, it is preserved unchanged here.
var_career_loc_reg_w3           <- c("cloc_beijing_w3", "cloc_shanghai_w3", "cloc_gzsz_w3", "cloc_tjcq_w3",
                                     "cloc_hkmc_w3", "cloc_taiwan_w3", "cloc_dom_w2", "cloc_for_w3")
var_career_loc_reg_w123         <- c("cloc_beijing", "cloc_shanghai", "cloc_gzsz", "cloc_tjcq",
                                     "cloc_hkmc", "cloc_taiwan", "cloc_dom", "cloc_for")

var_stock_invest_reg_w1         <- c("stock_participation_w1")
var_stock_invest_reg_w2         <- c("stock_participation_w2")
var_stock_invest_reg_w3         <- c("stock_participation_w3")
var_stock_invest_fig_w123       <- c("stock_participation")


# 	F: Demographics, background characteristics, and fundamental preferences

## demographics
var_demog_reg                   <- c("gender", "birth_year", "ethnicity_han", "birthplace_coastal", "residence_coastal", "hukou_urban",
                                     "religion_religious", "ccp_member", "university_elite", "hs_track_science", "department_ssh",
                                     "siblings_total", "father_edu_hsabove", "father_ccp", "mother_edu_hsabove", "mother_ccp",
                                     "hh_income", "domestic_english_atleast4", "foreign_english_yes", "travel_hktaiwan", "travel_foreign_yes")
var_demog_reg_personal          <- c("gender", "birth_year", "ethnicity_han", "birthplace_coastal", "residence_coastal",
                                     "hukou_urban", "religion_religious", "ccp_member")
var_demog_reg_education         <- c("university_elite", "hs_track_science", "department_ssh")
var_demog_reg_english           <- c("domestic_english_atleast4", "foreign_english_yes")
var_demog_reg_travel            <- c("travel_hktaiwan", "travel_foreign_yes")
var_demog_reg_household         <- c("siblings_total", "father_edu_hsabove", "work_father_govt", "father_ccp",
                                     "mother_edu_hsabove", "work_mother_govt", "mother_ccp", "hh_income")

var_demog_reg_imbalance         <- c("residence_coastal", "hs_track_science", "father_ccp", "mother_ccp", "gift_amount", "risk_preference_ce")

## fundamental preferences
var_preference_risk             <- c("willing_risk", "risk_preference_ce", "risk_lottery_choice")
var_preference_time              <- c("willing_future", "procrastinate")
var_preference_altruism         <- c("willing_goodcauses", "donate_amount")
var_preference_reciprocity      <- c("willing_returnfavor", "best_intentions", "gift_amount",
                                     "willing_punish_you", "willing_punish_others", "willing_revenge")


# 	G: List experiment
var_listexp_all                 <- c("list_exp_direct", "list_trust_direct_count", "list_trust_direct_yes", "list_trust_veiled_count")
var_listexp_reg                 <- c("list_count_trust")


# 	Y: Overall effects
var_overall_w1 <- c("info_freq_website_for_w1", "az_belief_media_value_w1", "az_belief_media_trust_w1", "bias_domestic_w1",
                    "az_belief_media_justif_w1", "bias_foreign_r_w1", "news_perccor_cen_w1", "news_perccor_unc_w1",
                    "protest_pcheard_china_w1", "protest_pcheard_foreign_w1", "az_knowledge_meta_w1", "az_belief_econ_perf_cn_r_w1",
                    "az_belief_econ_conf_cn_w1", "az_belief_instchange_w1", "az_belief_trust_govt_r_w1", "az_belief_trust_foreign_w1",
                    "az_belief_evalgovt_r_w1", "importance_live_in_demo_r_w1", "az_belief_willing_w1", "az_var_socialinteract_w1",
                    "az_var_polparticipation_w1", "plan_grad_gradschool_dom_w1", "plan_grad_foreignmaster_w1", "plan_grad_foreignphd_w1",
                    "plan_grad_military_w1", "plan_grad_work_w1", "cloc_beijing_w1", "cloc_shanghai_w1", "cloc_gzsz_w1",
                    "cloc_tjcq_w1", "cloc_hkmc_w1", "cloc_taiwan_w1", "cloc_dom_w2", "cloc_for_w1", "stock_participation_w1")
# NOTE: "cloc_dom_w2" inside the wave-1 list is preserved exactly as in the
# original Stata code (likely a typo for "cloc_dom_w1", but left unchanged
# per instructions to translate faithfully without altering content).

var_overall_w2 <- c("info_freq_website_for_w2", "az_belief_media_value_w2", "az_belief_media_trust_w2", "bias_domestic_w2",
                    "az_belief_media_justif_w2", "bias_foreign_r_w2", "news_perccor_cen_w2", "news_perccor_unc_w2",
                    "protest_pcheard_china_w2", "protest_pcheard_foreign_w2", "az_knowledge_meta_w2", "az_belief_econ_perf_cn_r_w2",
                    "az_belief_econ_conf_cn_w2", "az_belief_instchange_w2", "az_belief_trust_govt_r_w2", "az_belief_trust_foreign_w2",
                    "az_belief_evalgovt_r_w2", "importance_live_in_demo_r_w2", "az_belief_willing_w2", "az_var_socialinteract_w2",
                    "az_var_polparticipation_w2", "plan_grad_gradschool_dom_w2", "plan_grad_foreignmaster_w2", "plan_grad_foreignphd_w2",
                    "plan_grad_military_w2", "plan_grad_work_w2", "cloc_beijing_w2", "cloc_shanghai_w2", "cloc_gzsz_w2",
                    "cloc_tjcq_w2", "cloc_hkmc_w2", "cloc_taiwan_w2", "cloc_dom_w2", "cloc_for_w2", "stock_participation_w2")

var_overall_w3 <- c("info_freq_website_for_w3", "az_belief_media_value_w3", "az_belief_media_trust_w3", "bias_domestic_w3",
                    "az_belief_media_justif_w3", "bias_foreign_r_w3", "news_perccor_cen_w3", "news_perccor_unc_w3",
                    "protest_pcheard_china_w3", "protest_pcheard_foreign_w3", "az_knowledge_meta_w3", "az_belief_econ_perf_cn_r_w3",
                    "az_belief_econ_conf_cn_w3", "az_belief_instchange_w3", "az_belief_trust_govt_r_w3", "az_belief_trust_foreign_w3",
                    "az_belief_evalgovt_r_w3", "importance_live_in_demo_r_w3", "az_belief_willing_w3", "az_var_socialinteract_w3",
                    "az_var_polparticipation_w3", "plan_grad_gradschool_dom_w3", "plan_grad_foreignmaster_w3", "plan_grad_foreignphd_w3",
                    "plan_grad_military_w3", "plan_grad_work_w3", "cloc_beijing_w3", "cloc_shanghai_w3", "cloc_gzsz_w3",
                    "cloc_tjcq_w3", "cloc_hkmc_w3", "cloc_taiwan_w3", "cloc_dom_w3", "cloc_for_w3", "stock_participation_w3")



# Prepare data: merge additional variables --------------------------------

# Each join below is keyed on responseID_wave1, matching the original Stata `merge 1:1`.
vpn_date    <- zap_labels(read_dta(file.path(data_path, "vpn_date_adoption.dta")))
vpn_active  <- zap_labels(read_dta(file.path(data_path, "vpn_browsing_active_user.dta")))
vpn_purch   <- zap_labels(read_dta(file.path(data_path, "vpn_purchase.dta")))

# Stata's missing value `.` becomes R's `NA`.
# tidyr's replace_na() is equivalent of `replace var = 0 if
# var == .`.
cy <- cy |>
  left_join(vpn_date,   by = "responseID_wave1") |>
  left_join(vpn_active, by = "responseID_wave1") |>
  mutate(active_user = replace_na(active_user, 0)) |>
  left_join(vpn_purch,  by = "responseID_wave1")


# Prepare data: generate additional variables -----------------------------


# wave 1 ------------------------------------------------------------

cy <- cy %>%
  mutate(
    news_c_stock_rise            = as.integer(news_stock_rise == 1),
    news_c_exchange_depreciate   = as.integer(news_exchange_depreciate == 1),
    news_c_train_brazil_peru     = as.integer(news_train_brazil_peru == 0),
    news_c_army_reduce           = as.integer(news_army_reduce == 1),
    news_c_taiwan_election       = as.integer(news_taiwan_election == 0),
    news_c_china_us_network_coop = as.integer(news_china_us_network_coop == 0),
    news_c_lijiacheng_china      = as.integer(news_lijiacheng_china == 0),
    news_c_nanjing_memory_list   = as.integer(news_nanjing_memory_list == 1)
  )

## generate count variable: # news answered correctly
cy <- cy %>%
  mutate(
    news_totalcorrect = news_c_stock_rise + news_c_exchange_depreciate + news_c_train_brazil_peru +
      news_c_army_reduce + news_c_taiwan_election + news_c_china_us_network_coop +
      news_c_lijiacheng_china + news_c_nanjing_memory_list
  )

## generate percentage correct
cy <- cy %>% mutate(news_perccorrect = news_totalcorrect / 8)

## correct count in percentage: censored news
cy <- cy %>% mutate(news_perccor_cen = (news_c_taiwan_election + news_c_lijiacheng_china) / 2)

## correct count in percentage: uncensored news
cy <- cy %>%
  mutate(
    news_perccor_unc = (news_c_stock_rise + news_c_exchange_depreciate + news_c_train_brazil_peru +
                          news_c_army_reduce + news_c_china_us_network_coop + news_c_nanjing_memory_list) / 6
  )

## generate count variable: # people heard before
cy <- cy %>%
  mutate(
    people_totalheard = people_puzhiqiang_w1 + people_lizehou_w1 + people_huangzhifeng_w1 +
      people_chenguangcheng_w1 + people_lixiaolin_w1 + people_renzhiqiang_w1 +
      people_maoyushi_w1 + people_honghuang_w1 + people_liuqiangdong_w1,
    people_pcheard = people_totalheard / 9
  )

## generate count variable: # protests heard before
cy <- cy %>%
  mutate(
    protest_totalheard = protest_2014_europe_square_w1 + protest_2014_sun_flower_w1 +
      protest_2010_arabic_spring_w1 + protest_2014_crimea_vote_w1 + protest_2012_hk_curriculum_w1 +
      protest_2010_catal_indep_w1 + protest_2014_umbrella_w1,
    protest_pcheard_total = protest_totalheard / 9
    # only 7 protests but divided by 9?
  )

## generate percentage by category: full panel version
cy <- cy %>%
  mutate(
    protest_pcheard_china = (protest_2014_sun_flower_w1 + protest_2012_hk_curriculum_w1 + protest_2014_umbrella_w1) / 3,
    protest_pcheard_foreign = (protest_2014_europe_square_w1 + protest_2010_arabic_spring_w1 +
                                 protest_2014_crimea_vote_w1 + protest_2010_catal_indep_w1) / 4
  )

## political participation: recode into ever participated
cy <- cy %>%
  mutate(
    participate_stuunion_ever = as.integer(participate_stuunion < 3),
    participate_tuanwei_ever  = as.integer(participate_tuanwei < 3)
  )

## generate distance from neutrality
bias_vars <- c("us_pos_us", "us_pos_cn", "cn_pos_us", "cn_pos_cn",
               "us_neg_us", "us_neg_cn", "cn_neg_us", "cn_neg_cn")

for (v in bias_vars) {
  src <- paste0("bias_", v, "_media_w1")
  out <- paste0("distneutral_", v)
  cy[[out]] <- abs(cy[[src]] - 4)
}

## separate censorship and bias
for (v in bias_vars) {
  src <- paste0("bias_", v, "_media_w1")
  out <- paste0("bias_", v, "_cens")
  cy[[out]] <- as.integer(cy[[src]] == 1)
}

## generate indicator for top categories: bias_for
for_vars <- c("bias_for_govt_policy", "bias_for_firm_interest", "bias_for_media_pref", "bias_for_reader_demand")
for (v in for_vars) {
  src <- paste0(v, "_w1")
  out <- paste0(v, "_t1")
  cy[[out]] <- as.integer(cy[[src]] == 1)
}

## generate indicator for top categories: bias_dom
dom_vars <- c("bias_dom_govt_policy", "bias_dom_firm_interest", "bias_dom_media_pref", "bias_dom_reader_demand")
for (v in dom_vars) {
  src <- paste0(v, "_w1")
  out <- paste0(v, "_t1")
  cy[[out]] <- as.integer(cy[[src]] == 1)
}

## generate relative weights for evaluating government performance
cy <- cy %>%
  mutate(
    evalgovt_total_w1 = evalgovt_election_w1 + evalgovt_economy_w1 + evalgovt_equality_w1 +
      evalgovt_ruleoflaw_w1 + evalgovt_human_rights_w1 + evalgovt_freedom_speech_w1 +
      evalgovt_global_power_w1 + evalgovt_fair_history_w1
  )

evalgovt_vars <- c("evalgovt_election", "evalgovt_economy", "evalgovt_equality", "evalgovt_ruleoflaw",
                   "evalgovt_human_rights", "evalgovt_freedom_speech", "evalgovt_global_power",
                   "evalgovt_fair_history")
for (v in evalgovt_vars) {
  src <- paste0(v, "_w1")
  out <- paste0("r", v)
  cy[[out]] <- cy[[src]] / cy$evalgovt_total_w1
}

## generate dummies for plan_graduation
cy <- cy %>%
  mutate(
    plan_grad_gradschool_dom = as.integer(plan_graduation_w1 == 1),
    plan_grad_foreignmaster  = as.integer(plan_graduation_w1 == 2),
    plan_grad_foreignphd     = as.integer(plan_graduation_w1 == 3),
    plan_grad_military       = as.integer(plan_graduation_w1 == 4),
    plan_grad_work           = as.integer(plan_graduation_w1 == 5),
    plan_grad_dontknow       = as.integer(plan_graduation_w1 == 6)
  )

## generate top1 and top3 dummies for career preferences
# Use rowSums(is.na(...)) == 8 to check "all 8 source vars are NA for this row"
work_top3_stems <- c("national_civil", "local_civil", "military", "chinese_private",
                     "for_firm", "soe", "institutional", "entrepreneur")
work_top3_w1_cols <- paste0("work_top3_", work_top3_stems, "_w1")

cy <- cy %>%
  mutate(work_top3_missing = as.integer(rowSums(is.na(across(all_of(work_top3_w1_cols)))) == length(work_top3_w1_cols)))

for (v in work_top3_stems) {
  src <- paste0("work_top3_", v, "_w1")
  cy[[paste0("cp_t1_", v)]] <- ifelse(cy$work_top3_missing == 0, as.integer(cy[[src]] == 1), NA_integer_)
  cy[[paste0("cp_t3_", v)]] <- ifelse(cy$work_top3_missing == 0, as.integer(!is.na(cy[[src]])), NA_integer_)
}
cy$work_top3_missing <- NULL  # drop work_top3_missing

## generate dummies for career location preferences
cy <- cy %>%
  mutate(
    cloc_beijing  = as.integer(place_top_w1 == 1),
    cloc_shanghai = as.integer(place_top_w1 == 2),
    cloc_gzsz     = as.integer(place_top_w1 == 3 | place_top_w1 == 4),
    cloc_tjcq     = as.integer(place_top_w1 == 5 | place_top_w1 == 6),
    cloc_hkmc     = as.integer(place_top_w1 == 7 | place_top_w1 == 8),
    cloc_taiwan   = as.integer(place_top_w1 == 9),
    cloc_dom      = as.integer(place_top_w1 == 10),
    cloc_for      = as.integer(place_top_w1 == 11)
  )

## add subscript: wave 1
rename_w1 <- c("news_totalcorrect", "news_perccorrect", "news_perccor_cen", "news_perccor_unc",
               "people_totalheard", "people_pcheard", "protest_totalheard", "protest_pcheard_total",
               "protest_pcheard_china", "protest_pcheard_foreign", "distneutral_us_pos_us", "distneutral_us_pos_cn",
               "distneutral_cn_pos_us", "distneutral_cn_pos_cn", "distneutral_us_neg_us", "distneutral_us_neg_cn",
               "distneutral_cn_neg_us", "distneutral_cn_neg_cn", "bias_us_pos_us_cens", "bias_us_pos_cn_cens",
               "bias_cn_pos_us_cens", "bias_cn_pos_cn_cens", "bias_us_neg_us_cens", "bias_us_neg_cn_cens",
               "bias_cn_neg_us_cens", "bias_cn_neg_cn_cens", "bias_for_govt_policy_t1", "bias_for_firm_interest_t1",
               "bias_for_media_pref_t1", "bias_for_reader_demand_t1", "bias_dom_govt_policy_t1",
               "bias_dom_firm_interest_t1", "bias_dom_media_pref_t1", "bias_dom_reader_demand_t1",
               "revalgovt_election", "revalgovt_economy", "revalgovt_equality", "revalgovt_ruleoflaw",
               "revalgovt_human_rights", "revalgovt_freedom_speech", "revalgovt_global_power", "revalgovt_fair_history",
               "plan_grad_gradschool_dom", "plan_grad_foreignmaster", "plan_grad_foreignphd", "plan_grad_military",
               "plan_grad_work", "plan_grad_dontknow", "cp_t1_national_civil", "cp_t3_national_civil",
               "cp_t1_local_civil", "cp_t3_local_civil", "cp_t1_military", "cp_t3_military", "cp_t1_chinese_private",
               "cp_t3_chinese_private", "cp_t1_for_firm", "cp_t3_for_firm", "cp_t1_soe", "cp_t3_soe",
               "cp_t1_institutional", "cp_t3_institutional", "cp_t1_entrepreneur", "cp_t3_entrepreneur",
               "cloc_beijing", "cloc_shanghai", "cloc_gzsz", "cloc_tjcq", "cloc_hkmc", "cloc_taiwan", "cloc_dom", "cloc_for")

cy <- cy %>% rename_with(.fn = ~ paste0(.x, "_w1"), .cols = all_of(rename_w1))


# 	wave 2 ------------------------------------------------------------

## foreach news, indicator correct or not
cy <- cy %>%
  mutate(
    news_c_topinequality  = as.integer(news_topinequality == 0),
    news_c_cpicensorship  = as.integer(news_cpicensorship == 1),
    news_c_laborunrest    = as.integer(news_laborunrest == 1),
    news_c_waterpollution = as.integer(news_waterpollution == 0),
    news_c_stockcrash     = as.integer(news_stockcrash == 1),
    news_c_taiwanelection = as.integer(news_taiwanelection == 0),
    news_c_applefbi       = as.integer(news_applefbi == 0),
    news_c_tenyearshk     = as.integer(news_tenyearshk == 1),
    news_c_panamapapers   = as.integer(news_panamapapers == 1),
    news_c_yihehotel      = as.integer(news_yihehotel == 1),
    news_c_economistcensor = as.integer(news_economistcensor == 1)
  )

## generate count variable: # news answered correctly
cy <- cy %>%
  mutate(
    news_totalcorrect = news_c_topinequality + news_c_cpicensorship + news_c_laborunrest +
      news_c_waterpollution + news_c_stockcrash + news_c_taiwanelection + news_c_applefbi +
      news_c_tenyearshk + news_c_panamapapers + news_c_yihehotel + news_c_economistcensor
  )

## generate percentage correct
cy <- cy %>% mutate(news_perccorrect = news_totalcorrect / 11)

## correct count in percentage: quiz questions
cy <- cy %>%
  mutate(news_perccor_qui = (news_c_topinequality + news_c_cpicensorship + news_c_laborunrest + news_c_waterpollution) / 4)

## correct count in percentage: censored news
cy <- cy %>%
  mutate(news_perccor_cen = (news_c_stockcrash + news_c_panamapapers + news_c_tenyearshk + news_c_economistcensor) / 4)

## correct count in percentage: uncensored news
cy <- cy %>%
  mutate(news_perccor_unc = (news_c_taiwanelection + news_c_applefbi + news_c_yihehotel) / 3)

## generate count variable: # people heard before
cy <- cy %>%
  mutate(
    people_totalheard = people_puzhiqiang_w2 + people_lizehou_w2 + people_huangzhifeng_w2 +
      people_chenguangcheng_w2 + people_lixiaolin_w2 + people_renzhiqiang_w2 +
      people_maoyushi_w2 + people_honghuang_w2 + people_liuqiangdong_w2,
    people_pcheard = people_totalheard / 9
  )

## percentage of names heard by category
cy <- cy %>%
  mutate(
    people_perchd_censored = (people_lixiaolin_w2 + people_huangzhifeng_w2 + people_chenguangcheng_w2 + people_lizehou_w2) / 4,
    people_perchd_uncensor = (people_maoyushi_w2 + people_liuqiangdong_w2 + people_honghuang_w2) / 3,
    people_perchd_unctocen = (people_puzhiqiang_w2 + people_renzhiqiang_w2) / 2
  )

## generate count variable: # events heard before
cy <- cy %>%
  mutate(
    protest_totalheard = protest_2014_europe_square_w2 + protest_2014_sun_flower_w2 +
      protest_2010_arabic_spring_w2 + protest_2014_crimea_vote_w2 + protest_2012_hk_curriculum_w2 +
      protest_2010_catal_indep_w2 + protest_2014_umbrella_w2 + protest_2016_mongkok_riot_w2,
    protest_pcheard_total = protest_totalheard / 8
  )

## generate percentage by category: full panel version
cy <- cy %>%
  mutate(
    protest_pcheard_china = (protest_2014_sun_flower_w2 + protest_2012_hk_curriculum_w2 + protest_2014_umbrella_w2) / 3,
    protest_pcheard_foreign = (protest_2014_europe_square_w2 + protest_2010_arabic_spring_w2 +
                                 protest_2014_crimea_vote_w2 + protest_2010_catal_indep_w2) / 4
  )

## generate distance from neutrality
for (v in bias_vars) {
  src <- paste0("bias_", v, "_media_w2")
  out <- paste0("distneutral_", v)
  cy[[out]] <- abs(cy[[src]] - 4)
}

## separate censorship and bias
for (v in bias_vars) {
  src <- paste0("bias_", v, "_media_w2")
  out <- paste0("bias_", v, "_cens")
  cy[[out]] <- as.integer(cy[[src]] == 1)
}

## generate indicator for top categories: bias_for
for (v in for_vars) {
  src <- paste0(v, "_w2")
  out <- paste0(v, "_t1")
  cy[[out]] <- as.integer(cy[[src]] == 1)
}

## generate indicator for top categories: bias_dom
for (v in dom_vars) {
  src <- paste0(v, "_w2")
  out <- paste0(v, "_t1")
  cy[[out]] <- as.integer(cy[[src]] == 1)
}

## generate relative weights for evaluating government performance
cy <- cy %>%
  mutate(
    evalgovt_total_w2 = evalgovt_election_w2 + evalgovt_economy_w2 + evalgovt_equality_w2 +
      evalgovt_ruleoflaw_w2 + evalgovt_human_rights_w2 + evalgovt_freedom_speech_w2 +
      evalgovt_global_power_w2 + evalgovt_fair_history_w2
  )

for (v in evalgovt_vars) {
  src <- paste0(v, "_w2")
  out <- paste0("r", v)
  cy[[out]] <- cy[[src]] / cy$evalgovt_total_w2
}

## generate dummies for plan_graduation
cy <- cy %>%
  mutate(
    plan_grad_gradschool_dom = as.integer(plan_graduation_w2 == 1),
    plan_grad_foreignmaster  = as.integer(plan_graduation_w2 == 2),
    plan_grad_foreignphd     = as.integer(plan_graduation_w2 == 3),
    plan_grad_military       = as.integer(plan_graduation_w2 == 4),
    plan_grad_work           = as.integer(plan_graduation_w2 == 5),
    plan_grad_dontknow       = as.integer(plan_graduation_w2 == 6)
  )

## generate top1 and top3 dummies for career preferences
work_top3_w2_cols <- paste0("work_top3_", work_top3_stems, "_w2")
cy <- cy %>%
  mutate(work_top3_missing = as.integer(rowSums(is.na(across(all_of(work_top3_w2_cols)))) == length(work_top3_w2_cols)))

for (v in work_top3_stems) {
  src <- paste0("work_top3_", v, "_w2")
  cy[[paste0("cp_t1_", v)]] <- ifelse(cy$work_top3_missing == 0, as.integer(cy[[src]] == 1), NA_integer_)
  cy[[paste0("cp_t3_", v)]] <- ifelse(cy$work_top3_missing == 0, as.integer(!is.na(cy[[src]])), NA_integer_)
}
cy$work_top3_missing <- NULL

## generate dummies for career location preferences
cy <- cy %>%
  mutate(
    cloc_beijing  = as.integer(place_top_w2 == 1),
    cloc_shanghai = as.integer(place_top_w2 == 2),
    cloc_gzsz     = as.integer(place_top_w2 == 3 | place_top_w2 == 4),
    cloc_tjcq     = as.integer(place_top_w2 == 5 | place_top_w2 == 6),
    cloc_hkmc     = as.integer(place_top_w2 == 7 | place_top_w2 == 8),
    cloc_taiwan   = as.integer(place_top_w2 == 9),
    cloc_dom      = as.integer(place_top_w2 == 10),
    cloc_for      = as.integer(place_top_w2 == 11)
  )

## add subscript: wave 2
rename_w2 <- c("news_totalcorrect", "news_perccorrect", "news_perccor_cen", "news_perccor_unc",
               "people_totalheard", "people_pcheard", "people_perchd_censored", "people_perchd_uncensor",
               "people_perchd_unctocen", "protest_totalheard", "protest_pcheard_total", "protest_pcheard_china",
               "protest_pcheard_foreign", "distneutral_us_pos_us", "distneutral_us_pos_cn", "distneutral_cn_pos_us",
               "distneutral_cn_pos_cn", "distneutral_us_neg_us", "distneutral_us_neg_cn", "distneutral_cn_neg_us",
               "distneutral_cn_neg_cn", "bias_us_pos_us_cens", "bias_us_pos_cn_cens", "bias_cn_pos_us_cens",
               "bias_cn_pos_cn_cens", "bias_us_neg_us_cens", "bias_us_neg_cn_cens", "bias_cn_neg_us_cens",
               "bias_cn_neg_cn_cens", "bias_for_govt_policy_t1", "bias_for_firm_interest_t1", "bias_for_media_pref_t1",
               "bias_for_reader_demand_t1", "bias_dom_govt_policy_t1", "bias_dom_firm_interest_t1",
               "bias_dom_media_pref_t1", "bias_dom_reader_demand_t1", "revalgovt_election", "revalgovt_economy",
               "revalgovt_equality", "revalgovt_ruleoflaw", "revalgovt_human_rights", "revalgovt_freedom_speech",
               "revalgovt_global_power", "revalgovt_fair_history", "plan_grad_gradschool_dom", "plan_grad_foreignmaster",
               "plan_grad_foreignphd", "plan_grad_military", "plan_grad_work", "plan_grad_dontknow",
               "cp_t1_national_civil", "cp_t3_national_civil", "cp_t1_local_civil", "cp_t3_local_civil", "cp_t1_military",
               "cp_t3_military", "cp_t1_chinese_private", "cp_t3_chinese_private", "cp_t1_for_firm", "cp_t3_for_firm",
               "cp_t1_soe", "cp_t3_soe", "cp_t1_institutional", "cp_t3_institutional", "cp_t1_entrepreneur",
               "cp_t3_entrepreneur", "cloc_beijing", "cloc_shanghai", "cloc_gzsz", "cloc_tjcq", "cloc_hkmc",
               "cloc_taiwan", "cloc_dom", "cloc_for")

cy <- cy %>% rename_with(.fn = ~ paste0(.x, "_w2"), .cols = all_of(rename_w2))


# 	wave 3 ------------------------------------------------------------

## foreach news, indicator correct or not
cy <- cy %>%
  mutate(
    news_c_coalprod       = as.integer(news_coalprod == 0),
    news_c_trumpchina     = as.integer(news_trumpchina == 1),
    news_c_xiaojianhua    = as.integer(news_xiaojianhua == 0),
    news_c_northkoreacoal = as.integer(news_northkoreacoal == 1),
    news_c_birdflu        = as.integer(news_birdflu == 1),
    news_c_ethiopiatrain  = as.integer(news_ethiopiatrain == 0),
    news_c_foreignreserve = as.integer(news_foreignreserve == 0),
    news_c_xijiangcar     = as.integer(news_xijiangcar == 0),
    news_c_chinanorway    = as.integer(news_chinanorway == 1),
    news_c_womenrights    = as.integer(news_womenrights == 0),
    news_c_hkceelection   = as.integer(news_hkceelection == 0)
  )

## generate count variable: # news answered correctly
cy <- cy %>%
  mutate(
    news_totalcorrect = news_c_coalprod + news_c_trumpchina + news_c_xiaojianhua + news_c_northkoreacoal +
      news_c_birdflu + news_c_ethiopiatrain + news_c_foreignreserve + news_c_xijiangcar +
      news_c_chinanorway + news_c_womenrights + news_c_hkceelection
  )

## generate percentage correct
cy <- cy %>% mutate(news_perccorrect = news_totalcorrect / 11)

## correct count in percentage: censored news
cy <- cy %>%
  mutate(
    news_perccor_cen = (news_c_coalprod + news_c_trumpchina + news_c_xiaojianhua + news_c_xijiangcar +
                          news_c_chinanorway + news_c_womenrights + news_c_hkceelection) / 7
  )

## correct count in percentage: uncensored news
cy <- cy %>%
  mutate(news_perccor_unc = (news_c_northkoreacoal + news_c_birdflu + news_c_ethiopiatrain + news_c_foreignreserve) / 4)

## generate count variable: # events heard before
cy <- cy %>%
  mutate(
    protest_totalheard = protest_2014_europe_square_w3 + protest_2014_sun_flower_w3 +
      protest_2010_arabic_spring_w3 + protest_2014_crimea_vote_w3 + protest_2012_hk_curriculum_w3 +
      protest_2010_catal_indep_w3 + protest_2014_umbrella_w3 + protest_2016_mongkok_riot_w3 +
      protest_2017_women_march_w3,
    protest_pcheard_total = protest_totalheard / 9
  )

## generate percentage by category: full panel version
cy <- cy %>%
  mutate(
    protest_pcheard_china = (protest_2014_sun_flower_w3 + protest_2012_hk_curriculum_w3 + protest_2014_umbrella_w3) / 3,
    protest_pcheard_foreign = (protest_2014_europe_square_w3 + protest_2010_arabic_spring_w3 +
                                 protest_2014_crimea_vote_w3 + protest_2010_catal_indep_w3) / 4
  )

## generate indicator for top categories: bias_for
for (v in for_vars) {
  src <- paste0(v, "_w3")
  out <- paste0(v, "_t1")
  cy[[out]] <- as.integer(cy[[src]] == 1)
}

## generate indicator for top categories: bias_dom
for (v in dom_vars) {
  src <- paste0(v, "_w3")
  out <- paste0(v, "_t1")
  cy[[out]] <- as.integer(cy[[src]] == 1)
}

## generate indicators for regression
cy <- cy %>%
  mutate(
    vpn_purchase_yes = as.integer(vpn_purchase < 6),
    vpn_purchase_wmt = as.integer((vpn_purchase_yes == 1) & (vpn_purchase != 1))
  )

cy <- cy %>%
  mutate(
    vpn_purchase_premium = case_when(
      vpn_purchase == 6 ~ 0,
      vpn_purchase == 1 ~ 1,
      vpn_purchase == 4 ~ 2,
      vpn_purchase == 3 ~ 3,
      vpn_purchase == 5 ~ 4,
      TRUE ~ NA_real_
    )
  )

## generate dummies for plan_graduation
cy <- cy %>%
  mutate(
    plan_grad_gradschool_dom = as.integer(plan_graduation_w3 == 1),
    plan_grad_foreignmaster  = as.integer(plan_graduation_w3 == 2),
    plan_grad_foreignphd     = as.integer(plan_graduation_w3 == 3),
    plan_grad_military       = as.integer(plan_graduation_w3 == 4),
    plan_grad_work           = as.integer(plan_graduation_w3 == 5),
    plan_grad_dontknow       = as.integer(plan_graduation_w3 == 6)
  )

## generate top1 and top3 dummies for career preferences
work_top3_w3_cols <- paste0("work_top3_", work_top3_stems, "_w3")
cy <- cy %>%
  mutate(work_top3_missing = as.integer(rowSums(is.na(across(all_of(work_top3_w3_cols)))) == length(work_top3_w3_cols)))

for (v in work_top3_stems) {
  src <- paste0("work_top3_", v, "_w3")
  cy[[paste0("cp_t1_", v)]] <- ifelse(cy$work_top3_missing == 0, as.integer(cy[[src]] == 1), NA_integer_)
  cy[[paste0("cp_t3_", v)]] <- ifelse(cy$work_top3_missing == 0, as.integer(!is.na(cy[[src]])), NA_integer_)
}
cy$work_top3_missing <- NULL

## generate dummies for career location preferences
cy <- cy %>%
  mutate(
    cloc_beijing  = as.integer(place_top_w3 == 1),
    cloc_shanghai = as.integer(place_top_w3 == 2),
    cloc_gzsz     = as.integer(place_top_w3 == 3 | place_top_w3 == 4),
    cloc_tjcq     = as.integer(place_top_w3 == 5 | place_top_w3 == 6),
    cloc_hkmc     = as.integer(place_top_w3 == 7 | place_top_w3 == 8),
    cloc_taiwan   = as.integer(place_top_w3 == 9),
    cloc_dom      = as.integer(place_top_w3 == 10),
    cloc_for      = as.integer(place_top_w3 == 11)
  )

## add subscript: wave 3
rename_w3 <- c("news_totalcorrect", "news_perccorrect", "news_perccor_cen", "news_perccor_unc",
               "protest_totalheard", "protest_pcheard_total", "protest_pcheard_china", "protest_pcheard_foreign",
               "bias_for_govt_policy_t1", "bias_for_firm_interest_t1", "bias_for_media_pref_t1",
               "bias_for_reader_demand_t1", "bias_dom_govt_policy_t1", "bias_dom_firm_interest_t1",
               "bias_dom_media_pref_t1", "bias_dom_reader_demand_t1", "plan_grad_gradschool_dom",
               "plan_grad_foreignmaster", "plan_grad_foreignphd", "plan_grad_military", "plan_grad_work",
               "plan_grad_dontknow", "cp_t1_national_civil", "cp_t3_national_civil", "cp_t1_local_civil",
               "cp_t3_local_civil", "cp_t1_military", "cp_t3_military", "cp_t1_chinese_private", "cp_t3_chinese_private",
               "cp_t1_for_firm", "cp_t3_for_firm", "cp_t1_soe", "cp_t3_soe", "cp_t1_institutional",
               "cp_t3_institutional", "cp_t1_entrepreneur", "cp_t3_entrepreneur", "cloc_beijing", "cloc_shanghai",
               "cloc_gzsz", "cloc_tjcq", "cloc_hkmc", "cloc_taiwan", "cloc_dom", "cloc_for")

cy <- cy %>% rename_with(.fn = ~ paste0(.x, "_w3"), .cols = all_of(rename_w3))

## generate percentage of total news quiz
cy <- cy %>% mutate(news_perccor_cen_all = (news_perccor_cen_w2 + news_perccor_cen_w3) / 2)


# 	generate certainty equivalent from Falk et al. 2015 staircase risk preference

cy <- cy %>%
  mutate(
    risk_preference_ce = case_when(
      risk_preference_1==1 & risk_preference_17==1 & risk_preference_25==1 & risk_preference_29==1 & risk_preference_31==1 ~ 32,
      risk_preference_1==1 & risk_preference_17==1 & risk_preference_25==1 & risk_preference_29==1 & risk_preference_31==2 ~ 31,
      risk_preference_1==1 & risk_preference_17==1 & risk_preference_25==1 & risk_preference_29==2 & risk_preference_30==1 ~ 30,
      risk_preference_1==1 & risk_preference_17==1 & risk_preference_25==1 & risk_preference_29==2 & risk_preference_30==2 ~ 29,
      risk_preference_1==1 & risk_preference_17==1 & risk_preference_25==2 & risk_preference_26==1 & risk_preference_27==1 ~ 28,
      risk_preference_1==1 & risk_preference_17==1 & risk_preference_25==2 & risk_preference_26==1 & risk_preference_27==2 ~ 27,
      risk_preference_1==1 & risk_preference_17==1 & risk_preference_25==2 & risk_preference_26==2 & risk_preference_28==1 ~ 26,
      risk_preference_1==1 & risk_preference_17==1 & risk_preference_25==2 & risk_preference_26==2 & risk_preference_28==2 ~ 25,
      risk_preference_1==1 & risk_preference_17==2 & risk_preference_18==1 & risk_preference_22==1 & risk_preference_23==1 ~ 24,
      risk_preference_1==1 & risk_preference_17==2 & risk_preference_18==1 & risk_preference_22==1 & risk_preference_23==2 ~ 23,
      risk_preference_1==1 & risk_preference_17==2 & risk_preference_18==1 & risk_preference_22==2 & risk_preference_24==1 ~ 22,
      risk_preference_1==1 & risk_preference_17==2 & risk_preference_18==1 & risk_preference_22==2 & risk_preference_24==2 ~ 21,
      risk_preference_1==1 & risk_preference_17==2 & risk_preference_18==2 & risk_preference_19==1 & risk_preference_20==1 ~ 20,
      risk_preference_1==1 & risk_preference_17==2 & risk_preference_18==2 & risk_preference_19==1 & risk_preference_20==2 ~ 19,
      risk_preference_1==1 & risk_preference_17==2 & risk_preference_18==2 & risk_preference_19==2 & risk_preference_21==1 ~ 18,
      risk_preference_1==1 & risk_preference_17==2 & risk_preference_18==2 & risk_preference_19==2 & risk_preference_21==2 ~ 17,
      risk_preference_1==2 & risk_preference_2==1 & risk_preference_10==1 & risk_preference_14==1 & risk_preference_15==1 ~ 16,
      risk_preference_1==2 & risk_preference_2==1 & risk_preference_10==1 & risk_preference_14==1 & risk_preference_15==2 ~ 15,
      risk_preference_1==2 & risk_preference_2==1 & risk_preference_10==1 & risk_preference_14==2 & risk_preference_16==1 ~ 14,
      risk_preference_1==2 & risk_preference_2==1 & risk_preference_10==1 & risk_preference_14==2 & risk_preference_16==2 ~ 13,
      risk_preference_1==2 & risk_preference_2==1 & risk_preference_10==2 & risk_preference_11==1 & risk_preference_13==1 ~ 12,
      risk_preference_1==2 & risk_preference_2==1 & risk_preference_10==2 & risk_preference_11==1 & risk_preference_13==2 ~ 11,
      risk_preference_1==2 & risk_preference_2==1 & risk_preference_10==2 & risk_preference_11==2 & risk_preference_12==1 ~ 10,
      risk_preference_1==2 & risk_preference_2==1 & risk_preference_10==2 & risk_preference_11==2 & risk_preference_12==2 ~ 9,
      risk_preference_1==2 & risk_preference_2==2 & risk_preference_3==1 & risk_preference_4==1 & risk_preference_5==1 ~ 8,
      risk_preference_1==2 & risk_preference_2==2 & risk_preference_3==1 & risk_preference_4==1 & risk_preference_5==2 ~ 7,
      risk_preference_1==2 & risk_preference_2==2 & risk_preference_3==1 & risk_preference_4==2 & risk_preference_6==1 ~ 6,
      risk_preference_1==2 & risk_preference_2==2 & risk_preference_3==1 & risk_preference_4==2 & risk_preference_6==2 ~ 5,
      risk_preference_1==2 & risk_preference_2==2 & risk_preference_3==2 & risk_preference_7==1 & risk_preference_8==1 ~ 4,
      risk_preference_1==2 & risk_preference_2==2 & risk_preference_3==2 & risk_preference_7==1 & risk_preference_8==2 ~ 3,
      risk_preference_1==2 & risk_preference_2==2 & risk_preference_3==2 & risk_preference_7==2 & risk_preference_9==1 ~ 2,
      risk_preference_1==2 & risk_preference_2==2 & risk_preference_3==2 & risk_preference_7==2 & risk_preference_9==2 ~ 1,
      TRUE ~ NA_real_
    )
  )


# 	generate time preference

cy <- cy %>%
  mutate(
    time_preference_fe = case_when(
      time_preference_1==1 & time_preference_17==1 & time_preference_18==1 & time_preference_22==1 & time_preference_23==1 ~ 32,
      time_preference_1==1 & time_preference_17==1 & time_preference_18==1 & time_preference_22==1 & time_preference_23==2 ~ 31,
      time_preference_1==1 & time_preference_17==1 & time_preference_18==1 & time_preference_22==2 & time_preference_24==1 ~ 30,
      time_preference_1==1 & time_preference_17==1 & time_preference_18==1 & time_preference_22==2 & time_preference_24==2 ~ 29,
      time_preference_1==1 & time_preference_17==1 & time_preference_18==2 & time_preference_19==1 & time_preference_20==1 ~ 28,
      time_preference_1==1 & time_preference_17==1 & time_preference_18==2 & time_preference_19==1 & time_preference_20==2 ~ 27,
      time_preference_1==1 & time_preference_17==1 & time_preference_18==2 & time_preference_19==2 & time_preference_21==1 ~ 26,
      time_preference_1==1 & time_preference_17==1 & time_preference_18==2 & time_preference_19==2 & time_preference_21==2 ~ 25,
      time_preference_1==1 & time_preference_17==2 & time_preference_25==1 & time_preference_29==1 & time_preference_31==1 ~ 24,
      time_preference_1==1 & time_preference_17==2 & time_preference_25==1 & time_preference_29==1 & time_preference_31==2 ~ 23,
      time_preference_1==1 & time_preference_17==2 & time_preference_25==1 & time_preference_29==2 & time_preference_30==1 ~ 22,
      time_preference_1==1 & time_preference_17==2 & time_preference_25==1 & time_preference_29==2 & time_preference_30==2 ~ 21,
      time_preference_1==1 & time_preference_17==2 & time_preference_25==2 & time_preference_26==1 & time_preference_28==1 ~ 20,
      time_preference_1==1 & time_preference_17==2 & time_preference_25==2 & time_preference_26==1 & time_preference_28==2 ~ 19,
      time_preference_1==1 & time_preference_17==2 & time_preference_25==2 & time_preference_26==2 & time_preference_27==1 ~ 18,
      time_preference_1==1 & time_preference_17==2 & time_preference_25==2 & time_preference_26==2 & time_preference_27==2 ~ 17,
      time_preference_1==2 & time_preference_2==1 & time_preference_10==1 & time_preference_14==1 & time_preference_16==1 ~ 16,
      time_preference_1==2 & time_preference_2==1 & time_preference_10==1 & time_preference_14==1 & time_preference_16==2 ~ 15,
      time_preference_1==2 & time_preference_2==1 & time_preference_10==1 & time_preference_14==2 & time_preference_15==1 ~ 14,
      time_preference_1==2 & time_preference_2==1 & time_preference_10==1 & time_preference_14==2 & time_preference_15==2 ~ 13,
      time_preference_1==2 & time_preference_2==1 & time_preference_10==2 & time_preference_11==1 & time_preference_13==1 ~ 12,
      time_preference_1==2 & time_preference_2==1 & time_preference_10==2 & time_preference_11==1 & time_preference_13==2 ~ 11,
      time_preference_1==2 & time_preference_2==1 & time_preference_10==2 & time_preference_11==2 & time_preference_12==1 ~ 10,
      time_preference_1==2 & time_preference_2==1 & time_preference_10==2 & time_preference_11==2 & time_preference_12==2 ~ 9,
      time_preference_1==2 & time_preference_2==2 & time_preference_3==1 & time_preference_7==1 & time_preference_8==1 ~ 8,
      time_preference_1==2 & time_preference_2==2 & time_preference_3==1 & time_preference_7==1 & time_preference_8==2 ~ 7,
      time_preference_1==2 & time_preference_2==2 & time_preference_3==1 & time_preference_7==2 & time_preference_9==1 ~ 6,
      time_preference_1==2 & time_preference_2==2 & time_preference_3==1 & time_preference_7==2 & time_preference_9==2 ~ 5,
      time_preference_1==2 & time_preference_2==2 & time_preference_3==2 & time_preference_4==1 & time_preference_6==1 ~ 4,
      time_preference_1==2 & time_preference_2==2 & time_preference_3==2 & time_preference_4==1 & time_preference_6==2 ~ 3,
      time_preference_1==2 & time_preference_2==2 & time_preference_3==2 & time_preference_4==2 & time_preference_5==1 ~ 2,
      time_preference_1==2 & time_preference_2==2 & time_preference_3==2 & time_preference_4==2 & time_preference_5==2 ~ 1,
      TRUE ~ NA_real_
    )
  )


# 	demographic and background characteristics ------------------------

cy <- cy %>%
  mutate(
    ethnicity_han       = as.integer(ethnicity == 1),
    hukou_urban         = as.integer(hukou == 1),
    siblings_total      = number_bro_younger + number_bro_older + number_sis_younger + number_sis_older,
    father_edu_hsabove  = as.integer(father_education > 3),
    mother_edu_hsabove  = as.integer(mother_education > 3),
    religion_religious  = as.integer(religion > 1),
    hs_track_science    = as.integer(hs_track == 1)
  )

## indicator of university
cy <- cy %>% mutate(university_elite = as.integer(university == 1))

## recode birth_year
cy <- cy %>%
  mutate(
    birth_year = case_match(
      birth_year,
      1 ~ 1990, 2 ~ 1991, 3 ~ 1992, 4 ~ 1993, 5 ~ 1994, 6 ~ 1995, 7 ~ 1996, 8 ~ 1997,
      9 ~ 1998, 10 ~ 1999, 11 ~ 2000, 12 ~ 2001, 13 ~ 2002, 14 ~ 2003, 15 ~ 2004, 16 ~ 2005,
      .default = birth_year  # leaves unmatched/missing values unchanged by default; dentical to Stata
    )
  )

## indicator of coastal provinces
coastal_provinces <- c("Beijing", "Fujian", "Guangdong", "Hainan", "Hebei", "Jiangsu",
                       "Shandong", "Shanghai", "Tianjin", "Zhejiang", "Non-mainland")

cy <- cy %>%
  mutate(
    birthplace_coastal = as.integer(birthplace_province %in% coastal_provinces),
    residence_coastal   = as.integer(residence_province %in% coastal_provinces)
  )

## indicator of english credentials
cy <- cy %>%
  mutate(
    domestic_english_atleast4 = as.integer(english_qual_domestic > 1),
    foreign_english_yes        = as.integer(english_qual_foreign < 3)
  )

## indicator of parents' work sector
cy <- cy %>%
  mutate(
    work_father_govt = as.integer(father_work <= 4),
    work_mother_govt = as.integer(mother_work <= 4)
  )

## indicator of travel experience
cy <- cy %>% mutate(travel_foreign_yes = as.integer(travel_foreign > 1))


# 	treatment related indicators ---------------------------------------

## indicators and dummies

cy <- cy %>%
  mutate(
    treatment_control = as.integer(treatment_newsletter == 0 & treatment_vpn == 0 & vpn_current_paid_user == 0),
    treatment_vpnonly  = as.integer(treatment_newsletter == 0 & treatment_vpn == 1 & vpn_current_paid_user == 0),
    treatment_nlonly   = as.integer(treatment_newsletter == 1 & treatment_vpn == 0 & vpn_current_paid_user == 0),
    treatment_vpnnl    = as.integer(treatment_newsletter == 1 & treatment_vpn == 1 & vpn_current_paid_user == 0),
    treatment_user     = as.integer(vpn_current_paid_user == 1)
  )

cy <- cy %>%
  mutate(
    treatment_master = case_when(
      treatment_control == 1 ~ 1,
      treatment_vpnonly == 1 ~ 2,
      treatment_nlonly  == 1 ~ 3,
      treatment_vpnnl   == 1 ~ 4,
      treatment_user    == 1 ~ 5,
      TRUE ~ NA_real_
    )
  )

cy <- cy %>%
  mutate(
    treatment_master_lbl = factor(
      treatment_master,
      levels = 1:5,
      labels = c("Group-C", "Group-A", "Group-CE", "Group-AE", "Existing users")
    )
  )

cy <- cy %>%
  mutate(
    treatment_main = case_when(
      treatment_control == 1 | treatment_nlonly == 1 ~ 1,
      treatment_vpnonly == 1 ~ 2,
      treatment_vpnnl   == 1 ~ 3,
      treatment_user    == 1 ~ 4,
      TRUE ~ NA_real_
    ),
    treatment_main_lbl = factor(
      treatment_main,
      levels = 1:4,
      labels = c("Control", "Access", "Access + Encour.", "Existing users")
    )
  )

## treatment interactions
cy <- cy %>% mutate(treatvpnXnl = treatment_vpn * treatment_newsletter)


# 	vpn adoption status -------------------------------------------------
cy <- cy %>% mutate(vpn_adopted = as.integer(!is.na(vpn_date_adoption)))
cy <- cy %>% mutate(vpn_adopted = ifelse(!is.na(vpn_date_adoption), 1L, NA_integer_))

cy <- cy %>%
  mutate(
    vpn_adopted = ifelse(
      is.na(vpn_adopted) & (treatment_master == 2 | treatment_master == 4),
      0,
      vpn_adopted
    )
  )


# 	roommates' vpn usage -------------------------------------------------

## generate social learning primitives
cy <- cy %>%
  mutate(
    vpn_roommate_existing = vpn_roommate_usage_w1,
    vpn_roommate_new_w2   = vpn_roommate_usage_w2 - vpn_roommate_usage_w1,
    vpn_roommate_new_w3   = vpn_roommate_usage_w3 - vpn_roommate_usage_w1
  )

## top code at 2
cy <- cy %>%
  mutate(
    vpn_roommate_existing = ifelse(vpn_roommate_existing > 2, 2, vpn_roommate_existing),
    vpn_roommate_new_w2   = ifelse(vpn_roommate_new_w2 < 0, 0, vpn_roommate_new_w2),
    vpn_roommate_new_w2   = ifelse(vpn_roommate_new_w2 > 2, 2, vpn_roommate_new_w2),
    vpn_roommate_new_w3   = ifelse(vpn_roommate_new_w3 < 0, 0, vpn_roommate_new_w3),
    vpn_roommate_new_w3   = ifelse(vpn_roommate_new_w3 > 2, 2, vpn_roommate_new_w3)
  )

## finalize primitives
cy <- cy %>%
  mutate(
    soclearning_ownaccess  = as.integer(treatment_master >= 4),
    soclearning_rm_new_w2  = vpn_roommate_new_w2,
    soclearning_rm_new_w3  = vpn_roommate_new_w3,
    soclearning_rm_existing = vpn_roommate_existing,
    soclearning_ownXnew_w2 = soclearning_ownaccess * soclearning_rm_new_w2,
    soclearning_ownXnew_w3 = soclearning_ownaccess * soclearning_rm_new_w3
  )


# 	list experiment -------------------------------------------------------

## list experiment: indicator for provision of veil
cy <- cy %>% mutate(list_exp_veiled = 1 - list_exp_direct)

## list experiment: calculated overall count for direct group

list_exp_stems <- c("trust")

for (cc in list_exp_stems) {
  direct_count_col <- paste0("list_", cc, "_direct_count")
  direct_yes_col   <- paste0("list_", cc, "_direct_yes")
  veiled_count_col <- paste0("list_", cc, "_veiled_count")

  countd_direct_col <- paste0("list_countd_", cc, "_direct")
  count_col         <- paste0("list_count_", cc)

  cy[[countd_direct_col]] <- cy[[direct_count_col]] + cy[[direct_yes_col]]

  cy[[count_col]] <- ifelse(
    cy$list_exp_direct == 1,
    cy[[countd_direct_col]],
    ifelse(cy$list_exp_direct == 0, cy[[veiled_count_col]], NA_real_)
  )
}

### Prepare data: generate az-scores
#Uses the andersonz() function defined
#previously, passing the corresponding R character vector (e.g. var_x, defined
#in variable_lists.R) as the `varlist` argument, and assigning the result back
#to `cy`
#
#Each call below follows the pattern: cy <- andersonz(cy, varlist = var_x,
#generate = "az_y")


# 	A: Beliefs and attitudes regarding media -----------------------------

## A.2: Purchase of censorship circumvention tools
cy <- andersonz(cy, varlist = var_vpn_purchase_reg_w3, generate = "az_var_vpnpurchase")

## A.3: Valuation of access to foreign media outlets
cy <- andersonz(cy, varlist = var_media_valuation_reg_w1, generate = "az_belief_media_value_w1")
cy <- andersonz(cy, varlist = var_media_valuation_reg_w2, generate = "az_belief_media_value_w2")
cy <- andersonz(cy, varlist = var_media_valuation_reg_w3, generate = "az_belief_media_value_w3")

## A.4: Trust in media outlets
cy <- andersonz(cy, varlist = var_media_trust_reg_w1, generate = "az_belief_media_trust_w1")
cy <- andersonz(cy, varlist = var_media_trust_reg_w2, generate = "az_belief_media_trust_w2")
cy <- andersonz(cy, varlist = var_media_trust_reg_w3, generate = "az_belief_media_trust_w3")

## A.5: Calibration of news outlets' level of censorship and biases
cy <- andersonz(cy, varlist = var_percmediabias_reg_ce_cn_w2, generate = "az_belief_media_cens_cn")
cy <- andersonz(cy, varlist = var_percmediabias_reg_ce_us_w2, generate = "az_belief_media_cens_us")
cy <- andersonz(cy, varlist = var_percmediabias_reg_di_cn_w2, generate = "az_belief_media_bias_cn")
cy <- andersonz(cy, varlist = var_percmediabias_reg_di_us_w2, generate = "az_belief_media_bias_us")

## A.6: Justification of media censorship
cy <- andersonz(cy, varlist = var_censor_justif_reg_w1, generate = "az_belief_media_justif_w1")
cy <- andersonz(cy, varlist = var_censor_justif_reg_w2, generate = "az_belief_media_justif_w2")
cy <- andersonz(cy, varlist = var_censor_justif_reg_w3, generate = "az_belief_media_justif_w3")
cy <- andersonz(cy, varlist = var_censor_justif_reg_fp_w1, generate = "az_belief_media_justif_pf_w1")
cy <- andersonz(cy, varlist = var_censor_justif_reg_fp_w2, generate = "az_belief_media_justif_pf_w2")
cy <- andersonz(cy, varlist = var_censor_justif_reg_fp_w3, generate = "az_belief_media_justif_pf_w3")

## A.7: Belief regarding drivers of media censorship
cy <- andersonz(cy, varlist = var_censor_driver_reg_dom_w2, generate = "az_belief_media_dri_dom_w2")
cy <- andersonz(cy, varlist = var_censor_driver_reg_for_w2, generate = "az_belief_media_dri_for_w2")


# 	B: Knowledge -----------------------------------------------------------

## B.1: Current news events covered in the demand treatment
cy <- andersonz(cy, varlist = var_knowledge_news_reg_quiz, generate = "az_knowledge_news_quiz")

## B.2: Current news events not covered in the demand treatment
cy <- andersonz(cy, varlist = var_knowledge_news_reg_cen_w2, generate = "az_knowledge_news_censored")
cy <- andersonz(cy, varlist = var_knowledge_news_reg_unc_w2, generate = "az_knowledge_news_uncensor")
cy <- andersonz(cy, varlist = var_knowledge_news_reg_cen_w1, generate = "az_knowledge_news_censored_w1")
cy <- andersonz(cy, varlist = var_knowledge_news_reg_unc_w1, generate = "az_knowledge_news_uncensor_w1")

## B.3: Awareness of notable figures
cy <- andersonz(cy, varlist = var_knowledge_people_reg_tocens, generate = "az_knowledge_people_tocens")
cy <- andersonz(cy, varlist = var_knowledge_people_reg_censor, generate = "az_knowledge_people_censor")
cy <- andersonz(cy, varlist = var_knowledge_people_reg_uncens, generate = "az_knowledge_people_uncens")
cy <- andersonz(cy, varlist = var_knowledge_people_censor_w1, generate = "az_knowledge_people_censor_w1")
cy <- andersonz(cy, varlist = var_knowledge_people_uncens_w1, generate = "az_knowledge_people_uncens_w1")

## B.4: Awareness of protest events
cy <- andersonz(cy, varlist = var_knowledge_prot_reg_chi_w2, generate = "az_knowledge_protest_china")
cy <- andersonz(cy, varlist = var_knowledge_prot_reg_for_w2, generate = "az_knowledge_protest_forei")
cy <- andersonz(cy, varlist = var_knowledge_prot_reg_chi_w1, generate = "az_knowledge_protest_china_w1")

## B.5: Meta-knowledge
cy <- andersonz(cy, varlist = var_knowledge_meta_reg_w1, generate = "az_knowledge_meta_w1")
cy <- andersonz(cy, varlist = var_knowledge_meta_reg_w2, generate = "az_knowledge_meta_w2")
cy <- andersonz(cy, varlist = var_knowledge_meta_reg_w3, generate = "az_knowledge_meta_w3")


# 	C: Economic beliefs ------------------------------------------------------

## C.1: Belief on economic performance in China
cy <- andersonz(cy, varlist = var_econ_guess_reg_cn_perf_w1, generate = "az_belief_econ_perf_cn_w1")
cy <- andersonz(cy, varlist = var_econ_guess_reg_cn_perf_w2, generate = "az_belief_econ_perf_cn_w2")
cy <- andersonz(cy, varlist = var_econ_guess_reg_cn_perf_w3, generate = "az_belief_econ_perf_cn_w3")

## C.2: Confidence on guesses regarding economic performance in China
cy <- andersonz(cy, varlist = var_econ_guess_reg_cn_conf_w1, generate = "az_belief_econ_conf_cn_w1")
cy <- andersonz(cy, varlist = var_econ_guess_reg_cn_conf_w2, generate = "az_belief_econ_conf_cn_w2")
cy <- andersonz(cy, varlist = var_econ_guess_reg_cn_conf_w3, generate = "az_belief_econ_conf_cn_w3")

## C.3: Belief on economic performance in the US
cy <- andersonz(cy, varlist = var_econ_guess_reg_us_perf_w2, generate = "az_belief_econ_perf_us_w2")
cy <- andersonz(cy, varlist = var_econ_guess_reg_us_perf_w3, generate = "az_belief_econ_perf_us_w3")

## C.4: Confidence on guesses regarding economic performance in the US
cy <- andersonz(cy, varlist = var_econ_guess_reg_us_conf_w2, generate = "az_belief_econ_conf_us_w2")
cy <- andersonz(cy, varlist = var_econ_guess_reg_us_conf_w3, generate = "az_belief_econ_conf_us_w3")


# 	D: Political attitudes -----------------------------------------------------

## D.1: Demand for institutional change
cy <- andersonz(cy, varlist = var_demand_change_reg_w1, generate = "az_belief_instchange_w1")
cy <- andersonz(cy, varlist = var_demand_change_reg_w2, generate = "az_belief_instchange_w2")
cy <- andersonz(cy, varlist = var_demand_change_reg_w3, generate = "az_belief_instchange_w3")

## D.2: Trust in institutions
cy <- andersonz(cy, varlist = var_trust_inst_reg_govt_w1, generate = "az_belief_trust_govt_w1")
cy <- andersonz(cy, varlist = var_trust_inst_reg_govt_w2, generate = "az_belief_trust_govt_w2")
cy <- andersonz(cy, varlist = var_trust_inst_reg_govt_w3, generate = "az_belief_trust_govt_w3")
cy <- andersonz(cy, varlist = var_trust_inst_reg_foreign_w1, generate = "az_belief_trust_foreign_w1")
cy <- andersonz(cy, varlist = var_trust_inst_reg_foreign_w2, generate = "az_belief_trust_foreign_w2")
cy <- andersonz(cy, varlist = var_trust_inst_reg_foreign_w3, generate = "az_belief_trust_foreign_w3")
cy <- andersonz(cy, varlist = var_trust_inst_reg_copo_w2, generate = "az_belief_trust_copo_w2")

## D.3: Evaluation of government's performance
cy <- andersonz(cy, varlist = var_eval_govt_reg_w1, generate = "az_belief_evalgovt_w1")
cy <- andersonz(cy, varlist = var_eval_govt_reg_w2, generate = "az_belief_evalgovt_w2")
cy <- andersonz(cy, varlist = var_eval_govt_reg_w3, generate = "az_belief_evalgovt_w3")

## D.4: Performance evaluation criteria
cy <- andersonz(cy, varlist = var_eval_criteria_reg_w2, generate = "az_belief_evalcrit_w2")

## D.5: Evaluation of severity of socioeconomic issues
cy <- andersonz(cy, varlist = var_severity_reg_w2, generate = "az_belief_severity_w2")

## D.6: Evaluation of democracy and human rights protection in China
cy <- andersonz(cy, varlist = var_democracy_reg_fp_w1, generate = "az_belief_democracy_fp_w1")
cy <- andersonz(cy, varlist = var_democracy_reg_fp_w2, generate = "az_belief_democracy_fp_w2")
cy <- andersonz(cy, varlist = var_democracy_reg_fp_w3, generate = "az_belief_democracy_fp_w3")
cy <- andersonz(cy, varlist = var_democracy_reg_w2, generate = "az_belief_democracy_w2")

## D.7: Justification of controversial policies and issues
cy <- andersonz(cy, varlist = var_contro_justi_reg_policy_w2, generate = "az_belief_justify_policy_w2")
cy <- andersonz(cy, varlist = var_contro_justi_reg_liberal_w2, generate = "az_belief_justify_liberal_w2")

## D.8: Willingness to act
cy <- andersonz(cy, varlist = var_willing_fight_reg_w1, generate = "az_belief_willing_w1")
cy <- andersonz(cy, varlist = var_willing_fight_reg_w2, generate = "az_belief_willing_w2")
cy <- andersonz(cy, varlist = var_willing_fight_reg_w3, generate = "az_belief_willing_w3")

## D.9: Interest in politics and economics
cy <- andersonz(cy, varlist = var_interest_reg_w2, generate = "az_belief_interest_w2")


# 	E: behaviors ---------------------------------------------------------------

## E.1: Social interactions
cy <- andersonz(cy, varlist = var_socialinteract_reg_w1, generate = "az_var_socialinteract_w1")
cy <- andersonz(cy, varlist = var_socialinteract_reg_w2, generate = "az_var_socialinteract_w2")
cy <- andersonz(cy, varlist = var_socialinteract_reg_w3, generate = "az_var_socialinteract_w3")

## E.2: Participation
cy <- andersonz(cy, varlist = var_polparticipation_reg_pf_w1, generate = "az_var_polparticipation_w1")
cy <- andersonz(cy, varlist = var_polparticipation_reg_pf_w2, generate = "az_var_polparticipation_w2")
cy <- andersonz(cy, varlist = var_polparticipation_reg_pf_w3, generate = "az_var_polparticipation_w3")


# 	F: Demographics, background characteristics, and fundamental preferences ---

## F.1: Personal characteristics
cy <- andersonz(cy, varlist = var_demog_reg_personal, generate = "az_demographics_personal")

## F.2: Educational background
cy <- andersonz(cy, varlist = var_demog_reg_education, generate = "az_demographics_education")

## F.3: English ability and oversea travel experiences
cy <- andersonz(cy, varlist = var_demog_reg_english, generate = "az_demographics_english")
cy <- andersonz(cy, varlist = var_demog_reg_travel, generate = "az_demographics_travel")

## F.4: Household characteristics
cy <- andersonz(cy, varlist = var_demog_reg_household, generate = "az_demographics_household")

## F.5: Fundamental preferences
cy <- andersonz(cy, varlist = var_preference_risk, generate = "az_preference_risk")
cy <- andersonz(cy, varlist = var_preference_time, generate = "az_preference_time")
cy <- andersonz(cy, varlist = var_preference_altruism, generate = "az_preference_altruism")
cy <- andersonz(cy, varlist = var_preference_reciprocity, generate = "az_preference_reciprocity")


# 	X: Overall effect index ----------------------------------------------------

## flip variables to generate overall impact index
for (i in 1:3) {
  bias_for_col      <- paste0("bias_foreign_w", i)
  econ_perf_col     <- paste0("az_belief_econ_perf_cn_w", i)
  trust_govt_col    <- paste0("az_belief_trust_govt_w", i)
  evalgovt_col      <- paste0("az_belief_evalgovt_w", i)
  demo_col          <- paste0("importance_live_in_demo_w", i)

  cy[[paste0("bias_foreign_r_w", i)]]              <- 10 - cy[[bias_for_col]]
  cy[[paste0("az_belief_econ_perf_cn_r_w", i)]]     <- 0 - cy[[econ_perf_col]]
  cy[[paste0("az_belief_trust_govt_r_w", i)]]       <- 0 - cy[[trust_govt_col]]
  cy[[paste0("az_belief_evalgovt_r_w", i)]]         <- 0 - cy[[evalgovt_col]]
  cy[[paste0("importance_live_in_demo_r_w", i)]]    <- 10 - cy[[demo_col]]
}

cy <- andersonz(cy, varlist = var_overall_w1, generate = "az_overall_w1")
cy <- andersonz(cy, varlist = var_overall_w2, generate = "az_overall_w2")
cy <- andersonz(cy, varlist = var_overall_w3, generate = "az_overall_w3")

## five main categories: wave 1
cy <- andersonz(cy,
                varlist = c("info_freq_website_for_w1", "az_belief_media_value_w1", "az_belief_media_trust_w1",
                            "az_belief_media_justif_pf_w1", "bias_domestic_w1", "bias_foreign_r_w1"),
                generate = "az_overall_a_w1")

cy <- andersonz(cy,
                varlist = c("news_perccor_cen_w1", "news_perccor_unc_w1", "protest_pcheard_china_w1",
                            "protest_pcheard_foreign_w1", "az_knowledge_meta_w1"),
                generate = "az_overall_b_w1")

cy <- andersonz(cy,
                varlist = c("az_belief_econ_perf_cn_r_w1", "az_belief_econ_conf_cn_w1"),
                generate = "az_overall_c_w1")

cy <- andersonz(cy,
                varlist = c("az_belief_instchange_w1", "az_belief_trust_govt_r_w1", "az_belief_trust_foreign_w1",
                            "az_belief_evalgovt_r_w1", "importance_live_in_demo_r_w1", "az_belief_willing_w1"),
                generate = "az_overall_d_w1")

cy <- andersonz(cy,
                varlist = c("az_var_socialinteract_w1", "az_var_polparticipation_w1", "plan_grad_foreignmaster_w1",
                            "cloc_for_w1", "stock_participation_w1"),
                generate = "az_overall_e_w1")

## five main categories: wave 2
cy <- andersonz(cy,
                varlist = c("info_freq_website_for_w2", "az_belief_media_value_w2", "az_belief_media_trust_w2",
                            "az_belief_media_justif_pf_w2", "bias_domestic_w2", "bias_foreign_r_w2"),
                generate = "az_overall_a_w2")

cy <- andersonz(cy,
                varlist = c("news_perccor_cen_w2", "news_perccor_unc_w2", "protest_pcheard_china_w2",
                            "protest_pcheard_foreign_w2", "az_knowledge_meta_w2"),
                generate = "az_overall_b_w2")

cy <- andersonz(cy,
                varlist = c("az_belief_econ_perf_cn_r_w2", "az_belief_econ_conf_cn_w2"),
                generate = "az_overall_c_w2")

cy <- andersonz(cy,
                varlist = c("az_belief_instchange_w2", "az_belief_trust_govt_r_w2", "az_belief_trust_foreign_w2",
                            "az_belief_evalgovt_r_w2", "importance_live_in_demo_r_w2", "az_belief_willing_w2"),
                generate = "az_overall_d_w2")

cy <- andersonz(cy,
                varlist = c("az_var_socialinteract_w2", "az_var_polparticipation_w2", "plan_grad_foreignmaster_w2",
                            "cloc_for_w2", "stock_participation_w2"),
                generate = "az_overall_e_w2")

## five main categories: wave 3
cy <- andersonz(cy,
                varlist = c("info_freq_website_for_w3", "az_belief_media_value_w3", "az_belief_media_trust_w3",
                            "az_belief_media_justif_pf_w3", "bias_domestic_w3", "bias_foreign_r_w3"),
                generate = "az_overall_a_w3")

cy <- andersonz(cy,
                varlist = c("news_perccor_cen_w3", "news_perccor_unc_w3", "protest_pcheard_china_w3",
                            "protest_pcheard_foreign_w3", "az_knowledge_meta_w3"),
                generate = "az_overall_b_w3")

cy <- andersonz(cy,
                varlist = c("az_belief_econ_perf_cn_r_w3", "az_belief_econ_conf_cn_w3"),
                generate = "az_overall_c_w3")

cy <- andersonz(cy,
                varlist = c("az_belief_instchange_w3", "az_belief_trust_govt_r_w3", "az_belief_trust_foreign_w3",
                            "az_belief_evalgovt_r_w3", "importance_live_in_demo_r_w3", "az_belief_willing_w3"),
                generate = "az_overall_d_w3")

cy <- andersonz(cy,
                varlist = c("az_var_socialinteract_w3", "az_var_polparticipation_w3", "plan_grad_foreignmaster_w3",
                            "cloc_for_w3", "stock_participation_w3"),
                generate = "az_overall_e_w3")

## mega z-scores: wave 3
cy <- andersonz(cy,
                varlist = c("az_overall_a_w3", "az_overall_b_w3", "az_overall_c_w3", "az_overall_d_w3", "az_overall_e_w3"),
                generate = "az_overall_all_w3")

### Prepare data: generate heterogeneity cuts

cy <- cy %>%
  mutate(
    # demographics (category F)
    h_gender                       = as.integer(gender == 1),
    h_birth_year                   = as.integer(birth_year < 1996),
    h_ethnicity_han                 = as.integer(ethnicity_han == 1),
    h_birthplace_coastal            = as.integer(birthplace_coastal == 1),
    h_residence_coastal             = as.integer(residence_coastal == 1),
    h_hukou_urban                   = as.integer(hukou_urban == 1),
    h_religion_religious            = as.integer(religion_religious == 1),
    h_ccp_member                    = as.integer(ccp_member == 1),
    h_university_elite               = as.integer(university_elite == 1),
    h_hs_track_science               = as.integer(hs_track_science == 1),
    h_department_ssh                 = as.integer(department_ssh == 1),
    h_domestic_english_atleast4      = as.integer(domestic_english_atleast4 == 1),
    h_foreign_english_yes            = as.integer(foreign_english_yes == 1),
    h_travel_hktaiwan                = as.integer(travel_hktaiwan == 1),
    h_travel_foreign_yes             = as.integer(travel_foreign_yes == 1),
    h_siblings_total                 = as.integer(siblings_total > 0),
    h_father_edu_hsabove             = as.integer(father_edu_hsabove == 1),
    h_work_father_govt               = as.integer(work_father_govt == 1),
    h_father_ccp                     = as.integer(father_ccp == 1),
    h_mother_edu_hsabove             = as.integer(mother_edu_hsabove == 1),
    h_work_mother_govt               = as.integer(work_mother_govt == 1),
    h_mother_ccp                     = as.integer(mother_ccp == 1),
    h_hh_income                      = as.integer(hh_income >= 75000),
    h_az_preference_risk             = as.integer(az_preference_risk > 0),
    h_az_preference_time             = as.integer(az_preference_time > 0),
    h_az_preference_altruism         = as.integer(az_preference_altruism > 0),
    h_az_preference_reciprocity      = as.integer(az_preference_reciprocity > 0),

    # baseline outcomes
    h_az_overall_a_w1                = as.integer(az_overall_a_w1 > 0),
    h_az_overall_b_w1                = as.integer(az_overall_b_w1 > 0),
    h_az_overall_c_w1                = as.integer(az_overall_c_w1 > 0),
    h_az_overall_d_w1                = as.integer(az_overall_d_w1 > 0),
    h_az_overall_e_w1                = as.integer(az_overall_e_w1 > 0),
    h_az_belief_media_value_w1       = as.integer(az_belief_media_value_w1 > 0),
    h_az_belief_media_trust_w1       = as.integer(az_belief_media_trust_w1 > 0),
    h_az_knowledge_news_cens_w1      = as.integer(az_knowledge_news_censored_w1 > 0),
    h_az_knowledge_news_unce_w1      = as.integer(az_knowledge_news_uncensor_w1 > 0),
    h_az_knowledge_pp_censor_w1      = as.integer(az_knowledge_people_censor_w1 > 0),
    h_az_knowledge_pp_uncens_w1      = as.integer(az_knowledge_people_uncens_w1 > 0),
    h_az_knowledge_pr_china_w1       = as.integer(az_knowledge_protest_china_w1 > 0),
    h_az_belief_trust_govt_w1        = as.integer(az_belief_trust_govt_w1 > 0)
  )


### Prepare data: generate variable labels

var_labels <- c(
  # A. Beliefs and attitudes regarding media
  l_wtp_vpn                       = "WTP for uncensored Internet access ($/month)",
  l_added_value_foreign_media     = "Value added of foreign media access",
  l_az_belief_media_value         = "Valuation of access to foreign media outlets",

  l_trust_media_dom_state         = "Distrust in domestic state-owned media",
  l_trust_media_dom_private       = "Distrust in domestic privately-owned media",
  l_trust_media_foreign           = "Trust in foreign media",
  l_az_belief_media_trust         = "Trust in non-domestic media outlets",

  l_bias_domestic                 = "Degree of censorship on domestic news outlets",
  l_bias_foreign                  = "Degree of censorship on foreign news outlets",
  l_bias_dom_govt_policy_t1       = "Domestic cens. driven by govt. policies",
  l_bias_dom_firm_interest_t1     = "Domestic cens. driven by corp. interest",
  l_bias_dom_media_pref_t1        = "Domestic cens. driven by media\u2019s ideology",
  l_bias_dom_reader_demand_t1     = "Domestic cens. driven by readers\u2019 demand",
  l_bias_for_govt_policy_t1       = "Foreign cens. driven by govt. policies",
  l_bias_for_firm_interest_t1     = "Foreign cens. driven by corp. interest",
  l_bias_for_media_pref_t1        = "Foreign cens. driven by media\u2019s ideology",
  l_bias_for_reader_demand_t1     = "Foreign cens. driven by readers\u2019 demand",

  l_distneutral_cn_neg_cn         = "Bias: Chinese media on neg. news in China",
  l_distneutral_cn_pos_cn         = "Bias: Chinese media on pos. news in China",
  l_distneutral_us_neg_cn         = "Bias: Chinese media on neg. news in US",
  l_distneutral_us_pos_cn         = "Bias: Chinese media on pos. news in US",
  l_distneutral_cn_neg_us         = "Bias: US media on neg. news in China",
  l_distneutral_cn_pos_us         = "Bias: US media on pos. news in China",
  l_distneutral_us_neg_us         = "Bias: US media on neg. news in US",
  l_distneutral_us_pos_us         = "Bias: US media on pos. news in US",

  l_bias_cn_neg_cn_cens           = "Censorship: Chinese media on neg. news in China",
  l_bias_cn_pos_cn_cens           = "Censorship: Chinese media on pos. news in China",
  l_bias_us_neg_cn_cens           = "Censorship: Chinese media on neg. news in US",
  l_bias_us_pos_cn_cens           = "Censorship: Chinese media on pos. news in US",
  l_bias_cn_neg_us_cens           = "Censorship: US media on neg. news in China",
  l_bias_cn_pos_us_cens           = "Censorship: US media on pos. news in China",
  l_bias_us_neg_us_cens           = "Censorship: US media on neg. news in US",
  l_bias_us_pos_us_cens           = "Censorship: US media on pos. news in US",

  l_censor_just_dom_economic      = "Unjustified: censoring economic news",
  l_censor_just_dom_political     = "Unjustified: censoring political news",
  l_censor_just_dom_social        = "Unjustified: censoring social news",
  l_censor_just_for               = "Unjustified: censoring foreign news",
  l_censor_just_porn              = "Unjustified: censoring pornography",
  l_az_belief_media_justif_pf     = "Censorship is unjustified",

  # B. Knowledge
  l_panamapapers                  = "Panama Papers",
  l_tenyearshk                    = "HK independence",
  l_stockcrash                    = "2016 stock mkt crash",
  l_economistcensor               = "Censoring Economist",
  l_coalprod                      = "Steel prod. & pollution",
  l_trumpchina                    = "Trump trademark in China",
  l_xiaojianhua                   = "Jianhua Xiao kidnap",
  l_xijiangcar                    = "Tracking Xinjiang cars",
  l_chinanorway                   = "China-Norway relations",
  l_womenrights                   = "Women rights activisits",
  l_hkceelection                  = "HK CE election",
  l_news_perccor_cen               = "% quizzes answered correctly: sensitive",
  l_news_perccor_unc               = "% quizzes answered correctly: non-sensitive",

  l_people_puzhiqiang              = "Aware of Zhiqiang Pu",
  l_people_lizehou                 = "Aware of Zehou Li",
  l_people_huangzhifeng            = "Aware of Joshua Wong",
  l_people_chenguangcheng          = "Aware of Guangcheng Cheng",
  l_people_lixiaolin               = "Aware of Xiaolin Li",
  l_people_renzhiqiang             = "Aware of Zhiqiang Ren",
  l_people_maoyushi                = "Aware of Yushi Mao",
  l_people_honghuang               = "Aware of Huang Hong",
  l_people_liuqiangdong            = "Aware of Qiangdong Liu",
  l_people_jialequn                = "Aware of Lequn Jia",

  l_protest_2014_europe_square    = "Aware of 2014 Ukrainian Euromaidan Revolution",
  l_protest_2014_sun_flower       = "Aware of 2014 Taiwan Sunflower Stud. Movement",
  l_protest_2010_arabic_spring    = "Aware of 2010 Arab Spring",
  l_protest_2014_crimea_vote      = "Aware of 2014 Crimean Status Referendum",
  l_protest_2012_hk_curriculum    = "Aware of 2012 HK Anti-National Curr. Movement",
  l_protest_2010_catal_indep      = "Aware of 2010 Catalonian Indep. Movement",
  l_protest_2014_umbrella         = "Aware of 2014 HK Umbrella Revolution",
  l_protest_2011_tmrw_parade      = "Aware of 2011 Tomorrow Revolution",
  l_protest_2016_mongkok_riot     = "Aware of 2016 HK Mong Kok Revolution",
  l_protest_pcheard_china         = "Awareness of protests in Greater China",
  l_protest_pcheard_foreign       = "Awareness of foreign protests",

  l_familiar_china_issues_self     = "Informedness of issues in China",
  l_familiar_china_others          = "Greater informedness than peers",
  l_az_knowledge_meta               = "Self-assessment of knowledge level",

  # C. Economic beliefs
  l_guess_gdp_growth_china         = "Guess on GDP growth rate in China",
  l_guess_stock_index_sh           = "Guess on year-end SSCI",
  l_az_belief_econ_perf_cn         = "Optimistic belief of Chinese economy",
  l_az_belief_econ_conf_cn         = "Confidence of guesses on Chinese economy",

  l_guess_gdp_growth_china_con     = "Confidence of China GDP guess",
  l_guess_stock_index_sh_con       = "Confidence of SSCI guess",

  l_guess_gdp_growth_us            = "Guess on GDP growth rate in US",
  l_guess_stock_index_dj           = "Guess on year-end DJI",

  l_guess_gdp_growth_us_con        = "Confidence of US GDP guess",
  l_guess_stock_index_dj_con       = "Confidence of DJI guess",

  # D. Political attitudes
  l_inst_change_econ               = "Economic institution needs changes",
  l_inst_change_poli               = "Political institution needs changes",
  l_az_belief_instchange            = "Demand for institutional change",

  l_trust_central_govt             = "Trust in central govt. of China",
  l_trust_provincial_govt          = "Trust in provincial govt. of China",
  l_trust_local_govt                = "Trust in local govt. of China",
  l_az_belief_trust_govt            = "Trust in Chinese govt.",
  l_trust_japan_govt                = "Trust in central govt. of Japan",
  l_trust_us_govt                   = "Trust in federal govt. of US",
  l_az_belief_trust_foreign         = "Trust in foreign govt.",
  l_trust_court                     = "Trust in court",
  l_trust_police                    = "Trust in police",
  l_trust_ngo                       = "Trust in NGOs",
  l_trust_financial_domestic        = "Trust in domestic financial inst.",
  l_trust_financial_foreign         = "Trust in foreign financial inst.",

  l_eval_govt_economic              = "Satisfaction of economic dev.",
  l_eval_govt_dom_politics          = "Satisfaction of domestic politics",
  l_eval_govt_for_relations         = "Satisfaction of diplomatic affairs",
  l_az_belief_evalgovt              = "Satisfaction of govt. performance",

  l_revalgovt_election              = "Eval. importance: universal suffrage",
  l_revalgovt_economy               = "Eval. importance: economic dev.",
  l_revalgovt_equality               = "Eval. importance: income and wealth equality",
  l_revalgovt_ruleoflaw              = "Eval. importance: rule of law",
  l_revalgovt_human_rights           = "Eval. importance: civil and human rights",
  l_revalgovt_freedom_speech         = "Eval. importance: freedom of speech",
  l_revalgovt_global_power           = "Eval. importance: intl. affairs",
  l_revalgovt_fair_history           = "Eval. importance: handle history fairly",

  l_severity_welfare                 = "Severity: social security and welfare",
  l_severity_employment              = "Severity: employments",
  l_severity_pollution                = "Severity: environmental pollution",
  l_severity_inequality              = "Severity: wealth inequality",
  l_severity_corruption              = "Severity: govt. corruption",
  l_severity_dscrm_minority          = "Severity: minority discrimination",

  l_importance_live_in_demo          = "Living in democracy is not important",
  l_china_interest_group             = "China cares interest for masses",
  l_china_rate_democracy             = "Level of democracy in China",
  l_china_rate_humanrights           = "Level of human rights protection",
  l_az_belief_democracy_fp           = "Living in democracy is not important",

  l_justify_minority_policy          = "Justified: minority policies",
  l_justify_reduce_pollution         = "Justified: prod. cut to reduce pollution",
  l_justify_hukou                    = "Justified: migration restrictions",
  l_justify_one_child                = "Justified: one-child policy",
  l_justify_gaokao                   = "Justified: college admission policies",
  l_justify_hongkong_policy          = "Justified: policy towards HK",
  l_justify_taiwan_policy            = "Justified: policy towards Taiwan",
  l_justify_transgene                = "Justified: transgenetic food",
  l_justify_receive_refugee          = "Justified: refusal of DPRK refugees",
  l_justify_soe_privatize            = "Justified: privatization of SOEs",
  l_justify_homo_marriage            = "Justified: legal. of homosexual marriages",
  l_justify_legal_prostitute         = "Justified: legal. of prostitution",
  l_justify_abortion                 = "Justified: abortion",
  l_justify_soft_drug                = "Justified: soft drugs usage",
  l_justify_violence_stability       = "Justified: govt. use of violence",

  l_willing_against_illi_govt        = "Willing to battle illegal govt. acts",
  l_willing_report_mis               = "Willing to report govt. misconduct",
  l_willing_protect_weak             = "Willing to stand up for the weak",
  l_az_belief_willing                 = "Willingness to act",

  l_interest_economic                = "Interest in economics",
  l_interest_politics                = "Interest in politics",

  l_proud_being_chinese              = "Proud of being Chinese",

  l_fear_critic_govt_self            = "Fear to criticize govt. in public",

  # E. Behaviors
  l_info_domestic_website            = "Important info source: domestic websites",
  l_info_foreign_website             = "Important info source: foreign websites",
  l_info_social_media_dom            = "Important info source: domestic social media",
  l_info_word_of_mouth               = "Important info source: word of mouth",
  l_info_social_media_for            = "Important info source: foreign social media",

  l_info_freq_website_for            = "Frequently visit foreign websites for info",

  l_participate_social_protest       = "Participated in social protests",
  l_participate_plan_vote            = "Plan to vote in PCR election",
  l_participate_complain_school      = "Filed complaints to school",
  l_participate_ngo                  = "Participated in NGO activities",
  l_az_var_polparticipation           = "Participation behaviors",

  l_frequency_talk_politic           = "Frequency of discussing poli. with friends",
  l_frequency_persuade_friends       = "Frequency of persuading others",
  l_az_var_socialinteract             = "Social interaction in politics",

  l_plan_grad_gradschool_dom         = "Plan: grad. school in China",
  l_plan_grad_foreignmaster          = "Plan: master degree abroad",
  l_plan_grad_foreignphd             = "Plan: PhD degree abroad",
  l_plan_grad_military               = "Plan: military in China",
  l_plan_grad_work                   = "Plan: work right away",

  l_cp_t1_national_civil             = "Sector pref.: national civil service",
  l_cp_t1_local_civil                = "Sector pref.: local civil service",
  l_cp_t1_military                   = "Sector pref.: military",
  l_cp_t1_chinese_private            = "Sector pref.: private firm in China",
  l_cp_t1_for_firm                   = "Sector pref.: private firm in China",
  l_cp_t1_soe                        = "Sector pref.: SOEs",
  l_cp_t1_institutional              = "Sector pref.: inst. organizations",
  l_cp_t1_entrepreneur               = "Sector pref.: entrepreneurship",

  l_cloc_beijing                     = "Location pref.: Beijing",
  l_cloc_shanghai                    = "Location pref.: Shanghai",
  l_cloc_gzsz                        = "Location pref.: Shenzhen/Guangzhou",
  l_cloc_tjcq                        = "Location pref.: tier 2 cities in central",
  l_cloc_hkmc                        = "Location pref.: HK and Macau",
  l_cloc_taiwan                      = "Location pref.: Taiwan",
  l_cloc_dom                         = "Location pref.: other cities in China",
  l_cloc_for                         = "Location pref.: foreign cities",
  l_stock_participation              = "Currently invested in Chinese stock mkt."
)


### The changes start

####################################
##  Code for Table A.13           ##
####################################

## keep wave 3 non-existing users only
cy <- cy %>% filter(panelmerged_wave3 == 1)
cy <- cy %>% filter(treatment_user != 1)

## transform into dummy indicators: uncensored belief - above median

## transform into dummy indicators: uncensored belief - above median

median_above_vars_w3 <- c(
  "info_foreign_website_w3", "info_freq_website_for_w3",
  "az_belief_media_value_w3", "az_belief_media_trust_w3",
  "bias_domestic_w3", "az_belief_media_justif_w3",
  "news_perccor_cen_w3", "news_perccor_unc_w3",
  "protest_pcheard_china_w3", "protest_pcheard_foreign_w3",
  "az_knowledge_meta_w3", "az_belief_econ_conf_cn_w3",
  "az_belief_econ_perf_us_w3", "az_belief_econ_conf_us_w3",
  "az_belief_instchange_w3", "az_belief_trust_foreign_w3",
  "az_belief_willing_w3", "frequency_talk_politic_w3",
  "frequency_persuade_friends_w3"
)

for (Y in median_above_vars_w3) {
  m_col <- paste0(Y, "_m")
  p_col <- paste0(Y, "_p")

  cutoff <- median(cy[[Y]], na.rm = TRUE)

  cy[[m_col]] <- cutoff

  # Exact Stata replication:
  # Stata numeric missing values are greater than finite values,
  # so (Y > cutoff) evaluates to 1 when Y is missing.
  cy[[p_col]] <- ifelse(
    is.na(cy[[Y]]),
    1L,
    as.integer(cy[[Y]] > cutoff)
  )
}

median_above_vars_w1 <- c(
  "info_foreign_website_w1", "info_freq_website_for_w1",
  "az_belief_media_value_w1", "az_belief_media_trust_w1",
  "bias_domestic_w1", "az_belief_media_justif_w1",
  "news_perccor_cen_w1", "news_perccor_unc_w1",
  "protest_pcheard_china_w1", "protest_pcheard_foreign_w1",
  "az_knowledge_meta_w1", "az_belief_econ_conf_cn_w1",
  "az_belief_instchange_w1", "az_belief_trust_foreign_w1",
  "az_belief_willing_w1", "frequency_talk_politic_w1",
  "frequency_persuade_friends_w1"
)

for (Y in median_above_vars_w1) {
  m_col <- paste0(Y, "_m")
  p_col <- paste0(Y, "_p")

  cutoff <- median(cy[[Y]], na.rm = TRUE)

  cy[[m_col]] <- cutoff

  # Exact Stata replication:
  # missing source values become 1 under (Y > cutoff).
  cy[[p_col]] <- ifelse(
    is.na(cy[[Y]]),
    1L,
    as.integer(cy[[Y]] > cutoff)
  )
}

## transform into dummy indicators: uncensored belief - below median

median_below_vars_w3 <- c(
  "bias_foreign_w3", "importance_live_in_demo_w3",
  "az_belief_econ_perf_cn_w3", "az_belief_trust_govt_w3",
  "az_belief_evalgovt_w3"
)

for (Y in median_below_vars_w3) {
  m_col <- paste0(Y, "_m")
  p_col <- paste0(Y, "_p")

  cutoff <- median(cy[[Y]], na.rm = TRUE)

  cy[[m_col]] <- cutoff

  # Exact Stata replication:
  # Stata missing values are larger than finite values,
  # so (Y < cutoff) evaluates to 0 when Y is missing.
  cy[[p_col]] <- ifelse(
    is.na(cy[[Y]]),
    0L,
    as.integer(cy[[Y]] < cutoff)
  )
}

median_below_vars_w1 <- c(
  "bias_foreign_w1", "importance_live_in_demo_w1",
  "az_belief_econ_perf_cn_w1", "az_belief_trust_govt_w1",
  "az_belief_evalgovt_w1"
)

for (Y in median_below_vars_w1) {
  m_col <- paste0(Y, "_m")
  p_col <- paste0(Y, "_p")

  cutoff <- median(cy[[Y]], na.rm = TRUE)

  cy[[m_col]] <- cutoff

  cy[[p_col]] <- ifelse(
    is.na(cy[[Y]]),
    0L,
    as.integer(cy[[Y]] < cutoff)
  )
}

## already dummy, keep as it is
cy <- cy %>% rename(
  par_complain_school_w1 = participate_complain_school_w1,
  par_complain_school_w3 = participate_complain_school_w3
)

copy_as_is_w3 <- c(
  "vpn_purchase_wmt_record", "vpn_purchase_yes", "bias_dom_govt_policy_t1_w3",
  "protest_2011_tmrw_parade_w3", "participate_social_protest_w3", "participate_plan_vote_w3",
  "par_complain_school_w3", "plan_grad_foreignmaster_w3", "cp_t3_for_firm_w3", "cloc_for_w3"
)
for (Y in copy_as_is_w3) {
  cy[[paste0(Y, "_p")]] <- cy[[Y]]
}

copy_as_is_w1 <- c(
  "bias_dom_govt_policy_t1_w1", "protest_2011_tmrw_parade_w1", "participate_social_protest_w1",
  "participate_plan_vote_w1", "par_complain_school_w1", "plan_grad_foreignmaster_w1",
  "cp_t3_for_firm_w1", "cloc_for_w1"
)
for (Y in copy_as_is_w1) {
  cy[[paste0(Y, "_p")]] <- cy[[Y]]
}

## already dummy, flip sign
flip_sign_w3 <- c("bias_for_govt_policy_t1_w3", "stock_participation_w3")
for (Y in flip_sign_w3) {
  cy[[paste0(Y, "_p")]] <- 1 - cy[[Y]]
}

flip_sign_w1 <- c("bias_for_govt_policy_t1_w1", "stock_participation_w1")
for (Y in flip_sign_w1) {
  cy[[paste0(Y, "_p")]] <- 1 - cy[[Y]]
}


## final variable subset (Stata `keep`)
cy <- cy %>%
  select(
    treatment_main, active_user, treatment_control,
    treatment_vpnonly, treatment_nlonly, treatment_vpnnl,
    info_foreign_website_w3_p,        # A.1.2 Ranked high: foreign websites
    info_freq_website_for_w3_p,       # A.1.6 Freq. of visiting foreign websites for info.
    az_belief_media_value_w3_p,       # A.3 Valuation of access to foreign media outlets
    az_belief_media_trust_w3_p,       # A.4 Trust in non-domestic media outlets
    bias_domestic_w3_p,                # A.5.1 Degree of censorship on domestic news outlets
    bias_foreign_w3_p,                 # A.5.2 Degree of censorship on foreign news outlets
    az_belief_media_justif_w3_p,      # A.6 Censorship unjustified
    bias_dom_govt_policy_t1_w3_p,     # A.7.1 Domestic cens. driven by govt. policies
    bias_for_govt_policy_t1_w3_p,     # A.7.2 Foreign cens. driven by govt. policies
    info_foreign_website_w1_p,        # A.1.2 Ranked high: foreign websites
    info_freq_website_for_w1_p,       # A.1.6 Freq. of visiting foreign websites for info.
    az_belief_media_value_w1_p,       # A.3 Valuation of access to foreign media outlets
    az_belief_media_trust_w1_p,       # A.4 Trust in non-domestic media outlets
    bias_domestic_w1_p,                # A.5.1 Degree of censorship on domestic news outlets
    bias_foreign_w1_p,                 # A.5.2 Degree of censorship on foreign news outlets
    az_belief_media_justif_w1_p,      # A.6 Censorship unjustified
    bias_dom_govt_policy_t1_w1_p,     # A.7.1 Domestic cens. driven by govt. policies
    bias_for_govt_policy_t1_w1_p,     # A.7.2 Foreign cens. driven by govt. policies
    vpn_purchase_wmt_record_p,        # A.2.1 Purchase discounted tool we offered
    vpn_purchase_yes_p                # A.2.2 Purchase any tool
  )

## save ChenYang2019.rds in /output
saveRDS(
  cy %>% haven::zap_labels() %>% as_tibble(),
  file = file.path(output_dir, "ChenYang2019.rds")
)

# ==============================
# ANALYSES, not used for this paper but from the Chen and Yang, 2019, also translated into R
# since the local environment's cy database is now the filtered version, need to reload Chen and Yang data to run the following codes
# ==============================

# Figure 2 --------------------------------------------------------------------

# figure2_vars <- c("wtp_vpn")
#
# for (Y in figure2_vars) {
#
#   figure_dat <- cy %>%
#     filter(
#       panelmerged_wave2 == 1,
#       panelmerged_wave3 == 1
#     )
#
#
#   w1 <- paste0(Y, "_w1")
#   w2 <- paste0(Y, "_w2")
#   w3 <- paste0(Y, "_w3")
#
#   figure_dat <- figure_dat %>%
#     group_by(treatment_main) %>%
#     summarise(
#
#       !!w1 := mean(.data[[w1]], na.rm = TRUE),
#       !!w2 := mean(.data[[w2]], na.rm = TRUE),
#       !!w3 := mean(.data[[w3]], na.rm = TRUE),
#
#       sd_w1 = sd(.data[[w1]], na.rm = TRUE),
#       sd_w2 = sd(.data[[w2]], na.rm = TRUE),
#       sd_w3 = sd(.data[[w3]], na.rm = TRUE),
#
#       n_w1 = sum(!is.na(.data[[w1]])),
#       n_w2 = sum(!is.na(.data[[w2]])),
#       n_w3 = sum(!is.na(.data[[w3]])),
#
#       .groups = "drop"
#     )
#
#
# # qt(.95, df) gives the identical one-sided critical value of Stata's invttail(df,.05)
#
#   figure_dat <- figure_dat %>%
#     mutate(
#
#       h_w1 = .data[[w1]] + qt(.95, n_w1 - 1) * sd_w1 / sqrt(n_w1),
#       l_w1 = .data[[w1]] - qt(.95, n_w1 - 1) * sd_w1 / sqrt(n_w1),
#
#       h_w2 = .data[[w2]] + qt(.95, n_w2 - 1) * sd_w2 / sqrt(n_w2),
#       l_w2 = .data[[w2]] - qt(.95, n_w2 - 1) * sd_w2 / sqrt(n_w2),
#
#       h_w3 = .data[[w3]] + qt(.95, n_w3 - 1) * sd_w3 / sqrt(n_w3),
#       l_w3 = .data[[w3]] - qt(.95, n_w3 - 1) * sd_w3 / sqrt(n_w3)
#     ) %>%
#
#     filter(!is.na(treatment_main))
#
# # pivot_longer() produces the same long structure as Stata reshape long
#
#   figure_long <- bind_rows(
#
#     figure_dat %>%
#       transmute(
#         treatment_main,
#         wave = 1,
#         value = .data[[w1]],
#         lower = l_w1,
#         upper = h_w1
#       ),
#
#     figure_dat %>%
#       transmute(
#         treatment_main,
#         wave = 2,
#         value = .data[[w2]],
#         lower = l_w2,
#         upper = h_w2
#       ),
#
#     figure_dat %>%
#       transmute(
#         treatment_main,
#         wave = 4,
#         value = .data[[w3]],
#         lower = l_w3,
#         upper = h_w3
#       )
#
#   )
#
#   figure_long <- figure_long %>%
#     mutate(
#
#       treatment = factor(
#         treatment_main,
#         levels = c(1, 2, 3, 4),
#         labels = c(
#           "Control",
#           "Access",
#           "Access + Encour.",
#           "Existing users"
#         )
#       )
#
#     )
#
# # The colors below are identical to the Stata graph while following on ggplot2 conventions.
#
#   p <- ggplot(
#     figure_long,
#     aes(
#       x = wave,
#       y = value,
#       colour = treatment,
#       group = treatment
#     )
#   ) +
#
#     geom_line(linewidth = 1) +
#
#     geom_point(size = 2.5) +
#
#     geom_errorbar(
#       aes(
#         ymin = lower,
#         ymax = upper
#       ),
#       width = .08
#     ) +
#
#     geom_vline(
#       xintercept = 1.2,
#       linetype = "dashed"
#     ) +
#
#     scale_x_continuous(
#       breaks = c(1, 2, 4),
#       labels = c(
#         "Nov. 2015",
#         "Apr. 2016",
#         "May 2017"
#       ),
#       limits = c(.8, 4.2)
#     ) +
#
#     scale_colour_manual(
#       values = c(
#         "grey60",
#         "grey20",
#         "firebrick3",
#         "dodgerblue3"
#       )
#     ) +
#
#     labs(
#       title = variable_labels[[Y]],
#       x = NULL,
#       y = NULL,
#       colour = NULL
#     ) +
#
#     theme_bw() +
#
#     theme(
#       panel.grid = element_blank(),
#       legend.position = "bottom"
#     )
#
#   ggsave(
#     filename = paste0("figure_panel_", Y, "_s.pdf"),
#     plot = p,
#     width = 6,
#     height = 4
#   )
#
# }
#
# # Figure 3 --------------------------------------------------------------------
#
# figure3_dat <- cy %>%
#   group_by(treatment_main) %>%
#   summarise(
#
#     # Stata: collapse (mean)
#     m_vpn_purchase_yes = mean(vpn_purchase_yes, na.rm = TRUE),
#     m_vpn_purchase_wmt = mean(vpn_purchase_wmt_record, na.rm = TRUE),
#
#     # Stata: collapse (sd)
#     s_vpn_purchase_yes = sd(vpn_purchase_yes, na.rm = TRUE),
#     s_vpn_purchase_wmt = sd(vpn_purchase_wmt_record, na.rm = TRUE),
#
#     # Stata: collapse (count)
#     n_vpn_purchase_yes = sum(!is.na(vpn_purchase_yes)),
#     n_vpn_purchase_wmt = sum(!is.na(vpn_purchase_wmt_record)),
#
#     .groups = "drop"
#   ) %>%
#
#   mutate(
#
# # qt(.95, df) returns the same critical t value as Stata's invttail(df,.05)
#
#     h_vpn_purchase_yes =
#       m_vpn_purchase_yes +
#       qt(.95, n_vpn_purchase_yes - 1) *
#       s_vpn_purchase_yes /
#       sqrt(n_vpn_purchase_yes),
#
#     l_vpn_purchase_yes =
#       m_vpn_purchase_yes -
#       qt(.95, n_vpn_purchase_yes - 1) *
#       s_vpn_purchase_yes /
#       sqrt(n_vpn_purchase_yes),
#
#     h_vpn_purchase_wmt =
#       m_vpn_purchase_wmt +
#       qt(.95, n_vpn_purchase_wmt - 1) *
#       s_vpn_purchase_wmt /
#       sqrt(n_vpn_purchase_wmt),
#
#     l_vpn_purchase_wmt =
#       m_vpn_purchase_wmt -
#       qt(.95, n_vpn_purchase_wmt - 1) *
#       s_vpn_purchase_wmt /
#       sqrt(n_vpn_purchase_wmt),
#
#     # Stata benchmark variables
#     benchmark_1 = 0.65,
#     benchmark_2 = 0.47,
#     benchmark_x = treatment_main - 0.5,
#
#     treatment = factor(
#       treatment_main,
#       levels = c(1, 2, 3, 4),
#       labels = c(
#         "Control",
#         "Access",
#         "Access + Encour.",
#         "Existing users"
#       )
#     )
#   )
#
#
# # ggplot version of the Stata twoway graph.
#
# figure3_plot <-
#
#   ggplot(figure3_dat, aes(x = treatment_main)) +
#
#   ## Self-reported purchase (background bar)
#   geom_col(
#     aes(y = m_vpn_purchase_yes),
#     fill = "grey75",
#     alpha = 0.45,
#     width = 0.75
#   ) +
#
#   ## Walmart purchase (foreground bar)
#   geom_col(
#     aes(
#       y = m_vpn_purchase_wmt,
#       fill = treatment
#     ),
#     alpha = 0.75,
#     width = 0.75
#   ) +
#
#   ## Confidence intervals for self-reported purchase
#   geom_errorbar(
#     aes(
#       ymin = l_vpn_purchase_yes,
#       ymax = h_vpn_purchase_yes
#     ),
#     width = .12
#   ) +
#
#   ## Benchmark lines
#   geom_segment(
#     aes(
#       x = 2,
#       xend = 4,
#       y = benchmark_1,
#       yend = benchmark_1
#     ),
#     linetype = "dashed",
#     colour = "indianred3",
#     linewidth = 0.7
#   ) +
#
#   geom_segment(
#     aes(
#       x = 1,
#       xend = 3,
#       y = benchmark_2,
#       yend = benchmark_2
#     ),
#     linetype = "dashed",
#     colour = "indianred3",
#     linewidth = 0.7
#   ) +
#
#   scale_fill_manual(
#     values = c(
#       "Control" = "grey70",
#       "Access" = "grey70",
#       "Access + Encour." = "firebrick3",
#       "Existing users" = "dodgerblue3"
#     )
#   ) +
#
#   scale_x_continuous(
#     breaks = 1:4,
#     labels = levels(figure3_dat$treatment)
#   ) +
#
#   scale_y_continuous(
#     limits = c(0, 1),
#     breaks = seq(0, 1, 0.2),
#     labels = scales::percent_format(accuracy = 1)
#   ) +
#
#   labs(
#     x = NULL,
#     y = "% purchased any circumvention tool after April 2017"
#   ) +
#
#   theme_bw() +
#
#   theme(
#     legend.position = "none",
#     panel.grid = element_blank()
#   )
#
# ggsave(
#   filename = file.path(outreg_dir, "figure_vpn_purchase_yes_2.pdf"),
#   plot = figure3_plot,
#   width = 5,
#   height = 3
# )
#
# # Table 1 ---------------------------------------------------------------------
#
# # Rather than using Stata's preserve/append, the code computes one row per variable and `bind_rows()` them
#
# panelmerged_wave3_temp <-
#   cy %>%
#   filter(panelmerged_wave3 == 1) %>%
#   mutate(wave3 = 1) %>%
#   select(-treatment_master)
#
# table1_dat <-
#   bind_rows(
#     cy,
#     panelmerged_wave3_temp
#   ) %>%
#   mutate(
#     wave3 = replace_na(wave3, 0)
#   )
#
# table1_vars <-
#   c(
#     var_demog_reg_personal,
#     "az_demographics_education",
#     "az_demographics_english",
#     "az_demographics_travel",
#     "az_demographics_household",
#     "az_preference_risk",
#     "az_preference_time",
#     "az_preference_altruism",
#     "az_preference_reciprocity",
#     "panelmerged_wave3",
#     "treatment_control",
#     "treatment_vpnonly",
#     "treatment_nlonly",
#     "treatment_vpnnl",
#     "treatment_user"
#   )
#
#
# # Calculate one summary row for each variable.
# # This replaces Stata's postfile/post workflow.
#
# summstats_main <-
#
#   purrr::map_dfr(table1_vars, function(X) {
#
#     x <- table1_dat[[X]]
#
#     ## ----------------------------
#     ## Overall Wave 1
#     ## ----------------------------
#
#     w1 <- table1_dat %>%
#       filter(wave3 == 0)
#
#     ## ----------------------------
#     ## Wave 1 respondents who later
#     ## completed Wave 3
#     ## ----------------------------
#
#     w3 <- table1_dat %>%
#       filter(
#         panelmerged_wave3 == 1,
#         wave3 == 0
#       )
#
#     ## Treatment-specific summaries
#
#     treat_summary <-
#
#       map_dfr(
#         1:5,
#         function(tt) {
#
#           tmp <-
#             table1_dat %>%
#             filter(
#               panelmerged_wave3 == 1,
#               wave3 == 0,
#               treatment_master == tt
#             )
#
#           tibble(
#
#             treatment = tt,
#
#             mean = mean(tmp[[X]], na.rm = TRUE),
#             sd = sd(tmp[[X]], na.rm = TRUE),
#             n = sum(!is.na(tmp[[X]]))
#
#           )
#
#         }
#
#       )
#
#     ## ----------------------------
#     ## ANOVA across experimental arms
#     ## (exclude existing users,
#     ## exactly as in Stata)
#     ## ----------------------------
#
#     anova_dat <-
#
#       table1_dat %>%
#       filter(
#         panelmerged_wave3 == 1,
#         wave3 == 0,
#         treatment_master != 5
#       )
#
#     fit <-
#
#       aov(
#         reformulate(
#           "treatment_master",
#           response = X
#         ),
#         data = anova_dat
#       )
#
#     anova_tbl <- summary(fit)[[1]]
#
#     ## ----------------------------
#     ## Attrition test
#     ## Equivalent to:
#     ##
#     ## ttest X, by(wave3)
#     ## ----------------------------
#
#     attrition_test <-
#
#       t.test(
#         reformulate("wave3", response = X),
#         data = table1_dat
#       )
#
#     tibble(
#
#       variable = X,
#
#       mu_w1 = mean(w1[[X]], na.rm = TRUE),
#       sd_w1 = sd(w1[[X]], na.rm = TRUE),
#       N_w1 = sum(!is.na(w1[[X]])),
#
#       mu_w3 = mean(w3[[X]], na.rm = TRUE),
#       sd_w3 = sd(w3[[X]], na.rm = TRUE),
#       N_w3 = sum(!is.na(w3[[X]])),
#
#       attrition_pvalue = attrition_test$p.value,
#
#       mu_c = treat_summary$mean[1],
#       sd_c = treat_summary$sd[1],
#       N_c = treat_summary$n[1],
#
#       mu_a = treat_summary$mean[2],
#       sd_a = treat_summary$sd[2],
#       N_a = treat_summary$n[2],
#
#       mu_ce = treat_summary$mean[3],
#       sd_ce = treat_summary$sd[3],
#       N_ce = treat_summary$n[3],
#
#       mu_ae = treat_summary$mean[4],
#       sd_ae = treat_summary$sd[4],
#       N_ae = treat_summary$n[4],
#
#       mu_ex = treat_summary$mean[5],
#       sd_ex = treat_summary$sd[5],
#       N_ex = treat_summary$n[5],
#
#       anova_fstat = anova_tbl$`F value`[1],
#       anova_pvalue = anova_tbl$`Pr(>F)`[1]
#
#     )
#
#   })
#
# write.csv(
#   summstats_main,
#   file.path(outreg_dir, "summstats_main.csv"),
#   row.names = FALSE
# )
#
# # Additional packages for table 3 and table 4
# library(purrr)
# library(broom)
# library(sandwich)
# library(lmtest)
# library(AER)
# library(numDeriv)
# library(readr)
#
# # Table 3 ---------------------------------------------------------------------
#
# table3_outcomes <- c("a", "b", "c", "d", "e")
#
# robust_lm_results <- function(model, keep_terms) {
#
#   vc <- sandwich::vcovHC(model, type = "HC1")
#
#   lmtest::coeftest(model, vcov. = vc) %>%
#     broom::tidy() %>%
#     filter(term %in% keep_terms) %>%
#     transmute(
#       term,
#       estimate,
#       std_error = std.error,
#       statistic,
#       p_value = p.value
#     )
# }
#
# table3_panel_a <- purrr::map_dfr(table3_outcomes, function(Y) {
#
#   dv <- paste0("az_overall_", Y, "_w3")
#
#   table3_dat <- cy %>%
#     filter(panelmerged_wave3 == 1)
#
#   # Stata: reg ..., r
#   # R: lm() plus HC1 robust standard errors.
#   model <- lm(
#     reformulate(
#       c("treatment_vpnonly", "treatment_nlonly", "treatment_vpnnl"),
#       response = dv
#     ),
#     data = table3_dat %>% filter(treatment_user != 1)
#   )
#
#   reg_out <- robust_lm_results(
#     model,
#     keep_terms = c(
#       "treatment_vpnonly",
#       "treatment_nlonly",
#       "treatment_vpnnl"
#     )
#   )
#
#   # Stata addstat() summaries.
#   stats <- tibble(
#     outcome = Y,
#     dv = dv,
#     mean_dv_all = mean(table3_dat[[dv]][table3_dat$treatment_user != 1], na.rm = TRUE),
#     sd_dv_all = sd(table3_dat[[dv]][table3_dat$treatment_user != 1], na.rm = TRUE),
#     mean_dv_control = mean(table3_dat[[dv]][table3_dat$treatment_control == 1], na.rm = TRUE),
#     sd_dv_control = sd(table3_dat[[dv]][table3_dat$treatment_control == 1], na.rm = TRUE),
#     mean_dv_user = mean(table3_dat[[dv]][table3_dat$treatment_user == 1], na.rm = TRUE),
#     sd_dv_user = sd(table3_dat[[dv]][table3_dat$treatment_user == 1], na.rm = TRUE)
#   )
#
#   reg_out %>%
#     mutate(
#       outcome = Y,
#       dv = dv,
#       specification = "Raw diff"
#     ) %>%
#     left_join(stats, by = c("outcome", "dv"))
# })
#
# write_csv(
#   table3_panel_a,
#   file.path(outreg_dir, "outregs_maintreatmenteffects_panel_a.csv")
# )
#
#
# # Panel A addition: Lee bounds --------------------------------------------
#
# table3_lee_bounds <- purrr::map_dfr(table3_outcomes, function(Y) {
#
#   dv <- paste0("az_overall_", Y, "_w3")
#
#   cy %>%
#     mutate(weights = 1) %>%
#     filter(treatment_control == 1 | treatment_vpnnl == 1) %>%
#     trimbound(
#       y = dv,
#       d = "treatment_vpnnl",
#       wgt = "weights"
#     )
# })
#
# write_csv(
#   table3_lee_bounds,
#   file.path(outreg_dir, "outregs_maintreatmenteffects_lee_bounds.csv")
# )
#
#
# # Panel B: two-stage estimates --------------------------------------------
#
# table3_panel_b <- purrr::map_dfr(table3_outcomes, function(Y) {
#
#   dv <- paste0("az_overall_", Y, "_w3")
#
#   table3_dat <- cy %>%
#     filter(
#       panelmerged_wave3 == 1,
#       treatment_user != 1
#     )
#
# # R's AER::ivreg() implements the same 2SLS formula as Stata's
# # ivregress 2sls y (active_user = treatment_vpnonly treatment_nlonly treatment_vpnnl), first
#   model <- AER::ivreg(
#     as.formula(
#       paste0(
#         dv,
#         " ~ active_user | treatment_vpnonly + treatment_nlonly + treatment_vpnnl"
#       )
#     ),
#     data = table3_dat
#   )
#
#   vc <- sandwich::vcovHC(model, type = "HC1")
#
#   broom::tidy(lmtest::coeftest(model, vcov. = vc)) %>%
#     mutate(
#       outcome = Y,
#       dv = dv,
#       specification = "2SLS"
#     )
# })
#
# write_csv(
#   table3_panel_b,
#   file.path(outreg_dir, "outregs_maintreatmenteffects_panel_b.csv")
# )
#
# # Table 4 ---------------------------------------------------------------------
#
# # Helper: social learning estimates -----------------------------------------
#
# social_learning_stats <- function(data, Y, rm_new_var, ownXnew_var) {
#
#   reg_dat <- data %>%
#     filter(
#       vpn_roommate_existing == 0,
#       .data[[rm_new_var]] < 2
#     )
#
#   model <- lm(
#     reformulate(
#       c("soclearning_ownaccess", rm_new_var, ownXnew_var),
#       response = Y
#     ),
#     data = reg_dat
#   )
#
#   reg_out <- robust_lm_results(
#     model,
#     keep_terms = c(
#       "soclearning_ownaccess",
#       rm_new_var,
#       ownXnew_var
#     )
#   )
#
#   # Stata: local mu_100, mu_000, mu_101, mu_001, etc.
#   mu_100 <- data %>%
#     filter(treatment_master >= 4, vpn_roommate_existing == 0, .data[[rm_new_var]] == 0) %>%
#     summarise(x = mean(.data[[Y]], na.rm = TRUE)) %>%
#     pull(x)
#
#   mu_000 <- data %>%
#     filter(treatment_master < 4, vpn_roommate_existing == 0, .data[[rm_new_var]] == 0) %>%
#     summarise(x = mean(.data[[Y]], na.rm = TRUE)) %>%
#     pull(x)
#
#   mu_101 <- data %>%
#     filter(treatment_master >= 4, vpn_roommate_existing == 0, .data[[rm_new_var]] == 1) %>%
#     summarise(x = mean(.data[[Y]], na.rm = TRUE)) %>%
#     pull(x)
#
#   mu_001 <- data %>%
#     filter(treatment_master < 4, vpn_roommate_existing == 0, .data[[rm_new_var]] == 1) %>%
#     summarise(x = mean(.data[[Y]], na.rm = TRUE)) %>%
#     pull(x)
#
#   mu_002 <- data %>%
#     filter(treatment_master < 4, vpn_roommate_existing == 0, .data[[rm_new_var]] == 2) %>%
#     summarise(x = mean(.data[[Y]], na.rm = TRUE)) %>%
#     pull(x)
#
#   mu_102 <- data %>%
#     filter(treatment_master >= 4, vpn_roommate_existing == 0, .data[[rm_new_var]] == 2) %>%
#     summarise(x = mean(.data[[Y]], na.rm = TRUE)) %>%
#     pull(x)
#
#   p_all <- mu_100 - mu_000
#   q_new_1_0 <- (mu_001 - mu_000) / mu_100
#   q_new_1_1 <- (mu_101 - mu_100) / mu_100
#
#   pred_mu_102 <- mu_100 + (1 - (1 - mu_100 * q_new_1_1)^2)
#   pred_mu_002 <- mu_000 + (1 - (1 - mu_100 * q_new_1_0)^2)
#
#   reg_out %>%
#     mutate(
#       outcome = Y,
#       p_all = p_all,
#       q_new_1_0 = q_new_1_0,
#       q_new_1_1 = q_new_1_1,
#       pred_mu_002 = pred_mu_002,
#       mu_002 = mu_002,
#       pred_mu_102 = pred_mu_102,
#       mu_102 = mu_102,
#       specification = "Raw Diff"
#     )
# }
#
#
# # Main Table 4 regressions ---------------------------------------------------
#
# table4_w2 <- purrr::map_dfr(
#   c("news_c_panamapapers"),
#   ~ social_learning_stats(
#     data = cy,
#     Y = .x,
#     rm_new_var = "soclearning_rm_new_w2",
#     ownXnew_var = "soclearning_ownXnew_w2"
#   )
# )
#
# table4_w3 <- purrr::map_dfr(
#   c("news_c_coalprod", "news_c_hkceelection", "news_perccor_cen_all"),
#   ~ social_learning_stats(
#     data = cy,
#     Y = .x,
#     rm_new_var = "soclearning_rm_new_w3",
#     ownXnew_var = "soclearning_ownXnew_w3"
#   )
# )
#
# table4_main <- bind_rows(table4_w2, table4_w3)
#
# write_csv(
#   table4_main,
#   file.path(outreg_dir, "outregs_wave3sociallearning.csv")
# )
#
#
# # Standard errors of nonlinear estimates -------------------------------------
#
# nlcom_social_learning <- function(data, Y, rm_new_var) {
#
#   tmp <- data %>%
#     mutate(
#       sc_reg_100 = as.integer(treatment_master >= 4 & vpn_roommate_existing == 0 & .data[[rm_new_var]] == 0),
#       sc_reg_000 = as.integer(treatment_master <  4 & vpn_roommate_existing == 0 & .data[[rm_new_var]] == 0),
#       sc_reg_101 = as.integer(treatment_master >= 4 & vpn_roommate_existing == 0 & .data[[rm_new_var]] == 1),
#       sc_reg_001 = as.integer(treatment_master <  4 & vpn_roommate_existing == 0 & .data[[rm_new_var]] == 1)
#     ) %>%
#     filter(
#       vpn_roommate_existing == 0,
#       .data[[rm_new_var]] < 2
#     )
#
#   # Stata: reg Y sc_reg_100 sc_reg_000 sc_reg_101 sc_reg_001, nocons
#   model <- lm(
#     reformulate(
#       c("sc_reg_100", "sc_reg_000", "sc_reg_101", "sc_reg_001"),
#       response = Y,
#       intercept = FALSE
#     ),
#     data = tmp
#   )
#
#   b <- coef(model)
#   vc <- sandwich::vcovHC(model, type = "HC1")
#
#   # Stata nlcom:
#   # q_new_1_0, q_new_1_1, pred_mu_002, pred_mu_102.
#   f <- function(beta) {
#
#     names(beta) <- names(b)
#
#     q_new_1_0 <- (beta["sc_reg_001"] - beta["sc_reg_000"]) / beta["sc_reg_100"]
#     q_new_1_1 <- (beta["sc_reg_101"] - beta["sc_reg_100"]) / beta["sc_reg_100"]
#
#     pred_mu_002 <- beta["sc_reg_000"] +
#       (1 - (1 - beta["sc_reg_100"] * q_new_1_0)^2)
#
#     pred_mu_102 <- beta["sc_reg_100"] +
#       (1 - (1 - beta["sc_reg_100"] * q_new_1_1)^2)
#
#     c(
#       q_new_1_0 = q_new_1_0,
#       q_new_1_1 = q_new_1_1,
#       pred_mu_002 = pred_mu_002,
#       pred_mu_102 = pred_mu_102
#     )
#   }
#
#   est <- f(b)
#
#   # Delta-method SEs.
#   grad <- numDeriv::jacobian(f, b)
#   se <- sqrt(diag(grad %*% vc %*% t(grad)))
#
#   tibble(
#     outcome = Y,
#     term = names(est),
#     estimate = as.numeric(est),
#     std_error = as.numeric(se)
#   )
# }
#
# table4_nlcom_w2 <- purrr::map_dfr(
#   c("news_c_panamapapers"),
#   ~ nlcom_social_learning(
#     data = cy,
#     Y = .x,
#     rm_new_var = "vpn_roommate_new_w2"
#   )
# )
#
# table4_nlcom_w3 <- purrr::map_dfr(
#   c("news_c_coalprod", "news_c_hkceelection", "news_perccor_cen_all"),
#   ~ nlcom_social_learning(
#     data = cy,
#     Y = .x,
#     rm_new_var = "vpn_roommate_new_w3"
#   )
# )
#
# table4_nlcom <- bind_rows(table4_nlcom_w2, table4_nlcom_w3)
#
# write_csv(
#   table4_nlcom,
#   file.path(outreg_dir, "outregs_wave3sociallearning_nlcom.csv")
# )
#
# # Figure A.8, A.11, A.12, A.13 ---------------------------------------------
#
# library(ggplot2)
# library(purrr)
# library(broom)
# library(sandwich)
# library(lmtest)
# library(patchwork)
#
# fig_dotplot_vars <- c(
#   var_media_valuation_reg_w3, "az_belief_media_value_w3",
#   var_media_trust_reg_w3, "az_belief_media_trust_w3",
#   var_censor_level_reg_w3,
#   var_censor_justif_reg_w3, "az_belief_media_justif_w3",
#   var_censor_driver_reg_dom_w3, var_censor_driver_reg_for_w3,
#   var_knowledge_news_reg_cen_w3, "news_perccor_cen_w3",
#   var_knowledge_news_reg_unc_w3, "news_perccor_unc_w3",
#   var_knowledge_prot_reg_chi_w3, "protest_pcheard_china_w3",
#   var_knowledge_prot_reg_for_w3, "protest_pcheard_foreign_w3",
#   var_knowledge_prot_reg_fak_w3,
#   var_knowledge_meta_reg_w3, "az_knowledge_meta_w3",
#   var_econ_guess_reg_cn_perf_w3, "az_belief_econ_perf_cn_w3",
#   var_econ_guess_reg_cn_conf_w3, "az_belief_econ_conf_cn_w3",
#   var_econ_guess_reg_us_perf_w3, "az_belief_econ_perf_us_w3",
#   var_econ_guess_reg_us_conf_w3, "az_belief_econ_conf_us_w3",
#   var_demand_change_reg_w3, "az_belief_instchange_w3",
#   var_trust_inst_reg_govt_w3, "az_belief_trust_govt_w3",
#   var_trust_inst_reg_foreign_w3, "az_belief_trust_foreign_w3",
#   var_eval_govt_reg_w3, "az_belief_evalgovt_w3",
#   var_democracy_reg_w3,
#   var_willing_fight_reg_w3, "az_belief_willing_w3",
#   var_vpn_purchase_reg_w3,
#   var_info_ranking_reg_w3,
#   var_info_freq_reg_w3,
#   var_socialinteract_reg_w3,
#   var_polparticipation_reg_w3,
#   var_stock_invest_reg_w3,
#   var_planaftergrad_reg_w3,
#   var_career_sector_reg_w3,
#   var_career_loc_reg_w3
# )
#
# robust_coef <- function(model, term) {
#   vc <- sandwich::vcovHC(model, type = "HC1")
#   lmtest::coeftest(model, vcov. = vc) |>
#     broom::tidy() |>
#     filter(term == !!term)
# }
#
# fig_dotplot_dat <- map2_dfr(fig_dotplot_vars, seq_along(fig_dotplot_vars), function(Y, v) {
#
#   # Stata first standardizes Y using the regression sample from:
#   # reg Y treatment_vpnonly treatment_vpnnl, r
#   sample_fit <- lm(
#     reformulate(c("treatment_vpnonly", "treatment_vpnnl"), response = Y),
#     data = cy
#   )
#
#   sample_rows <- as.integer(names(residuals(sample_fit)))
#   mu <- mean(cy[[Y]][sample_rows], na.rm = TRUE)
#   sig <- sd(cy[[Y]][sample_rows], na.rm = TRUE)
#
#   tmp <- cy %>%
#     mutate(eb = (.data[[Y]] - mu) / sig)
#
#   # Access and Access + Encouragement, excluding current paid users.
#   fit_main <- lm(
#     eb ~ treatment_vpnonly + treatment_vpnnl,
#     data = tmp %>% filter(vpn_current_paid_user != 1)
#   )
#
#   vc_main <- sandwich::vcovHC(fit_main, type = "HC1")
#   ct_main <- lmtest::coeftest(fit_main, vcov. = vc_main) |>
#     broom::tidy()
#
#   rcoef <- ct_main %>% filter(term == "treatment_vpnnl") %>% pull(estimate)
#   rvpn  <- ct_main %>% filter(term == "treatment_vpnonly") %>% pull(estimate)
#   rse   <- ct_main %>% filter(term == "treatment_vpnnl") %>% pull(std.error)
#
#   # Existing users relative to control.
#   fit_user <- lm(
#     eb ~ treatment_user,
#     data = tmp %>% filter(treatment_control == 1 | treatment_user == 1)
#   )
#
#   vc_user <- sandwich::vcovHC(fit_user, type = "HC1")
#   ruser <- lmtest::coeftest(fit_user, vcov. = vc_user) |>
#     broom::tidy() |>
#     filter(term == "treatment_user") |>
#     pull(estimate)
#
#   tibble(
#     variable = Y,
#     vnum_raw = v,
#     rcoef = rcoef,
#     rvpn = rvpn,
#     rse = rse,
#     rcil = rcoef - 1.9 * rse,
#     rcih = rcoef + 1.9 * rse,
#     ruser = ruser
#   )
# })
#
# # Stata manually inserts blank rows by shifting vnum.
# fig_dotplot_dat <- fig_dotplot_dat %>%
#   mutate(
#     vnum = vnum_raw,
#     vnum = ifelse(between(vnum, 110, 117), vnum + 35, vnum),
#     vnum = ifelse(between(vnum, 102, 109), vnum + 34, vnum),
#     vnum = ifelse(between(vnum, 97, 101), vnum + 33, vnum),
#     vnum = ifelse(between(vnum, 96, 96), vnum + 32, vnum),
#     vnum = ifelse(between(vnum, 93, 95), vnum + 31, vnum),
#     vnum = ifelse(between(vnum, 91, 92), vnum + 30, vnum),
#     vnum = ifelse(between(vnum, 90, 90), vnum + 29, vnum),
#     vnum = ifelse(between(vnum, 85, 89), vnum + 28, vnum),
#     vnum = ifelse(between(vnum, 83, 84), vnum + 27, vnum),
#     vnum = ifelse(between(vnum, 79, 82), vnum + 25, vnum),
#     vnum = ifelse(between(vnum, 78, 78), vnum + 24, vnum),
#     vnum = ifelse(between(vnum, 74, 77), vnum + 23, vnum),
#     vnum = ifelse(between(vnum, 71, 73), vnum + 22, vnum),
#     vnum = ifelse(between(vnum, 67, 70), vnum + 21, vnum),
#     vnum = ifelse(between(vnum, 64, 66), vnum + 20, vnum),
#     vnum = ifelse(between(vnum, 61, 63), vnum + 18, vnum),
#     vnum = ifelse(between(vnum, 58, 60), vnum + 17, vnum),
#     vnum = ifelse(between(vnum, 55, 57), vnum + 16, vnum),
#     vnum = ifelse(between(vnum, 52, 54), vnum + 15, vnum),
#     vnum = ifelse(between(vnum, 49, 51), vnum + 13, vnum),
#     vnum = ifelse(between(vnum, 48, 48), vnum + 12, vnum),
#     vnum = ifelse(between(vnum, 42, 47), vnum + 11, vnum),
#     vnum = ifelse(between(vnum, 37, 41), vnum + 10, vnum),
#     vnum = ifelse(between(vnum, 32, 36), vnum + 9, vnum),
#     vnum = ifelse(between(vnum, 24, 31), vnum + 8, vnum),
#     vnum = ifelse(between(vnum, 20, 23), vnum + 6, vnum),
#     vnum = ifelse(between(vnum, 16, 19), vnum + 5, vnum),
#     vnum = ifelse(between(vnum, 10, 15), vnum + 4, vnum),
#     vnum = ifelse(between(vnum, 9, 9), vnum + 3, vnum),
#     vnum = ifelse(between(vnum, 8, 8), vnum + 2, vnum),
#     vnum = ifelse(between(vnum, 4, 7), vnum + 1, vnum)
#   )
#
# fig_dotplot_labels <- c(
#   `1` = "Willingness to pay for circumvention tool",
#   `2` = "Value added of foreign media access",
#   `3` = "z-score: valuation of access to foreign media outlets",
#   `5` = "Distrust in domestic state-owned media",
#   `6` = "Distrust in domestic privately-owned media",
#   `7` = "Trust in foreign media",
#   `8` = "z-score: trust in non-domestic media outlets",
#   `10` = "Degree of censorship on domestic news outlets",
#   `12` = "Degree of censorship on foreign news outlets",
#   `14` = "Unjustified: censoring economic news",
#   `15` = "Unjustified: censoring political news",
#   `16` = "Unjustified: censoring social news",
#   `17` = "Unjustified: censoring foreign news",
#   `18` = "Unjustified: censoring pornography",
#   `19` = "z-score: censorship unjustified",
#   `21` = "Domestic cens. driven by govt. policies",
#   `22` = "Domestic cens. driven by corp. interest",
#   `23` = "Domestic cens. driven by media’s ideology",
#   `24` = "Domestic cens. driven by readers’ demand",
#   `26` = "Foreign cens. driven by govt. policies",
#   `27` = "Foreign cens. driven by corp. interest",
#   `28` = "Foreign cens. driven by media’s ideology",
#   `29` = "Foreign cens. driven by readers’ demand",
#   `32` = "Steel production reduction reaches target",
#   `33` = "Trump registered trademarks in China",
#   `34` = "Jianhua Xiao kidnapped in Hong Kong",
#   `35` = "Xinjiang installed GPS on all automobiles",
#   `36` = "China and Norway re-normalize ties",
#   `37` = "Feminist groups fight women's rights",
#   `38` = "Carrie Lam becomes HK Chief Executive",
#   `39` = "% quizzes answered correctly: poli. sensitive news",
#   `41` = "China stops importing coal from North Korea",
#   `42` = "H7N9 influenza epidemic",
#   `43` = "Transnational railway in Ethiopia",
#   `44` = "Foreign reserves fall below threshold",
#   `45` = "% quizzes answered correctly: nonsensitive news",
#   `47` = "2012 HK Anti-National Curr. Movement",
#   `48` = "2014 HK Umbrella Revolution",
#   `49` = "2016 HK Mong Kok Revolution",
#   `50` = "2014 Taiwan Sunflower Stud. Movement",
#   `51` = "% protests in Greater China heard of",
#   `53` = "2014 Ukrainian Euromaidan Revolution",
#   `54` = "2010 Arab Spring",
#   `55` = "2014 Crimean Status Referendum",
#   `56` = "2010 Catalonian Indep. Movement",
#   `57` = "2017 Women's March",
#   `58` = "% foreign protests heard of",
#   `60` = "2011 Tomorrow Revolution [fake]",
#   `62` = "Informedness of issues in China",
#   `63` = "Greater informedness than peers",
#   `64` = "z-score: self-assessment of knowledge level",
#   `67` = "Guess on GDP growth rate in 2016 China",
#   `68` = "Guess on SSCI by end of 2016",
#   `69` = "z-score: optimistic belief of Chinese economy",
#   `71` = "Confidence of China GDP guess",
#   `72` = "Confidence of SSCI guess",
#   `73` = "z-score: confidence of guesses on Chinese economy",
#   `75` = "Guess on GDP growth rate in 2016 US",
#   `76` = "Guess on DJI by end of 2016",
#   `77` = "z-score: optimistic belief of US economy",
#   `79` = "Confidence of US GDP guess",
#   `80` = "Confidence of DJI guess",
#   `81` = "z-score: confidence of guesses on US economy",
#   `84` = "Economic system needs changes",
#   `85` = "Political system needs changes",
#   `86` = "z-score: demand for institutional change",
#   `88` = "Trust in central govt. of China",
#   `89` = "Trust in provincial govt. of China",
#   `90` = "Trust in local govt. of China",
#   `91` = "z-score: trust in Chinese govt.",
#   `93` = "Trust in central govt. of Japan",
#   `94` = "Trust in federal govt. of US",
#   `95` = "z-score: trust in foreign govt.",
#   `97` = "Satisfaction of economic dev.",
#   `98` = "Satisfaction of domestic politics",
#   `99` = "Satisfaction of diplomatic affairs",
#   `100` = "z-score: satisfaction of govt’s performance",
#   `102` = "Living in democracy is not important",
#   `104` = "Willing to battle illegal govt. acts",
#   `105` = "Willing to report govt. misconduct",
#   `106` = "Willing to stand up for the weak",
#   `107` = "z-score: willingness to act",
#   `110` = "Purchase discounted tool we offered",
#   `111` = "Purchase any tool",
#   `113` = "Ranked high: domestic websites",
#   `114` = "Ranked high: foreign websites",
#   `115` = "Ranked high: domestic social media",
#   `116` = "Ranked high: foreign social media",
#   `117` = "Ranked high: word of mouth",
#   `119` = "Frequency of visiting foreign websites for info.",
#   `121` = "Frequency of discussing poli. with friends",
#   `122` = "Frequency of persuading others",
#   `124` = "Protests concerning social issues",
#   `125` = "Plan to vote for local PCR",
#   `126` = "Complain to school authorities",
#   `128` = "Currently invested in Chinese stock mkt.",
#   `130` = "Plan: grad. school in China",
#   `131` = "Plan: master degree abroad",
#   `132` = "Plan: PhD degree abroad",
#   `133` = "Plan: military in China",
#   `134` = "Plan: work right away",
#   `136` = "Sector pref.: national civil service",
#   `137` = "Sector pref.: local civil service",
#   `138` = "Sector pref.: military",
#   `139` = "Sector pref.: private firm in China",
#   `140` = "Sector pref.: foreign firm in China",
#   `141` = "Sector pref.: SOEs",
#   `142` = "Sector pref.: inst. organizations",
#   `143` = "Sector pref.: entrepreneurship",
#   `145` = "Location pref.: Beijing",
#   `146` = "Location pref.: Shanghai",
#   `147` = "Location pref.: Guangzhou and Shenzhen",
#   `148` = "Location pref.: tier 2 cities in central",
#   `149` = "Location pref.: other cities in China",
#   `150` = "Location pref.: HK and Macau",
#   `151` = "Location pref.: Taiwan",
#   `152` = "Location pref.: foreign cities"
# )
#
# fig_dotplot_long <- fig_dotplot_dat %>%
#   transmute(
#     vnum,
#     label = fig_dotplot_labels[as.character(vnum)],
#     `Access` = rvpn,
#     `Access + Encour.` = rcoef,
#     `Existing users` = ruser,
#     lower = rcil,
#     upper = rcih
#   ) %>%
#   pivot_longer(
#     cols = c("Access", "Access + Encour.", "Existing users"),
#     names_to = "group",
#     values_to = "estimate"
#   )
#
# figure_dotplot_master_w3 <- ggplot(fig_dotplot_long, aes(x = estimate, y = vnum)) +
#   geom_vline(xintercept = 0) +
#   geom_segment(
#     data = fig_dotplot_dat,
#     aes(x = rcil, xend = rcih, y = vnum, yend = vnum),
#     inherit.aes = FALSE
#   ) +
#   geom_point(aes(shape = group, colour = group), size = 1.6) +
#   geom_hline(yintercept = c(4, 9, 13, 20, 46, 61, 70, 74, 78, 87, 96, 101, 103, 112, 120, 123, 127, 129, 135),
#              linetype = "longdash", linewidth = .2) +
#   geom_hline(yintercept = c(11, 25, 40, 52, 59, 92, 118, 144),
#              linetype = "dashed", linewidth = .2) +
#   geom_hline(yintercept = c(31, 66, 83, 109),
#              linewidth = .35) +
#   scale_y_reverse(
#     breaks = as.integer(names(fig_dotplot_labels)),
#     labels = unname(fig_dotplot_labels)
#   ) +
#   labs(
#     x = "Standardized means (Control = 0)",
#     y = NULL,
#     colour = NULL,
#     shape = NULL
#   ) +
#   theme_bw() +
#   theme(
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.position = "bottom",
#     axis.text.y = element_text(size = 5)
#   )
#
# ggsave(
#   file.path(outreg_dir, "figure_dotplot_master_w3.pdf"),
#   figure_dotplot_master_w3,
#   width = 4,
#   height = 20
# )
#
# # Figure A.10, A.14, A.15, A.16 -------------------------------------------
#
# fig_dotplot_cce_vars <- c(
#   var_media_valuation_reg_w3, "az_belief_media_value_w3",
#   var_media_trust_reg_w3, "az_belief_media_trust_w3",
#   var_censor_level_reg_w3,
#   var_censor_justif_reg_w3, "az_belief_media_justif_w3",
#   var_censor_driver_reg_dom_w3, var_censor_driver_reg_for_w3,
#   var_knowledge_news_reg_cen_w3, "news_perccor_cen_w3",
#   var_knowledge_news_reg_unc_w3, "news_perccor_unc_w3",
#   var_knowledge_prot_reg_chi_w3, "protest_pcheard_china_w3",
#   var_knowledge_prot_reg_for_w3, "protest_pcheard_foreign_w3",
#   var_knowledge_prot_reg_fak_w3,
#   var_knowledge_meta_reg_w3, "az_knowledge_meta_w3",
#   var_econ_guess_reg_cn_perf_w3, "az_belief_econ_perf_cn_w3",
#   var_econ_guess_reg_cn_conf_w3, "az_belief_econ_conf_cn_w3",
#   var_econ_guess_reg_us_perf_w3, "az_belief_econ_perf_us_w3",
#   var_econ_guess_reg_us_conf_w3, "az_belief_econ_conf_us_w3",
#   var_demand_change_reg_w3, "az_belief_instchange_w3",
#   var_trust_inst_reg_govt_w3, "az_belief_trust_govt_w3",
#   var_trust_inst_reg_foreign_w3, "az_belief_trust_foreign_w3",
#   var_eval_govt_reg_w3, "az_belief_evalgovt_w3",
#   var_democracy_reg_w3,
#   var_willing_fight_reg_w3, "az_belief_willing_w3",
#   var_vpn_purchase_reg_w3,
#   var_info_ranking_reg_w3,
#   var_info_freq_reg_w3,
#   var_socialinteract_reg_w3,
#   var_polparticipation_reg_w3,
#   var_stock_invest_reg_w3,
#   var_planaftergrad_reg_w3,
#   var_career_sector_reg_w3,
#   var_career_loc_reg_w3
# )
#
# fig_dotplot_cce_dat <- purrr::map2_dfr(
#   fig_dotplot_cce_vars,
#   seq_along(fig_dotplot_cce_vars),
#   function(Y, v) {
#
#     # Stata:
#     # qui xi: reg `Y' treatment_vpnonly treatment_vpnnl, r
#     # su `Y' if e(sample)
#     #
#     # R:
#     # Use the complete regression sample from the first regression to
#     # standardize each outcome.
#     sample_fit <- lm(
#       reformulate(
#         c("treatment_vpnonly", "treatment_vpnnl"),
#         response = Y
#       ),
#       data = cy
#     )
#
#     sample_rows <- as.integer(names(residuals(sample_fit)))
#
#     eb_mean <- mean(cy[[Y]][sample_rows], na.rm = TRUE)
#     eb_sd   <- sd(cy[[Y]][sample_rows], na.rm = TRUE)
#
#     tmp <- cy %>%
#       mutate(
#         eb = (.data[[Y]] - eb_mean) / eb_sd
#       )
#
#     # Stata:
#     # reg eb treatment_nlonly if treatment_control == 1 | treatment_nlonly == 1, r
#     #
#     # This estimates the newsletter-only / control-plus-encouragement contrast.
#     fit_nl <- lm(
#       eb ~ treatment_nlonly,
#       data = tmp %>%
#         filter(
#           treatment_control == 1 |
#             treatment_nlonly == 1
#         )
#     )
#
#     vc_nl <- sandwich::vcovHC(fit_nl, type = "HC1")
#
#     nl_coef <- lmtest::coeftest(fit_nl, vcov. = vc_nl) %>%
#       broom::tidy() %>%
#       filter(term == "treatment_nlonly")
#
#     tibble(
#       variable = Y,
#       vnum_raw = v,
#       rnl = nl_coef$estimate,
#       rse = nl_coef$std.error,
#       rcil = nl_coef$estimate - 1.9 * nl_coef$std.error,
#       rcih = nl_coef$estimate + 1.9 * nl_coef$std.error
#     )
#   }
# )
#
# # Same vnum spacing as the Stata code.
# fig_dotplot_cce_dat <- fig_dotplot_cce_dat %>%
#   mutate(
#     vnum = vnum_raw,
#     vnum = ifelse(between(vnum, 110, 117), vnum + 35, vnum),
#     vnum = ifelse(between(vnum, 102, 109), vnum + 34, vnum),
#     vnum = ifelse(between(vnum, 97, 101), vnum + 33, vnum),
#     vnum = ifelse(between(vnum, 96, 96), vnum + 32, vnum),
#     vnum = ifelse(between(vnum, 93, 95), vnum + 31, vnum),
#     vnum = ifelse(between(vnum, 91, 92), vnum + 30, vnum),
#     vnum = ifelse(between(vnum, 90, 90), vnum + 29, vnum),
#     vnum = ifelse(between(vnum, 85, 89), vnum + 28, vnum),
#     vnum = ifelse(between(vnum, 83, 84), vnum + 27, vnum),
#     vnum = ifelse(between(vnum, 79, 82), vnum + 25, vnum),
#     vnum = ifelse(between(vnum, 78, 78), vnum + 24, vnum),
#     vnum = ifelse(between(vnum, 74, 77), vnum + 23, vnum),
#     vnum = ifelse(between(vnum, 71, 73), vnum + 22, vnum),
#     vnum = ifelse(between(vnum, 67, 70), vnum + 21, vnum),
#     vnum = ifelse(between(vnum, 64, 66), vnum + 20, vnum),
#     vnum = ifelse(between(vnum, 61, 63), vnum + 18, vnum),
#     vnum = ifelse(between(vnum, 58, 60), vnum + 17, vnum),
#     vnum = ifelse(between(vnum, 55, 57), vnum + 16, vnum),
#     vnum = ifelse(between(vnum, 52, 54), vnum + 15, vnum),
#     vnum = ifelse(between(vnum, 49, 51), vnum + 13, vnum),
#     vnum = ifelse(between(vnum, 48, 48), vnum + 12, vnum),
#     vnum = ifelse(between(vnum, 42, 47), vnum + 11, vnum),
#     vnum = ifelse(between(vnum, 37, 41), vnum + 10, vnum),
#     vnum = ifelse(between(vnum, 32, 36), vnum + 9, vnum),
#     vnum = ifelse(between(vnum, 24, 31), vnum + 8, vnum),
#     vnum = ifelse(between(vnum, 20, 23), vnum + 6, vnum),
#     vnum = ifelse(between(vnum, 16, 19), vnum + 5, vnum),
#     vnum = ifelse(between(vnum, 10, 15), vnum + 4, vnum),
#     vnum = ifelse(between(vnum, 9, 9), vnum + 3, vnum),
#     vnum = ifelse(between(vnum, 8, 8), vnum + 2, vnum),
#     vnum = ifelse(between(vnum, 4, 7), vnum + 1, vnum),
#     label = fig_dotplot_labels[as.character(vnum)]
#   )
#
# figure_dotplot_master_w3_cce <- ggplot(fig_dotplot_cce_dat, aes(x = rnl, y = vnum)) +
#   geom_vline(xintercept = 0) +
#   geom_segment(
#     aes(
#       x = rcil,
#       xend = rcih,
#       y = vnum,
#       yend = vnum
#     )
#   ) +
#   geom_point(size = 1.8, shape = 15) +
#   geom_hline(
#     yintercept = c(4, 9, 13, 20, 46, 61, 70, 74, 78, 87, 96, 101, 103,
#                    112, 120, 123, 127, 129, 135),
#     linetype = "longdash",
#     linewidth = .2
#   ) +
#   geom_hline(
#     yintercept = c(11, 25, 40, 52, 59, 92, 118, 144),
#     linetype = "dashed",
#     linewidth = .2
#   ) +
#   geom_hline(
#     yintercept = c(31, 66, 83, 109),
#     linewidth = .35
#   ) +
#   scale_y_reverse(
#     breaks = as.integer(names(fig_dotplot_labels)),
#     labels = unname(fig_dotplot_labels)
#   ) +
#   labs(
#     x = "Standardized means (Control = 0)",
#     y = NULL
#   ) +
#   theme_bw() +
#   theme(
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.position = "none",
#     axis.text.y = element_text(size = 5)
#   )
#
# ggsave(
#   file.path(outreg_dir, "figure_dotplot_master_w3_cce.pdf"),
#   figure_dotplot_master_w3_cce,
#   width = 4,
#   height = 20
# )
#
# # Figure A.17 -----------------------------------------------------------------
#
# figure_a17_dat <- cy %>%
#   filter(panelmerged_wave3 == 1)
#
# heterogeneity_vars <- c(
#   "gender",
#   "birth_year",
#   "residence_coastal",
#   "hukou_urban",
#   "university_elite",
#   "hs_track_science",
#   "department_ssh",
#   "domestic_english_atleast4",
#   "foreign_english_yes",
#   "travel_hktaiwan",
#   "travel_foreign_yes",
#   "father_edu_hsabove",
#   "work_father_govt",
#   "father_ccp",
#   "mother_edu_hsabove",
#   "work_mother_govt",
#   "mother_ccp",
#   "hh_income",
#   "az_preference_risk",
#   "az_preference_time",
#   "az_preference_altruism",
#   "az_preference_reciprocity",
#   "az_overall_a_w1",
#   "az_overall_b_w1",
#   "az_overall_c_w1",
#   "az_overall_d_w1",
#   "az_overall_e_w1"
# )
#
# heterogeneity_labels <- c(
#   `0`  = "All subjects",
#   `2`  = "Female vs. male",
#   `5`  = "Lower class vs. upper class",
#   `8`  = "Non-coastal vs. coastal",
#   `11` = "Rural vs. urban",
#   `14` = "2nd tier vs. elite univ.",
#   `17` = "Humanities vs. science track",
#   `20` = "Sc/Eng vs. SocS/Hum major",
#   `23` = "Not passed vs. at least Eng Level 4",
#   `26` = "Not taken vs. taken TOEFL/IELTS",
#   `29` = "Not been vs. been to HK/TW",
#   `32` = "Not been vs. been abroad",
#   `35` = "Father below vs. above hs",
#   `38` = "Father not work vs. work for govt.",
#   `41` = "Father not vs. is CCP member",
#   `44` = "Mother below vs. above hs",
#   `47` = "Mother not work vs. work for govt.",
#   `50` = "Mother not vs. is CCP member",
#   `53` = "HH income < vs. > median",
#   `56` = "Risk pref. < vs. > median",
#   `59` = "Time pref. < vs. > median",
#   `62` = "Altruism < vs. > median",
#   `65` = "Recipro. < vs. > median",
#   `68` = "(A) media-related < vs. > median",
#   `71` = "(B) knowledge < vs. > median",
#   `74` = "(C) economic beliefs < vs. > median",
#   `77` = "(D) pol. attitudes < vs. > median",
#   `80` = "(E) behaviors < vs. > median"
# )
#
# get_group_ae_coef <- function(data, Y, i, subset_expr) {
#
#   dat <- data %>%
#     filter({{ subset_expr }})
#
#   fit <- lm(
#     reformulate("treatment_vpnnl", response = Y),
#     data = dat
#   )
#
#   vc <- sandwich::vcovHC(fit, type = "HC1")
#
#   out <- lmtest::coeftest(fit, vcov. = vc) %>%
#     broom::tidy() %>%
#     filter(term == "treatment_vpnnl")
#
#   tibble(
#     i = i,
#     b = out$estimate,
#     se = out$std.error,
#     u = out$estimate + 1.68 * out$std.error,
#     l = out$estimate - 1.68 * out$std.error
#   )
# }
#
# figure_a17_results <- purrr::map_dfr(
#   c("az_overall_all_w3"),
#   function(Y) {
#
#     # Stata: all subjects, i = 0.
#     all_subjects <- get_group_ae_coef(
#       data = figure_a17_dat,
#       Y = Y,
#       i = 0,
#       subset_expr = treatment_user != 1
#     )
#
#     # Stata starts subgroup labels at i = 2 and then jumps by 3:
#     # i, i + 1, then i + 2 before the next subgroup pair.
#     subgroup_results <- purrr::imap_dfr(
#       heterogeneity_vars,
#       function(h, idx) {
#
#         base_i <- 2 + (idx - 1) * 3
#         h_var <- paste0("h_", h)
#
#         bind_rows(
#           get_group_ae_coef(
#             data = figure_a17_dat,
#             Y = Y,
#             i = base_i,
#             subset_expr = treatment_user != 1 & .data[[h_var]] == 0
#           ),
#           get_group_ae_coef(
#             data = figure_a17_dat,
#             Y = Y,
#             i = base_i + 1,
#             subset_expr = treatment_user != 1 & .data[[h_var]] == 1
#           )
#         )
#       }
#     )
#
#     bind_rows(all_subjects, subgroup_results) %>%
#       mutate(outcome = Y)
#   }
# )
#
# b_all <- figure_a17_results %>%
#   filter(i == 0) %>%
#   pull(b)
#
# figure_a17_plot <- ggplot(figure_a17_results, aes(x = i, y = b)) +
#   geom_hline(
#     yintercept = b_all,
#     linetype = "dashed"
#   ) +
#   geom_vline(
#     xintercept = 1,
#     linewidth = .4
#   ) +
#   geom_vline(
#     xintercept = c(13, 22, 34, 55, 67),
#     linewidth = .35
#   ) +
#   geom_vline(
#     xintercept = c(4, 7, 10, 16, 19, 25, 28, 31, 37, 40, 43, 46,
#                    49, 52, 58, 61, 64, 70, 73, 76, 79),
#     linetype = "longdash",
#     linewidth = .25
#   ) +
#   geom_errorbar(
#     aes(ymin = l, ymax = u),
#     width = .15
#   ) +
#   geom_point(
#     data = figure_a17_results %>% filter(i > 0),
#     size = 2
#   ) +
#   geom_point(
#     data = figure_a17_results %>% filter(i == 0),
#     size = 3,
#     shape = 15
#   ) +
#   scale_x_continuous(
#     limits = c(-.5, 81),
#     breaks = as.integer(names(heterogeneity_labels)),
#     labels = unname(heterogeneity_labels)
#   ) +
#   labs(
#     x = NULL,
#     y = "Coefficient on Group-AE effect"
#   ) +
#   theme_bw() +
#   theme(
#     legend.position = "none",
#     panel.grid = element_blank(),
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)
#   )
#
# ggsave(
#   file.path(outreg_dir, "figure_heterogeneity_az_overall_all_w3.pdf"),
#   figure_a17_plot,
#   width = 20,
#   height = 7
# )
#
# # Figure A.18 -----------------------------------------------------------------
#
# library(minpack.lm)
# library(ggplot2)
# library(dplyr)
# library(purrr)
# library(patchwork)
#
# estimate_social_learning_nls <- function(data, Y, wave) {
#
#   rm_new_var <- paste0("soclearning_rm_new_w", wave)
#   y_var <- paste0("news_c_", Y)
#
#   dat <- data %>%
#     transmute(
#       y = .data[[y_var]],
#       soclearning_ownaccess,
#       soclearning_rm_new = .data[[rm_new_var]],
#       soclearning_rm_existing
#     ) %>%
#     filter(
#       soclearning_rm_existing == 0,
#       soclearning_rm_new < 2
#     )
#
#   # Stata:
#   # nl (news_c_Y = {alpha} + ownaccess*{p} + ...)
#   #
#   # R:
#   # nlsLM() is more stable than base nls() with the same starting values.
#   fit <- minpack.lm::nlsLM(
#     y ~ alpha +
#       soclearning_ownaccess * p +
#       (
#         1 -
#           (
#             (alpha + p) *
#               (
#                 1 -
#                   (1 - soclearning_ownaccess) * q0 -
#                   soclearning_ownaccess * q1
#               ) +
#               (1 - alpha - p)
#           ) ^ soclearning_rm_new
#       ),
#     data = dat,
#     start = list(
#       alpha = 0.5,
#       p = 0.25,
#       q0 = 0.2,
#       q1 = 0.1
#     )
#   )
#
#   coefs <- coef(fit)
#
#   label <- variable_labels[[paste0("l_", Y)]]
#
#   if (is.null(label)) {
#     label <- Y
#   }
#
#   tibble(
#     news = label,
#     alpha = coefs[["alpha"]],
#     p = coefs[["p"]],
#     q0 = coefs[["q0"]],
#     q1 = coefs[["q1"]]
#   )
# }
#
# slest_w2 <- map_dfr(
#   c("panamapapers", "tenyearshk", "stockcrash", "economistcensor"),
#   ~ estimate_social_learning_nls(cy, .x, wave = 2)
# )
#
# slest_w3 <- map_dfr(
#   c("coalprod", "trumpchina", "xiaojianhua", "xijiangcar",
#     "chinanorway", "womenrights", "hkceelection"),
#   ~ estimate_social_learning_nls(cy, .x, wave = 3)
# )
#
# slest_compiled <- bind_rows(slest_w2, slest_w3) %>%
#   arrange(p) %>%
#   mutate(
#     i = row_number(),
#     p_l = p - 0.003,
#     p_h = p + 0.003,
#     p_u = p + 0.009,
#     p_v = p + 0.012,
#     label_p = case_when(
#       i == 3  ~ p_u,
#       i == 5  ~ p_h,
#       i == 9  ~ p_l,
#       i == 10 ~ p_v,
#       i == 11 ~ p_h,
#       TRUE ~ p
#     ),
#     label_q = ifelse(i == 1, q1, q0)
#   )
#
# figure_sociallearning_est <- ggplot(slest_compiled, aes(x = p)) +
#   geom_segment(
#     aes(xend = p, y = q0, yend = q1),
#     linetype = "dashed"
#   ) +
#   geom_point(aes(y = q0), shape = 21, size = 3) +
#   geom_point(aes(y = q1), shape = 23, size = 2.5) +
#   geom_text(
#     aes(x = label_p, y = label_q, label = news),
#     angle = 90,
#     hjust = 0,
#     size = 2.5,
#     colour = "grey40"
#   ) +
#   scale_y_continuous(
#     limits = c(-0.1, 0.45),
#     breaks = seq(-0.1, 0.4, 0.1)
#   ) +
#   scale_x_continuous(
#     limits = c(0.05, 0.35),
#     breaks = seq(0.05, 0.35, 0.05)
#   ) +
#   labs(
#     x = "Direct learning rate: p",
#     y = "Social transmission rates: q"
#   ) +
#   theme_bw() +
#   theme(
#     panel.grid = element_blank(),
#     legend.position = "none"
#   )
#
# ggsave(
#   file.path(outreg_dir, "figure_sociallearning_est.pdf"),
#   figure_sociallearning_est,
#   width = 6,
#   height = 4
# )
#
# # Figure A.19 -----------------------------------------------------------------
#
# time_labels <- c(
#   news_time_lastclick = "Time (sec) spent on news quizzes",
#   news_time_totalclick = "# clicks on news quizzes",
#   pppr_time_lastclick = "Time (sec) spent on notable figures & protests",
#   pppr_time_totalclick = "# clicks on notable figures & protests"
# )
#
# top_bottom_code_time <- function(data) {
#   data %>%
#     mutate(
#       news_time_lastclick = ifelse(news_time_lastclick > 300, 300, news_time_lastclick),
#       news_time_totalclick = ifelse(news_time_totalclick > 20, 20, news_time_totalclick),
#       pppr_time_lastclick = ifelse(pppr_time_lastclick > 200, 200, pppr_time_lastclick),
#       pppr_time_lastclick = ifelse(pppr_time_lastclick < 20, 20, pppr_time_lastclick),
#       pppr_time_totalclick = ifelse(pppr_time_totalclick > 30, 30, pppr_time_totalclick)
#     )
# }
#
# make_time_boxplot <- function(data, g, suffix) {
#
#   plot_dat <- data %>%
#     filter(!is.na(treatment_main)) %>%
#     mutate(
#       treatment_main_lbl = factor(
#         treatment_main,
#         levels = 1:4,
#         labels = c("Control", "Access", "Access + Encour.", "Existing users")
#       )
#     )
#
#   p <- ggplot(plot_dat, aes(x = treatment_main_lbl, y = .data[[g]])) +
#     geom_boxplot() +
#     labs(
#       x = NULL,
#       y = time_labels[[g]]
#     ) +
#     theme_bw() +
#     theme(
#       panel.grid = element_blank(),
#       axis.text.x = element_text(angle = 30, hjust = 1)
#     )
#
#   ggsave(
#     file.path(outreg_dir, paste0("figure_time_", g, "_", suffix, ".pdf")),
#     p,
#     width = 5,
#     height = 4
#   )
#
#   p
# }
#
# # Descriptive summaries by treatment arm.
# figure_a19_summaries <- cy %>%
#   filter(!is.na(treatment_main)) %>%
#   group_by(treatment_main) %>%
#   summarise(
#     across(
#       c(
#         news_time_lastclick,
#         news_time_totalclick,
#         pppr_time_lastclick,
#         pppr_time_totalclick
#       ),
#       list(
#         mean = ~ mean(.x, na.rm = TRUE),
#         sd = ~ sd(.x, na.rm = TRUE),
#         p50 = ~ median(.x, na.rm = TRUE),
#         min = ~ min(.x, na.rm = TRUE),
#         max = ~ max(.x, na.rm = TRUE)
#       )
#     ),
#     .groups = "drop"
#   )
#
# # Everyone.
# time_everyone_dat <- cy %>%
#   top_bottom_code_time()
#
# time_everyone_plots <- map(
#   names(time_labels),
#   ~ make_time_boxplot(time_everyone_dat, .x, "everyone")
# )
#
# # Conditional on answering news quiz correctly.
# time_news_correct_dat <- cy %>%
#   top_bottom_code_time() %>%
#   filter(news_totalcorrect_w2 >= 6)
#
# time_news_correct_plots <- map(
#   c("news_time_lastclick", "news_time_totalclick"),
#   ~ make_time_boxplot(time_news_correct_dat, .x, "correct")
# )
#
# # Conditional on hearing enough people/protest items.
# time_pppr_correct_dat <- cy %>%
#   top_bottom_code_time() %>%
#   filter(
#     people_totalheard_w2 >= 3,
#     protest_totalheard_w2 >= 3
#   )
#
# time_pppr_correct_plots <- map(
#   c("pppr_time_lastclick", "pppr_time_totalclick"),
#   ~ make_time_boxplot(time_pppr_correct_dat, .x, "correct")
# )
#
# figure_time_news_everyone <-
#   time_everyone_plots[[1]] +
#   time_everyone_plots[[2]] +
#   plot_annotation(title = "Panel A: among all participants")
#
# ggsave(
#   file.path(outreg_dir, "figure_time_news_everyone.pdf"),
#   figure_time_news_everyone,
#   width = 13,
#   height = 5
# )
#
# figure_time_news_correct <-
#   time_news_correct_plots[[1]] +
#   time_news_correct_plots[[2]] +
#   plot_annotation(title = "Panel B: among those who answered > half questions correctly")
#
# ggsave(
#   file.path(outreg_dir, "figure_time_news_correct.pdf"),
#   figure_time_news_correct,
#   width = 13,
#   height = 5
# )
#
# # Table A.1 -------------------------------------------------------------------
#
# summstats_vars <- c(
#   var_demog_reg_personal,
#   "az_demographics_personal",
#   var_demog_reg_education,
#   "az_demographics_education",
#   var_demog_reg_english,
#   "az_demographics_english",
#   var_demog_reg_travel,
#   "az_demographics_travel",
#   var_demog_reg_household,
#   "az_demographics_household",
#   var_preference_risk,
#   "az_preference_risk",
#   var_preference_time,
#   "az_preference_time",
#   var_preference_altruism,
#   "az_preference_altruism",
#   var_preference_reciprocity,
#   "az_preference_reciprocity",
#   "treatment_control",
#   "treatment_vpnonly",
#   "treatment_nlonly",
#   "treatment_vpnnl",
#   "treatment_user"
# )
#
# summstats_w1 <- purrr::map_dfr(summstats_vars, function(X) {
#
#   treat_stats <- purrr::map_dfr(1:5, function(t) {
#
#     tmp <- cy %>%
#       filter(treatment_master == t)
#
#     tibble(
#       treatment_master = t,
#       mu = mean(tmp[[X]], na.rm = TRUE),
#       sd = sd(tmp[[X]], na.rm = TRUE),
#       N = sum(!is.na(tmp[[X]]))
#     )
#   })
#
#   anova_dat <- cy %>%
#     filter(treatment_master != 5)
#
#   fit <- aov(
#     reformulate(
#       "factor(treatment_master)",
#       response = X
#     ),
#     data = anova_dat
#   )
#
#   anova_tbl <- summary(fit)[[1]]
#
#   tibble(
#     variable = X,
#
#     mu_0 = mean(cy[[X]], na.rm = TRUE),
#     sd_0 = sd(cy[[X]], na.rm = TRUE),
#     N_0 = sum(!is.na(cy[[X]])),
#
#     mu_1 = treat_stats$mu[treat_stats$treatment_master == 1],
#     sd_1 = treat_stats$sd[treat_stats$treatment_master == 1],
#     N_1 = treat_stats$N[treat_stats$treatment_master == 1],
#
#     mu_2 = treat_stats$mu[treat_stats$treatment_master == 2],
#     sd_2 = treat_stats$sd[treat_stats$treatment_master == 2],
#     N_2 = treat_stats$N[treat_stats$treatment_master == 2],
#
#     mu_3 = treat_stats$mu[treat_stats$treatment_master == 3],
#     sd_3 = treat_stats$sd[treat_stats$treatment_master == 3],
#     N_3 = treat_stats$N[treat_stats$treatment_master == 3],
#
#     mu_4 = treat_stats$mu[treat_stats$treatment_master == 4],
#     sd_4 = treat_stats$sd[treat_stats$treatment_master == 4],
#     N_4 = treat_stats$N[treat_stats$treatment_master == 4],
#
#     mu_5 = treat_stats$mu[treat_stats$treatment_master == 5],
#     sd_5 = treat_stats$sd[treat_stats$treatment_master == 5],
#     N_5 = treat_stats$N[treat_stats$treatment_master == 5],
#
#     f_stat = anova_tbl$`F value`[1],
#     p_val = anova_tbl$`Pr(>F)`[1]
#   )
# })
#
# write.csv(
#   summstats_w1,
#   file.path(outreg_dir, "summstats_w1.csv"),
#   row.names = FALSE
# )
#
#
# # Table A.2 -------------------------------------------------------------------
#
# cy_w2 <- cy %>%
#   filter(panelmerged_wave2 == 1)
#
# summstats_w2 <- purrr::map_dfr(summstats_vars, function(X) {
#
#   treat_stats <- purrr::map_dfr(1:5, function(t) {
#
#     tmp <- cy_w2 %>%
#       filter(treatment_master == t)
#
#     tibble(
#       treatment_master = t,
#       mu = mean(tmp[[X]], na.rm = TRUE),
#       sd = sd(tmp[[X]], na.rm = TRUE),
#       N = sum(!is.na(tmp[[X]]))
#     )
#   })
#
#   anova_dat <- cy_w2 %>%
#     filter(treatment_master != 5)
#
#   fit <- aov(
#     reformulate(
#       "factor(treatment_master)",
#       response = X
#     ),
#     data = anova_dat
#   )
#
#   anova_tbl <- summary(fit)[[1]]
#
#   tibble(
#     variable = X,
#
#     mu_0 = mean(cy_w2[[X]], na.rm = TRUE),
#     sd_0 = sd(cy_w2[[X]], na.rm = TRUE),
#     N_0 = sum(!is.na(cy_w2[[X]])),
#
#     mu_1 = treat_stats$mu[treat_stats$treatment_master == 1],
#     sd_1 = treat_stats$sd[treat_stats$treatment_master == 1],
#     N_1 = treat_stats$N[treat_stats$treatment_master == 1],
#
#     mu_2 = treat_stats$mu[treat_stats$treatment_master == 2],
#     sd_2 = treat_stats$sd[treat_stats$treatment_master == 2],
#     N_2 = treat_stats$N[treat_stats$treatment_master == 2],
#
#     mu_3 = treat_stats$mu[treat_stats$treatment_master == 3],
#     sd_3 = treat_stats$sd[treat_stats$treatment_master == 3],
#     N_3 = treat_stats$N[treat_stats$treatment_master == 3],
#
#     mu_4 = treat_stats$mu[treat_stats$treatment_master == 4],
#     sd_4 = treat_stats$sd[treat_stats$treatment_master == 4],
#     N_4 = treat_stats$N[treat_stats$treatment_master == 4],
#
#     mu_5 = treat_stats$mu[treat_stats$treatment_master == 5],
#     sd_5 = treat_stats$sd[treat_stats$treatment_master == 5],
#     N_5 = treat_stats$N[treat_stats$treatment_master == 5],
#
#     f_stat = anova_tbl$`F value`[1],
#     p_val = anova_tbl$`Pr(>F)`[1]
#   )
# })
#
# write.csv(
#   summstats_w2,
#   file.path(outreg_dir, "summstats_w2.csv"),
#   row.names = FALSE
# )
#
# # Table A.3 -------------------------------------------------------------------
#
# cy_w3 <- cy %>%
#   filter(panelmerged_wave3 == 1)
#
# summstats_w3 <- purrr::map_dfr(summstats_vars, function(X) {
#
#   treat_stats <- purrr::map_dfr(1:5, function(t) {
#
#     tmp <- cy_w3 %>%
#       filter(treatment_master == t)
#
#     tibble(
#       treatment_master = t,
#       mu = mean(tmp[[X]], na.rm = TRUE),
#       sd = sd(tmp[[X]], na.rm = TRUE),
#       N = sum(!is.na(tmp[[X]]))
#     )
#   })
#
#   anova_dat <- cy_w3 %>%
#     filter(treatment_master != 5)
#
#   fit <- aov(
#     reformulate(
#       "factor(treatment_master)",
#       response = X
#     ),
#     data = anova_dat
#   )
#
#   anova_tbl <- summary(fit)[[1]]
#
#   tibble(
#     variable = X,
#
#     mu_0 = mean(cy_w3[[X]], na.rm = TRUE),
#     sd_0 = sd(cy_w3[[X]], na.rm = TRUE),
#     N_0 = sum(!is.na(cy_w3[[X]])),
#
#     mu_1 = treat_stats$mu[treat_stats$treatment_master == 1],
#     sd_1 = treat_stats$sd[treat_stats$treatment_master == 1],
#     N_1 = treat_stats$N[treat_stats$treatment_master == 1],
#
#     mu_2 = treat_stats$mu[treat_stats$treatment_master == 2],
#     sd_2 = treat_stats$sd[treat_stats$treatment_master == 2],
#     N_2 = treat_stats$N[treat_stats$treatment_master == 2],
#
#     mu_3 = treat_stats$mu[treat_stats$treatment_master == 3],
#     sd_3 = treat_stats$sd[treat_stats$treatment_master == 3],
#     N_3 = treat_stats$N[treat_stats$treatment_master == 3],
#
#     mu_4 = treat_stats$mu[treat_stats$treatment_master == 4],
#     sd_4 = treat_stats$sd[treat_stats$treatment_master == 4],
#     N_4 = treat_stats$N[treat_stats$treatment_master == 4],
#
#     mu_5 = treat_stats$mu[treat_stats$treatment_master == 5],
#     sd_5 = treat_stats$sd[treat_stats$treatment_master == 5],
#     N_5 = treat_stats$N[treat_stats$treatment_master == 5],
#
#     f_stat = anova_tbl$`F value`[1],
#     p_val = anova_tbl$`Pr(>F)`[1]
#   )
# })
#
# write.csv(
#   summstats_w3,
#   file.path(outreg_dir, "summstats_w3.csv"),
#   row.names = FALSE
# )
#
#
# # Table A.4 -------------------------------------------------------------------
#
# panelmerged_wave2_temp <- cy %>%
#   filter(panelmerged_wave2 == 1) %>%
#   mutate(wave2 = 1) %>%
#   select(-treatment_master)
#
# panelmerged_wave3_temp <- cy %>%
#   filter(panelmerged_wave3 == 1) %>%
#   mutate(wave3 = 1) %>%
#   select(-treatment_master)
#
# attrition_dat <- bind_rows(
#   cy,
#   panelmerged_wave2_temp
# ) %>%
#   mutate(
#     wave2 = replace_na(wave2, 0)
#   ) %>%
#   bind_rows(
#     panelmerged_wave3_temp
#   ) %>%
#   mutate(
#     wave3 = ifelse(is.na(wave3) & wave2 != 1, 0, wave3)
#   )
#
# attrition_vars <- c(
#   var_info_ranking_reg_w1,
#   var_info_freq_reg_w1,
#   var_media_valuation_reg_w1,
#   "az_belief_media_value_w1",
#   var_media_trust_reg_w1,
#   "az_belief_media_trust_w1",
#   var_censor_level_reg_w1,
#   var_censor_justif_reg_w1,
#   "az_belief_media_justif_w1",
#   var_censor_driver_reg_dom_w1,
#   var_censor_driver_reg_for_w1,
#   var_knowledge_news_reg_cen_w1,
#   "news_perccor_cen_w1",
#   var_knowledge_news_reg_unc_w1,
#   "news_perccor_unc_w1",
#   var_knowledge_prot_reg_chi_w1,
#   "protest_pcheard_china_w1",
#   var_knowledge_prot_reg_for_w1,
#   "protest_pcheard_foreign_w1",
#   var_knowledge_prot_reg_fak_w1,
#   var_knowledge_meta_reg_w1,
#   "az_knowledge_meta_w1",
#   var_econ_guess_reg_cn_perf_w1,
#   "az_belief_econ_perf_cn_w1",
#   var_econ_guess_reg_cn_conf_w1,
#   "az_belief_econ_conf_cn_w1",
#   var_demand_change_reg_w1,
#   "az_belief_instchange_w1",
#   var_trust_inst_reg_govt_w1,
#   "az_belief_trust_govt_w1",
#   var_trust_inst_reg_foreign_w1,
#   "az_belief_trust_foreign_w1",
#   var_eval_govt_reg_w1,
#   "az_belief_evalgovt_w1",
#   var_democracy_reg_fp_w1,
#   "az_belief_democracy_fp_w1",
#   var_willing_fight_reg_w1,
#   "az_belief_willing_w1",
#   var_socialinteract_reg_w1,
#   "az_var_socialinteract_w1",
#   var_polparticipation_reg_pf_w1,
#   "az_var_polparticipation_w1",
#   var_stock_invest_reg_w1,
#   var_planaftergrad_reg_w1,
#   var_career_sector_reg_w1,
#   var_career_loc_reg_w1,
#   var_demog_reg_personal,
#   "az_demographics_personal",
#   var_demog_reg_education,
#   "az_demographics_education",
#   var_demog_reg_english,
#   "az_demographics_english",
#   var_demog_reg_travel,
#   "az_demographics_travel",
#   var_demog_reg_household,
#   "az_demographics_household",
#   var_preference_risk,
#   "az_preference_risk",
#   var_preference_time,
#   "az_preference_time",
#   var_preference_altruism,
#   "az_preference_altruism",
#   var_preference_reciprocity,
#   "az_preference_reciprocity",
#   "treatment_control",
#   "treatment_vpnonly",
#   "treatment_nlonly",
#   "treatment_vpnnl",
#   "treatment_user"
# )
#
# attrition_appendix <- purrr::map_dfr(attrition_vars, function(X) {
#
#   w1 <- attrition_dat %>%
#     filter(wave2 == 0, wave3 == 0)
#
#   w2 <- attrition_dat %>%
#     filter(
#       panelmerged_wave2 == 1,
#       wave2 == 0,
#       wave3 == 0
#     )
#
#   w3 <- attrition_dat %>%
#     filter(
#       panelmerged_wave3 == 1,
#       wave2 == 0,
#       wave3 == 0
#     )
#
#   attrition_w2_test <- t.test(
#     reformulate("wave2", response = X),
#     data = attrition_dat
#   )
#
#   attrition_w3_test <- t.test(
#     reformulate("wave3", response = X),
#     data = attrition_dat
#   )
#
#   tibble(
#     variable = X,
#
#     mu_w1 = mean(w1[[X]], na.rm = TRUE),
#     sd_w1 = sd(w1[[X]], na.rm = TRUE),
#     N_w1 = sum(!is.na(w1[[X]])),
#
#     mu_w2 = mean(w2[[X]], na.rm = TRUE),
#     sd_w2 = sd(w2[[X]], na.rm = TRUE),
#     N_w2 = sum(!is.na(w2[[X]])),
#
#     attrition_w2_pvalue = attrition_w2_test$p.value,
#
#     mu_w3 = mean(w3[[X]], na.rm = TRUE),
#     sd_w3 = sd(w3[[X]], na.rm = TRUE),
#     N_w3 = sum(!is.na(w3[[X]])),
#
#     attrition_w3_pvalue = attrition_w3_test$p.value
#   )
# })
#
# write.csv(
#   attrition_appendix,
#   file.path(outreg_dir, "attrition_appendix.csv"),
#   row.names = FALSE
# )
#
# # Table A.5 -------------------------------------------------------------------
#
# attrition_outcomes <- c(
#   "az_overall_a_w1",
#   "az_overall_b_w1",
#   "az_overall_c_w1",
#   "az_overall_d_w1",
#   "az_overall_e_w1"
# )
#
# table_a5_dat <- cy %>%
#   mutate(
#     attrited_wave3 = 1 - panelmerged_wave3
#   )
#
# table_a5 <- purrr::map_dfr(attrition_outcomes, function(a) {
#
#   # Stata generated interactions:
#   # gen a_i = a * (treatment_master == i)
#   #
#   # Only interactions 2--5 enter the regressions because treatment 1 is omitted.
#   table_a5_dat <- table_a5_dat %>%
#     mutate(
#       a_2 = .data[[a]] * as.integer(treatment_master == 2),
#       a_3 = .data[[a]] * as.integer(treatment_master == 3),
#       a_4 = .data[[a]] * as.integer(treatment_master == 4),
#       a_5 = .data[[a]] * as.integer(treatment_master == 5)
#     )
#
#   fit <- lm(
#     attrited_wave3 ~
#       .data[[a]] +
#       treatment_vpnonly +
#       treatment_nlonly +
#       treatment_vpnnl +
#       treatment_user +
#       a_2 + a_3 + a_4 + a_5,
#     data = table_a5_dat
#   )
#
#   vc <- sandwich::vcovHC(fit, type = "HC1")
#
#   lmtest::coeftest(fit, vcov. = vc) %>%
#     broom::tidy() %>%
#     mutate(outcome = a)
# })
#
# write.csv(
#   table_a5,
#   file.path(outreg_dir, "outregs_attrition_prediction.csv"),
#   row.names = FALSE
# )
#
# # Table A.9 -------------------------------------------------------------------
#
# cy_w3 <- cy %>%
#   filter(panelmerged_wave3 == 1)
#
# table_a9_vars_w3_suffix <- c(
#   var_info_ranking_fig_w23,
#   var_info_freq_fig_w123,
#   var_vpn_purchase_reg_w3,
#   var_media_valuation_fig_w123,
#   var_media_trust_fig_w123,
#   var_censor_level_fig_w123,
#   var_censor_justif_fig_w123,
#   var_censor_justif_fig_w23,
#   var_censor_driver_fig_w123,
#   var_knowledge_news_reg_cen_w3,
#   var_knowledge_news_reg_unc_w3,
#   var_knowledge_news_fig_w123,
#   var_knowledge_prot_fig_w123,
#   var_knowledge_prot_fig_w23,
#   var_knowledge_prot_fig_w3,
#   var_knowledge_meta_fig_w123,
#   var_demand_change_fig_w123,
#   var_trust_inst_fig_w123,
#   var_eval_govt_fig_w123,
#   var_democracy_fig_w123,
#   var_willing_fight_fig_123,
#   var_socialinteract_fig_w123,
#   var_polparticipation_fig_w123,
#   var_stock_invest_fig_w123,
#   var_planaftergrad_fig_w123,
#   var_career_sector_fig_w123,
#   var_career_loc_reg_w123
# )
#
# table_a9 <- purrr::map_dfr(table_a9_vars_w3_suffix, function(Y) {
#
#   dv_w3 <- paste0(Y, "_w3")
#   dv_w1 <- paste0(Y, "_w1")
#
#   if (!dv_w3 %in% names(cy_w3)) {
#     return(tibble())
#   }
#
#   out <- list()
#
#   # Raw Diff ---------------------------------------------------------------
#
#   fit_raw <- lm(
#     reformulate(
#       c("treatment_vpnonly", "treatment_nlonly", "treatment_vpnnl"),
#       response = dv_w3
#     ),
#     data = cy_w3 %>%
#       filter(treatment_user != 1)
#   )
#
#   vc_raw <- sandwich::vcovHC(fit_raw, type = "HC1")
#
#   out[["raw"]] <- lmtest::coeftest(fit_raw, vcov. = vc_raw) %>%
#     broom::tidy() %>%
#     filter(
#       term %in% c(
#         "treatment_vpnonly",
#         "treatment_nlonly",
#         "treatment_vpnnl"
#       )
#     ) %>%
#     mutate(
#       variable = Y,
#       outcome = dv_w3,
#       specification = "Raw Diff",
#       mean_dv_all = mean(cy_w3[[dv_w3]][cy_w3$treatment_user != 1], na.rm = TRUE),
#       sd_dv_all = sd(cy_w3[[dv_w3]][cy_w3$treatment_user != 1], na.rm = TRUE),
#       mean_dv_control = mean(cy_w3[[dv_w3]][cy_w3$treatment_control == 1], na.rm = TRUE),
#       sd_dv_control = sd(cy_w3[[dv_w3]][cy_w3$treatment_control == 1], na.rm = TRUE),
#       mean_dv_user = mean(cy_w3[[dv_w3]][cy_w3$treatment_user == 1], na.rm = TRUE),
#       sd_dv_user = sd(cy_w3[[dv_w3]][cy_w3$treatment_user == 1], na.rm = TRUE),
#       p_value_treatment_vpnnl = p.value[term == "treatment_vpnnl"][1]
#     )
#
#   # Control ---------------------------------------------------------------
#
#   fit_control <- lm(
#     reformulate(
#       c(
#         "treatment_vpnonly",
#         "treatment_nlonly",
#         "treatment_vpnnl",
#         var_demog_reg_imbalance
#       ),
#       response = dv_w3
#     ),
#     data = cy_w3 %>%
#       filter(treatment_user != 1)
#   )
#
#   vc_control <- sandwich::vcovHC(fit_control, type = "HC1")
#
#   out[["control"]] <- lmtest::coeftest(fit_control, vcov. = vc_control) %>%
#     broom::tidy() %>%
#     filter(
#       term %in% c(
#         "treatment_vpnonly",
#         "treatment_nlonly",
#         "treatment_vpnnl"
#       )
#     ) %>%
#     mutate(
#       variable = Y,
#       outcome = dv_w3,
#       specification = "Control",
#       mean_dv_all = mean(cy_w3[[dv_w3]][cy_w3$treatment_user != 1], na.rm = TRUE),
#       sd_dv_all = sd(cy_w3[[dv_w3]][cy_w3$treatment_user != 1], na.rm = TRUE),
#       mean_dv_control = mean(cy_w3[[dv_w3]][cy_w3$treatment_control == 1], na.rm = TRUE),
#       sd_dv_control = sd(cy_w3[[dv_w3]][cy_w3$treatment_control == 1], na.rm = TRUE),
#       mean_dv_user = mean(cy_w3[[dv_w3]][cy_w3$treatment_user == 1], na.rm = TRUE),
#       sd_dv_user = sd(cy_w3[[dv_w3]][cy_w3$treatment_user == 1], na.rm = TRUE),
#       p_value_treatment_vpnnl = p.value[term == "treatment_vpnnl"][1]
#     )
#
#   # Panel -----------------------------------------------------------------
#
#   if (dv_w1 %in% names(cy_w3)) {
#
#     fit_panel <- lm(
#       reformulate(
#         c(
#           "treatment_vpnonly",
#           "treatment_nlonly",
#           "treatment_vpnnl",
#           dv_w1
#         ),
#         response = dv_w3
#       ),
#       data = cy_w3 %>%
#         filter(treatment_user != 1)
#     )
#
#     vc_panel <- sandwich::vcovHC(fit_panel, type = "HC1")
#
#     out[["panel"]] <- lmtest::coeftest(fit_panel, vcov. = vc_panel) %>%
#       broom::tidy() %>%
#       filter(
#         term %in% c(
#           "treatment_vpnonly",
#           "treatment_nlonly",
#           "treatment_vpnnl"
#         )
#       ) %>%
#       mutate(
#         variable = Y,
#         outcome = dv_w3,
#         specification = "Panel",
#         mean_dv_all = mean(cy_w3[[dv_w3]][cy_w3$treatment_user != 1], na.rm = TRUE),
#         sd_dv_all = sd(cy_w3[[dv_w3]][cy_w3$treatment_user != 1], na.rm = TRUE),
#         mean_dv_control = mean(cy_w3[[dv_w3]][cy_w3$treatment_control == 1], na.rm = TRUE),
#         sd_dv_control = sd(cy_w3[[dv_w3]][cy_w3$treatment_control == 1], na.rm = TRUE),
#         mean_dv_user = mean(cy_w3[[dv_w3]][cy_w3$treatment_user == 1], na.rm = TRUE),
#         sd_dv_user = sd(cy_w3[[dv_w3]][cy_w3$treatment_user == 1], na.rm = TRUE),
#         p_value_treatment_vpnnl = p.value[term == "treatment_vpnnl"][1]
#       )
#   }
#
#   bind_rows(out)
# })
#
#
# # Outcomes already named with _w3 --------------------------------------------
#
# table_a9_vars_already_w3 <- c(
#   var_knowledge_news_reg_cen_w3,
#   var_knowledge_news_reg_unc_w3
# )
#
# table_a9_extra_w3 <- purrr::map_dfr(table_a9_vars_already_w3, function(Y) {
#
#   if (!Y %in% names(cy_w3)) {
#     return(tibble())
#   }
#
#   out <- list()
#
#   fit_raw <- lm(
#     reformulate(
#       c("treatment_vpnonly", "treatment_nlonly", "treatment_vpnnl"),
#       response = Y
#     ),
#     data = cy_w3 %>%
#       filter(treatment_user != 1)
#   )
#
#   vc_raw <- sandwich::vcovHC(fit_raw, type = "HC1")
#
#   out[["raw"]] <- lmtest::coeftest(fit_raw, vcov. = vc_raw) %>%
#     broom::tidy() %>%
#     filter(
#       term %in% c(
#         "treatment_vpnonly",
#         "treatment_nlonly",
#         "treatment_vpnnl"
#       )
#     ) %>%
#     mutate(
#       variable = Y,
#       outcome = Y,
#       specification = "Raw Diff",
#       mean_dv_all = mean(cy_w3[[Y]][cy_w3$treatment_user != 1], na.rm = TRUE),
#       sd_dv_all = sd(cy_w3[[Y]][cy_w3$treatment_user != 1], na.rm = TRUE),
#       mean_dv_control = mean(cy_w3[[Y]][cy_w3$treatment_control == 1], na.rm = TRUE),
#       sd_dv_control = sd(cy_w3[[Y]][cy_w3$treatment_control == 1], na.rm = TRUE),
#       mean_dv_user = mean(cy_w3[[Y]][cy_w3$treatment_user == 1], na.rm = TRUE),
#       sd_dv_user = sd(cy_w3[[Y]][cy_w3$treatment_user == 1], na.rm = TRUE),
#       p_value_treatment_vpnnl = p.value[term == "treatment_vpnnl"][1]
#     )
#
#   fit_control <- lm(
#     reformulate(
#       c(
#         "treatment_vpnonly",
#         "treatment_nlonly",
#         "treatment_vpnnl",
#         var_demog_reg_imbalance
#       ),
#       response = Y
#     ),
#     data = cy_w3 %>%
#       filter(treatment_user != 1)
#   )
#
#   vc_control <- sandwich::vcovHC(fit_control, type = "HC1")
#
#   out[["control"]] <- lmtest::coeftest(fit_control, vcov. = vc_control) %>%
#     broom::tidy() %>%
#     filter(
#       term %in% c(
#         "treatment_vpnonly",
#         "treatment_nlonly",
#         "treatment_vpnnl"
#       )
#     ) %>%
#     mutate(
#       variable = Y,
#       outcome = Y,
#       specification = "Control",
#       mean_dv_all = mean(cy_w3[[Y]][cy_w3$treatment_user != 1], na.rm = TRUE),
#       sd_dv_all = sd(cy_w3[[Y]][cy_w3$treatment_user != 1], na.rm = TRUE),
#       mean_dv_control = mean(cy_w3[[Y]][cy_w3$treatment_control == 1], na.rm = TRUE),
#       sd_dv_control = sd(cy_w3[[Y]][cy_w3$treatment_control == 1], na.rm = TRUE),
#       mean_dv_user = mean(cy_w3[[Y]][cy_w3$treatment_user == 1], na.rm = TRUE),
#       sd_dv_user = sd(cy_w3[[Y]][cy_w3$treatment_user == 1], na.rm = TRUE),
#       p_value_treatment_vpnnl = p.value[term == "treatment_vpnnl"][1]
#     )
#
#   bind_rows(out)
# })
#
#
# # Economic guess outcomes ------------------------------------------------------
#
# table_a9_vars_econ <- c(
#   var_econ_guess_fig_w123,
#   var_econ_guess_fig_w23
# )
#
# table_a9_econ <- purrr::map_dfr(table_a9_vars_econ, function(Y) {
#
#   dv_w3 <- paste0(Y, "_w3")
#   dv_w1 <- paste0(Y, "_w1")
#
#   if (!dv_w3 %in% names(cy_w3)) {
#     return(tibble())
#   }
#
#   out <- list()
#
#   fit_raw <- lm(
#     reformulate(
#       c("treatment_vpnonly", "treatment_nlonly", "treatment_vpnnl"),
#       response = dv_w3
#     ),
#     data = cy_w3 %>%
#       filter(treatment_user != 1)
#   )
#
#   vc_raw <- sandwich::vcovHC(fit_raw, type = "HC1")
#
#   out[["raw"]] <- lmtest::coeftest(fit_raw, vcov. = vc_raw) %>%
#     broom::tidy() %>%
#     filter(
#       term %in% c(
#         "treatment_vpnonly",
#         "treatment_nlonly",
#         "treatment_vpnnl"
#       )
#     ) %>%
#     mutate(
#       variable = Y,
#       outcome = dv_w3,
#       specification = "Raw Diff",
#       mean_dv_all = mean(cy_w3[[dv_w3]][cy_w3$treatment_user != 1], na.rm = TRUE),
#       sd_dv_all = sd(cy_w3[[dv_w3]][cy_w3$treatment_user != 1], na.rm = TRUE),
#       mean_dv_control = mean(cy_w3[[dv_w3]][cy_w3$treatment_control == 1], na.rm = TRUE),
#       sd_dv_control = sd(cy_w3[[dv_w3]][cy_w3$treatment_control == 1], na.rm = TRUE),
#       mean_dv_user = mean(cy_w3[[dv_w3]][cy_w3$treatment_user == 1], na.rm = TRUE),
#       sd_dv_user = sd(cy_w3[[dv_w3]][cy_w3$treatment_user == 1], na.rm = TRUE),
#       p_value_treatment_vpnnl = p.value[term == "treatment_vpnnl"][1],
#       decimal_places = 8
#     )
#
#   fit_control <- lm(
#     reformulate(
#       c(
#         "treatment_vpnonly",
#         "treatment_nlonly",
#         "treatment_vpnnl",
#         var_demog_reg_imbalance
#       ),
#       response = dv_w3
#     ),
#     data = cy_w3 %>%
#       filter(treatment_user != 1)
#   )
#
#   vc_control <- sandwich::vcovHC(fit_control, type = "HC1")
#
#   out[["control"]] <- lmtest::coeftest(fit_control, vcov. = vc_control) %>%
#     broom::tidy() %>%
#     filter(
#       term %in% c(
#         "treatment_vpnonly",
#         "treatment_nlonly",
#         "treatment_vpnnl"
#       )
#     ) %>%
#     mutate(
#       variable = Y,
#       outcome = dv_w3,
#       specification = "Control",
#       mean_dv_all = mean(cy_w3[[dv_w3]][cy_w3$treatment_user != 1], na.rm = TRUE),
#       sd_dv_all = sd(cy_w3[[dv_w3]][cy_w3$treatment_user != 1], na.rm = TRUE),
#       mean_dv_control = mean(cy_w3[[dv_w3]][cy_w3$treatment_control == 1], na.rm = TRUE),
#       sd_dv_control = sd(cy_w3[[dv_w3]][cy_w3$treatment_control == 1], na.rm = TRUE),
#       mean_dv_user = mean(cy_w3[[dv_w3]][cy_w3$treatment_user == 1], na.rm = TRUE),
#       sd_dv_user = sd(cy_w3[[dv_w3]][cy_w3$treatment_user == 1], na.rm = TRUE),
#       p_value_treatment_vpnnl = p.value[term == "treatment_vpnnl"][1],
#       decimal_places = 8
#     )
#
#   if (dv_w1 %in% names(cy_w3)) {
#
#     fit_panel <- lm(
#       reformulate(
#         c(
#           "treatment_vpnonly",
#           "treatment_nlonly",
#           "treatment_vpnnl",
#           dv_w1
#         ),
#         response = dv_w3
#       ),
#       data = cy_w3 %>%
#         filter(treatment_user != 1)
#     )
#
#     vc_panel <- sandwich::vcovHC(fit_panel, type = "HC1")
#
#     out[["panel"]] <- lmtest::coeftest(fit_panel, vcov. = vc_panel) %>%
#       broom::tidy() %>%
#       filter(
#         term %in% c(
#           "treatment_vpnonly",
#           "treatment_nlonly",
#           "treatment_vpnnl"
#         )
#       ) %>%
#       mutate(
#         variable = Y,
#         outcome = dv_w3,
#         specification = "Panel",
#         mean_dv_all = mean(cy_w3[[dv_w3]][cy_w3$treatment_user != 1], na.rm = TRUE),
#         sd_dv_all = sd(cy_w3[[dv_w3]][cy_w3$treatment_user != 1], na.rm = TRUE),
#         mean_dv_control = mean(cy_w3[[dv_w3]][cy_w3$treatment_control == 1], na.rm = TRUE),
#         sd_dv_control = sd(cy_w3[[dv_w3]][cy_w3$treatment_control == 1], na.rm = TRUE),
#         mean_dv_user = mean(cy_w3[[dv_w3]][cy_w3$treatment_user == 1], na.rm = TRUE),
#         sd_dv_user = sd(cy_w3[[dv_w3]][cy_w3$treatment_user == 1], na.rm = TRUE),
#         p_value_treatment_vpnnl = p.value[term == "treatment_vpnnl"][1],
#         decimal_places = 8
#       )
#   }
#
#   bind_rows(out)
# })
#
#
# table_a9 <- bind_rows(
#   table_a9,
#   table_a9_extra_w3,
#   table_a9_econ
# )
#
# write.csv(
#   table_a9,
#   file.path(outreg_dir, "outregs_wave3_panelregression.csv"),
#   row.names = FALSE
# )
#
# # Table A.11 ------------------------------------------------------------------
#
# cy_w3 <- cy %>%
#   filter(panelmerged_wave3 == 1)
#
# table_a11 <- purrr::map_dfr(c("a", "b", "c", "d", "e"), function(Y) {
#
#   dv_w3 <- paste0("az_overall_", Y, "_w3")
#   dv_w1 <- paste0("az_overall_", Y, "_w1")
#
#   out <- list()
#
#   # Raw diff ---------------------------------------------------------------
#
#   fit_raw <- lm(
#     reformulate(
#       c("treatment_vpnonly", "treatment_nlonly", "treatment_vpnnl"),
#       response = dv_w3
#     ),
#     data = cy_w3 %>%
#       filter(treatment_user != 1)
#   )
#
#   vc_raw <- sandwich::vcovHC(fit_raw, type = "HC1")
#
#   out[["raw"]] <- lmtest::coeftest(fit_raw, vcov. = vc_raw) %>%
#     broom::tidy() %>%
#     filter(
#       term %in% c(
#         "treatment_vpnonly",
#         "treatment_nlonly",
#         "treatment_vpnnl"
#       )
#     ) %>%
#     mutate(
#       outcome = Y,
#       dv = dv_w3,
#       specification = "Raw diff",
#       mean_dv_all = mean(cy_w3[[dv_w3]][cy_w3$treatment_user != 1], na.rm = TRUE),
#       sd_dv_all = sd(cy_w3[[dv_w3]][cy_w3$treatment_user != 1], na.rm = TRUE),
#       mean_dv_control = mean(cy_w3[[dv_w3]][cy_w3$treatment_control == 1], na.rm = TRUE),
#       sd_dv_control = sd(cy_w3[[dv_w3]][cy_w3$treatment_control == 1], na.rm = TRUE),
#       mean_dv_user = mean(cy_w3[[dv_w3]][cy_w3$treatment_user == 1], na.rm = TRUE),
#       sd_dv_user = sd(cy_w3[[dv_w3]][cy_w3$treatment_user == 1], na.rm = TRUE),
#       p_value_treatment_vpnnl = p.value[term == "treatment_vpnnl"][1]
#     )
#
#   # Control for imbalance --------------------------------------------------
#
#   fit_control <- lm(
#     reformulate(
#       c(
#         "treatment_vpnonly",
#         "treatment_nlonly",
#         "treatment_vpnnl",
#         var_demog_reg_imbalance
#       ),
#       response = dv_w3
#     ),
#     data = cy_w3 %>%
#       filter(treatment_user != 1)
#   )
#
#   vc_control <- sandwich::vcovHC(fit_control, type = "HC1")
#
#   out[["control_imbalance"]] <- lmtest::coeftest(fit_control, vcov. = vc_control) %>%
#     broom::tidy() %>%
#     filter(
#       term %in% c(
#         "treatment_vpnonly",
#         "treatment_nlonly",
#         "treatment_vpnnl"
#       )
#     ) %>%
#     mutate(
#       outcome = Y,
#       dv = dv_w3,
#       specification = "Control for imbalance",
#       mean_dv_all = mean(cy_w3[[dv_w3]][cy_w3$treatment_user != 1], na.rm = TRUE),
#       sd_dv_all = sd(cy_w3[[dv_w3]][cy_w3$treatment_user != 1], na.rm = TRUE),
#       mean_dv_control = mean(cy_w3[[dv_w3]][cy_w3$treatment_control == 1], na.rm = TRUE),
#       sd_dv_control = sd(cy_w3[[dv_w3]][cy_w3$treatment_control == 1], na.rm = TRUE),
#       mean_dv_user = mean(cy_w3[[dv_w3]][cy_w3$treatment_user == 1], na.rm = TRUE),
#       sd_dv_user = sd(cy_w3[[dv_w3]][cy_w3$treatment_user == 1], na.rm = TRUE),
#       p_value_treatment_vpnnl = p.value[term == "treatment_vpnnl"][1]
#     )
#
#   # Control for baseline ---------------------------------------------------
#
#   fit_baseline <- lm(
#     reformulate(
#       c(
#         "treatment_vpnonly",
#         "treatment_nlonly",
#         "treatment_vpnnl",
#         dv_w1
#       ),
#       response = dv_w3
#     ),
#     data = cy_w3 %>%
#       filter(treatment_user != 1)
#   )
#
#   vc_baseline <- sandwich::vcovHC(fit_baseline, type = "HC1")
#
#   out[["control_baseline"]] <- lmtest::coeftest(fit_baseline, vcov. = vc_baseline) %>%
#     broom::tidy() %>%
#     filter(
#       term %in% c(
#         "treatment_vpnonly",
#         "treatment_nlonly",
#         "treatment_vpnnl"
#       )
#     ) %>%
#     mutate(
#       outcome = Y,
#       dv = dv_w3,
#       specification = "Control for baseline",
#       mean_dv_all = mean(cy_w3[[dv_w3]][cy_w3$treatment_user != 1], na.rm = TRUE),
#       sd_dv_all = sd(cy_w3[[dv_w3]][cy_w3$treatment_user != 1], na.rm = TRUE),
#       mean_dv_control = mean(cy_w3[[dv_w3]][cy_w3$treatment_control == 1], na.rm = TRUE),
#       sd_dv_control = sd(cy_w3[[dv_w3]][cy_w3$treatment_control == 1], na.rm = TRUE),
#       mean_dv_user = mean(cy_w3[[dv_w3]][cy_w3$treatment_user == 1], na.rm = TRUE),
#       sd_dv_user = sd(cy_w3[[dv_w3]][cy_w3$treatment_user == 1], na.rm = TRUE),
#       p_value_treatment_vpnnl = p.value[term == "treatment_vpnnl"][1]
#     )
#
#   # Drop politically correct respondents ----------------------------------
#
#   fit_drop <- lm(
#     reformulate(
#       c("treatment_vpnonly", "treatment_nlonly", "treatment_vpnnl"),
#       response = dv_w3
#     ),
#     data = cy_w3 %>%
#       filter(
#         treatment_user != 1,
#         trust_central_govt_w3 != 10
#       )
#   )
#
#   vc_drop <- sandwich::vcovHC(fit_drop, type = "HC1")
#
#   out[["drop_pol_correct"]] <- lmtest::coeftest(fit_drop, vcov. = vc_drop) %>%
#     broom::tidy() %>%
#     filter(
#       term %in% c(
#         "treatment_vpnonly",
#         "treatment_nlonly",
#         "treatment_vpnnl"
#       )
#     ) %>%
#     mutate(
#       outcome = Y,
#       dv = dv_w3,
#       specification = "Drop pol correct",
#       mean_dv_all = mean(
#         cy_w3[[dv_w3]][cy_w3$treatment_user != 1 & cy_w3$trust_central_govt_w3 != 10],
#         na.rm = TRUE
#       ),
#       sd_dv_all = sd(
#         cy_w3[[dv_w3]][cy_w3$treatment_user != 1 & cy_w3$trust_central_govt_w3 != 10],
#         na.rm = TRUE
#       ),
#       mean_dv_control = mean(
#         cy_w3[[dv_w3]][cy_w3$treatment_control == 1 & cy_w3$trust_central_govt_w3 != 10],
#         na.rm = TRUE
#       ),
#       sd_dv_control = sd(
#         cy_w3[[dv_w3]][cy_w3$treatment_control == 1 & cy_w3$trust_central_govt_w3 != 10],
#         na.rm = TRUE
#       ),
#       mean_dv_user = mean(
#         cy_w3[[dv_w3]][cy_w3$treatment_user == 1 & cy_w3$trust_central_govt_w3 != 10],
#         na.rm = TRUE
#       ),
#       sd_dv_user = sd(
#         cy_w3[[dv_w3]][cy_w3$treatment_user == 1 & cy_w3$trust_central_govt_w3 != 10],
#         na.rm = TRUE
#       ),
#       p_value_treatment_vpnnl = p.value[term == "treatment_vpnnl"][1]
#     )
#
#   bind_rows(out)
# })
#
# write.csv(
#   table_a11,
#   file.path(outreg_dir, "outregs_maintreatmenteffects_robustness.csv"),
#   row.names = FALSE
# )
#
# # Table A.12 ------------------------------------------------------------------
#
# set.seed(12345)
#
# table_a12_panel_vars <- c(
#   "info_foreign_website",
#   "info_freq_website_for",
#   "az_belief_media_value",
#   "az_belief_media_trust",
#   "bias_domestic",
#   "bias_foreign_r",
#   "az_belief_media_justif",
#   "bias_dom_govt_policy_t1",
#   "bias_for_govt_policy_t1",
#   "news_perccor_cen",
#   "news_perccor_unc",
#   "protest_pcheard_china",
#   "protest_pcheard_foreign",
#   "protest_2011_tmrw_parade",
#   "az_knowledge_meta",
#   "az_belief_econ_perf_cn_r",
#   "az_belief_econ_conf_cn",
#   "az_belief_instchange",
#   "az_belief_trust_govt_r",
#   "az_belief_trust_foreign",
#   "az_belief_evalgovt_r",
#   "importance_live_in_demo_r",
#   "az_belief_willing",
#   "frequency_talk_politic",
#   "frequency_persuade_friends",
#   "participate_social_protest",
#   "participate_plan_vote",
#   "participate_complain_school",
#   "stock_participation",
#   "plan_grad_foreignmaster",
#   "cp_t3_for_firm",
#   "cloc_for",
#   "az_overall"
# )
#
# table_a12_nonpanel_vars <- c(
#   "vpn_purchase_wmt_record",
#   "vpn_purchase_yes",
#   "az_belief_econ_perf_us_w3",
#   "az_belief_econ_conf_us_w3"
# )
#
# # Version A: ranked among all Wave 3 participants -----------------------------
#
# table_a12_all_dat <- cy %>%
#   filter(panelmerged_wave3 == 1)
#
# outregs_quantilemovement_wave3_all_panel <-
#   purrr::map_dfr(table_a12_panel_vars, function(Y) {
#
#     y_w1 <- paste0(Y, "_w1")
#     y_w3 <- paste0(Y, "_w3")
#
#     tmp <- table_a12_all_dat %>%
#       mutate(
#         y_w1_n = .data[[y_w1]] + runif(n()) / 2,
#         y_w3_n = .data[[y_w3]] + runif(n()) / 2,
#         y_p1 = ntile(y_w1_n, 100),
#         y_p3 = ntile(y_w3_n, 100)
#       )
#
#     p1 <- tmp %>%
#       filter(treatment_vpnnl == 1) %>%
#       summarise(p1 = median(y_p1, na.rm = TRUE)) %>%
#       pull(p1)
#
#     p3 <- tmp %>%
#       filter(treatment_vpnnl == 1) %>%
#       summarise(p3 = median(y_p3, na.rm = TRUE)) %>%
#       pull(p3)
#
#     tibble(
#       variable = Y,
#       p1 = p1,
#       p3 = p3,
#       change = p3 - p1
#     )
#   })
#
# outregs_quantilemovement_wave3_all_nonpanel <-
#   purrr::map_dfr(table_a12_nonpanel_vars, function(Y) {
#
#     tmp <- table_a12_all_dat %>%
#       mutate(
#         y_n = .data[[Y]] + runif(n()) / 2,
#         y_p = ntile(y_n, 100)
#       )
#
#     p1 <- tmp %>%
#       filter(
#         treatment_control == 1 |
#           treatment_vpnonly == 1 |
#           treatment_nlonly == 1
#       ) %>%
#       summarise(p1 = median(y_p, na.rm = TRUE)) %>%
#       pull(p1)
#
#     p3 <- tmp %>%
#       filter(treatment_vpnnl == 1) %>%
#       summarise(p3 = median(y_p, na.rm = TRUE)) %>%
#       pull(p3)
#
#     tibble(
#       variable = Y,
#       p1 = p1,
#       p3 = p3,
#       change = p3 - p1
#     )
#   })
#
# outregs_quantilemovement_wave3_all <- bind_rows(
#   outregs_quantilemovement_wave3_all_panel,
#   outregs_quantilemovement_wave3_all_nonpanel
# )
#
# write.csv(
#   outregs_quantilemovement_wave3_all,
#   file.path(outreg_dir, "outregs_quantilemovement_wave3_all.csv"),
#   row.names = FALSE
# )
#
#
# # Version B: ranked among non-existing users ----------------------------------
#
# table_a12_nonexisting_dat <- cy %>%
#   filter(
#     panelmerged_wave3 == 1,
#     treatment_user != 1
#   )
#
# outregs_quantilemovement_wave3_nonexistinguser_panel <-
#   purrr::map_dfr(table_a12_panel_vars, function(Y) {
#
#     y_w1 <- paste0(Y, "_w1")
#     y_w3 <- paste0(Y, "_w3")
#
#     tmp <- table_a12_nonexisting_dat %>%
#       mutate(
#         y_w1_n = .data[[y_w1]] + runif(n()) / 500,
#         y_w3_n = .data[[y_w3]] + runif(n()) / 500,
#         y_p1 = ntile(y_w1_n, 100),
#         y_p3 = ntile(y_w3_n, 100)
#       )
#
#     p1 <- tmp %>%
#       filter(treatment_vpnnl == 1) %>%
#       summarise(p1 = median(y_p1, na.rm = TRUE)) %>%
#       pull(p1)
#
#     p3 <- tmp %>%
#       filter(treatment_vpnnl == 1) %>%
#       summarise(p3 = median(y_p3, na.rm = TRUE)) %>%
#       pull(p3)
#
#     tibble(
#       variable = Y,
#       p1 = p1,
#       p3 = p3,
#       change = p3 - p1
#     )
#   })
#
# outregs_quantilemovement_wave3_nonexistinguser_nonpanel <-
#   purrr::map_dfr(table_a12_nonpanel_vars, function(Y) {
#
#     tmp <- table_a12_nonexisting_dat %>%
#       mutate(
#         y_n = .data[[Y]] + runif(n()) / 500,
#         y_p = ntile(y_n, 100)
#       )
#
#     p1 <- tmp %>%
#       filter(
#         treatment_control == 1 |
#           treatment_vpnonly == 1 |
#           treatment_nlonly == 1
#       ) %>%
#       summarise(p1 = median(y_p, na.rm = TRUE)) %>%
#       pull(p1)
#
#     p3 <- tmp %>%
#       filter(treatment_vpnnl == 1) %>%
#       summarise(p3 = median(y_p, na.rm = TRUE)) %>%
#       pull(p3)
#
#     tibble(
#       variable = Y,
#       p1 = p1,
#       p3 = p3,
#       change = p3 - p1
#     )
#   })
#
# outregs_quantilemovement_wave3_nonexistinguser <- bind_rows(
#   outregs_quantilemovement_wave3_nonexistinguser_panel,
#   outregs_quantilemovement_wave3_nonexistinguser_nonpanel
# )
#
# write.csv(
#   outregs_quantilemovement_wave3_nonexistinguser,
#   file.path(outreg_dir, "outregs_quantilemovement_wave3_nonexistinguser.csv"),
#   row.names = FALSE
# )
#
# # Table A.13 ------------------------------------------------------------------
#
# table_a13_dat <- cy %>%
#   filter(
#     panelmerged_wave3 == 1,
#     treatment_user != 1
#   )
#
# # Transform into dummy indicators: above median -------------------------------
#
# above_median_w3 <- c(
#   "info_foreign_website_w3",
#   "info_freq_website_for_w3",
#   "az_belief_media_value_w3",
#   "az_belief_media_trust_w3",
#   "bias_domestic_w3",
#   "az_belief_media_justif_w3",
#   "news_perccor_cen_w3",
#   "news_perccor_unc_w3",
#   "protest_pcheard_china_w3",
#   "protest_pcheard_foreign_w3",
#   "az_knowledge_meta_w3",
#   "az_belief_econ_conf_cn_w3",
#   "az_belief_econ_perf_us_w3",
#   "az_belief_econ_conf_us_w3",
#   "az_belief_instchange_w3",
#   "az_belief_trust_foreign_w3",
#   "az_belief_willing_w3",
#   "frequency_talk_politic_w3",
#   "frequency_persuade_friends_w3"
# )
#
# above_median_w1 <- c(
#   "info_foreign_website_w1",
#   "info_freq_website_for_w1",
#   "az_belief_media_value_w1",
#   "az_belief_media_trust_w1",
#   "bias_domestic_w1",
#   "az_belief_media_justif_w1",
#   "news_perccor_cen_w1",
#   "news_perccor_unc_w1",
#   "protest_pcheard_china_w1",
#   "protest_pcheard_foreign_w1",
#   "az_knowledge_meta_w1",
#   "az_belief_econ_conf_cn_w1",
#   "az_belief_instchange_w1",
#   "az_belief_trust_foreign_w1",
#   "az_belief_willing_w1",
#   "frequency_talk_politic_w1",
#   "frequency_persuade_friends_w1"
# )
#
# for (Y in c(above_median_w3, above_median_w1)) {
#
#   median_y <- median(table_a13_dat[[Y]], na.rm = TRUE)
#
#   table_a13_dat <- table_a13_dat %>%
#     mutate(
#       !!paste0(Y, "_p") := as.integer(.data[[Y]] > median_y)
#     )
# }
#
# # Transform into dummy indicators: below median -------------------------------
#
# below_median_w3 <- c(
#   "bias_foreign_w3",
#   "importance_live_in_demo_w3",
#   "az_belief_econ_perf_cn_w3",
#   "az_belief_trust_govt_w3",
#   "az_belief_evalgovt_w3"
# )
#
# below_median_w1 <- c(
#   "bias_foreign_w1",
#   "importance_live_in_demo_w1",
#   "az_belief_econ_perf_cn_w1",
#   "az_belief_trust_govt_w1",
#   "az_belief_evalgovt_w1"
# )
#
# for (Y in c(below_median_w3, below_median_w1)) {
#
#   median_y <- median(table_a13_dat[[Y]], na.rm = TRUE)
#
#   table_a13_dat <- table_a13_dat %>%
#     mutate(
#       !!paste0(Y, "_p") := as.integer(.data[[Y]] < median_y)
#     )
# }
#
# # Already dummy, keep as is ----------------------------------------------------
#
# table_a13_dat <- table_a13_dat %>%
#   rename(
#     par_complain_school_w1 = participate_complain_school_w1,
#     par_complain_school_w3 = participate_complain_school_w3
#   )
#
# dummy_keep_w3 <- c(
#   "vpn_purchase_wmt_record",
#   "vpn_purchase_yes",
#   "bias_dom_govt_policy_t1_w3",
#   "protest_2011_tmrw_parade_w3",
#   "participate_social_protest_w3",
#   "participate_plan_vote_w3",
#   "par_complain_school_w3",
#   "plan_grad_foreignmaster_w3",
#   "cp_t3_for_firm_w3",
#   "cloc_for_w3"
# )
#
# dummy_keep_w1 <- c(
#   "bias_dom_govt_policy_t1_w1",
#   "protest_2011_tmrw_parade_w1",
#   "participate_social_protest_w1",
#   "participate_plan_vote_w1",
#   "par_complain_school_w1",
#   "plan_grad_foreignmaster_w1",
#   "cp_t3_for_firm_w1",
#   "cloc_for_w1"
# )
#
# for (Y in c(dummy_keep_w3, dummy_keep_w1)) {
#
#   table_a13_dat <- table_a13_dat %>%
#     mutate(
#       !!paste0(Y, "_p") := .data[[Y]]
#     )
# }
#
# # Already dummy, flip sign -----------------------------------------------------
#
# dummy_flip_w3 <- c(
#   "bias_for_govt_policy_t1_w3",
#   "stock_participation_w3"
# )
#
# dummy_flip_w1 <- c(
#   "bias_for_govt_policy_t1_w1",
#   "stock_participation_w1"
# )
#
# for (Y in c(dummy_flip_w3, dummy_flip_w1)) {
#
#   table_a13_dat <- table_a13_dat %>%
#     mutate(
#       !!paste0(Y, "_p") := 1 - .data[[Y]]
#     )
# }
#
# # Persuasion rates for panel variables ----------------------------------------
#
# table_a13_panel_vars <- c(
#   "info_foreign_website",
#   "info_freq_website_for",
#   "az_belief_media_value",
#   "az_belief_media_trust",
#   "bias_domestic",
#   "az_belief_media_justif",
#   "news_perccor_cen",
#   "news_perccor_unc",
#   "protest_pcheard_china",
#   "protest_pcheard_foreign",
#   "az_knowledge_meta",
#   "az_belief_econ_conf_cn",
#   "az_belief_instchange",
#   "az_belief_trust_foreign",
#   "importance_live_in_demo",
#   "az_belief_willing",
#   "frequency_talk_politic",
#   "frequency_persuade_friends",
#   "bias_foreign",
#   "az_belief_econ_perf_cn",
#   "az_belief_trust_govt",
#   "az_belief_evalgovt",
#   "bias_dom_govt_policy_t1",
#   "protest_2011_tmrw_parade",
#   "participate_social_protest",
#   "participate_plan_vote",
#   "par_complain_school",
#   "plan_grad_foreignmaster",
#   "cp_t3_for_firm",
#   "cloc_for",
#   "bias_for_govt_policy_t1",
#   "stock_participation"
# )
#
# outregs_persuasionrates_wave3_panel <- purrr::map_dfr(
#   table_a13_panel_vars,
#   function(Y) {
#
#     y_w3_p <- paste0(Y, "_w3_p")
#     y_w1_p <- paste0(Y, "_w1_p")
#
#     fit <- lm(
#       reformulate(
#         c(
#           "treatment_vpnonly",
#           "treatment_nlonly",
#           "treatment_vpnnl"
#         ),
#         response = y_w3_p
#       ),
#       data = table_a13_dat
#     )
#
#     vc <- sandwich::vcovHC(fit, type = "HC1")
#
#     reg_out <- lmtest::coeftest(fit, vcov. = vc) %>%
#       broom::tidy()
#
#     b_vpnnl <- reg_out %>%
#       filter(term == "treatment_vpnnl") %>%
#       pull(estimate)
#
#     se_vpnnl <- reg_out %>%
#       filter(term == "treatment_vpnnl") %>%
#       pull(std.error)
#
#     df_r <- df.residual(fit)
#
#     p_treatvpnXnl <- 2 * pt(
#       abs(b_vpnnl / se_vpnnl),
#       df = df_r,
#       lower.tail = FALSE
#     )
#
#     baseline_mean <- table_a13_dat %>%
#       filter(treatment_vpnnl == 1) %>%
#       summarise(x = mean(.data[[y_w1_p]], na.rm = TRUE)) %>%
#       pull(x)
#
#     persuasion <- b_vpnnl / (0.646 * (1 - baseline_mean))
#
#     reg_out %>%
#       filter(term == "treatment_vpnnl") %>%
#       mutate(
#         variable = Y,
#         outcome = y_w3_p,
#         specification = "Raw Diff",
#         mean_dv_all = mean(table_a13_dat[[y_w3_p]], na.rm = TRUE),
#         sd_dv_all = sd(table_a13_dat[[y_w3_p]], na.rm = TRUE),
#         mean_dv_control = mean(
#           table_a13_dat[[y_w3_p]][table_a13_dat$treatment_control == 1],
#           na.rm = TRUE
#         ),
#         p_value_treatment_vpnnl = p_treatvpnXnl,
#         persuasion = persuasion
#       )
#   }
# )
#
# # Persuasion rates for non-panel variables ------------------------------------
#
# table_a13_nonpanel_vars <- c(
#   "vpn_purchase_wmt_record",
#   "vpn_purchase_yes",
#   "az_belief_econ_perf_us_w3",
#   "az_belief_econ_conf_us_w3"
# )
#
# outregs_persuasionrates_wave3_nonpanel <- purrr::map_dfr(
#   table_a13_nonpanel_vars,
#   function(Y) {
#
#     y_p <- paste0(Y, "_p")
#
#     fit <- lm(
#       reformulate(
#         c(
#           "treatment_vpnonly",
#           "treatment_nlonly",
#           "treatment_vpnnl"
#         ),
#         response = y_p
#       ),
#       data = table_a13_dat
#     )
#
#     vc <- sandwich::vcovHC(fit, type = "HC1")
#
#     reg_out <- lmtest::coeftest(fit, vcov. = vc) %>%
#       broom::tidy()
#
#     b_vpnnl <- reg_out %>%
#       filter(term == "treatment_vpnnl") %>%
#       pull(estimate)
#
#     se_vpnnl <- reg_out %>%
#       filter(term == "treatment_vpnnl") %>%
#       pull(std.error)
#
#     b_cons <- reg_out %>%
#       filter(term == "(Intercept)") %>%
#       pull(estimate)
#
#     df_r <- df.residual(fit)
#
#     p_treatvpnXnl <- 2 * pt(
#       abs(b_vpnnl / se_vpnnl),
#       df = df_r,
#       lower.tail = FALSE
#     )
#
#     persuasion <- b_vpnnl / (0.646 * (1 - b_cons))
#
#     reg_out %>%
#       filter(term == "treatment_vpnnl") %>%
#       mutate(
#         variable = Y,
#         outcome = y_p,
#         specification = "Raw Diff",
#         mean_dv_all = mean(table_a13_dat[[y_p]], na.rm = TRUE),
#         sd_dv_all = sd(table_a13_dat[[y_p]], na.rm = TRUE),
#         mean_dv_control = mean(
#           table_a13_dat[[y_p]][table_a13_dat$treatment_control == 1],
#           na.rm = TRUE
#         ),
#         p_value_treatment_vpnnl = p_treatvpnXnl,
#         persuasion = persuasion
#       )
#   }
# )
#
# outregs_persuasionrates_wave3 <- bind_rows(
#   outregs_persuasionrates_wave3_panel,
#   outregs_persuasionrates_wave3_nonpanel
# )
#
# write.csv(
#   outregs_persuasionrates_wave3,
#   file.path(outreg_dir, "outregs_persuasionrates_wave3.csv"),
#   row.names = FALSE
# )
#
# # Table A.14 ------------------------------------------------------------------
#
# table_a14_dat <- cy %>%
#   filter(
#     panelmerged_wave3 == 1,
#     treatment_user != 1
#   )
#
# table_a14_hvars <- c(
#   "gender", "birth_year", "residence_coastal", "hukou_urban",
#   "university_elite", "hs_track_science", "department_ssh",
#   "domestic_english_atleast4", "foreign_english_yes",
#   "travel_hktaiwan", "travel_foreign_yes", "father_edu_hsabove",
#   "work_father_govt", "father_ccp", "mother_edu_hsabove",
#   "work_mother_govt", "mother_ccp", "hh_income",
#   "az_preference_risk", "az_preference_time",
#   "az_preference_altruism", "az_preference_reciprocity",
#   "az_overall_a_w1", "az_overall_b_w1", "az_overall_c_w1",
#   "az_overall_d_w1", "az_overall_e_w1"
# )
#
# for (h in table_a14_hvars) {
#   table_a14_dat <- table_a14_dat %>%
#     mutate(
#       !!paste0("aeX", h) := treatment_vpnnl * .data[[paste0("h_", h)]]
#     )
# }
#
# table_a14_rhs <- c(
#   "treatment_vpnnl",
#   paste0("h_", table_a14_hvars),
#   paste0("aeX", table_a14_hvars)
# )
#
# table_a14 <- purrr::map_dfr(
#   c(
#     "az_overall_a_w3",
#     "az_overall_b_w3",
#     "az_overall_c_w3",
#     "az_overall_d_w3",
#     "az_overall_e_w3"
#   ),
#   function(Y) {
#
#     fit <- lm(
#       reformulate(
#         table_a14_rhs,
#         response = Y
#       ),
#       data = table_a14_dat
#     )
#
#     broom::tidy(fit) %>%
#       mutate(outcome = Y)
#   }
# )
#
# write.csv(
#   table_a14,
#   file.path(outreg_dir, "outregs_heterogeneity_all.csv"),
#   row.names = FALSE
# )
#
#
# # Table A.15 ------------------------------------------------------------------
#
# table_a15_w2 <- purrr::map_dfr(var_knowledge_news_reg_cen_w2, function(Y) {
#
#   fit <- lm(
#     reformulate(
#       c(
#         "soclearning_ownaccess",
#         "soclearning_rm_new_w2",
#         "soclearning_ownXnew_w2"
#       ),
#       response = Y
#     ),
#     data = cy %>%
#       filter(
#         vpn_roommate_existing == 0,
#         soclearning_rm_new_w2 < 2
#       )
#   )
#
#   reg_out <- broom::tidy(fit)
#
#   mu_100 <- mean(
#     cy[[Y]][
#       cy$treatment_master >= 4 &
#         cy$vpn_roommate_existing == 0 &
#         cy$vpn_roommate_new_w2 == 0
#     ],
#     na.rm = TRUE
#   )
#
#   mu_000 <- mean(
#     cy[[Y]][
#       cy$treatment_master < 4 &
#         cy$vpn_roommate_existing == 0 &
#         cy$vpn_roommate_new_w2 == 0
#     ],
#     na.rm = TRUE
#   )
#
#   p_all <- mu_100 - mu_000
#
#   mu_101 <- mean(
#     cy[[Y]][
#       cy$treatment_master >= 4 &
#         cy$vpn_roommate_existing == 0 &
#         cy$vpn_roommate_new_w2 == 1
#     ],
#     na.rm = TRUE
#   )
#
#   mu_001 <- mean(
#     cy[[Y]][
#       cy$treatment_master < 4 &
#         cy$vpn_roommate_existing == 0 &
#         cy$vpn_roommate_new_w2 == 1
#     ],
#     na.rm = TRUE
#   )
#
#   q_new_1_0 <- (mu_001 - mu_000) / mu_100
#   q_new_1_1 <- (mu_101 - mu_100) / mu_100
#
#   pred_mu_102 <- mu_100 + (1 - (1 - mu_100 * q_new_1_1)^2)
#   pred_mu_002 <- mu_000 + (1 - (1 - mu_100 * q_new_1_0)^2)
#
#   mu_002 <- mean(
#     cy[[Y]][
#       cy$treatment_master < 4 &
#         cy$vpn_roommate_existing == 0 &
#         cy$vpn_roommate_new_w2 == 2
#     ],
#     na.rm = TRUE
#   )
#
#   mu_102 <- mean(
#     cy[[Y]][
#       cy$treatment_master >= 4 &
#         cy$vpn_roommate_existing == 0 &
#         cy$vpn_roommate_new_w2 == 2
#     ],
#     na.rm = TRUE
#   )
#
#   reg_out %>%
#     filter(
#       term %in% c(
#         "soclearning_ownaccess",
#         "soclearning_rm_new_w2",
#         "soclearning_ownXnew_w2"
#       )
#     ) %>%
#     mutate(
#       outcome = Y,
#       specification = "Raw Diff",
#       p_all = p_all,
#       q_new_1_0 = q_new_1_0,
#       q_new_1_1 = q_new_1_1,
#       pred_mu_002 = pred_mu_002,
#       mu_002 = mu_002,
#       pred_mu_102 = pred_mu_102,
#       mu_102 = mu_102
#     )
# })
#
# table_a15_w3 <- purrr::map_dfr(
#   c(var_knowledge_news_reg_cen_w3, "news_perccor_cen_all"),
#   function(Y) {
#
#     fit <- lm(
#       reformulate(
#         c(
#           "soclearning_ownaccess",
#           "soclearning_rm_new_w3",
#           "soclearning_ownXnew_w3"
#         ),
#         response = Y
#       ),
#       data = cy %>%
#         filter(
#           vpn_roommate_existing == 0,
#           soclearning_rm_new_w3 < 2
#         )
#     )
#
#     reg_out <- broom::tidy(fit)
#
#     mu_100 <- mean(
#       cy[[Y]][
#         cy$treatment_master >= 4 &
#           cy$vpn_roommate_existing == 0 &
#           cy$vpn_roommate_new_w3 == 0
#       ],
#       na.rm = TRUE
#     )
#
#     mu_000 <- mean(
#       cy[[Y]][
#         cy$treatment_master < 4 &
#           cy$vpn_roommate_existing == 0 &
#           cy$vpn_roommate_new_w3 == 0
#       ],
#       na.rm = TRUE
#     )
#
#     p_all <- mu_100 - mu_000
#
#     mu_101 <- mean(
#       cy[[Y]][
#         cy$treatment_master >= 4 &
#           cy$vpn_roommate_existing == 0 &
#           cy$vpn_roommate_new_w3 == 1
#       ],
#       na.rm = TRUE
#     )
#
#     mu_001 <- mean(
#       cy[[Y]][
#         cy$treatment_master < 4 &
#           cy$vpn_roommate_existing == 0 &
#           cy$vpn_roommate_new_w3 == 1
#       ],
#       na.rm = TRUE
#     )
#
#     q_new_1_0 <- (mu_001 - mu_000) / mu_100
#     q_new_1_1 <- (mu_101 - mu_100) / mu_100
#
#     pred_mu_102 <- mu_100 + (1 - (1 - mu_100 * q_new_1_1)^2)
#     pred_mu_002 <- mu_000 + (1 - (1 - mu_100 * q_new_1_0)^2)
#
#     mu_002 <- mean(
#       cy[[Y]][
#         cy$treatment_master < 4 &
#           cy$vpn_roommate_existing == 0 &
#           cy$vpn_roommate_new_w3 == 2
#       ],
#       na.rm = TRUE
#     )
#
#     mu_102 <- mean(
#       cy[[Y]][
#         cy$treatment_master >= 4 &
#           cy$vpn_roommate_existing == 0 &
#           cy$vpn_roommate_new_w3 == 2
#       ],
#       na.rm = TRUE
#     )
#
#     reg_out %>%
#       filter(
#         term %in% c(
#           "soclearning_ownaccess",
#           "soclearning_rm_new_w3",
#           "soclearning_ownXnew_w3"
#         )
#       ) %>%
#       mutate(
#         outcome = Y,
#         specification = "Raw Diff",
#         p_all = p_all,
#         q_new_1_0 = q_new_1_0,
#         q_new_1_1 = q_new_1_1,
#         pred_mu_002 = pred_mu_002,
#         mu_002 = mu_002,
#         pred_mu_102 = pred_mu_102,
#         mu_102 = mu_102
#       )
#   }
# )
#
# outregs_wave3sociallearning_full <- bind_rows(
#   table_a15_w2,
#   table_a15_w3
# )
#
# write.csv(
#   outregs_wave3sociallearning_full,
#   file.path(outreg_dir, "outregs_wave3sociallearning_full.csv"),
#   row.names = FALSE
# )
#
#
# # Table A.15: standard errors of nonlinear estimates --------------------------
#
# table_a15_nlcom_w2_dat <- cy %>%
#   mutate(
#     sc_reg_100 = as.integer(
#       treatment_master >= 4 &
#         vpn_roommate_existing == 0 &
#         vpn_roommate_new_w2 == 0
#     ),
#     sc_reg_000 = as.integer(
#       treatment_master < 4 &
#         vpn_roommate_existing == 0 &
#         vpn_roommate_new_w2 == 0
#     ),
#     sc_reg_101 = as.integer(
#       treatment_master >= 4 &
#         vpn_roommate_existing == 0 &
#         vpn_roommate_new_w2 == 1
#     ),
#     sc_reg_001 = as.integer(
#       treatment_master < 4 &
#         vpn_roommate_existing == 0 &
#         vpn_roommate_new_w2 == 1
#     )
#   )
#
# table_a15_nlcom_w2 <- purrr::map_dfr(var_knowledge_news_reg_cen_w2, function(Y) {
#
#   fit <- lm(
#     reformulate(
#       c("sc_reg_100", "sc_reg_000", "sc_reg_101", "sc_reg_001"),
#       response = Y,
#       intercept = FALSE
#     ),
#     data = table_a15_nlcom_w2_dat %>%
#       filter(
#         vpn_roommate_existing == 0,
#         vpn_roommate_new_w2 < 2
#       )
#   )
#
#   b <- coef(fit)
#   vc <- vcov(fit)
#
#   nl_fun <- function(beta) {
#
#     names(beta) <- names(b)
#
#     q_new_1_0 <-
#       (beta["sc_reg_001"] - beta["sc_reg_000"]) /
#       beta["sc_reg_100"]
#
#     q_new_1_1 <-
#       (beta["sc_reg_101"] - beta["sc_reg_100"]) /
#       beta["sc_reg_100"]
#
#     pred_mu_002 <-
#       beta["sc_reg_000"] +
#       (
#         1 -
#           (
#             1 -
#               beta["sc_reg_100"] *
#               (
#                 (beta["sc_reg_001"] - beta["sc_reg_000"]) /
#                   beta["sc_reg_100"]
#               )
#           )^2
#       )
#
#     pred_mu_102 <-
#       beta["sc_reg_100"] +
#       (
#         1 -
#           (
#             1 -
#               beta["sc_reg_100"] *
#               (
#                 (beta["sc_reg_101"] - beta["sc_reg_100"]) /
#                   beta["sc_reg_100"]
#               )
#           )^2
#       )
#
#     c(
#       q_new_1_0 = q_new_1_0,
#       q_new_1_1 = q_new_1_1,
#       pred_mu_002 = pred_mu_002,
#       pred_mu_102 = pred_mu_102
#     )
#   }
#
#   estimates <- nl_fun(b)
#
#   grad <- numDeriv::jacobian(nl_fun, b)
#
#   ses <- sqrt(
#     diag(
#       grad %*% vc %*% t(grad)
#     )
#   )
#
#   tibble(
#     outcome = Y,
#     term = names(estimates),
#     estimate = as.numeric(estimates),
#     std_error = as.numeric(ses)
#   )
# })
#
# table_a15_nlcom_w3_dat <- cy %>%
#   mutate(
#     sc_reg_100 = as.integer(
#       treatment_master >= 4 &
#         vpn_roommate_existing == 0 &
#         vpn_roommate_new_w3 == 0
#     ),
#     sc_reg_000 = as.integer(
#       treatment_master < 4 &
#         vpn_roommate_existing == 0 &
#         vpn_roommate_new_w3 == 0
#     ),
#     sc_reg_101 = as.integer(
#       treatment_master >= 4 &
#         vpn_roommate_existing == 0 &
#         vpn_roommate_new_w3 == 1
#     ),
#     sc_reg_001 = as.integer(
#       treatment_master < 4 &
#         vpn_roommate_existing == 0 &
#         vpn_roommate_new_w3 == 1
#     )
#   )
#
# table_a15_nlcom_w3 <- purrr::map_dfr(
#   c(var_knowledge_news_reg_cen_w3, "news_perccor_cen_all"),
#   function(Y) {
#
#     fit <- lm(
#       reformulate(
#         c("sc_reg_100", "sc_reg_000", "sc_reg_101", "sc_reg_001"),
#         response = Y,
#         intercept = FALSE
#       ),
#       data = table_a15_nlcom_w3_dat %>%
#         filter(
#           vpn_roommate_existing == 0,
#           vpn_roommate_new_w3 < 2
#         )
#     )
#
#     b <- coef(fit)
#     vc <- vcov(fit)
#
#     nl_fun <- function(beta) {
#
#       names(beta) <- names(b)
#
#       q_new_1_0 <-
#         (beta["sc_reg_001"] - beta["sc_reg_000"]) /
#         beta["sc_reg_100"]
#
#       q_new_1_1 <-
#         (beta["sc_reg_101"] - beta["sc_reg_100"]) /
#         beta["sc_reg_100"]
#
#       pred_mu_002 <-
#         beta["sc_reg_000"] +
#         (
#           1 -
#             (
#               1 -
#                 beta["sc_reg_100"] *
#                 (
#                   (beta["sc_reg_001"] - beta["sc_reg_000"]) /
#                     beta["sc_reg_100"]
#                 )
#             )^2
#         )
#
#       pred_mu_102 <-
#         beta["sc_reg_100"] +
#         (
#           1 -
#             (
#               1 -
#                 beta["sc_reg_100"] *
#                 (
#                   (beta["sc_reg_101"] - beta["sc_reg_100"]) /
#                     beta["sc_reg_100"]
#                 )
#             )^2
#         )
#
#       c(
#         q_new_1_0 = q_new_1_0,
#         q_new_1_1 = q_new_1_1,
#         pred_mu_002 = pred_mu_002,
#         pred_mu_102 = pred_mu_102
#       )
#     }
#
#     estimates <- nl_fun(b)
#
#     grad <- numDeriv::jacobian(nl_fun, b)
#
#     ses <- sqrt(
#       diag(
#         grad %*% vc %*% t(grad)
#       )
#     )
#
#     tibble(
#       outcome = Y,
#       term = names(estimates),
#       estimate = as.numeric(estimates),
#       std_error = as.numeric(ses)
#     )
#   }
# )
#
# table_a15_nlcom <- bind_rows(
#   table_a15_nlcom_w2,
#   table_a15_nlcom_w3
# )
#
# write.csv(
#   table_a15_nlcom,
#   file.path(outreg_dir, "outregs_wave3sociallearning_full_nlcom.csv"),
#   row.names = FALSE
# )
#
#
# # Table A.15: bootstrap block from Stata, kept commented out ------------------
# # The original Stata code comments out this section. The R translation below is
# # also commented out to preserve the original workflow.
#
# # set.seed(12345)
# #
# # bootstrap_source <- cy %>%
# #   select(
# #     treatment_master,
# #     vpn_roommate_existing,
# #     vpn_roommate_new_w2,
# #     vpn_roommate_new_w3,
# #     all_of(var_knowledge_news_reg_cen_w2),
# #     all_of(var_knowledge_news_reg_cen_w3),
# #     news_perccor_cen_all
# #   ) %>%
# #   filter(!is.na(news_perccor_cen_all))
# #
# # sociallearning_pred_bootstrap_full <- purrr::map_dfr(1:1000, function(i) {
# #
# #   boot_dat <- bootstrap_source %>%
# #     slice_sample(
# #       n = n(),
# #       replace = TRUE
# #     )
# #
# #   boot_w2 <- purrr::map_dfr(var_knowledge_news_reg_cen_w2, function(Y) {
# #
# #     mu_100 <- mean(
# #       boot_dat[[Y]][
# #         boot_dat$treatment_master >= 4 &
# #           boot_dat$vpn_roommate_existing == 0 &
# #           boot_dat$vpn_roommate_new_w2 == 0
# #       ],
# #       na.rm = TRUE
# #     )
# #
# #     mu_000 <- mean(
# #       boot_dat[[Y]][
# #         boot_dat$treatment_master < 4 &
# #           boot_dat$vpn_roommate_existing == 0 &
# #           boot_dat$vpn_roommate_new_w2 == 0
# #       ],
# #       na.rm = TRUE
# #     )
# #
# #     mu_101 <- mean(
# #       boot_dat[[Y]][
# #         boot_dat$treatment_master >= 4 &
# #           boot_dat$vpn_roommate_existing == 0 &
# #           boot_dat$vpn_roommate_new_w2 == 1
# #       ],
# #       na.rm = TRUE
# #     )
# #
# #     mu_001 <- mean(
# #       boot_dat[[Y]][
# #         boot_dat$treatment_master < 4 &
# #           boot_dat$vpn_roommate_existing == 0 &
# #           boot_dat$vpn_roommate_new_w2 == 1
# #       ],
# #       na.rm = TRUE
# #     )
# #
# #     q_new_1_0 <- (mu_001 - mu_000) / mu_100
# #     q_new_1_1 <- (mu_101 - mu_100) / mu_100
# #
# #     tibble(
# #       i = i,
# #       y = Y,
# #       pred_mu_102 = mu_100 + (1 - (1 - mu_100 * q_new_1_1)^2),
# #       pred_mu_002 = mu_000 + (1 - (1 - mu_100 * q_new_1_0)^2)
# #     )
# #   })
# #
# #   boot_w3 <- purrr::map_dfr(
# #     c(var_knowledge_news_reg_cen_w3, "news_perccor_cen_all"),
# #     function(Y) {
# #
# #       mu_100 <- mean(
# #         boot_dat[[Y]][
# #           boot_dat$treatment_master >= 4 &
# #             boot_dat$vpn_roommate_existing == 0 &
# #             boot_dat$vpn_roommate_new_w3 == 0
# #         ],
# #         na.rm = TRUE
# #       )
# #
# #       mu_000 <- mean(
# #         boot_dat[[Y]][
# #           boot_dat$treatment_master < 4 &
# #             boot_dat$vpn_roommate_existing == 0 &
# #             boot_dat$vpn_roommate_new_w3 == 0
# #         ],
# #         na.rm = TRUE
# #       )
# #
# #       mu_101 <- mean(
# #         boot_dat[[Y]][
# #           boot_dat$treatment_master >= 4 &
# #             boot_dat$vpn_roommate_existing == 0 &
# #             boot_dat$vpn_roommate_new_w3 == 1
# #         ],
# #         na.rm = TRUE
# #       )
# #
# #       mu_001 <- mean(
# #         boot_dat[[Y]][
# #           boot_dat$treatment_master < 4 &
# #             boot_dat$vpn_roommate_existing == 0 &
# #             boot_dat$vpn_roommate_new_w3 == 1
# #         ],
# #         na.rm = TRUE
# #       )
# #
# #       q_new_1_0 <- (mu_001 - mu_000) / mu_100
# #       q_new_1_1 <- (mu_101 - mu_100) / mu_100
# #
# #       tibble(
# #         i = i,
# #         y = Y,
# #         pred_mu_102 = mu_100 + (1 - (1 - mu_100 * q_new_1_1)^2),
# #         pred_mu_002 = mu_000 + (1 - (1 - mu_100 * q_new_1_0)^2)
# #       )
# #     }
# #   )
# #
# #   bind_rows(boot_w2, boot_w3)
# # })
# #
# # sociallearning_pred_bootstrap_full <- sociallearning_pred_bootstrap_full %>%
# #   group_by(y) %>%
# #   mutate(
# #     pred_mu_102_bse = sd(pred_mu_102, na.rm = TRUE),
# #     pred_mu_002_bse = sd(pred_mu_002, na.rm = TRUE)
# #   ) %>%
# #   ungroup() %>%
# #   select(
# #     y,
# #     pred_mu_102_bse,
# #     pred_mu_002_bse
# #   ) %>%
# #   distinct()
# #
# # write.csv(
# #   sociallearning_pred_bootstrap_full,
# #   file.path(outreg_dir, "sociallearning_pred_bootstrap_full.csv"),
# #   row.names = FALSE
# # )
#
# # Table A.16 ------------------------------------------------------------------
#
# cy_w3 <- cy %>%
#   filter(panelmerged_wave3 == 1)
#
# table_a16 <- purrr::map_dfr(
#   c(
#     "az_overall_a_w3",
#     "az_overall_b_w3",
#     "az_overall_c_w3",
#     "az_overall_d_w3",
#     "az_overall_e_w3"
#   ),
#   function(Y) {
#
#     fit <- lm(
#       reformulate(
#         c(
#           "soclearning_ownaccess",
#           "soclearning_rm_new_w3",
#           "soclearning_ownXnew_w3"
#         ),
#         response = Y
#       ),
#       data = cy_w3 %>%
#         filter(
#           vpn_roommate_existing == 0,
#           soclearning_rm_new_w3 < 2
#         )
#     )
#
#     reg_out <- broom::tidy(fit)
#
#     mu_100 <- mean(
#       cy_w3[[Y]][
#         cy_w3$treatment_master >= 4 &
#           cy_w3$vpn_roommate_existing == 0 &
#           cy_w3$vpn_roommate_new_w3 == 0
#       ],
#       na.rm = TRUE
#     )
#
#     mu_000 <- mean(
#       cy_w3[[Y]][
#         cy_w3$treatment_master < 4 &
#           cy_w3$vpn_roommate_existing == 0 &
#           cy_w3$vpn_roommate_new_w3 == 0
#       ],
#       na.rm = TRUE
#     )
#
#     p_all <- mu_100 - mu_000
#
#     mu_101 <- mean(
#       cy_w3[[Y]][
#         cy_w3$treatment_master >= 4 &
#           cy_w3$vpn_roommate_existing == 0 &
#           cy_w3$vpn_roommate_new_w3 == 1
#       ],
#       na.rm = TRUE
#     )
#
#     mu_001 <- mean(
#       cy_w3[[Y]][
#         cy_w3$treatment_master < 4 &
#           cy_w3$vpn_roommate_existing == 0 &
#           cy_w3$vpn_roommate_new_w3 == 1
#       ],
#       na.rm = TRUE
#     )
#
#     q_new_1_0 <- (mu_001 - mu_000) / mu_100
#     q_new_1_1 <- (mu_101 - mu_100) / mu_100
#
#     pred_mu_102 <- mu_100 + (1 - (1 - mu_100 * q_new_1_1)^2)
#     pred_mu_002 <- mu_000 + (1 - (1 - mu_100 * q_new_1_0)^2)
#
#     mu_002 <- mean(
#       cy_w3[[Y]][
#         cy_w3$treatment_master < 4 &
#           cy_w3$vpn_roommate_existing == 0 &
#           cy_w3$vpn_roommate_new_w3 == 2
#       ],
#       na.rm = TRUE
#     )
#
#     mu_102 <- mean(
#       cy_w3[[Y]][
#         cy_w3$treatment_master >= 4 &
#           cy_w3$vpn_roommate_existing == 0 &
#           cy_w3$vpn_roommate_new_w3 == 2
#       ],
#       na.rm = TRUE
#     )
#
#     reg_out %>%
#       filter(
#         term %in% c(
#           "soclearning_ownaccess",
#           "soclearning_rm_new_w3",
#           "soclearning_ownXnew_w3"
#         )
#       ) %>%
#       mutate(
#         outcome = Y,
#         specification = "Raw Diff",
#         p_all = p_all,
#         q_new_1_0 = q_new_1_0,
#         q_new_1_1 = q_new_1_1,
#         pred_mu_002 = pred_mu_002,
#         mu_002 = mu_002,
#         pred_mu_102 = pred_mu_102,
#         mu_102 = mu_102
#       )
#   }
# )
#
# write.csv(
#   table_a16,
#   file.path(outreg_dir, "outregs_wave3sociallearning_additional.csv"),
#   row.names = FALSE
# )
#
# # Table A.17 ------------------------------------------------------------------
#
# cy_w2 <- cy %>%
#   filter(panelmerged_wave2 == 1)
#
# table_a17_vars <- c(
#   var_info_ranking_reg_w2,
#   var_info_freq_reg_w2,
#   var_media_valuation_reg_w2,
#   "az_belief_media_value_w2",
#   var_media_trust_reg_w2,
#   "az_belief_media_trust_w2",
#   var_censor_level_reg_w2,
#   var_censor_justif_reg_w2,
#   "az_belief_media_justif_w2",
#   var_censor_driver_reg_dom_w2,
#   var_censor_driver_reg_for_w2,
#   var_percmediabias_reg_ce_cn_w2,
#   "az_belief_media_cens_cn",
#   var_percmediabias_reg_ce_us_w2,
#   "az_belief_media_cens_us",
#   var_percmediabias_reg_di_cn_w2,
#   "az_belief_media_bias_cn",
#   var_percmediabias_reg_di_us_w2,
#   "az_belief_media_bias_us",
#   var_knowledge_news_reg_quiz,
#   "news_perccor_qui",
#   var_knowledge_news_reg_cen_w2,
#   "news_perccor_cen_w2",
#   var_knowledge_news_reg_unc_w2,
#   "news_perccor_unc_w2",
#   var_knowledge_prot_reg_chi_w2,
#   "protest_pcheard_china_w2",
#   var_knowledge_prot_reg_for_w2,
#   "protest_pcheard_foreign_w2",
#   var_knowledge_prot_reg_fak_w2,
#   var_knowledge_people_reg_tocens,
#   "people_perchd_unctocen",
#   var_knowledge_people_reg_censor,
#   "people_perchd_censored",
#   var_knowledge_people_reg_uncens,
#   "people_perchd_uncensor",
#   var_knowledge_people_reg_fake,
#   var_knowledge_meta_reg_w2,
#   "az_knowledge_meta_w2",
#   var_econ_guess_reg_cn_perf_w2,
#   "az_belief_econ_perf_cn_w2",
#   var_econ_guess_reg_cn_conf_w2,
#   "az_belief_econ_conf_cn_w2",
#   var_econ_guess_reg_us_perf_w2,
#   "az_belief_econ_perf_us_w2",
#   var_econ_guess_reg_us_conf_w2,
#   "az_belief_econ_conf_us_w2",
#   var_demand_change_reg_w2,
#   "az_belief_instchange_w2",
#   var_trust_inst_reg_govt_w2,
#   "az_belief_trust_govt_w2",
#   var_trust_inst_reg_copo_w2,
#   "az_belief_trust_copo_w2",
#   "trust_financial_domestic_w2",
#   var_trust_inst_reg_foreign_w2,
#   "az_belief_trust_foreign_w2",
#   "trust_financial_foreign_w2",
#   "trust_ngo_w2",
#   var_eval_govt_reg_w2,
#   "az_belief_evalgovt_w2",
#   var_eval_criteria_reg_w2,
#   var_severity_reg_w2,
#   "az_belief_severity_w2",
#   var_democracy_reg_w2,
#   "az_belief_democracy_w2",
#   var_contro_justi_reg_policy_w2,
#   "az_belief_justify_policy_w2",
#   var_contro_justi_reg_liberal_w2,
#   "az_belief_justify_liberal_w2",
#   var_willing_fight_reg_w2,
#   "az_belief_willing_w2",
#   var_interest_reg_w2,
#   "az_belief_interest_w2",
#   var_patriotism_reg_w2,
#   var_fear_critgovt_reg_w2,
#   var_socialinteract_reg_w2,
#   var_polparticipation_reg_w2,
#   var_stock_invest_reg_w2,
#   var_planaftergrad_reg_w2,
#   var_career_sector_reg_w2,
#   var_career_loc_reg_w2
# )
#
# table_a17 <- purrr::map_dfr(table_a17_vars, function(Y) {
#
#   if (!Y %in% names(cy_w2)) {
#     return(tibble())
#   }
#
#   fit <- lm(
#     reformulate(
#       c(
#         "treatment_vpnonly",
#         "treatment_nlonly",
#         "treatment_vpnnl"
#       ),
#       response = Y
#     ),
#     data = cy_w2 %>%
#       filter(treatment_user != 1)
#   )
#
#   vc <- sandwich::vcovHC(fit, type = "HC1")
#
#   reg_out <- lmtest::coeftest(fit, vcov. = vc) %>%
#     broom::tidy()
#
#   reg_out %>%
#     filter(
#       term %in% c(
#         "treatment_vpnonly",
#         "treatment_nlonly",
#         "treatment_vpnnl"
#       )
#     ) %>%
#     mutate(
#       outcome = Y,
#       specification = "Raw Diff",
#       mean_dv_all = mean(
#         cy_w2[[Y]][cy_w2$treatment_user != 1],
#         na.rm = TRUE
#       ),
#       sd_dv_all = sd(
#         cy_w2[[Y]][cy_w2$treatment_user != 1],
#         na.rm = TRUE
#       ),
#       mean_dv_control = mean(
#         cy_w2[[Y]][cy_w2$treatment_control == 1],
#         na.rm = TRUE
#       ),
#       sd_dv_control = sd(
#         cy_w2[[Y]][cy_w2$treatment_control == 1],
#         na.rm = TRUE
#       ),
#       mean_dv_user = mean(
#         cy_w2[[Y]][cy_w2$treatment_user == 1],
#         na.rm = TRUE
#       ),
#       sd_dv_user = sd(
#         cy_w2[[Y]][cy_w2$treatment_user == 1],
#         na.rm = TRUE
#       ),
#       p_value_treatment_vpnnl =
#         p.value[term == "treatment_vpnnl"][1],
#       decimal_places = 3
#     )
# })
#
# table_a17_econ <- purrr::map_dfr(
#   c(
#     var_econ_guess_reg_cn_perf_w2,
#     var_econ_guess_reg_us_perf_w2
#   ),
#   function(Y) {
#
#     if (!Y %in% names(cy_w2)) {
#       return(tibble())
#     }
#
#     fit <- lm(
#       reformulate(
#         c(
#           "treatment_vpnonly",
#           "treatment_nlonly",
#           "treatment_vpnnl"
#         ),
#         response = Y
#       ),
#       data = cy_w2 %>%
#         filter(treatment_user != 1)
#     )
#
#     vc <- sandwich::vcovHC(fit, type = "HC1")
#
#     reg_out <- lmtest::coeftest(fit, vcov. = vc) %>%
#       broom::tidy()
#
#     reg_out %>%
#       filter(
#         term %in% c(
#           "treatment_vpnonly",
#           "treatment_nlonly",
#           "treatment_vpnnl"
#         )
#       ) %>%
#       mutate(
#         outcome = Y,
#         specification = "Raw Diff",
#         mean_dv_all = mean(
#           cy_w2[[Y]][cy_w2$treatment_user != 1],
#           na.rm = TRUE
#         ),
#         sd_dv_all = sd(
#           cy_w2[[Y]][cy_w2$treatment_user != 1],
#           na.rm = TRUE
#         ),
#         mean_dv_control = mean(
#           cy_w2[[Y]][cy_w2$treatment_control == 1],
#           na.rm = TRUE
#         ),
#         sd_dv_control = sd(
#           cy_w2[[Y]][cy_w2$treatment_control == 1],
#           na.rm = TRUE
#         ),
#         mean_dv_user = mean(
#           cy_w2[[Y]][cy_w2$treatment_user == 1],
#           na.rm = TRUE
#         ),
#         sd_dv_user = sd(
#           cy_w2[[Y]][cy_w2$treatment_user == 1],
#           na.rm = TRUE
#         ),
#         p_value_treatment_vpnnl =
#           p.value[term == "treatment_vpnnl"][1],
#         decimal_places = 8
#       )
#   }
# )
#
# outregs_wave2_maineffects <- bind_rows(
#   table_a17,
#   table_a17_econ
# )
#
# write.csv(
#   outregs_wave2_maineffects,
#   file.path(outreg_dir, "outregs_wave2_maineffects.csv"),
#   row.names = FALSE
# )
