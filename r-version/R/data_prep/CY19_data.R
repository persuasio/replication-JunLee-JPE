# ******************************************************************************
# SCRIPT:   CY19_data.R
# PURPOSE:  Creates a dataset for Table 2
#           Translation of CY19_data.do
#
# NOTE: Run from repo root. Requires panelsurvey_raw.dta and auxiliary .dta
#       files in data/ChenYang2019/. ChenYang2019.dta is already pre-built;
#       this script is provided for transparency and is commented out in run.R.
# ******************************************************************************

library(haven)
library(dplyr)
library(tidyr)

data_path <- here::here("data/ChenYang2019")

# ==============================================================================
# HELPER: Anderson (2008) GLS z-score
# Replicates the Stata `andersonz` program
# ==============================================================================

andersonz <- function(df, vars) {
  # Step 1: standardize each variable (effect size relative to full-sample SD)
  effect_sizes <- lapply(vars, function(v) {
    x <- df[[v]]
    sd_x <- sd(x, na.rm = TRUE)
    mean_x <- mean(x, na.rm = TRUE)
    (x - mean_x) / sd_x
  })
  names(effect_sizes) <- paste0(vars, "e")
  es_df <- as.data.frame(effect_sizes)

  # Step 2: GLS weighting matrix (inverse of correlation matrix)
  complete_rows <- complete.cases(es_df)
  R <- cor(es_df[complete_rows, ], use = "complete.obs")
  R_inv <- solve(R)
  weights <- rowSums(R_inv)           # sum of each row of inverse = GLS weight

  # Step 3: weighted sum per observation (skip NAs, scale by non-missing weight)
  n_vars <- length(vars)
  result <- rep(0, nrow(df))
  sample_weight <- rep(0, nrow(df))

  for (i in seq_along(vars)) {
    obs <- !is.na(es_df[[i]])
    result[obs]        <- result[obs]        + es_df[[i]][obs] * weights[i]
    sample_weight[obs] <- sample_weight[obs] + weights[i]
  }

  raw <- ifelse(sample_weight > 0, result / sample_weight, NA)

  # Step 4: re-standardize
  (raw - mean(raw, na.rm = TRUE)) / sd(raw, na.rm = TRUE)
}

# ==============================================================================
# LOAD DATA
# ==============================================================================

cy <- zap_labels(read_dta(file.path(data_path, "panelsurvey_raw.dta")))

# Merge auxiliary files
vpn_date    <- zap_labels(read_dta(file.path(data_path, "vpn_date_adoption.dta")))
vpn_active  <- zap_labels(read_dta(file.path(data_path, "vpn_browsing_active_user.dta")))
vpn_purch   <- zap_labels(read_dta(file.path(data_path, "vpn_purchase.dta")))

cy <- cy |>
  left_join(vpn_date,   by = "responseID_wave1") |>
  left_join(vpn_active, by = "responseID_wave1") |>
  mutate(active_user = replace_na(active_user, 0)) |>
  left_join(vpn_purch,  by = "responseID_wave1")

# ==============================================================================
# WAVE 1: generate derived variables
# ==============================================================================

cy <- cy |> mutate(
  # News correct/incorrect indicators
  news_c_stock_rise            = if_else(!is.na(news_stock_rise),            as.integer(news_stock_rise == 1),            NA_integer_),
  news_c_exchange_depreciate   = if_else(!is.na(news_exchange_depreciate),   as.integer(news_exchange_depreciate == 1),   NA_integer_),
  news_c_train_brazil_peru     = if_else(!is.na(news_train_brazil_peru),     as.integer(news_train_brazil_peru == 0),     NA_integer_),
  news_c_army_reduce           = if_else(!is.na(news_army_reduce),           as.integer(news_army_reduce == 1),           NA_integer_),
  news_c_taiwan_election       = if_else(!is.na(news_taiwan_election),       as.integer(news_taiwan_election == 0),       NA_integer_),
  news_c_china_us_network_coop = if_else(!is.na(news_china_us_network_coop), as.integer(news_china_us_network_coop == 0), NA_integer_),
  news_c_lijiacheng_china      = if_else(!is.na(news_lijiacheng_china),      as.integer(news_lijiacheng_china == 0),      NA_integer_),
  news_c_nanjing_memory_list   = if_else(!is.na(news_nanjing_memory_list),   as.integer(news_nanjing_memory_list == 1),   NA_integer_),

  news_totalcorrect = rowSums(across(c(news_c_stock_rise, news_c_exchange_depreciate,
                                       news_c_train_brazil_peru, news_c_army_reduce, news_c_taiwan_election,
                                       news_c_china_us_network_coop, news_c_lijiacheng_china, news_c_nanjing_memory_list)),
                              na.rm = FALSE),
  news_perccorrect   = news_totalcorrect / 8,
  news_perccor_cen   = (news_c_taiwan_election + news_c_lijiacheng_china) / 2,
  news_perccor_unc   = (news_c_stock_rise + news_c_exchange_depreciate + news_c_train_brazil_peru +
                          news_c_army_reduce + news_c_china_us_network_coop + news_c_nanjing_memory_list) / 6,

  # People heard
  people_totalheard = rowSums(across(c(people_puzhiqiang_w1, people_lizehou_w1, people_huangzhifeng_w1,
                                       people_chenguangcheng_w1, people_lixiaolin_w1, people_renzhiqiang_w1,
                                       people_maoyushi_w1, people_honghuang_w1, people_liuqiangdong_w1)), na.rm = FALSE),
  people_pcheard = people_totalheard / 9,

  # Protests heard
  protest_totalheard  = rowSums(across(c(protest_2014_europe_square_w1, protest_2014_sun_flower_w1,
                                         protest_2010_arabic_spring_w1, protest_2014_crimea_vote_w1,
                                         protest_2012_hk_curriculum_w1, protest_2010_catal_indep_w1,
                                         protest_2014_umbrella_w1)), na.rm = FALSE),
  protest_pcheard_total   = protest_totalheard / 9,
  protest_pcheard_china   = (protest_2014_sun_flower_w1 + protest_2012_hk_curriculum_w1 +
                               protest_2014_umbrella_w1) / 3,
  protest_pcheard_foreign = (protest_2014_europe_square_w1 + protest_2010_arabic_spring_w1 +
                               protest_2014_crimea_vote_w1 + protest_2010_catal_indep_w1) / 4,

  # Political participation recodes
  participate_stuunion_ever = if_else(!is.na(participate_stuunion), as.integer(participate_stuunion < 3), NA_integer_),
  participate_tuanwei_ever  = if_else(!is.na(participate_tuanwei),  as.integer(participate_tuanwei < 3),  NA_integer_),
)

# Distance from neutrality and censorship indicators (wave 1)
media_vars <- c("us_pos_us","us_pos_cn","cn_pos_us","cn_pos_cn",
                "us_neg_us","us_neg_cn","cn_neg_us","cn_neg_cn")

for (v in media_vars) {
  src <- paste0("bias_", v, "_media_w1")
  cy[[paste0("distneutral_", v)]] <- if_else(!is.na(cy[[src]]), abs(cy[[src]] - 4), NA_real_)
  cy[[paste0("bias_", v, "_cens")]] <- if_else(!is.na(cy[[src]]), as.integer(cy[[src]] == 1), NA_integer_)
}

# Top-category indicators: bias_for and bias_dom (wave 1)
for (v in c("bias_for_govt_policy","bias_for_firm_interest","bias_for_media_pref","bias_for_reader_demand",
            "bias_dom_govt_policy","bias_dom_firm_interest","bias_dom_media_pref","bias_dom_reader_demand")) {
  src <- paste0(v, "_w1")
  cy[[paste0(v, "_t1")]] <- if_else(!is.na(cy[[src]]), as.integer(cy[[src]] == 1), NA_integer_)
}

# Government evaluation relative weights (wave 1)
evalgovt_vars_w1 <- c("evalgovt_election_w1","evalgovt_economy_w1","evalgovt_equality_w1",
                      "evalgovt_ruleoflaw_w1","evalgovt_human_rights_w1","evalgovt_freedom_speech_w1",
                      "evalgovt_global_power_w1","evalgovt_fair_history_w1")
cy <- cy |>
  mutate(evalgovt_total_w1 = rowSums(across(all_of(evalgovt_vars_w1)), na.rm = FALSE))
for (v in c("election","economy","equality","ruleoflaw","human_rights","freedom_speech","global_power","fair_history")) {
  cy[[paste0("revalgovt_", v)]] <- cy[[paste0("evalgovt_", v, "_w1")]] / cy$evalgovt_total_w1
}

# Graduation plan dummies (wave 1)
cy <- cy |> mutate(
  plan_grad_gradschool_dom = if_else(!is.na(plan_graduation_w1), as.integer(plan_graduation_w1 == 1), NA_integer_),
  plan_grad_foreignmaster  = if_else(!is.na(plan_graduation_w1), as.integer(plan_graduation_w1 == 2), NA_integer_),
  plan_grad_foreignphd     = if_else(!is.na(plan_graduation_w1), as.integer(plan_graduation_w1 == 3), NA_integer_),
  plan_grad_military       = if_else(!is.na(plan_graduation_w1), as.integer(plan_graduation_w1 == 4), NA_integer_),
  plan_grad_work           = if_else(!is.na(plan_graduation_w1), as.integer(plan_graduation_w1 == 5), NA_integer_),
  plan_grad_dontknow       = if_else(!is.na(plan_graduation_w1), as.integer(plan_graduation_w1 == 6), NA_integer_),
)

# Career top-1 and top-3 dummies (wave 1)
career_vars <- c("national_civil","local_civil","military","chinese_private",
                 "for_firm","soe","institutional","entrepreneur")
cy <- cy |> mutate(
  work_top3_missing_w1 = rowSums(across(paste0("work_top3_", career_vars, "_w1"), ~ is.na(.))) == length(career_vars)
)
for (v in career_vars) {
  cy[[paste0("cp_t1_", v)]] <- if_else(!cy$work_top3_missing_w1,
                                       as.integer(cy[[paste0("work_top3_", v, "_w1")]] == 1), NA_integer_)
  cy[[paste0("cp_t3_", v)]] <- if_else(!cy$work_top3_missing_w1,
                                       as.integer(!is.na(cy[[paste0("work_top3_", v, "_w1")]])), NA_integer_)
}
cy <- cy |> select(-work_top3_missing_w1)

# Career location dummies (wave 1)
cy <- cy |> mutate(
  cloc_beijing = if_else(!is.na(place_top_w1), as.integer(place_top_w1 == 1), NA_integer_),
  cloc_shanghai= if_else(!is.na(place_top_w1), as.integer(place_top_w1 == 2), NA_integer_),
  cloc_gzsz    = if_else(!is.na(place_top_w1), as.integer(place_top_w1 %in% 3:4), NA_integer_),
  cloc_tjcq    = if_else(!is.na(place_top_w1), as.integer(place_top_w1 %in% 5:6), NA_integer_),
  cloc_hkmc    = if_else(!is.na(place_top_w1), as.integer(place_top_w1 %in% 7:8), NA_integer_),
  cloc_taiwan  = if_else(!is.na(place_top_w1), as.integer(place_top_w1 == 9), NA_integer_),
  cloc_dom     = if_else(!is.na(place_top_w1), as.integer(place_top_w1 == 10), NA_integer_),
  cloc_for     = if_else(!is.na(place_top_w1), as.integer(place_top_w1 == 11), NA_integer_),
)

# Rename wave-1 derived variables with _w1 suffix
w1_rename <- c("news_totalcorrect","news_perccorrect","news_perccor_cen","news_perccor_unc",
               "people_totalheard","people_pcheard","protest_totalheard","protest_pcheard_total",
               "protest_pcheard_china","protest_pcheard_foreign",
               paste0("distneutral_", media_vars), paste0("bias_", media_vars, "_cens"),
               paste0(c("bias_for_govt_policy","bias_for_firm_interest","bias_for_media_pref","bias_for_reader_demand",
                        "bias_dom_govt_policy","bias_dom_firm_interest","bias_dom_media_pref","bias_dom_reader_demand"), "_t1"),
               paste0("revalgovt_", c("election","economy","equality","ruleoflaw","human_rights","freedom_speech","global_power","fair_history")),
               "plan_grad_gradschool_dom","plan_grad_foreignmaster","plan_grad_foreignphd",
               "plan_grad_military","plan_grad_work","plan_grad_dontknow",
               paste0("cp_t1_", career_vars), paste0("cp_t3_", career_vars),
               "cloc_beijing","cloc_shanghai","cloc_gzsz","cloc_tjcq","cloc_hkmc","cloc_taiwan","cloc_dom","cloc_for")

for (v in w1_rename) {
  if (v %in% names(cy)) cy <- cy |> rename(!!paste0(v, "_w1") := !!v)
}

# ==============================================================================
# WAVE 2: generate derived variables (same structure as wave 1)
# ==============================================================================

cy <- cy |> mutate(
  news_c_topinequality  = if_else(!is.na(news_topinequality),  as.integer(news_topinequality == 0),  NA_integer_),
  news_c_cpicensorship  = if_else(!is.na(news_cpicensorship),  as.integer(news_cpicensorship == 1),  NA_integer_),
  news_c_laborunrest    = if_else(!is.na(news_laborunrest),    as.integer(news_laborunrest == 1),    NA_integer_),
  news_c_waterpollution = if_else(!is.na(news_waterpollution), as.integer(news_waterpollution == 0), NA_integer_),
  news_c_stockcrash     = if_else(!is.na(news_stockcrash),     as.integer(news_stockcrash == 1),     NA_integer_),
  news_c_taiwanelection = if_else(!is.na(news_taiwanelection), as.integer(news_taiwanelection == 0), NA_integer_),
  news_c_applefbi       = if_else(!is.na(news_applefbi),       as.integer(news_applefbi == 0),       NA_integer_),
  news_c_tenyearshk     = if_else(!is.na(news_tenyearshk),     as.integer(news_tenyearshk == 1),     NA_integer_),
  news_c_panamapapers   = if_else(!is.na(news_panamapapers),   as.integer(news_panamapapers == 1),   NA_integer_),
  news_c_yihehotel      = if_else(!is.na(news_yihehotel),      as.integer(news_yihehotel == 1),      NA_integer_),
  news_c_economistcensor= if_else(!is.na(news_economistcensor),as.integer(news_economistcensor == 1),NA_integer_),

  news_totalcorrect = rowSums(across(c(news_c_topinequality, news_c_cpicensorship, news_c_laborunrest,
                                       news_c_waterpollution, news_c_stockcrash, news_c_taiwanelection, news_c_applefbi,
                                       news_c_tenyearshk, news_c_panamapapers, news_c_yihehotel, news_c_economistcensor)), na.rm = FALSE),
  news_perccorrect  = news_totalcorrect / 11,
  news_perccor_qui  = (news_c_topinequality + news_c_cpicensorship + news_c_laborunrest + news_c_waterpollution) / 4,
  news_perccor_cen  = (news_c_stockcrash + news_c_panamapapers + news_c_tenyearshk + news_c_economistcensor) / 4,
  news_perccor_unc  = (news_c_taiwanelection + news_c_applefbi + news_c_yihehotel) / 3,

  people_totalheard = rowSums(across(c(people_puzhiqiang_w2, people_lizehou_w2, people_huangzhifeng_w2,
                                       people_chenguangcheng_w2, people_lixiaolin_w2, people_renzhiqiang_w2,
                                       people_maoyushi_w2, people_honghuang_w2, people_liuqiangdong_w2)), na.rm = FALSE),
  people_pcheard         = people_totalheard / 9,
  people_perchd_censored = (people_lixiaolin_w2 + people_huangzhifeng_w2 +
                              people_chenguangcheng_w2 + people_lizehou_w2) / 4,
  people_perchd_uncensor = (people_maoyushi_w2 + people_liuqiangdong_w2 + people_honghuang_w2) / 3,
  people_perchd_unctocen = (people_puzhiqiang_w2 + people_renzhiqiang_w2) / 2,

  protest_totalheard  = rowSums(across(c(protest_2014_europe_square_w2, protest_2014_sun_flower_w2,
                                         protest_2010_arabic_spring_w2, protest_2014_crimea_vote_w2, protest_2012_hk_curriculum_w2,
                                         protest_2010_catal_indep_w2, protest_2014_umbrella_w2, protest_2016_mongkok_riot_w2)), na.rm = FALSE),
  protest_pcheard_total   = protest_totalheard / 8,
  protest_pcheard_china   = (protest_2014_sun_flower_w2 + protest_2012_hk_curriculum_w2 +
                               protest_2014_umbrella_w2) / 3,
  protest_pcheard_foreign = (protest_2014_europe_square_w2 + protest_2010_arabic_spring_w2 +
                               protest_2014_crimea_vote_w2 + protest_2010_catal_indep_w2) / 4,
)

for (v in media_vars) {
  src <- paste0("bias_", v, "_media_w2")
  cy[[paste0("distneutral_", v)]] <- if_else(!is.na(cy[[src]]), abs(cy[[src]] - 4), NA_real_)
  cy[[paste0("bias_", v, "_cens")]] <- if_else(!is.na(cy[[src]]), as.integer(cy[[src]] == 1), NA_integer_)
}
for (v in c("bias_for_govt_policy","bias_for_firm_interest","bias_for_media_pref","bias_for_reader_demand",
            "bias_dom_govt_policy","bias_dom_firm_interest","bias_dom_media_pref","bias_dom_reader_demand")) {
  src <- paste0(v, "_w2")
  cy[[paste0(v, "_t1")]] <- if_else(!is.na(cy[[src]]), as.integer(cy[[src]] == 1), NA_integer_)
}

evalgovt_vars_w2 <- c("evalgovt_election_w2","evalgovt_economy_w2","evalgovt_equality_w2",
                      "evalgovt_ruleoflaw_w2","evalgovt_human_rights_w2","evalgovt_freedom_speech_w2",
                      "evalgovt_global_power_w2","evalgovt_fair_history_w2")
cy <- cy |>
  mutate(evalgovt_total_w2 = rowSums(across(all_of(evalgovt_vars_w2)), na.rm = FALSE))
for (v in c("election","economy","equality","ruleoflaw","human_rights","freedom_speech","global_power","fair_history")) {
  cy[[paste0("revalgovt_", v)]] <- cy[[paste0("evalgovt_", v, "_w2")]] / cy$evalgovt_total_w2
}

cy <- cy |> mutate(
  plan_grad_gradschool_dom = if_else(!is.na(plan_graduation_w2), as.integer(plan_graduation_w2 == 1), NA_integer_),
  plan_grad_foreignmaster  = if_else(!is.na(plan_graduation_w2), as.integer(plan_graduation_w2 == 2), NA_integer_),
  plan_grad_foreignphd     = if_else(!is.na(plan_graduation_w2), as.integer(plan_graduation_w2 == 3), NA_integer_),
  plan_grad_military       = if_else(!is.na(plan_graduation_w2), as.integer(plan_graduation_w2 == 4), NA_integer_),
  plan_grad_work           = if_else(!is.na(plan_graduation_w2), as.integer(plan_graduation_w2 == 5), NA_integer_),
  plan_grad_dontknow       = if_else(!is.na(plan_graduation_w2), as.integer(plan_graduation_w2 == 6), NA_integer_),
)

cy <- cy |> mutate(
  work_top3_missing_w2 = rowSums(across(paste0("work_top3_", career_vars, "_w2"), ~ is.na(.))) == length(career_vars)
)
for (v in career_vars) {
  cy[[paste0("cp_t1_", v)]] <- if_else(!cy$work_top3_missing_w2,
                                       as.integer(cy[[paste0("work_top3_", v, "_w2")]] == 1), NA_integer_)
  cy[[paste0("cp_t3_", v)]] <- if_else(!cy$work_top3_missing_w2,
                                       as.integer(!is.na(cy[[paste0("work_top3_", v, "_w2")]])), NA_integer_)
}
cy <- cy |> select(-work_top3_missing_w2)

cy <- cy |> mutate(
  cloc_beijing = if_else(!is.na(place_top_w2), as.integer(place_top_w2 == 1), NA_integer_),
  cloc_shanghai= if_else(!is.na(place_top_w2), as.integer(place_top_w2 == 2), NA_integer_),
  cloc_gzsz    = if_else(!is.na(place_top_w2), as.integer(place_top_w2 %in% 3:4), NA_integer_),
  cloc_tjcq    = if_else(!is.na(place_top_w2), as.integer(place_top_w2 %in% 5:6), NA_integer_),
  cloc_hkmc    = if_else(!is.na(place_top_w2), as.integer(place_top_w2 %in% 7:8), NA_integer_),
  cloc_taiwan  = if_else(!is.na(place_top_w2), as.integer(place_top_w2 == 9), NA_integer_),
  cloc_dom     = if_else(!is.na(place_top_w2), as.integer(place_top_w2 == 10), NA_integer_),
  cloc_for     = if_else(!is.na(place_top_w2), as.integer(place_top_w2 == 11), NA_integer_),
)

w2_rename <- c("news_totalcorrect","news_perccorrect","news_perccor_cen","news_perccor_unc",
               "news_perccor_qui","people_totalheard","people_pcheard",
               "people_perchd_censored","people_perchd_uncensor","people_perchd_unctocen",
               "protest_totalheard","protest_pcheard_total","protest_pcheard_china","protest_pcheard_foreign",
               paste0("distneutral_", media_vars), paste0("bias_", media_vars, "_cens"),
               paste0(c("bias_for_govt_policy","bias_for_firm_interest","bias_for_media_pref","bias_for_reader_demand",
                        "bias_dom_govt_policy","bias_dom_firm_interest","bias_dom_media_pref","bias_dom_reader_demand"), "_t1"),
               paste0("revalgovt_", c("election","economy","equality","ruleoflaw","human_rights","freedom_speech","global_power","fair_history")),
               "plan_grad_gradschool_dom","plan_grad_foreignmaster","plan_grad_foreignphd",
               "plan_grad_military","plan_grad_work","plan_grad_dontknow",
               paste0("cp_t1_", career_vars), paste0("cp_t3_", career_vars),
               "cloc_beijing","cloc_shanghai","cloc_gzsz","cloc_tjcq","cloc_hkmc","cloc_taiwan","cloc_dom","cloc_for")

for (v in w2_rename) {
  if (v %in% names(cy)) cy <- cy |> rename(!!paste0(v, "_w2") := !!v)
}

# ==============================================================================
# WAVE 3: generate derived variables
# ==============================================================================

cy <- cy |> mutate(
  news_c_coalprod    = if_else(!is.na(news_coalprod),    as.integer(news_coalprod == 0),    NA_integer_),
  news_c_trumpchina  = if_else(!is.na(news_trumpchina),  as.integer(news_trumpchina == 1),  NA_integer_),
  news_c_xiaojianhua = if_else(!is.na(news_xiaojianhua), as.integer(news_xiaojianhua == 0), NA_integer_),
  news_c_northkoreacoal = if_else(!is.na(news_northkoreacoal), as.integer(news_northkoreacoal == 1), NA_integer_),
  news_c_birdflu     = if_else(!is.na(news_birdflu),     as.integer(news_birdflu == 1),     NA_integer_),
  news_c_ethiopiatrain = if_else(!is.na(news_ethiopiatrain), as.integer(news_ethiopiatrain == 0), NA_integer_),
  news_c_foreignreserve = if_else(!is.na(news_foreignreserve), as.integer(news_foreignreserve == 0), NA_integer_),
  news_c_xijiangcar  = if_else(!is.na(news_xijiangcar),  as.integer(news_xijiangcar == 0),  NA_integer_),
  news_c_chinanorway = if_else(!is.na(news_chinanorway), as.integer(news_chinanorway == 1), NA_integer_),
  news_c_womenrights = if_else(!is.na(news_womenrights), as.integer(news_womenrights == 0), NA_integer_),
  news_c_hkceelection= if_else(!is.na(news_hkceelection),as.integer(news_hkceelection == 0),NA_integer_),

  news_totalcorrect = rowSums(across(c(news_c_coalprod, news_c_trumpchina, news_c_xiaojianhua,
                                       news_c_northkoreacoal, news_c_birdflu, news_c_ethiopiatrain, news_c_foreignreserve,
                                       news_c_xijiangcar, news_c_chinanorway, news_c_womenrights, news_c_hkceelection)), na.rm = FALSE),
  news_perccorrect  = news_totalcorrect / 11,
  news_perccor_cen  = (news_c_coalprod + news_c_trumpchina + news_c_xiaojianhua + news_c_xijiangcar +
                         news_c_chinanorway + news_c_womenrights + news_c_hkceelection) / 7,
  news_perccor_unc  = (news_c_northkoreacoal + news_c_birdflu + news_c_ethiopiatrain + news_c_foreignreserve) / 4,

  protest_totalheard  = rowSums(across(c(protest_2014_europe_square_w3, protest_2014_sun_flower_w3,
                                         protest_2010_arabic_spring_w3, protest_2014_crimea_vote_w3, protest_2012_hk_curriculum_w3,
                                         protest_2010_catal_indep_w3, protest_2014_umbrella_w3, protest_2016_mongkok_riot_w3,
                                         protest_2017_women_march_w3)), na.rm = FALSE),
  protest_pcheard_total   = protest_totalheard / 9,
  protest_pcheard_china   = (protest_2014_sun_flower_w3 + protest_2012_hk_curriculum_w3 +
                               protest_2014_umbrella_w3) / 3,
  protest_pcheard_foreign = (protest_2014_europe_square_w3 + protest_2010_arabic_spring_w3 +
                               protest_2014_crimea_vote_w3 + protest_2010_catal_indep_w3) / 4,
)

for (v in c("bias_for_govt_policy","bias_for_firm_interest","bias_for_media_pref","bias_for_reader_demand",
            "bias_dom_govt_policy","bias_dom_firm_interest","bias_dom_media_pref","bias_dom_reader_demand")) {
  src <- paste0(v, "_w3")
  cy[[paste0(v, "_t1")]] <- if_else(!is.na(cy[[src]]), as.integer(cy[[src]] == 1), NA_integer_)
}

cy <- cy |> mutate(
  vpn_purchase_yes     = if_else(!is.na(vpn_purchase), as.integer(vpn_purchase < 6),  NA_integer_),
  vpn_purchase_wmt     = if_else(!is.na(vpn_purchase), as.integer(vpn_purchase_yes == 1 & vpn_purchase != 1), NA_integer_),
  vpn_purchase_premium = case_when(
    vpn_purchase == 6 ~ 0L, vpn_purchase == 1 ~ 1L, vpn_purchase == 4 ~ 2L,
    vpn_purchase == 3 ~ 3L, vpn_purchase == 5 ~ 4L, TRUE ~ NA_integer_),

  plan_grad_gradschool_dom = if_else(!is.na(plan_graduation_w3), as.integer(plan_graduation_w3 == 1), NA_integer_),
  plan_grad_foreignmaster  = if_else(!is.na(plan_graduation_w3), as.integer(plan_graduation_w3 == 2), NA_integer_),
  plan_grad_foreignphd     = if_else(!is.na(plan_graduation_w3), as.integer(plan_graduation_w3 == 3), NA_integer_),
  plan_grad_military       = if_else(!is.na(plan_graduation_w3), as.integer(plan_graduation_w3 == 4), NA_integer_),
  plan_grad_work           = if_else(!is.na(plan_graduation_w3), as.integer(plan_graduation_w3 == 5), NA_integer_),
  plan_grad_dontknow       = if_else(!is.na(plan_graduation_w3), as.integer(plan_graduation_w3 == 6), NA_integer_),
)

cy <- cy |> mutate(
  work_top3_missing_w3 = rowSums(across(paste0("work_top3_", career_vars, "_w3"), ~ is.na(.))) == length(career_vars)
)
for (v in career_vars) {
  cy[[paste0("cp_t1_", v)]] <- if_else(!cy$work_top3_missing_w3,
                                       as.integer(cy[[paste0("work_top3_", v, "_w3")]] == 1), NA_integer_)
  cy[[paste0("cp_t3_", v)]] <- if_else(!cy$work_top3_missing_w3,
                                       as.integer(!is.na(cy[[paste0("work_top3_", v, "_w3")]])), NA_integer_)
}
cy <- cy |> select(-work_top3_missing_w3)

cy <- cy |> mutate(
  cloc_beijing = if_else(!is.na(place_top_w3), as.integer(place_top_w3 == 1), NA_integer_),
  cloc_shanghai= if_else(!is.na(place_top_w3), as.integer(place_top_w3 == 2), NA_integer_),
  cloc_gzsz    = if_else(!is.na(place_top_w3), as.integer(place_top_w3 %in% 3:4), NA_integer_),
  cloc_tjcq    = if_else(!is.na(place_top_w3), as.integer(place_top_w3 %in% 5:6), NA_integer_),
  cloc_hkmc    = if_else(!is.na(place_top_w3), as.integer(place_top_w3 %in% 7:8), NA_integer_),
  cloc_taiwan  = if_else(!is.na(place_top_w3), as.integer(place_top_w3 == 9), NA_integer_),
  cloc_dom     = if_else(!is.na(place_top_w3), as.integer(place_top_w3 == 10), NA_integer_),
  cloc_for     = if_else(!is.na(place_top_w3), as.integer(place_top_w3 == 11), NA_integer_),
)

w3_rename <- c("news_totalcorrect","news_perccorrect","news_perccor_cen","news_perccor_unc",
               "protest_totalheard","protest_pcheard_total","protest_pcheard_china","protest_pcheard_foreign",
               paste0(c("bias_for_govt_policy","bias_for_firm_interest","bias_for_media_pref","bias_for_reader_demand",
                        "bias_dom_govt_policy","bias_dom_firm_interest","bias_dom_media_pref","bias_dom_reader_demand"), "_t1"),
               "plan_grad_gradschool_dom","plan_grad_foreignmaster","plan_grad_foreignphd",
               "plan_grad_military","plan_grad_work","plan_grad_dontknow",
               paste0("cp_t1_", career_vars), paste0("cp_t3_", career_vars),
               "cloc_beijing","cloc_shanghai","cloc_gzsz","cloc_tjcq","cloc_hkmc","cloc_taiwan","cloc_dom","cloc_for")

for (v in w3_rename) {
  if (v %in% names(cy)) cy <- cy |> rename(!!paste0(v, "_w3") := !!v)
}

# Cross-wave average
cy <- cy |> mutate(
  news_perccor_cen_all = (news_perccor_cen_w2 + news_perccor_cen_w3) / 2
)

# ==============================================================================
# RISK PREFERENCE: Falk et al. (2015) staircase certainty equivalent
# ==============================================================================

# Build lookup table matching the 32 replace conditions
risk_ce_lookup <- function(r) {
  r1 <- r["risk_preference_1"]; r2 <- r["risk_preference_2"]
  r3 <- r["risk_preference_3"]; r4 <- r["risk_preference_4"]
  r5 <- r["risk_preference_5"]; r6 <- r["risk_preference_6"]
  r7 <- r["risk_preference_7"]; r8 <- r["risk_preference_8"]
  r9 <- r["risk_preference_9"]; r10<- r["risk_preference_10"]
  r11<- r["risk_preference_11"];r12<- r["risk_preference_12"]
  r13<- r["risk_preference_13"];r14<- r["risk_preference_14"]
  r15<- r["risk_preference_15"];r16<- r["risk_preference_16"]
  r17<- r["risk_preference_17"];r18<- r["risk_preference_18"]
  r19<- r["risk_preference_19"];r20<- r["risk_preference_20"]
  r21<- r["risk_preference_21"];r22<- r["risk_preference_22"]
  r23<- r["risk_preference_23"];r24<- r["risk_preference_24"]
  r25<- r["risk_preference_25"];r26<- r["risk_preference_26"]
  r27<- r["risk_preference_27"];r28<- r["risk_preference_28"]
  r29<- r["risk_preference_29"];r30<- r["risk_preference_30"]
  r31<- r["risk_preference_31"]
  if (anyNA(c(r1,r17))) return(NA_integer_)
  dplyr::case_when(
    r1==1&r17==1&r25==1&r29==1&r31==1 ~ 32L, r1==1&r17==1&r25==1&r29==1&r31==2 ~ 31L,
    r1==1&r17==1&r25==1&r29==2&r30==1 ~ 30L, r1==1&r17==1&r25==1&r29==2&r30==2 ~ 29L,
    r1==1&r17==1&r25==2&r26==1&r27==1 ~ 28L, r1==1&r17==1&r25==2&r26==1&r27==2 ~ 27L,
    r1==1&r17==1&r25==2&r26==2&r28==1 ~ 26L, r1==1&r17==1&r25==2&r26==2&r28==2 ~ 25L,
    r1==1&r17==2&r18==1&r22==1&r23==1 ~ 24L, r1==1&r17==2&r18==1&r22==1&r23==2 ~ 23L,
    r1==1&r17==2&r18==1&r22==2&r24==1 ~ 22L, r1==1&r17==2&r18==1&r22==2&r24==2 ~ 21L,
    r1==1&r17==2&r18==2&r19==1&r20==1 ~ 20L, r1==1&r17==2&r18==2&r19==1&r20==2 ~ 19L,
    r1==1&r17==2&r18==2&r19==2&r21==1 ~ 18L, r1==1&r17==2&r18==2&r19==2&r21==2 ~ 17L,
    r1==2&r2==1&r10==1&r14==1&r15==1  ~ 16L, r1==2&r2==1&r10==1&r14==1&r15==2  ~ 15L,
    r1==2&r2==1&r10==1&r14==2&r16==1  ~ 14L, r1==2&r2==1&r10==1&r14==2&r16==2  ~ 13L,
    r1==2&r2==1&r10==2&r11==1&r13==1  ~ 12L, r1==2&r2==1&r10==2&r11==1&r13==2  ~ 11L,
    r1==2&r2==1&r10==2&r11==2&r12==1  ~ 10L, r1==2&r2==1&r10==2&r11==2&r12==2  ~  9L,
    r1==2&r2==2&r3==1&r4==1&r5==1     ~  8L, r1==2&r2==2&r3==1&r4==1&r5==2     ~  7L,
    r1==2&r2==2&r3==1&r4==2&r6==1     ~  6L, r1==2&r2==2&r3==1&r4==2&r6==2     ~  5L,
    r1==2&r2==2&r3==2&r7==1&r8==1     ~  4L, r1==2&r2==2&r3==2&r7==1&r8==2     ~  3L,
    r1==2&r2==2&r3==2&r7==2&r9==1     ~  2L, r1==2&r2==2&r3==2&r7==2&r9==2     ~  1L,
    TRUE ~ NA_integer_)
}

risk_cols <- paste0("risk_preference_", 1:31)
cy$risk_preference_ce <- apply(cy[, risk_cols], 1, risk_ce_lookup)

# Time preference (same staircase structure — note Stata code has different branching order)
time_ce_lookup <- function(r) {
  # Mirrors the time_preference replace conditions in the .do file exactly
  r1<-r["time_preference_1"];r2<-r["time_preference_2"];r3<-r["time_preference_3"]
  r4<-r["time_preference_4"];r5<-r["time_preference_5"];r6<-r["time_preference_6"]
  r7<-r["time_preference_7"];r8<-r["time_preference_8"];r9<-r["time_preference_9"]
  r10<-r["time_preference_10"];r11<-r["time_preference_11"];r12<-r["time_preference_12"]
  r13<-r["time_preference_13"];r14<-r["time_preference_14"];r15<-r["time_preference_15"]
  r16<-r["time_preference_16"];r17<-r["time_preference_17"];r18<-r["time_preference_18"]
  r19<-r["time_preference_19"];r20<-r["time_preference_20"];r21<-r["time_preference_21"]
  r22<-r["time_preference_22"];r23<-r["time_preference_23"];r24<-r["time_preference_24"]
  r25<-r["time_preference_25"];r26<-r["time_preference_26"];r27<-r["time_preference_27"]
  r28<-r["time_preference_28"];r29<-r["time_preference_29"];r30<-r["time_preference_30"]
  r31<-r["time_preference_31"]
  if (anyNA(c(r1,r17))) return(NA_integer_)
  dplyr::case_when(
    r1==1&r17==1&r18==1&r22==1&r23==1 ~ 32L, r1==1&r17==1&r18==1&r22==1&r23==2 ~ 31L,
    r1==1&r17==1&r18==1&r22==2&r24==1 ~ 30L, r1==1&r17==1&r18==1&r22==2&r24==2 ~ 29L,
    r1==1&r17==1&r18==2&r19==1&r20==1 ~ 28L, r1==1&r17==1&r18==2&r19==1&r20==2 ~ 27L,
    r1==1&r17==1&r18==2&r19==2&r21==1 ~ 26L, r1==1&r17==1&r18==2&r19==2&r21==2 ~ 25L,
    r1==1&r17==2&r25==1&r29==1&r31==1 ~ 24L, r1==1&r17==2&r25==1&r29==1&r31==2 ~ 23L,
    r1==1&r17==2&r25==1&r29==2&r30==1 ~ 22L, r1==1&r17==2&r25==1&r29==2&r30==2 ~ 21L,
    r1==1&r17==2&r25==2&r26==1&r28==1 ~ 20L, r1==1&r17==2&r25==2&r26==1&r28==2 ~ 19L,
    r1==1&r17==2&r25==2&r26==2&r27==1 ~ 18L, r1==1&r17==2&r25==2&r26==2&r27==2 ~ 17L,
    r1==2&r2==1&r10==1&r14==1&r16==1  ~ 16L, r1==2&r2==1&r10==1&r14==1&r16==2  ~ 15L,
    r1==2&r2==1&r10==1&r14==2&r15==1  ~ 14L, r1==2&r2==1&r10==1&r14==2&r15==2  ~ 13L,
    r1==2&r2==1&r10==2&r11==1&r13==1  ~ 12L, r1==2&r2==1&r10==2&r11==1&r13==2  ~ 11L,
    r1==2&r2==1&r10==2&r11==2&r12==1  ~ 10L, r1==2&r2==1&r10==2&r11==2&r12==2  ~  9L,
    r1==2&r2==2&r3==1&r7==1&r8==1     ~  8L, r1==2&r2==2&r3==1&r7==1&r8==2     ~  7L,
    r1==2&r2==2&r3==1&r7==2&r9==1     ~  6L, r1==2&r2==2&r3==1&r7==2&r9==2     ~  5L,
    r1==2&r2==2&r3==2&r4==1&r6==1     ~  4L, r1==2&r2==2&r3==2&r4==1&r6==2     ~  3L,
    r1==2&r2==2&r3==2&r4==2&r5==1     ~  2L, r1==2&r2==2&r3==2&r4==2&r5==2     ~  1L,
    TRUE ~ NA_integer_)
}

time_cols <- paste0("time_preference_", 1:31)
cy$time_preference_fe <- apply(cy[, time_cols], 1, time_ce_lookup)

# ==============================================================================
# DEMOGRAPHICS
# ==============================================================================

coastal_provinces <- c("Beijing","Fujian","Guangdong","Hainan","Hebei","Jiangsu",
                       "Shandong","Shanghai","Tianjin","Zhejiang","Non-mainland")

cy <- cy |> mutate(
  ethnicity_han          = if_else(!is.na(ethnicity),          as.integer(ethnicity == 1),          NA_integer_),
  hukou_urban            = if_else(!is.na(hukou),              as.integer(hukou == 1),              NA_integer_),
  siblings_total         = number_bro_younger + number_bro_older + number_sis_younger + number_sis_older,
  father_edu_hsabove     = if_else(!is.na(father_education),   as.integer(father_education > 3),    NA_integer_),
  mother_edu_hsabove     = if_else(!is.na(mother_education),   as.integer(mother_education > 3),    NA_integer_),
  religion_religious     = if_else(!is.na(religion),           as.integer(religion > 1),            NA_integer_),
  hs_track_science       = if_else(!is.na(hs_track),           as.integer(hs_track == 1),           NA_integer_),
  university_elite       = if_else(!is.na(university),         as.integer(university == 1),         NA_integer_),
  birth_year             = recode(birth_year, `1`=1990L,`2`=1991L,`3`=1992L,`4`=1993L,`5`=1994L,
                                  `6`=1995L,`7`=1996L,`8`=1997L,`9`=1998L,`10`=1999L,
                                  `11`=2000L,`12`=2001L,`13`=2002L,`14`=2003L,`15`=2004L,`16`=2005L),
  birthplace_coastal     = as.integer(birthplace_province %in% coastal_provinces),
  residence_coastal      = as.integer(residence_province  %in% coastal_provinces),
  domestic_english_atleast4 = if_else(!is.na(english_qual_domestic), as.integer(english_qual_domestic > 1), NA_integer_),
  foreign_english_yes    = if_else(!is.na(english_qual_foreign),  as.integer(english_qual_foreign < 3),  NA_integer_),
  work_father_govt       = if_else(!is.na(father_work), as.integer(father_work <= 4), NA_integer_),
  work_mother_govt       = if_else(!is.na(mother_work), as.integer(mother_work <= 4), NA_integer_),
  travel_foreign_yes     = if_else(!is.na(travel_foreign), as.integer(travel_foreign > 1), NA_integer_),
)

# ==============================================================================
# TREATMENT INDICATORS
# ==============================================================================

cy <- cy |> mutate(
  treatment_control = as.integer(treatment_newsletter == 0 & treatment_vpn == 0 & vpn_current_paid_user == 0),
  treatment_vpnonly = as.integer(treatment_newsletter == 0 & treatment_vpn == 1 & vpn_current_paid_user == 0),
  treatment_nlonly  = as.integer(treatment_newsletter == 1 & treatment_vpn == 0 & vpn_current_paid_user == 0),
  treatment_vpnnl   = as.integer(treatment_newsletter == 1 & treatment_vpn == 1 & vpn_current_paid_user == 0),
  treatment_user    = as.integer(vpn_current_paid_user == 1),

  treatment_master = case_when(
    treatment_control == 1 ~ 1L,
    treatment_vpnonly == 1 ~ 2L,
    treatment_nlonly  == 1 ~ 3L,
    treatment_vpnnl   == 1 ~ 4L,
    treatment_user    == 1 ~ 5L,
    TRUE ~ NA_integer_),

  treatment_main = case_when(
    treatment_control == 1 | treatment_nlonly == 1 ~ 1L,
    treatment_vpnonly == 1 ~ 2L,
    treatment_vpnnl   == 1 ~ 3L,
    treatment_user    == 1 ~ 4L,
    TRUE ~ NA_integer_),

  treatvpnXnl = treatment_vpn * treatment_newsletter,
)

# ==============================================================================
# VPN ADOPTION & SOCIAL LEARNING
# ==============================================================================

cy <- cy |> mutate(
  vpn_adopted = case_when(
    !is.na(vpn_date_adoption) ~ 1L,
    is.na(vpn_date_adoption) & treatment_master %in% c(2L, 4L) ~ 0L,
    TRUE ~ NA_integer_),

  vpn_roommate_existing = pmin(replace_na(vpn_roommate_usage_w1, 0), 2),
  vpn_roommate_new_w2   = pmax(pmin(vpn_roommate_usage_w2 - vpn_roommate_usage_w1, 2), 0),
  vpn_roommate_new_w3   = pmax(pmin(vpn_roommate_usage_w3 - vpn_roommate_usage_w1, 2), 0),

  soclearning_ownaccess   = if_else(!is.na(treatment_master), as.integer(treatment_master >= 4), NA_integer_),
  soclearning_rm_new_w2   = vpn_roommate_new_w2,
  soclearning_rm_new_w3   = vpn_roommate_new_w3,
  soclearning_rm_existing = vpn_roommate_existing,
  soclearning_ownXnew_w2  = soclearning_ownaccess * soclearning_rm_new_w2,
  soclearning_ownXnew_w3  = soclearning_ownaccess * soclearning_rm_new_w3,
)

# ==============================================================================
# LIST EXPERIMENT
# ==============================================================================

cy <- cy |> mutate(
  list_exp_veiled = 1L - list_exp_direct,
  list_countd_trust_direct = list_trust_direct_count + list_trust_direct_yes,
  list_count_trust = if_else(list_exp_direct == 1, list_countd_trust_direct, list_trust_veiled_count),
)

# ==============================================================================
# ANDERSON Z-SCORES
# ==============================================================================

# Variable group definitions (mirrors Stata locals)
var_groups <- list(
  vpn_purchase_reg_w3       = c("vpn_purchase_wmt_record","vpn_purchase_yes"),
  media_valuation_reg_w1    = c("wtp_vpn_w1","added_value_foreign_media_w1"),
  media_valuation_reg_w2    = c("wtp_vpn_w2","added_value_foreign_media_w2"),
  media_valuation_reg_w3    = c("wtp_vpn_w3","added_value_foreign_media_w3"),
  media_trust_reg_w1        = c("trust_media_dom_state_w1","trust_media_dom_private_w1","trust_media_foreign_w1"),
  media_trust_reg_w2        = c("trust_media_dom_state_w2","trust_media_dom_private_w2","trust_media_foreign_w2"),
  media_trust_reg_w3        = c("trust_media_dom_state_w3","trust_media_dom_private_w3","trust_media_foreign_w3"),
  censor_justif_reg_fp_w1   = c("censor_just_dom_economic_w1","censor_just_dom_political_w1","censor_just_dom_social_w1","censor_just_for_w1"),
  censor_justif_reg_fp_w2   = c("censor_just_dom_economic_w2","censor_just_dom_political_w2","censor_just_dom_social_w2","censor_just_for_w2"),
  censor_justif_reg_fp_w3   = c("censor_just_dom_economic_w3","censor_just_dom_political_w3","censor_just_dom_social_w3","censor_just_for_w3"),
  knowledge_meta_reg_w1     = c("familiar_china_issues_self_w1","familiar_china_others_w1"),
  knowledge_meta_reg_w2     = c("familiar_china_issues_self_w2","familiar_china_others_w2"),
  knowledge_meta_reg_w3     = c("familiar_china_issues_self_w3","familiar_china_others_w3"),
  econ_guess_reg_cn_perf_w1 = c("guess_gdp_growth_china_w1","guess_stock_index_sh_w1"),
  econ_guess_reg_cn_perf_w2 = c("guess_gdp_growth_china_w2","guess_stock_index_sh_w2"),
  econ_guess_reg_cn_perf_w3 = c("guess_gdp_growth_china_w3","guess_stock_index_sh_w3"),
  econ_guess_reg_cn_conf_w1 = c("guess_gdp_growth_china_con_w1","guess_stock_index_sh_con_w1"),
  econ_guess_reg_cn_conf_w2 = c("guess_gdp_growth_china_con_w2","guess_stock_index_sh_con_w2"),
  econ_guess_reg_cn_conf_w3 = c("guess_gdp_growth_china_con_w3","guess_stock_index_sh_con_w3"),
  demand_change_reg_w1      = c("inst_change_econ_w1","inst_change_poli_w1"),
  demand_change_reg_w2      = c("inst_change_econ_w2","inst_change_poli_w2"),
  demand_change_reg_w3      = c("inst_change_econ_w3","inst_change_poli_w3"),
  trust_inst_reg_govt_w1    = c("trust_central_govt_w1","trust_provincial_govt_w1","trust_local_govt_w1"),
  trust_inst_reg_govt_w2    = c("trust_central_govt_w2","trust_provincial_govt_w2","trust_local_govt_w2"),
  trust_inst_reg_govt_w3    = c("trust_central_govt_w3","trust_provincial_govt_w3","trust_local_govt_w3"),
  trust_inst_reg_foreign_w1 = c("trust_japan_govt_w1","trust_us_govt_w1"),
  trust_inst_reg_foreign_w2 = c("trust_japan_govt_w2","trust_us_govt_w2"),
  trust_inst_reg_foreign_w3 = c("trust_japan_govt_w3","trust_us_govt_w3"),
  eval_govt_reg_w1          = c("eval_govt_economic_w1","eval_govt_dom_politics_w1","eval_govt_for_relations_w1"),
  eval_govt_reg_w2          = c("eval_govt_economic_w2","eval_govt_dom_politics_w2","eval_govt_for_relations_w2"),
  eval_govt_reg_w3          = c("eval_govt_economic_w3","eval_govt_dom_politics_w3","eval_govt_for_relations_w3"),
  democracy_reg_fp_w1       = c("importance_live_in_demo_w1"),
  democracy_reg_fp_w2       = c("importance_live_in_demo_w2"),
  democracy_reg_fp_w3       = c("importance_live_in_demo_w3"),
  willing_fight_reg_w1      = c("willing_against_illi_govt_w1","willing_report_mis_w1","willing_protect_weak_w1"),
  willing_fight_reg_w2      = c("willing_against_illi_govt_w2","willing_report_mis_w2","willing_protect_weak_w2"),
  willing_fight_reg_w3      = c("willing_against_illi_govt_w3","willing_report_mis_w3","willing_protect_weak_w3"),
  socialinteract_reg_w1     = c("frequency_talk_politic_w1","frequency_persuade_friends_w1"),
  socialinteract_reg_w2     = c("frequency_talk_politic_w2","frequency_persuade_friends_w2"),
  socialinteract_reg_w3     = c("frequency_talk_politic_w3","frequency_persuade_friends_w3"),
  polparticipation_reg_pf_w1= c("participate_social_protest_w1","participate_plan_vote_w1","participate_complain_school_w1"),
  polparticipation_reg_pf_w2= c("participate_social_protest_w2","participate_plan_vote_w2","participate_complain_school_w2"),
  polparticipation_reg_pf_w3= c("participate_social_protest_w3","participate_plan_vote_w3","participate_complain_school_w3")
)

az_map <- list(
  az_var_vpnpurchase      = "vpn_purchase_reg_w3",
  az_belief_media_value_w1= "media_valuation_reg_w1",
  az_belief_media_value_w2= "media_valuation_reg_w2",
  az_belief_media_value_w3= "media_valuation_reg_w3",
  az_belief_media_trust_w1= "media_trust_reg_w1",
  az_belief_media_trust_w2= "media_trust_reg_w2",
  az_belief_media_trust_w3= "media_trust_reg_w3",
  az_belief_media_justif_pf_w1 = "censor_justif_reg_fp_w1",
  az_belief_media_justif_pf_w2 = "censor_justif_reg_fp_w2",
  az_belief_media_justif_pf_w3 = "censor_justif_reg_fp_w3",
  az_knowledge_meta_w1    = "knowledge_meta_reg_w1",
  az_knowledge_meta_w2    = "knowledge_meta_reg_w2",
  az_knowledge_meta_w3    = "knowledge_meta_reg_w3",
  az_belief_econ_perf_cn_w1 = "econ_guess_reg_cn_perf_w1",
  az_belief_econ_perf_cn_w2 = "econ_guess_reg_cn_perf_w2",
  az_belief_econ_perf_cn_w3 = "econ_guess_reg_cn_perf_w3",
  az_belief_econ_conf_cn_w1 = "econ_guess_reg_cn_conf_w1",
  az_belief_econ_conf_cn_w2 = "econ_guess_reg_cn_conf_w2",
  az_belief_econ_conf_cn_w3 = "econ_guess_reg_cn_conf_w3",
  az_belief_instchange_w1 = "demand_change_reg_w1",
  az_belief_instchange_w2 = "demand_change_reg_w2",
  az_belief_instchange_w3 = "demand_change_reg_w3",
  az_belief_trust_govt_w1 = "trust_inst_reg_govt_w1",
  az_belief_trust_govt_w2 = "trust_inst_reg_govt_w2",
  az_belief_trust_govt_w3 = "trust_inst_reg_govt_w3",
  az_belief_trust_foreign_w1 = "trust_inst_reg_foreign_w1",
  az_belief_trust_foreign_w2 = "trust_inst_reg_foreign_w2",
  az_belief_trust_foreign_w3 = "trust_inst_reg_foreign_w3",
  az_belief_evalgovt_w1   = "eval_govt_reg_w1",
  az_belief_evalgovt_w2   = "eval_govt_reg_w2",
  az_belief_evalgovt_w3   = "eval_govt_reg_w3",
  az_belief_democracy_fp_w1 = "democracy_reg_fp_w1",
  az_belief_democracy_fp_w2 = "democracy_reg_fp_w2",
  az_belief_democracy_fp_w3 = "democracy_reg_fp_w3",
  az_belief_willing_w1    = "willing_fight_reg_w1",
  az_belief_willing_w2    = "willing_fight_reg_w2",
  az_belief_willing_w3    = "willing_fight_reg_w3",
  az_var_socialinteract_w1= "socialinteract_reg_w1",
  az_var_socialinteract_w2= "socialinteract_reg_w2",
  az_var_socialinteract_w3= "socialinteract_reg_w3",
  az_var_polparticipation_w1 = "polparticipation_reg_pf_w1",
  az_var_polparticipation_w2 = "polparticipation_reg_pf_w2",
  az_var_polparticipation_w3 = "polparticipation_reg_pf_w3"
)

for (az_name in names(az_map)) {
  grp <- var_groups[[az_map[[az_name]]]]
  cy[[az_name]] <- andersonz(cy, grp)
}

# ==============================================================================
# OVERALL EFFECT INDEX: flipped variables + mega z-scores
# ==============================================================================

for (i in 1:3) {
  cy[[paste0("bias_foreign_r_w", i)]]         <- 10 - cy[[paste0("bias_foreign_w", i)]]
  cy[[paste0("az_belief_econ_perf_cn_r_w", i)]] <- -cy[[paste0("az_belief_econ_perf_cn_w", i)]]
  cy[[paste0("az_belief_trust_govt_r_w", i)]]  <- -cy[[paste0("az_belief_trust_govt_w", i)]]
  cy[[paste0("az_belief_evalgovt_r_w", i)]]    <- -cy[[paste0("az_belief_evalgovt_w", i)]]
  cy[[paste0("importance_live_in_demo_r_w", i)]] <- 10 - cy[[paste0("importance_live_in_demo_w", i)]]
}

# Category sub-indices per wave
for (i in 1:3) {
  wi <- paste0("_w", i)
  cy[[paste0("az_overall_a", wi)]] <- andersonz(cy, c(
    paste0("info_freq_website_for", wi), paste0("az_belief_media_value", wi),
    paste0("az_belief_media_trust", wi), paste0("az_belief_media_justif_pf", wi),
    paste0("bias_domestic", wi), paste0("bias_foreign_r", wi)))
  cy[[paste0("az_overall_b", wi)]] <- andersonz(cy, c(
    paste0("news_perccor_cen", wi), paste0("news_perccor_unc", wi),
    paste0("protest_pcheard_china", wi), paste0("protest_pcheard_foreign", wi),
    paste0("az_knowledge_meta", wi)))
  cy[[paste0("az_overall_c", wi)]] <- andersonz(cy, c(
    paste0("az_belief_econ_perf_cn_r", wi), paste0("az_belief_econ_conf_cn", wi)))
  cy[[paste0("az_overall_d", wi)]] <- andersonz(cy, c(
    paste0("az_belief_instchange", wi), paste0("az_belief_trust_govt_r", wi),
    paste0("az_belief_trust_foreign", wi), paste0("az_belief_evalgovt_r", wi),
    paste0("importance_live_in_demo_r", wi), paste0("az_belief_willing", wi)))
  plan_var <- if (i == 1) "plan_grad_foreignmaster_w1" else paste0("plan_grad_foreignmaster_w", i)
  cy[[paste0("az_overall_e", wi)]] <- andersonz(cy, c(
    paste0("az_var_socialinteract", wi), paste0("az_var_polparticipation", wi),
    plan_var, paste0("cloc_for", wi), paste0("stock_participation", wi)))
}

# Mega z-score wave 3
cy$az_overall_all_w3 <- andersonz(cy, c("az_overall_a_w3","az_overall_b_w3",
                                        "az_overall_c_w3","az_overall_d_w3","az_overall_e_w3"))

# ==============================================================================
# HETEROGENEITY CUTS
# ==============================================================================

hetero_above <- c("gender","birth_year","ethnicity_han","birthplace_coastal","residence_coastal",
                  "hukou_urban","religion_religious","ccp_member","university_elite","hs_track_science",
                  "department_ssh","domestic_english_atleast4","foreign_english_yes","travel_hktaiwan",
                  "travel_foreign_yes","father_edu_hsabove","work_father_govt","father_ccp",
                  "mother_edu_hsabove","work_mother_govt","mother_ccp","siblings_total",
                  "az_preference_risk","az_preference_time","az_preference_altruism","az_preference_reciprocity",
                  "az_overall_a_w1","az_overall_b_w1","az_overall_c_w1","az_overall_d_w1","az_overall_e_w1",
                  "az_belief_media_value_w1","az_belief_media_trust_w1","az_knowledge_news_censored_w1",
                  "az_knowledge_news_uncensor_w1","az_knowledge_people_censor_w1","az_knowledge_people_uncens_w1",
                  "az_knowledge_protest_china_w1","az_belief_trust_govt_w1")

for (v in hetero_above) {
  ref <- if (v == "hh_income") 75000 else if (v == "siblings_total") 0 else 0
  op  <- if (v == "siblings_total") `>` else `==`
  if (v %in% names(cy)) {
    cy[[paste0("h_", v)]] <- if_else(!is.na(cy[[v]]),
                                     if (v == "hh_income") as.integer(cy[[v]] >= 75000)
                                     else if (v == "siblings_total") as.integer(cy[[v]] > 0)
                                     else as.integer(cy[[v]] > 0), NA_integer_)
  }
}

# ==============================================================================
# FINAL SUBSET: keep wave 3, non-existing users; create _p indicators for Table 2
# ==============================================================================

cy_out <- cy |>
  filter(panelmerged_wave3 == 1, treatment_user != 1)

# Above-median dummies
above_med_vars_w3 <- c("info_foreign_website_w3","info_freq_website_for_w3",
                       "az_belief_media_value_w3","az_belief_media_trust_w3","bias_domestic_w3",
                       "az_belief_media_justif_w3","news_perccor_cen_w3","news_perccor_unc_w3",
                       "protest_pcheard_china_w3","protest_pcheard_foreign_w3","az_knowledge_meta_w3",
                       "az_belief_econ_conf_cn_w3","az_belief_econ_perf_us_w3","az_belief_econ_conf_us_w3",
                       "az_belief_instchange_w3","az_belief_trust_foreign_w3","az_belief_willing_w3",
                       "frequency_talk_politic_w3","frequency_persuade_friends_w3")
above_med_vars_w1 <- c("info_foreign_website_w1","info_freq_website_for_w1",
                       "az_belief_media_value_w1","az_belief_media_trust_w1","bias_domestic_w1",
                       "az_belief_media_justif_w1","news_perccor_cen_w1","news_perccor_unc_w1",
                       "protest_pcheard_china_w1","protest_pcheard_foreign_w1","az_knowledge_meta_w1",
                       "az_belief_econ_conf_cn_w1","az_belief_instchange_w1","az_belief_trust_foreign_w1",
                       "az_belief_willing_w1","frequency_talk_politic_w1","frequency_persuade_friends_w1")

for (v in c(above_med_vars_w3, above_med_vars_w1)) {
  med <- median(cy_out[[v]], na.rm = TRUE)
  cy_out[[paste0(v, "_p")]] <- as.integer(cy_out[[v]] > med)
}

# Below-median dummies (flipped direction)
below_med_vars_w3 <- c("bias_foreign_w3","importance_live_in_demo_w3",
                       "az_belief_econ_perf_cn_w3","az_belief_trust_govt_w3","az_belief_evalgovt_w3")
below_med_vars_w1 <- c("bias_foreign_w1","importance_live_in_demo_w1",
                       "az_belief_econ_perf_cn_w1","az_belief_trust_govt_w1","az_belief_evalgovt_w1")

for (v in c(below_med_vars_w3, below_med_vars_w1)) {
  med <- median(cy_out[[v]], na.rm = TRUE)
  cy_out[[paste0(v, "_p")]] <- as.integer(cy_out[[v]] < med)
}

# Rename for complain_school
cy_out <- cy_out |>
  rename(par_complain_school_w1 = participate_complain_school_w1,
         par_complain_school_w3 = participate_complain_school_w3)

# Already-dummy variables: keep as-is
passthrough_w3 <- c("vpn_purchase_wmt_record","vpn_purchase_yes","bias_dom_govt_policy_t1_w3",
                    "protest_2011_tmrw_parade_w3","participate_social_protest_w3","participate_plan_vote_w3",
                    "par_complain_school_w3","plan_grad_foreignmaster_w3","cp_t3_for_firm_w3","cloc_for_w3")
passthrough_w1 <- c("bias_dom_govt_policy_t1_w1","protest_2011_tmrw_parade_w1",
                    "participate_social_protest_w1","participate_plan_vote_w1","par_complain_school_w1",
                    "plan_grad_foreignmaster_w1","cp_t3_for_firm_w1","cloc_for_w1")

for (v in c(passthrough_w3, passthrough_w1)) {
  cy_out[[paste0(v, "_p")]] <- cy_out[[v]]
}

# Flipped dummies
for (v in c("bias_for_govt_policy_t1_w3","stock_participation_w3",
            "bias_for_govt_policy_t1_w1","stock_participation_w1")) {
  cy_out[[paste0(v, "_p")]] <- 1L - cy_out[[v]]
}

# ==============================================================================
# SELECT FINAL COLUMNS (mirrors the `keep` statement in the .do file)
# ==============================================================================

keep_cols <- c(
  "treatment_main","active_user","treatment_control",
  "treatment_vpnonly","treatment_nlonly","treatment_vpnnl",
  "info_foreign_website_w3_p","info_freq_website_for_w3_p",
  "az_belief_media_value_w3_p","az_belief_media_trust_w3_p",
  "bias_domestic_w3_p","bias_foreign_w3_p","az_belief_media_justif_w3_p",
  "bias_dom_govt_policy_t1_w3_p","bias_for_govt_policy_t1_w3_p",
  "info_foreign_website_w1_p","info_freq_website_for_w1_p",
  "az_belief_media_value_w1_p","az_belief_media_trust_w1_p",
  "bias_domestic_w1_p","bias_foreign_w1_p","az_belief_media_justif_w1_p",
  "bias_dom_govt_policy_t1_w1_p","bias_for_govt_policy_t1_w1_p",
  "vpn_purchase_wmt_record_p","vpn_purchase_yes_p"
)

cy_final <- cy_out |> select(all_of(keep_cols))

# ==============================================================================
# SAVE
# ==============================================================================

saveRDS(cy_final, file.path(data_path, "ChenYang2019.rds"))
# Also save as .dta for cross-format verification if needed:
# write_dta(cy_final, file.path(data_path, "ChenYang2019.dta"))

cat("ChenYang2019.rds saved:", nrow(cy_final), "obs,", ncol(cy_final), "variables\n")
