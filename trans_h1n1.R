
# H1N1 circulates mostly in 2016 and 2018

rm(list = ls())

library(tidyverse)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(lme4)
library(sjPlot)
library(jtools)
library(merTools)
source('export.R')

all <- read.csv('transmission/PHIRST flu shedding and profile 2016-2018.csv') %>%
  dplyr::select(-X) %>% 
  filter(year == 2016 | year == 2018)

# MANIPULATING DATA ####

# household exclusions - those that include individuals that had episodes that 
# could not have shedding estimated
h1_hh_ex <- c("K048", "K036", "K043", "A236", "A258", "A210", "A212", "A255", 
              "K224", "A221", "A240")
h1 <- all %>% filter(!(hh_id %in% h1_hh_ex))

# read in hai titer data
hai_dat <- read.csv('data/PHIRST qry_flu_hai_2016-2018_2023-02-13.csv') 
# only want pre-season - first draw 
hai_dat <- hai_dat %>% 
  filter(draw == 1) %>% 
  dplyr::select(2,8:19) %>%
  rename(indid = ind_id)
# also make groups categories for titer level
hai_dat$h1n1_1_group <- NA 
hai_dat$h1n1_1_group <- ifelse(hai_dat$flua_h1n1pdm_gm < 40, "<40", hai_dat$h1n1_1_group)
hai_dat$h1n1_1_group <- ifelse(hai_dat$flua_h1n1pdm_gm >= 40, ">=40", hai_dat$h1n1_1_group)

h1 <- h1 %>% 
  left_join(hai_dat, by = "indid")

# initiate regression matrix data
indid_l <- c(unique(h1$indid))
h1_profile <- data.frame()
for (i in indid_l) {
  sub <- h1 %>% filter(indid == i)
  h1_profile <- rbind(h1_profile, sub[1, ])
}
rm(i, indid_l, hai_dat, h1_hh_ex, sub)

## add anysmokenow, hh_income, under5children, crowding, hh_cat (hh_size)

h1_matrix <- data.frame(indid = rep(c(h1_profile$indid), each = 290), date = rep(0:289, times = 1035),
                        year = rep(c(h1_profile$year), each = 290),
                        hh = rep(c(h1_profile$hh_id), each = 290), site = rep(c(h1_profile$site), each = 290), 
                        age = rep(c(h1_profile$age_at_consent), each = 290), age_cat = rep(c(h1_profile$age_cat), each = 290),
                        sex = rep(c(h1_profile$sex), each = 290), bmi = rep(c(h1_profile$bmicat), each = 290), 
                        hiv = rep(c(h1_profile$HIV), each = 290), cleanhiv = rep(c(h1_profile$cleanhiv), each = 290),
                        hh_s = rep(c(h1_profile$true_hh_size), each = 290),
                        hai_h1 = rep(c(h1_profile$flua_h1n1pdm_gm), each = 290), 
                        hai_cat = rep(c(h1_profile$h1n1_1_group), each = 290),
                        anysmokenow = rep(c(h1_profile$anysmokenow), each = 290),
                        hh_income = rep(c(h1_profile$hh_income), each = 290),
                        under5children = rep(c(h1_profile$under5children), each = 290),
                        crowding = rep(c(h1_profile$crowding), each = 290),
                        hh_cat = rep(c(h1_profile$hh_size), each = 290))

# adding start date and viral load information for the primary infections 
h1_clean <- h1 %>% filter(TYPE == "AH1" | TYPE == "Missing" | TYPE == "Negative")
h1_one <- data.frame()
indid_inf <- unique(h1_clean$indid_inf[which(!is.na(h1_clean$indid_inf))])
for (n in indid_inf) {
  sub <- h1_clean %>% filter(indid_inf == n)
  h1_one <- rbind(h1_one, sub[1, ])
}
rm(sub, indid_inf, n, h1_clean)

h1_matrix$h1_present <- NA
h1_matrix$prolif_load <- NA
h1_matrix$clear_load <- NA
new <- data.frame()
indid_l <- c(h1_one$indid_inf)
for (i in indid_l) {
  ind <- h1_one$indid[which(h1_one$indid_inf == i)]
  sub <- h1_matrix %>% filter(indid == ind)
  s <- h1_one$round_start[which(h1_one$indid_inf == i)]
  sub$h1_present[1:(s-1)] <- 0
  # one DAY BUFFER (because s is by row and days are off by h1_one)
  sub$h1_present[s] <- 1   
  sub$h1_present[(s+1):290] <- NA
  
  sub$prolif_load <- ifelse(sub$date <= h1_one$start_point[which(h1_one$indid_inf == i)], 0, sub$prolif_load)
  sub$clear_load <- ifelse(sub$date <= h1_one$start_point[which(h1_one$indid_inf == i)], 0, sub$clear_load)
  sub$prolif_load <- ifelse(sub$date > h1_one$start_point[which(h1_one$indid_inf == i)] &
                                 sub$date <= h1_one$center_point[which(h1_one$indid_inf == i)],
                               (sub$date - h1_one$start_point[which(h1_one$indid_inf == i)]) * 
                                 h1_one$prolif_slope[which(h1_one$indid_inf == i)],
                               sub$prolif_load)
  sub$prolif_load <- ifelse(sub$date > h1_one$center_point[which(h1_one$indid_inf == i)] &
                              sub$date <= h1_one$end_point[which(h1_one$indid_inf == i)], 0, sub$prolif_load)
  sub$clear_load <- ifelse(sub$date > h1_one$center_point[which(h1_one$indid_inf == i)] &
                                 sub$date <= h1_one$end_point[which(h1_one$indid_inf == i)],
                               h1_one$meanct[which(h1_one$indid_inf == i)] + 
                                 ((sub$date - h1_one$center_point[which(h1_one$indid_inf == i)]) * 
                                    h1_one$clear_slope[which(h1_one$indid_inf == i)]),
                               sub$clear_load)
  sub$clear_load <- ifelse(sub$date > h1_one$start_point[which(h1_one$indid_inf == i)] &
                             sub$date <= h1_one$center_point[which(h1_one$indid_inf == i)], 0, sub$clear_load)
  sub$prolif_load <- ifelse(sub$date >= h1_one$end_point[which(h1_one$indid_inf == i)], 0, sub$prolif_load)
  sub$clear_load <- ifelse(sub$date >= h1_one$end_point[which(h1_one$indid_inf == i)], 0, sub$clear_load)
  new <- rbind(new, sub)
}
# no infections recorded
indid_one <- unique(c(h1_one$indid))
noinf <- setdiff(unique(c(h1_matrix$indid)), indid_one)
h1_matrix$h1_present <- ifelse(h1_matrix$indid %in% noinf, 0, h1_matrix$h1_present)
h1_matrix$prolif_load <- ifelse(h1_matrix$indid %in% noinf, 0, h1_matrix$prolif_load)
h1_matrix$clear_load <- ifelse(h1_matrix$indid %in% noinf, 0, h1_matrix$clear_load)
noinf_df <- h1_matrix %>% filter(indid %in% noinf)

h1_matrix <- rbind(new, noinf_df) %>% 
  arrange(indid)
rm(i, s, sub, ind, indid_l, h1_one, new, indid_one, noinf, noinf_df)

# household load
h1_matrix$hh_prof_load <- NA
h1_matrix$hh_clear_load <- NA
hh_l <- c(unique(h1_matrix$hh))
for (h in hh_l) {
  for (d in 0:289){
    h1_matrix$hh_prof_load[which(h1_matrix$hh == h & h1_matrix$date == d)] <- 
      sum(h1_matrix$prolif_load[which(h1_matrix$hh == h & h1_matrix$date == d)])
    h1_matrix$hh_clear_load[which(h1_matrix$hh == h & h1_matrix$date == d)] <- 
      sum(h1_matrix$clear_load[which(h1_matrix$hh == h & h1_matrix$date == d)])
  }
}
rm(h, d, hh_l)
# remove individuals own household load 
h1_matrix <- h1_matrix %>%
  mutate(hh_prof_load = hh_prof_load - prolif_load) %>%
  mutate(hh_clear_load = hh_clear_load - clear_load) %>% 
  mutate(ind_lod_tot = prolif_load + clear_load)

# add proxy for community surveillance - number of infections present at the time
# split by study site and by year
h1_matrix$community <- 0
h1_matrix$community <- ifelse(h1_matrix$ind_lod_tot > 0, h1_matrix$community - 1, h1_matrix$community)
for (s in c("Agincourt", "Klerksdorp")) {
  for (d in 0:289) { 
    for (y in c(2016, 2018)) {
    h1_matrix$community[which(h1_matrix$date == d & h1_matrix$site == s & h1_matrix$year == y)] <- 
      (h1_matrix$community[which(h1_matrix$date == d & h1_matrix$site == s & h1_matrix$year == y)] + 
         length(which(h1_matrix$ind_lod_tot[which(h1_matrix$date == d & h1_matrix$site == s & h1_matrix$year == y)] > 0))) / 
      length(h1_matrix$indid[which(h1_matrix$date == d & h1_matrix$site == s & h1_matrix$year == y)])
    }
  }
}
rm(d, s, y)

# other observations to exclude - those after the infection point
# also get rid of extra observations in cases of reinfections
h1_regress <- h1_matrix %>% filter(!is.na(h1_present))
indid_l <- c(unique(h1_regress$indid))
new <- data.frame()
for (i in indid_l) {
  sub <- h1_regress %>% filter(indid == i)
  if (length(which(sub$h1_present == 1)) > 1) {
    maxd <- max(sub$date)
    new2 <- data.frame()
    for(d in 0:maxd){
      sub2 <- sub %>% filter(date == d)
      new2 <- rbind(new2, sub2[1,])
    }
    sub <- new2
  }
  new <- rbind(new, sub)
}
h1_regress <- new
rm(indid_l, new, sub, i, maxd, new2, sub2, d)

# moving window smoothing - 4 day window to span between samples 
library(zoo)
h1_regress$comsmooth <- rollapply(h1_regress$community, width = 4, function(...) {round(mean(...), digits = 3)}, partial = TRUE)
h1_regress$comsmooth <- h1_regress$comsmooth * 100


# no hai data - only exclude individual
h1_regress <- h1_regress %>% filter(!is.na(hai_h1))

# now remove for those data points that shouldnt be included based on year length
# 2016 only goes
h1_regress <- h1_regress %>% filter(!(year == 2016 & date > 181))

# remove co-infections and/or is not the first infection of the year
co_l <- c(h1$indid[which(h1$coinfect == 1 & h1$TYPE == "AH1")])
re_l <- c(h1$indid[which(h1$TYPE == "AH1" & h1$infcluster != 1)])
ex <- unique(c(co_l, re_l))
h1_regress <- h1_regress %>% filter(!(indid %in% ex))
rm(co_l, re_l, ex) 

# REGRESSION ##############################

#write.csv(h1_regress, 'transmission/h1 regress.csv')
h1_regress <- read.csv('transmission/h1 regress.csv')
h1_regress$age_cat <- relevel(as.factor(h1_regress$age_cat), ref = ">40")
h1_regress$year <- as.factor(h1_regress$year)
h1_regress <- h1_regress %>% mutate(hh_load = hh_prof_load + hh_clear_load)

# bmi not reported 
h1_regress <- h1_regress %>% filter(bmi != "")

length(which(h1_regress$h1_present == 1))
length(unique(h1_regress$hh))


# main model
h1_mod <- glmer(h1_present ~ comsmooth + hh_load + log(hai_h1, base = 10) + age_cat + 
                sex + cleanhiv + bmi + hh_cat + under5children + site + (1|year),
              data = h1_regress, family = poisson())
summary(h1_mod)


# hh cluster model
h1_modhh <- glmer(h1_present ~ comsmooth + hh_load + log(hai_h1, base = 10) + age_cat + 
                  sex + cleanhiv + bmi + hh_cat + under5children + site + (1|year) +
                  (1|hh),
                data = h1_regress, family = poisson())
summary(h1_modhh)

# age split models
h1_young <- h1_regress %>% filter(age_cat == "<5" | age_cat == "5-11" | age_cat == "12-18")
h1_modyoung <- glm(h1_present ~ comsmooth + hh_load + log(hai_h1, base = 10) + age_cat + 
                     sex + cleanhiv + bmi + hh_cat + under5children + site,
                data = h1_young, family = poisson())
summary(h1_modyoung)
exp(cbind("Odds ratio" = coef(h1_modyoung), confint.default(h1_modyoung, level = 0.95)))
h1_old <- h1_regress %>% filter(age_cat == "19-40" | age_cat == ">40")
h1_modold <- glm(h1_present ~ comsmooth + hh_load + log(hai_h1, base = 10) + age_cat + 
                   sex + cleanhiv + bmi + hh_cat + under5children + site,
                     data = h1_old, family = poisson())
summary(h1_modold)
exp(cbind("Odds ratio" = coef(h1_modold), confint.default(h1_modold, level = 0.95)))



## PLOTS ## NEED TO ADJUST
h1plot <- plot_model(h1_mod, show.values = TRUE, value.offset = 0.35, color = "black",
                     ci_method="wald",
                     order.terms = c(1:3, 4, 7, 5, 6, 8:11, 14, 13, 12, 15:18),
                     axis.labels = c("Site Klerkdorp", "Children under 5 in Household", 
                                     "Household Size 6-10 Members", "Household Size 3-5", 
                                     "BMI Obese", "BMI Overweight", "BMI Underweight", 
                                     "HIV Status Unknown", "PLWH CD4>=200",
                                     "PLWH CD4<200", "Sex Male", "Age 19-40", "Age 12-18", 
                                     "Age 5-11", "Age <5",
                                     "H1N1 Pre-Season HAI Titer", 
                                      "Household FOI", "Community FOI")) +
  theme(panel.background = element_rect(fill = '#e4f0f6')) + labs(title = "A(H1N1)pdm09") 
export_plot(h1plot, "trans h1", 3.5, 5.5)

smh1plot <- plot_model(h1_mod, show.values = TRUE, value.offset = 0.35, color = "black",
                     ci_method="wald",
                     rm.terms = c("sex [Male]", "cleanhiv [Positive CD4 <200]", 
                                  "cleanhiv [Positive CD4 >=200]", "cleanhiv [Unknown]",
                                  "site [Klerksdorp]", "bmi [Obese]", "bmi [Overweight]",
                                  "bmi [Underweight]", "hh_s"),
                     order.terms = c(1:3, 4, 7, 5, 6),
                     axis.labels = c("Age 19-40", "Age 12-18", "Age 5-11", "Age <5",
                                     "H1N1 Titer", 
                                     "Community FOI", "Household FOI")) +
  theme(panel.background = element_rect(fill = '#e4f0f6')) + labs(title = "A(H1N1)pdm09") 
export_plot(smh1plot, "small trans h1", 2.8, 2.8)

h1plothh <- plot_model(h1_modhh, show.values = TRUE, value.offset = 0.35, color = "black",
                       ci_method="wald",
                       order.terms = c(1:3, 4, 7, 5, 6, 8:11, 14, 13, 12, 15:18),
                       axis.labels = c("Site Klerkdorp", "Children under 5 in Household", 
                                       "Household Size 6-10 Members", "Household Size 3-5", 
                                       "BMI Obese", "BMI Overweight", "BMI Underweight", 
                                       "HIV Status Unknown", "PLWH CD4>=200",
                                       "PLWH CD4<200", "Sex Male", "Age 19-40", "Age 12-18", 
                                       "Age 5-11", "Age <5",
                                       "H1N1 Pre-Season HAI Titer", 
                                       "Household FOI", "Community FOI")) +
  theme(panel.background = element_rect(fill = '#e4f0f6')) + labs(title = "A(H1N1)pdm09") 
export_plot(h1plothh, "trans h1 cluster", 3.5, 5.5)



