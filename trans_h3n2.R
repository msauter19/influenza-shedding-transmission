rm(list = ls())

library(tidyverse)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(lme4)
library(sjPlot)
library(jtools)
library(merTools)

all <- read.csv('transmission/PHIRST flu shedding and profile 2016-2018.csv') %>%
  dplyr::select(-X) %>% 
  filter(year == 2017)

h3_hh_ex <- c("K139", "A112", "K133")
h3 <- all %>% filter(!(hh_id %in% h3_hh_ex))
rm(h3_hh_ex)

# read in hai titer data
hai_dat <- read.csv('data/PHIRST qry_flu_hai_2016-2018_2023-02-13.csv') 
# only want pre-season - first draw 
hai_dat <- hai_dat %>% 
  filter(draw == 1) %>% 
  dplyr::select(2,8:19) %>%
  rename(indid = ind_id)
# also make groups categories for titer level
hai_dat$h3_group <- NA 
hai_dat$h3_group <- ifelse(hai_dat$flua_h3n2_gm < 40, "<40", hai_dat$h3_group)
hai_dat$h3_group <- ifelse(hai_dat$flua_h3n2_gm >= 40, ">=40", hai_dat$h3_group)

h3 <- h3 %>% 
  left_join(hai_dat, by = "indid")

# initiate regression matrix data
indid_l <- c(unique(h3$indid))
h3_profile <- data.frame()
for (i in indid_l) {
  sub <- h3 %>% filter(indid == i)
  h3_profile <- rbind(h3_profile, sub[1, ])
}
rm(i, indid_l, hai_dat, sub)

h3_matrix <- data.frame(indid = rep(c(h3_profile$indid), each = 286), date = rep(0:285, times = 539),
                         year = rep(c(h3_profile$year), each = 286),
                         hh = rep(c(h3_profile$hh_id), each = 286), site = rep(c(h3_profile$site), each = 286), 
                         age = rep(c(h3_profile$age_at_consent), each = 286), age_cat = rep(c(h3_profile$age_cat), each = 286),
                         sex = rep(c(h3_profile$sex), each = 286), bmi = rep(c(h3_profile$bmicat), each = 286), 
                         hiv = rep(c(h3_profile$cleanhiv), each = 286), hh_s = rep(c(h3_profile$true_hh_size), each = 286),
                         hai_h3 = rep(c(h3_profile$flua_h3n2_gm), each = 286), 
                         hai_cat = rep(c(h3_profile$h3_group), each = 286),
                        anysmokenow = rep(c(h3_profile$anysmokenow), each = 286),
                        hh_income = rep(c(h3_profile$hh_income), each = 286),
                        under5children = rep(c(h3_profile$under5children), each = 286),
                        crowding = rep(c(h3_profile$crowding), each = 286),
                        hh_cat = rep(c(h3_profile$hh_size), each = 286))

# adding start date and viral load information for the primary infections 
h3_clean <- h3 %>% filter(TYPE == "AH3" | TYPE == "Missing" | TYPE == "Negative")
h3_one <- data.frame()
indid_inf <- unique(h3_clean$indid_inf[which(!is.na(h3_clean$indid_inf))])
for (n in indid_inf) {
  sub <- h3_clean %>% filter(indid_inf == n)
  h3_one <- rbind(h3_one, sub[1, ])
}
rm(sub, indid_inf, n, h3_clean)

h3_matrix$h3_present <- NA
h3_matrix$prolif_load <- NA
h3_matrix$clear_load <- NA
new <- data.frame()
indid_l <- c(h3_one$indid_inf)
for (i in indid_l) {
  ind <- h3_one$indid[which(h3_one$indid_inf == i)]
  sub <- h3_matrix %>% filter(indid == ind)
  s <- h3_one$round_start[which(h3_one$indid_inf == i)]
  sub$h3_present[1:(s-1)] <- 0
  # one DAY BUFFER (because s is by row and days are off by h3_one)
  sub$h3_present[s] <- 1   
  sub$h3_present[(s+1):286] <- NA
  
  sub$prolif_load <- ifelse(sub$date <= h3_one$start_point[which(h3_one$indid_inf == i)], 0, sub$prolif_load)
  sub$clear_load <- ifelse(sub$date <= h3_one$start_point[which(h3_one$indid_inf == i)], 0, sub$clear_load)
  sub$prolif_load <- ifelse(sub$date > h3_one$start_point[which(h3_one$indid_inf == i)] &
                               sub$date <= h3_one$center_point[which(h3_one$indid_inf == i)],
                             (sub$date - h3_one$start_point[which(h3_one$indid_inf == i)]) * 
                               h3_one$prolif_slope[which(h3_one$indid_inf == i)],
                             sub$prolif_load)
  sub$clear_load <- ifelse(sub$date > h3_one$start_point[which(h3_one$indid_inf == i)] &
                              sub$date <= h3_one$center_point[which(h3_one$indid_inf == i)],0, sub$clear_load)
  sub$clear_load <- ifelse(sub$date > h3_one$center_point[which(h3_one$indid_inf == i)] &
                               sub$date <= h3_one$end_point[which(h3_one$indid_inf == i)],
                             h3_one$meanct[which(h3_one$indid_inf == i)] + 
                               ((sub$date - h3_one$center_point[which(h3_one$indid_inf == i)]) * 
                                  h3_one$clear_slope[which(h3_one$indid_inf == i)]),
                             sub$clear_load)
  sub$prolif_load <- ifelse(sub$date > h3_one$center_point[which(h3_one$indid_inf == i)] &
                             sub$date <= h3_one$end_point[which(h3_one$indid_inf == i)], 0, sub$prolif_load)
  sub$prolif_load <- ifelse(sub$date >= h3_one$end_point[which(h3_one$indid_inf == i)], 0, sub$prolif_load)
  sub$clear_load <- ifelse(sub$date >= h3_one$end_point[which(h3_one$indid_inf == i)], 0, sub$clear_load)
  new <- rbind(new, sub)
}

# no infections recorded
indid_one <- unique(c(h3_one$indid))
noinf <- setdiff(unique(c(h3_matrix$indid)), indid_one)
h3_matrix$h3_present <- ifelse(h3_matrix$indid %in% noinf, 0, h3_matrix$h3_present)
h3_matrix$prolif_load <- ifelse(h3_matrix$indid %in% noinf, 0, h3_matrix$prolif_load)
h3_matrix$clear_load <- ifelse(h3_matrix$indid %in% noinf, 0, h3_matrix$clear_load)
noinf_df <- h3_matrix %>% filter(indid %in% noinf)

h3_matrix <- rbind(new, noinf_df) %>% 
  arrange(indid)
rm(i, s, sub, ind, indid_l, h3_one, new, indid_one, noinf, noinf_df)

# household load
h3_matrix$hh_prof_load <- NA
h3_matrix$hh_clear_load <- NA
hh_l <- c(unique(h3_matrix$hh))
for (h in hh_l) {
  for (d in 0:286){
    h3_matrix$hh_prof_load[which(h3_matrix$hh == h & h3_matrix$date == d)] <- 
      sum(h3_matrix$prolif_load[which(h3_matrix$hh == h & h3_matrix$date == d)])
    h3_matrix$hh_clear_load[which(h3_matrix$hh == h & h3_matrix$date == d)] <- 
      sum(h3_matrix$clear_load[which(h3_matrix$hh == h & h3_matrix$date == d)])
  }
}
rm(h, d, hh_l)
# remove individuals own household load 
h3_matrix <- h3_matrix %>%
  mutate(hh_prof_load = hh_prof_load - prolif_load) %>%
  mutate(hh_clear_load = hh_clear_load - clear_load) %>% 
  mutate(ind_lod_tot = prolif_load + clear_load)

# add proxy for community surveillance - number of infections present at the time
# split by study site and by year
h3_matrix$community <- 0
h3_matrix$community <- ifelse(h3_matrix$ind_lod_tot > 0, h3_matrix$community - 1, h3_matrix$community)
for (s in c("Agincourt", "Klerksdorp")) {
  for (d in 0:286) { 
    h3_matrix$community[which(h3_matrix$date == d & h3_matrix$site == s)] <- 
        (h3_matrix$community[which(h3_matrix$date == d & h3_matrix$site == s)] + 
           length(which(h3_matrix$ind_lod_tot[which(h3_matrix$date == d & h3_matrix$site == s)] > 0))) / 
        length(h3_matrix$indid[which(h3_matrix$date == d & h3_matrix$site == s)])
  }
}
rm(d, s)
# moving window smoothing - 3 day window to span between samples 
h3_matrix$comsmooth <- rollapply(h3_matrix$community, width = 4, function(...) {round(mean(...), digits = 3)}, partial = TRUE)

# transform smooth community 
h3_matrix$comsmooth <- h3_matrix$comsmooth * 100

# other observations to exclude - those after the infection point
# also get rid of extra observations in cases of reinfections
h3_regress <- h3_matrix %>% filter(!is.na(h3_present))
indid_l <- c(unique(h3_regress$indid))
new <- data.frame()
for (i in indid_l) {
  sub <- h3_regress %>% filter(indid == i)
  if (length(which(sub$h3_present == 1)) > 1) {
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
h3_regress <- new
rm(indid_l, new, sub, i, maxd, new2, sub2, d)

# no hai data - only exclude individual
h3_regress <- h3_regress %>% filter(!is.na(hai_h3))

# remove co-infections and/or is not the first infection of the year
co_l <- c(h3$indid[which(h3$coinfect == 1 & h3$TYPE == "AH3")])
re_l <- c(h3$indid[which(h3$TYPE == "AH3" & h3$infcluster != 1)])
ex <- unique(c(co_l, re_l))
h3_regress <- h3_regress %>% filter(!(indid %in% ex))
rm(co_l, re_l, ex)

# REGRESSION #####
#write.csv(h3_regress, 'transmission/h3 regress.csv')
h3_regress <- read.csv('transmission/h3 regress.csv')
h3_regress$age_cat <- relevel(as.factor(h3_regress$age_cat), ref = ">40")
h3_regress <- h3_regress %>% mutate(hh_load = hh_prof_load + hh_clear_load)

# hh income/bmi not reported 
h3_regress <- h3_regress %>% filter(bmi != "")

length(which(h3_regress$h3_present == 1))
length(unique(h3_regress$hh))

# main model
h3_mod <- glm(h3_present ~ comsmooth + hh_load + log(hai_h3, base = 10) + 
                age_cat + sex + hiv + bmi + hh_cat + under5children + site,
               data = h3_regress, family = poisson())
summary(h3_mod)

# hh cluster model
h3_modhh <- glmer(h3_present ~ comsmooth + hh_load + log(hai_h3, base = 10) + 
                age_cat + sex + hiv + bmi + hh_cat + under5children + site + 
                  (1|hh),
              data = h3_regress, family = poisson())
summary(h3_modhh)

# age split models
h3_young <- h3_regress %>% filter(age_cat == "<5" | age_cat == "5-11" | age_cat == "12-18")
h3_modyoung <- glm(h3_present ~ comsmooth + hh_load + log(hai_h3, base = 10) + 
                age_cat + sex + hiv + bmi + hh_cat + under5children + site,
              data = h3_young, family = poisson())
summary(h3_modyoung)
exp(cbind("Odds ratio" = coef(h3_modyoung), confint.default(h3_modyoung, level = 0.95)))

h3_old <- h3_regress %>% filter(age_cat == "19-40" | age_cat == ">40")
h3_modold <- glm(h3_present ~ comsmooth + hh_load + log(hai_h3, base = 10) + 
                     age_cat + sex + hiv + bmi + hh_cat + under5children + site,
                   data = h3_old, family = poisson())
summary(h3_modold)
exp(cbind("Odds ratio" = coef(h3_modold), confint.default(h3_modold, level = 0.95)))



## plots ## NEED TO ADJUST 

h3plot <- plot_model(h3_mod, show.values = TRUE, value.offset = 0.35, color = "black",
                     ci_method="wald",
                     order.terms = c(1:3, 4, 7, 5, 6, 8:11, 14, 13, 12, 15:18),
                     axis.labels = c("Site Klerkdorp", "Children under 5 in Household", 
                                     "Household Size 6-10 Members", "Household Size 3-5", 
                                     "BMI Obese", "BMI Overweight", "BMI Underweight", 
                                     "HIV Status Unknown", "PLWH CD4>=200",
                                     "PLWH CD4<200", "Sex Male", "Age 19-40", "Age 12-18", 
                                     "Age 5-11", "Age <5",
                                     "H3N2 Pre-Season HAI Titer", 
                                     "Household FOI", "Community FOI")) +
  theme(panel.background = element_rect(fill = '#bacfde')) + labs(title = "A(H3N2)") 
export_plot(h3plot, "trans h3", 3.5, 5.5)


h3plot <- plot_model(h3_mod, show.values = TRUE, value.offset = 0.35, color = "black",
                     ci_method="wald",
                     rm.terms = c("sex [Male]", "hiv [Positive CD4 <200]", 
                                  "hiv [Positive CD4 >=200]", "hiv [Unknown]",
                                  "site [Klerksdorp]", "bmi [Obese]", "bmi [Overweight]",
                                  "bmi [Underweight]", "hh_s"),
                     order.terms = c(1:3, 4, 7, 5, 6),
                     axis.labels = c("Age 19-40", "Age 12-18", "Age 5-11", "Age <5",
                                     "H3N2 Titer", 
                                     "Community FOI", "Household FOI")) +
  theme(panel.background = element_rect(fill = '#bacfde')) + labs(title = "A(H3N2)") 
export_plot(h3plot, "small trans h3", 2.8, 2.8)

h3plothh <- plot_model(h3_modhh, show.values = TRUE, value.offset = 0.35, color = "black",
                     ci_method="wald",
                     order.terms = c(1:3, 4, 7, 5, 6, 8:11, 14, 13, 12, 15:18),
                     axis.labels = c("Site Klerkdorp", "Children under 5 in Household", 
                                     "Household Size 6-10 Members", "Household Size 3-5", 
                                     "BMI Obese", "BMI Overweight", "BMI Underweight", 
                                     "HIV Status Unknown", "PLWH CD4>=200",
                                     "PLWH CD4<200", "Sex Male", "Age 19-40", "Age 12-18", 
                                     "Age 5-11", "Age <5",
                                     "H3N2 Pre-Season HAI Titer", 
                                     "Household FOI", "Community FOI")) +
  theme(panel.background = element_rect(fill = '#bacfde')) + labs(title = "A(H3N2)") 
export_plot(h3plothh, "trans h3 cluster", 3.5, 5.5)




