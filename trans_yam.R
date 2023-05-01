# victoria circulates mostly in 2016 and 2018

rm(list = ls())

library(tidyverse)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(lme4)
library(sjPlot)
library(jtools)
library(merTools)
library(zoo)

all <- read.csv('transmission/PHIRST flu shedding and profile 2016-2018.csv') %>%
  dplyr::select(-X) %>% 
  filter(year == 2017)

# MANIPULATING DATA ####

# household exclusions - those that include individuals that had episodes that 
# could not have shedding estimated
yam_hh_ex <- c("A119", "A136", "A114", "A141")
yam <- all %>% filter(!(hh_id %in% yam_hh_ex))

# read in hai titer data
hai_dat <- read.csv('data/PHIRST qry_flu_hai_2016-2018_2023-02-13.csv') 
# only want pre-season - first draw 
hai_dat <- hai_dat %>% 
  filter(draw == 1) %>% 
  dplyr::select(2,8:19) %>%
  rename(indid = ind_id)
# also make groups categories for titer level
hai_dat$yam_group <- NA 
hai_dat$yam_group <- ifelse(hai_dat$flub_yamagata_gm < 40, "<40", hai_dat$yam_group)
hai_dat$yam_group <- ifelse(hai_dat$flub_yamagata_gm >= 40, ">=40", hai_dat$yam_group)

yam <- yam %>% 
  left_join(hai_dat, by = "indid")

# initiate regression matrix data
indid_l <- c(unique(yam$indid))
yam_profile <- data.frame()
for (i in indid_l) {
  sub <- yam %>% filter(indid == i)
  yam_profile <- rbind(yam_profile, sub[1, ])
}
rm(i, indid_l, hai_dat, yam_hh_ex, sub)

yam_matrix <- data.frame(indid = rep(c(yam_profile$indid), each = 286), date = rep(0:285, times = 541),
                         year = rep(c(yam_profile$year), each = 286),
                         hh = rep(c(yam_profile$hh_id), each = 286), site = rep(c(yam_profile$site), each = 286), 
                         age = rep(c(yam_profile$age_at_consent), each = 286), age_cat = rep(c(yam_profile$age_cat), each = 286),
                         sex = rep(c(yam_profile$sex), each = 286), bmi = rep(c(yam_profile$bmicat), each = 286), 
                         hiv = rep(c(yam_profile$cleanhiv), each = 286), hh_s = rep(c(yam_profile$true_hh_size), each = 286),
                         hai_yam = rep(c(yam_profile$flub_yamagata_gm), each = 286), 
                         hai_cat = rep(c(yam_profile$yam_group), each = 286),
                         anysmokenow = rep(c(yam_profile$anysmokenow), each = 286),
                         hh_income = rep(c(yam_profile$hh_income), each = 286),
                         under5children = rep(c(yam_profile$under5children), each = 286),
                         crowding = rep(c(yam_profile$crowding), each = 286),
                         hh_cat = rep(c(yam_profile$hh_size), each = 286))

# adding start date and viral load information for the primary infections 
yam_clean <- yam %>% filter(TYPE == "BYam" | TYPE == "Missing" | TYPE == "Negative")
yam_one <- data.frame()
indid_inf <- unique(yam_clean$indid_inf[which(!is.na(yam_clean$indid_inf))])
for (n in indid_inf) {
  sub <- yam_clean %>% filter(indid_inf == n)
  yam_one <- rbind(yam_one, sub[1, ])
}
rm(sub, indid_inf, n, yam_clean)

yam_matrix$yam_present <- NA
yam_matrix$prolif_load <- NA
yam_matrix$clear_load <- NA
new <- data.frame()
indid_l <- c(yam_one$indid_inf)
for (i in indid_l) {
  ind <- yam_one$indid[which(yam_one$indid_inf == i)]
  sub <- yam_matrix %>% filter(indid == ind)
  s <- yam_one$round_start[which(yam_one$indid_inf == i)]
  sub$yam_present[1:(s-1)] <- 0
  # one DAY BUFFER (because s is by row and days are off by yam_one)
  sub$yam_present[s] <- 1   
  sub$yam_present[(s+1):286] <- NA
  
  sub$prolif_load <- ifelse(sub$date <= yam_one$start_point[which(yam_one$indid_inf == i)], 0, sub$prolif_load)
  sub$clear_load <- ifelse(sub$date <= yam_one$start_point[which(yam_one$indid_inf == i)], 0, sub$clear_load)
  sub$prolif_load <- ifelse(sub$date > yam_one$start_point[which(yam_one$indid_inf == i)] &
                               sub$date <= yam_one$center_point[which(yam_one$indid_inf == i)],
                             (sub$date - yam_one$start_point[which(yam_one$indid_inf == i)]) * 
                               yam_one$prolif_slope[which(yam_one$indid_inf == i)],
                             sub$prolif_load)
  sub$prolif_load <- ifelse(sub$date > yam_one$center_point[which(yam_one$indid_inf == i)] &
                             sub$date <= yam_one$end_point[which(yam_one$indid_inf == i)], 0, sub$prolif_load)
  sub$clear_load <- ifelse(sub$date > yam_one$center_point[which(yam_one$indid_inf == i)] &
                               sub$date <= yam_one$end_point[which(yam_one$indid_inf == i)],
                             yam_one$meanct[which(yam_one$indid_inf == i)] + 
                               ((sub$date - yam_one$center_point[which(yam_one$indid_inf == i)]) * 
                                  yam_one$clear_slope[which(yam_one$indid_inf == i)]),
                             sub$clear_load)
  sub$clear_load <- ifelse(sub$date > yam_one$start_point[which(yam_one$indid_inf == i)] &
                              sub$date <= yam_one$center_point[which(yam_one$indid_inf == i)], 0, sub$clear_load)
  sub$prolif_load <- ifelse(sub$date >= yam_one$end_point[which(yam_one$indid_inf == i)], 0, sub$prolif_load)
  sub$clear_load <- ifelse(sub$date >= yam_one$end_point[which(yam_one$indid_inf == i)], 0, sub$clear_load)
  new <- rbind(new, sub)
}
# no infections recorded
indid_one <- unique(c(yam_one$indid))
noinf <- setdiff(unique(c(yam_matrix$indid)), indid_one)
yam_matrix$yam_present <- ifelse(yam_matrix$indid %in% noinf, 0, yam_matrix$yam_present)
yam_matrix$prolif_load <- ifelse(yam_matrix$indid %in% noinf, 0, yam_matrix$prolif_load)
yam_matrix$clear_load <- ifelse(yam_matrix$indid %in% noinf, 0, yam_matrix$clear_load)
noinf_df <- yam_matrix %>% filter(indid %in% noinf)

yam_matrix <- rbind(new, noinf_df) %>% 
  arrange(indid)
rm(i, s, sub, ind, indid_l, yam_one, new, indid_one, noinf, noinf_df)

# household load
yam_matrix$hh_prof_load <- NA
yam_matrix$hh_clear_load <- NA
hh_l <- c(unique(yam_matrix$hh))
for (h in hh_l) {
  for (d in 0:286){
    yam_matrix$hh_prof_load[which(yam_matrix$hh == h & yam_matrix$date == d)] <- 
      sum(yam_matrix$prolif_load[which(yam_matrix$hh == h & yam_matrix$date == d)])
    yam_matrix$hh_clear_load[which(yam_matrix$hh == h & yam_matrix$date == d)] <- 
      sum(yam_matrix$clear_load[which(yam_matrix$hh == h & yam_matrix$date == d)])
  }
}
rm(h, d, hh_l)
# remove individuals own household load 
yam_matrix <- yam_matrix %>%
  mutate(hh_prof_load = hh_prof_load - prolif_load) %>%
  mutate(hh_clear_load = hh_clear_load - clear_load) %>% 
  mutate(ind_lod_tot = prolif_load + clear_load)

# add proxy for community surveillance - number of infections present at the time
# split by study site and by year
yam_matrix$community <- 0
yam_matrix$community <- ifelse(yam_matrix$ind_lod_tot > 0, yam_matrix$community - 1, yam_matrix$community)
for (s in c("Agincourt", "Klerksdorp")) {
  for (d in 0:286) { 
    yam_matrix$community[which(yam_matrix$date == d & yam_matrix$site == s)] <- 
        (yam_matrix$community[which(yam_matrix$date == d & yam_matrix$site == s)] + 
           length(which(yam_matrix$ind_lod_tot[which(yam_matrix$date == d & yam_matrix$site == s)] > 0))) / 
        length(yam_matrix$indid[which(yam_matrix$date == d & yam_matrix$site == s)])
  }
}
rm(d, s, y)
# moving window smoothing - 3 day window to span between samples 
yam_matrix$comsmooth <- rollapply(yam_matrix$community, width = 4, function(...) {round(mean(...), digits = 3)}, partial = TRUE)
# transform smooth community 
yam_matrix$comsmooth <- yam_matrix$comsmooth * 100

# other observations to exclude - those after the infection point
# also get rid of extra observations in cases of reinfections
yam_regress <- yam_matrix %>% filter(!is.na(yam_present))
indid_l <- c(unique(yam_regress$indid))
new <- data.frame()
for (i in indid_l) {
  sub <- yam_regress %>% filter(indid == i)
  if (length(which(sub$yam_present == 1)) > 1) {
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
yam_regress <- new
rm(indid_l, new, sub, i, maxd, new2, sub2, d)

# no hai data - only exclude individual
yam_regress <- yam_regress %>% filter(!is.na(hai_yam))

# remove co-infections and/or is not the first infection of the year
co_l <- c(yam$indid[which(yam$coinfect == 1 & yam$TYPE == "BYam")])
re_l <- c(yam$indid[which(yam$TYPE == "BYam" & yam$infcluster != 1)])
ex <- unique(c(co_l, re_l))
yam_regress <- yam_regress %>% filter(!(indid %in% ex))
rm(co_l, re_l, ex)

# REGRESSION #####
#write.csv(yam_regress, 'transmission/yam regress.csv')
yam_regress <- read.csv('transmission/yam regress.csv')
yam_regress$age_cat <- relevel(as.factor(yam_regress$age_cat), ref = ">40")
yam_regress <- yam_regress %>% mutate(hh_load = hh_prof_load + hh_clear_load)

# hh income/bmi not reported  
yam_regress <- yam_regress %>% filter(bmi != "")

length(which(yam_regress$yam_present == 1)) #41
length(unique(yam_regress$indid)) #432 - so 391 individuals without infection
length(unique(yam_regress$hh)) 

# main model
yam_mod <- glm(yam_present ~ comsmooth + hh_load + log(hai_yam, base = 10) + 
                 age_cat + sex + hiv + bmi + hh_cat + under5children + site,
               data = yam_regress, family = poisson())
summary(yam_mod)
exp(cbind("Odds ratio" = coef(yam_mod), confint.default(yam_mod, level = 0.95)))



# hh cluster model 
yam_modhh <- glmer(yam_present ~ comsmooth + hh_load + log(hai_yam, base = 10) + 
                 age_cat + sex + hiv + bmi + hh_cat + under5children + site + (1|hh),
               data = yam_regress, family = poisson())
summary(yam_modhh)


# age split models
yam_young <- yam_regress %>% filter(age_cat == "<5" | age_cat == "5-11" | age_cat == "12-18")
yam_modyoung <- glm(yam_present ~ comsmooth + hh_load + log(hai_yam, base = 10) + 
                 age_cat + sex + hiv + bmi + hh_cat + under5children + site,
               data = yam_young, family = poisson())
summary(yam_modyoung)
exp(cbind("Odds ratio" = coef(yam_modyoung), confint.default(yam_modyoung, level = 0.95)))

yam_old <- yam_regress %>% filter(age_cat == "19-40" | age_cat == ">40")
yam_modold <- glm(yam_present ~ comsmooth + hh_load + log(hai_yam, base = 10) + 
                      age_cat + sex + hiv + bmi + hh_cat + under5children + site,
                    data = yam_old, family = poisson())
summary(yam_modold)
exp(cbind("Odds ratio" = coef(yam_modold), confint.default(yam_modold, level = 0.95)))



# plots NEED TO ADJUST

yamplot <- plot_model(yam_mod, show.values = TRUE, value.offset = 0.35, color = "black",
                      ci_method="wald",
                      order.terms = c(1:3, 4, 7, 5, 6, 8:11, 14, 13, 12, 15:18),
                      axis.labels = c("Site Klerkdorp", "Children under 5 in Household", 
                                      "Household Size 6-10 Members", "Household Size 3-5", 
                                      "BMI Obese", "BMI Overweight", "BMI Underweight", 
                                      "HIV Status Unknown", "PLWH CD4>=200",
                                      "PLWH CD4<200", "Sex Male", "Age 19-40", "Age 12-18", 
                                      "Age 5-11", "Age <5",
                                      "Yamagata Pre-Season HAI Titer", 
                                      "Household FOI", "Community FOI")) +
  theme(panel.background = element_rect(fill = '#def1ce')) + labs(title = "B/Yamagata")
export_plot(yamplot, "trans yam", 3.5, 5.5)


yamplot <- plot_model(yam_mod, show.values = TRUE, value.offset = 0.35, color = "black",
                      ci_method="wald",
                      rm.terms = c("sex [Male]", "hiv [Positive CD4 <200]", 
                                   "hiv [Positive CD4 >=200]", "hiv [Unknown]",
                                   "site [Klerksdorp]", "bmi [Obese]", "bmi [Overweight]",
                                   "bmi [Underweight]", "hh_s"),
                      order.terms = c(1:3, 4, 7, 5, 6),
                      axis.labels = c("Age 19-40", "Age 12-18", "Age 5-11", "Age <5",
                                      "Yamagata Titer", 
                                      "Community FOI", "Household FOI")) +
  theme(panel.background = element_rect(fill = '#def1ce')) + labs(title = "B/Yamagata") 
export_plot(yamplot, "small trans yam", 2.8, 2.8)

yamplothh <- plot_model(yam_modhh, show.values = TRUE, value.offset = 0.35, color = "black",
                      ci_method="wald",
                      order.terms = c(1:3, 4, 7, 5, 6, 8:11, 14, 13, 12, 15:18),
                      axis.labels = c("Site Klerkdorp", "Children under 5 in Household", 
                                      "Household Size 6-10 Members", "Household Size 3-5", 
                                      "BMI Obese", "BMI Overweight", "BMI Underweight", 
                                      "HIV Status Unknown", "PLWH CD4>=200",
                                      "PLWH CD4<200", "Sex Male", "Age 19-40", "Age 12-18", 
                                      "Age 5-11", "Age <5",
                                      "Yamagata Pre-Season HAI Titer", 
                                      "Household FOI", "Community FOI")) +
  theme(panel.background = element_rect(fill = '#def1ce')) + labs(title = "B/Yamagata")
export_plot(yamplothh, "trans yam cluster", 3.5, 5.5)


