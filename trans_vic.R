
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
  filter(year == 2016 | year == 2018)

# MANIPULATING DATA ####

# household exclusions - those that include individuals that had episodes that 
# could not have shedding estimated
vic_hh_ex <- c("A025", "A026", "K011", "K040", "K042", "A002", "A016", "A027", 
               "A037", "A063", "K007", "K027", "A257", "A253", "A275", "A201", 
               "A206", "A217", "A236", "A241", "A248", "A250", "A251", "A252", 
               "A261", "A264", "A266")
vic <- all %>% filter(!(hh_id %in% vic_hh_ex))

# read in hai titer data
hai_dat <- read.csv('data/PHIRST qry_flu_hai_2016-2018_2023-02-13.csv') 
# only want pre-season - first draw 
hai_dat <- hai_dat %>% 
  filter(draw == 1) %>% 
  dplyr::select(2,8:19) %>%
  rename(indid = ind_id)
# also make groups categories for titer level
hai_dat$vic_group <- NA 
hai_dat$vic_group <- ifelse(hai_dat$flub_victoria_gm < 40, "<40", hai_dat$vic_group)
hai_dat$vic_group <- ifelse(hai_dat$flub_victoria_gm >= 40, ">=40", hai_dat$vic_group)

vic <- vic %>% 
  left_join(hai_dat, by = "indid")

# initiate regression matrix data
indid_l <- c(unique(vic$indid))
vic_profile <- data.frame()
for (i in indid_l) {
  sub <- vic %>% filter(indid == i)
  vic_profile <- rbind(vic_profile, sub[1, ])
}
rm(i, indid_l, hai_dat, vic_hh_ex, sub)

vic_matrix <- data.frame(indid = rep(c(vic_profile$indid), each = 290), date = rep(0:289, times = 962),
                        year = rep(c(vic_profile$year), each = 290),
                        hh = rep(c(vic_profile$hh_id), each = 290), site = rep(c(vic_profile$site), each = 290), 
                        age = rep(c(vic_profile$age_at_consent), each = 290), age_cat = rep(c(vic_profile$age_cat), each = 290),
                        sex = rep(c(vic_profile$sex), each = 290), bmi = rep(c(vic_profile$bmicat), each = 290), 
                        hiv = rep(c(vic_profile$cleanhiv), each = 290), hh_s = rep(c(vic_profile$true_hh_size), each = 290),
                        hai_vic = rep(c(vic_profile$flub_victoria_gm), each = 290), 
                        hai_cat = rep(c(vic_profile$vic_group), each = 290),
                        anysmokenow = rep(c(vic_profile$anysmokenow), each = 290),
                        hh_income = rep(c(vic_profile$hh_income), each = 290),
                        under5children = rep(c(vic_profile$under5children), each = 290),
                        crowding = rep(c(vic_profile$crowding), each = 290),
                        hh_cat = rep(c(vic_profile$hh_size), each = 290))

# adding start date and viral load information for the primary infections 
vic_clean <- vic %>% filter(TYPE == "BVic" | TYPE == "Missing" | TYPE == "Negative")
vic_one <- data.frame()
indid_inf <- unique(vic_clean$indid_inf[which(!is.na(vic_clean$indid_inf))])
for (n in indid_inf) {
  sub <- vic_clean %>% filter(indid_inf == n)
  vic_one <- rbind(vic_one, sub[1, ])
}
rm(sub, indid_inf, n, vic_clean)

vic_matrix$vic_present <- NA
vic_matrix$prolif_load <- NA
vic_matrix$clear_load <- NA
new <- data.frame()
indid_l <- c(vic_one$indid_inf)
for (i in indid_l) {
  ind <- vic_one$indid[which(vic_one$indid_inf == i)]
  sub <- vic_matrix %>% filter(indid == ind)
  s <- vic_one$round_start[which(vic_one$indid_inf == i)]
  sub$vic_present[1:(s-1)] <- 0
  # one DAY BUFFER (because s is by row and days are off by vic_one)
  sub$vic_present[s] <- 1   
  sub$vic_present[(s+1):290] <- NA
  
  sub$prolif_load <- ifelse(sub$date <= vic_one$start_point[which(vic_one$indid_inf == i)], 0, sub$prolif_load)
  sub$clear_load <- ifelse(sub$date <= vic_one$start_point[which(vic_one$indid_inf == i)], 0, sub$clear_load)
  sub$prolif_load <- ifelse(sub$date > vic_one$start_point[which(vic_one$indid_inf == i)] &
                               sub$date <= vic_one$center_point[which(vic_one$indid_inf == i)],
                             (sub$date - vic_one$start_point[which(vic_one$indid_inf == i)]) * 
                               vic_one$prolif_slope[which(vic_one$indid_inf == i)],
                             sub$prolif_load)
  sub$prolif_load <- ifelse(sub$date > vic_one$center_point[which(vic_one$indid_inf == i)] &
                             sub$date <= vic_one$end_point[which(vic_one$indid_inf == i)], 0, sub$prolif_load)
  sub$clear_load <- ifelse(sub$date > vic_one$center_point[which(vic_one$indid_inf == i)] &
                               sub$date <= vic_one$end_point[which(vic_one$indid_inf == i)],
                             vic_one$meanct[which(vic_one$indid_inf == i)] + 
                               ((sub$date - vic_one$center_point[which(vic_one$indid_inf == i)]) * 
                                  vic_one$clear_slope[which(vic_one$indid_inf == i)]),
                             sub$clear_load)
  sub$clear_load <- ifelse(sub$date > vic_one$start_point[which(vic_one$indid_inf == i)] &
                              sub$date <= vic_one$center_point[which(vic_one$indid_inf == i)], 0, sub$clear_load)
  sub$prolif_load <- ifelse(sub$date >= vic_one$end_point[which(vic_one$indid_inf == i)], 0, sub$prolif_load)
  sub$clear_load <- ifelse(sub$date >= vic_one$end_point[which(vic_one$indid_inf == i)], 0, sub$clear_load)
  new <- rbind(new, sub)
}
# no infections recorded
indid_one <- unique(c(vic_one$indid))
noinf <- setdiff(unique(c(vic_matrix$indid)), indid_one)
vic_matrix$vic_present <- ifelse(vic_matrix$indid %in% noinf, 0, vic_matrix$vic_present)
vic_matrix$prolif_load <- ifelse(vic_matrix$indid %in% noinf, 0, vic_matrix$prolif_load)
vic_matrix$clear_load <- ifelse(vic_matrix$indid %in% noinf, 0, vic_matrix$clear_load)
noinf_df <- vic_matrix %>% filter(indid %in% noinf)

vic_matrix <- rbind(new, noinf_df) %>% 
  arrange(indid)
rm(i, s, sub, ind, indid_l, vic_one, new, indid_one, noinf, noinf_df)

# household load
vic_matrix$hh_prof_load <- NA
vic_matrix$hh_clear_load <- NA
hh_l <- c(unique(vic_matrix$hh))
for (h in hh_l) {
  for (d in 0:289){
    vic_matrix$hh_prof_load[which(vic_matrix$hh == h & vic_matrix$date == d)] <- 
      sum(vic_matrix$prolif_load[which(vic_matrix$hh == h & vic_matrix$date == d)])
    vic_matrix$hh_clear_load[which(vic_matrix$hh == h & vic_matrix$date == d)] <- 
      sum(vic_matrix$clear_load[which(vic_matrix$hh == h & vic_matrix$date == d)])
  }
}
rm(h, d, hh_l)
# remove individuals own household load 
vic_matrix <- vic_matrix %>%
  mutate(hh_prof_load = hh_prof_load - prolif_load) %>%
  mutate(hh_clear_load = hh_clear_load - clear_load) %>% 
  mutate(ind_lod_tot = prolif_load + clear_load)

# add proxy for community surveillance - number of infections present at the time
# split by study site and by year
vic_matrix$community <- 0
vic_matrix$community <- ifelse(vic_matrix$ind_lod_tot > 0, vic_matrix$community - 1, vic_matrix$community)
for (s in c("Agincourt", "Klerksdorp")) {
  for (d in 0:289) { 
    for (y in c(2016, 2018)) {
      vic_matrix$community[which(vic_matrix$date == d & vic_matrix$site == s & vic_matrix$year == y)] <- 
        (vic_matrix$community[which(vic_matrix$date == d & vic_matrix$site == s & vic_matrix$year == y)] + 
           length(which(vic_matrix$ind_lod_tot[which(vic_matrix$date == d & vic_matrix$site == s & vic_matrix$year == y)] > 0))) / 
        length(vic_matrix$indid[which(vic_matrix$date == d & vic_matrix$site == s & vic_matrix$year == y)])
    }
  }
}
rm(d, s, y)
# moving window smoothing - 3 day window to span between samples 
vic_matrix$comsmooth <- rollapply(vic_matrix$community, width = 4, function(...) {round(mean(...), digits = 3)}, partial = TRUE)

# other observations to exclude - those after the infection point
# also get rid of extra observations in cases of reinfections
vic_regress <- vic_matrix %>% filter(!is.na(vic_present))
indid_l <- c(unique(vic_regress$indid))
new <- data.frame()
for (i in indid_l) {
  sub <- vic_regress %>% filter(indid == i)
  if (length(which(sub$vic_present == 1)) > 1) {
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
vic_regress <- new
rm(indid_l, new, sub, i, maxd, new2, sub2, d)

# no hai data - only exclude individual
vic_regress <- vic_regress %>% filter(!is.na(hai_vic))

# now remove for those data points that shouldnt be included based on year length
# 2016 only goes 
vic_regress <- vic_regress %>% filter(!(year == 2016 & date > 181))

# remove co-infections and/or is not the first infection of the year
co_l <- c(vic$indid[which(vic$coinfect == 1 & vic$TYPE == "BVic")])
re_l <- c(vic$indid[which(vic$TYPE == "BVic" & vic$infcluster != 1)])
ex <- unique(c(co_l, re_l))
vic_regress <- vic_regress %>% filter(!(indid %in% ex))
rm(co_l, re_l, ex)

# transform smooth community 
vic_regress$comsmooth <- vic_regress$comsmooth * 100

# REGRESSION #####
#write.csv(vic_regress, 'transmission/vic regress.csv')
vic_regress <- read.csv('transmission/vic regress.csv')
vic_regress$age_cat <- relevel(as.factor(vic_regress$age_cat), ref = ">40")
vic_regress$year <- as.factor(vic_regress$year)
vic_regress <- vic_regress %>% mutate(hh_load = hh_prof_load + hh_clear_load)

# hh income/bmi not reported
vic_regress <- vic_regress %>% filter(bmi != "")

length(which(vic_regress$vic_present == 1))
length(unique(vic_regress$hh))

vic_mod <- glmer(vic_present ~ comsmooth + hh_load + log(hai_vic, base = 10) +
                age_cat  + sex + hiv + bmi + hh_cat + under5children + site + (1|year),
              data = vic_regress, family = poisson())
summary(vic_mod)

# hh cluster model
vic_modhh <- glmer(vic_present ~ comsmooth + hh_load + log(hai_vic, base = 10) +
                   age_cat  + sex + hiv + bmi + hh_cat + under5children + site + 
                   (1|year) + (1|hh),
                 data = vic_regress, family = poisson())
summary(vic_modhh)


# age split models
vic_young <- vic_regress %>% filter(age_cat == "<5" | age_cat == "5-11" | age_cat == "12-18")
vic_modyoung <- glm(vic_present ~ comsmooth + hh_load + log(hai_vic, base = 10) +
                        age_cat  + sex + hiv + bmi + hh_cat + under5children + site,
                   data = vic_young, family = poisson())
summary(vic_modyoung)
exp(cbind("Odds ratio" = coef(vic_modyoung), confint.default(vic_modyoung, level = 0.95)))

vic_old <- vic_regress %>% filter(age_cat == "19-40" | age_cat == ">40")
vic_modold <- glm(vic_present ~ comsmooth + hh_load + log(hai_vic, base = 10) +
                      age_cat  + sex + hiv + bmi + hh_cat + under5children + site,
                 data = vic_old, family = poisson())
summary(vic_modold)
exp(cbind("Odds ratio" = coef(vic_modold), confint.default(vic_modold, level = 0.95)))


# plots
vicplot <- plot_model(vic_mod, show.values = TRUE, value.offset = 0.35, color = "black",
                      ci_method="wald",
                      order.terms = c(1:3, 4, 7, 5, 6, 8:11, 14, 13, 12, 15:18),
                      axis.labels = c("Site Klerkdorp", "Children under 5 in Household", 
                                      "Household Size 6-10 Members", "Household Size 3-5", 
                                      "BMI Obese", "BMI Overweight", "BMI Underweight", 
                                      "HIV Status Unknown", "PLWH CD4>=200",
                                      "PLWH CD4<200", "Sex Male", "Age 19-40", "Age 12-18", 
                                      "Age 5-11", "Age <5",
                                      "Victoria Pre-Season HAI Titer", 
                                      "Household FOI", "Community FOI")) +
  theme(panel.background = element_rect(fill = '#a7c6a5')) + labs(title = "B/Victoria") 
export_plot(vicplot, "trans vic", 3.5, 5.5)


vicplot <- plot_model(vic_mod, show.values = TRUE, value.offset = 0.35, color = "black",
                     ci_method="wald",
                     rm.terms = c("sex [Male]", "hiv [Positive CD4 <200]", 
                                  "hiv [Positive CD4 >=200]", "hiv [Unknown]",
                                  "site [Klerksdorp]", "bmi [Obese]", "bmi [Overweight]",
                                  "bmi [Underweight]", "hh_s"),
                     order.terms = c(1:3, 4, 7, 5, 6),
                     axis.labels = c("Age 19-40", "Age 12-18", "Age 5-11", "Age <5",
                                     "Victoria Titer", 
                                     "Community FOI", "Household FOI")) +
  theme(panel.background = element_rect(fill = '#a7c6a5')) + labs(title = "B/Victoria") 
export_plot(vicplot, "small trans vic", 2.8, 2.8)

vicplothh <- plot_model(vic_modhh, show.values = TRUE, value.offset = 0.35, color = "black",
                      ci_method="wald",
                      order.terms = c(1:3, 4, 7, 5, 6, 8:11, 14, 13, 12, 15:18),
                      axis.labels = c("Site Klerkdorp", "Children under 5 in Household", 
                                      "Household Size 6-10 Members", "Household Size 3-5", 
                                      "BMI Obese", "BMI Overweight", "BMI Underweight", 
                                      "HIV Status Unknown", "PLWH CD4>=200",
                                      "PLWH CD4<200", "Sex Male", "Age 19-40", "Age 12-18", 
                                      "Age 5-11", "Age <5",
                                      "Victoria Pre-Season HAI Titer", 
                                      "Household FOI", "Community FOI")) +
  theme(panel.background = element_rect(fill = '#a7c6a5')) + labs(title = "B/Victoria") 
export_plot(vicplothh, "trans vic cluster", 3.5, 5.5)




### SAMPLE SIZE ANALYSIS AGAINST YAMAGATA #####################################

# NEED 41 INFECTIONS AND 432 TOTAL INDIVIDUALS - OR 391 OTHER INDIVIDUALS W/O INFECTION

winfect <- c(unique(vic_regress$indid[which(vic_regress$vic_present == 1)]))
noinfect <- setdiff(c(unique(vic_regress$indid)), winfect)

c <- data.frame()
d <- NA
for (i in 1:100) {
  randinf <- sample(winfect, size=41)
  randno <- sample(noinfect, size=391)
  randdf <- vic_regress %>% filter(indid %in% randinf | indid %in% randno)
  sample_vic <- glm(vic_present ~ comsmooth + hh_load + log(hai_vic, base = 10) +
                        age_cat  + sex + hiv + bmi + hh_cat + under5children + site,
                      data = randdf, family = poisson())
  d <- c(d, coef(sample_vic)['log(hai_vic, base = 10)'])
  c <- rbind(c, confint.default(sample_vic, "log(hai_vic, base = 10)"))
}
c <- c %>% filter(!(is.na(`2.5 %`)))
exp(mean(c$`2.5 %`))
exp(mean(c$`97.5 %`))




