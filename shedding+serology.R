
rm(list = ls())

library(tidyverse)
library(dplyr)
library(lme4)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(lmtest)
library(nlme)
source('export.R')

# pre-season serology vs. shedding dynamics

# PREPPING AND VISUALIZING DATA ################################################
# read in shedding data 
shed_dat <- read.csv('transmission/PHIRST flu shedding and profile 2016-2018.csv') %>% 
  dplyr::select(-X) %>%
  mutate(tot_duration = meanwp + meanwc)

# comparing duration between the subtypes
#pos <- shed_dat %>% filter(TYPE != "Missing" & TYPE != "Negative" & !is.na(tot_duration))
#ggplot(data = pos) +
#  geom_boxplot(aes(y = tot_duration, x = TYPE, fill = TYPE)) +
#  labs(title = "Estimated Duration of Shedding", x = "Subtype", y = "Duration (Days)") + 
#  scale_x_discrete(labels = c("A(H1N1)pdm09", "A(H3N2)", "B/Victoria", "B/Yamagata")) + 
#  theme_light()+
#  scale_fill_manual(values = c(AH1 = "#482677FF", AH3 = "#2D708EFF", 
#                    BYam = "#FDE725FF", BVic = "#73D055FF"), guide = "none")

# read in hai titer data
hai_dat <- read.csv('data/qry_flu_hai_2016-2018_23Mar2023.csv') 
# only want pre-season - first draw 
hai_dat <- hai_dat %>% 
  filter(draw == 1) %>% 
  dplyr::select(2,8:19) %>%
  rename(indid = ind_id)

# also make groups categories for titer level
hai_dat$h1n1_1_group <- NA 
hai_dat$h1n1_1_group <- ifelse(hai_dat$flua_h1n1pdm_gm < 40, "<40", hai_dat$h1n1_1_group)
hai_dat$h1n1_1_group <- ifelse(hai_dat$flua_h1n1pdm_gm >= 40, ">=40", hai_dat$h1n1_1_group)

hai_dat$h3n2_1_group <- NA 
hai_dat$h3n2_1_group <- ifelse(hai_dat$flua_h3n2_gm < 40, "<40", hai_dat$h3n2_1_group)
hai_dat$h3n2_1_group <- ifelse(hai_dat$flua_h3n2_gm >= 40, ">=40", hai_dat$h3n2_1_group)

hai_dat$yam_1_group <- NA 
hai_dat$yam_1_group <- ifelse(hai_dat$flub_yamagata_gm < 40, "<40", hai_dat$yam_1_group)
hai_dat$yam_1_group <- ifelse(hai_dat$flub_yamagata_gm >= 40, ">=40", hai_dat$yam_1_group)

hai_dat$vic_1_group <- NA 
hai_dat$vic_1_group <- ifelse(hai_dat$flub_victoria_gm < 40, "<40", hai_dat$vic_1_group)
hai_dat$vic_1_group <- ifelse(hai_dat$flub_victoria_gm >= 40, ">=40", hai_dat$vic_1_group)

# combine
all_dat <- shed_dat %>% 
  left_join(hai_dat, by = "indid")

# start sample size 
# table(all_dat$TYPE)
# AH1      AH3     BVic     BYam Negative 
# 147      198      255       81

# only want first infection of season that is not a coinfection - then only need 
# one observation per infection
all_dat <- all_dat %>%
  filter(!(is.na(meanwc))) %>% # this gets rid of samples that could not have shedding estimated
  filter(!(is.na(flua_h1n1pdm_gm))) # this gets rid of based on HAI titer avail
all_dat <- distinct(all_dat, indid_inf, .keep_all = TRUE)
# table(all_dat$TYPE)
#  AH1  AH3 BVic BYam 
#  101  155  179   64 

# remove second infection and coinfections 
all_dat <- all_dat %>% 
  filter(coinfect == 0) %>% 
  filter(infcluster == 1)

# final sample size
# table(all_dat$TYPE)
#  AH1  AH3 BVic BYam 
# 80  130  155   43

# removing non-first infection (removes coinfections as well) - 21 H1, 13 H3, 24 Victoria, 20 Yamagata
# removing co-infections - 11 H3, 1 Yam (these H3 are all because H3 starts first for that 
# infection before H1)

# REGRESSION ##################################################################

# secondary model - check interaction / subset age for different titer levels
# partial effects plot --- variable against response --- use for titer 
# package for model iterative checking - mumin - function: dredge 

#AH1 = "#79b5d5", AH3 = "#196293",
#BYam = "#94d25c", BVic = "#24721f"

#library(MuMIn)
#options(na.action = "na.fail")

# h1n1 ######################################
h1n1 <- all_dat %>% filter(TYPE == "AH1") 
h1n1$age_cat <- relevel(as.factor(h1n1$age_cat), ref = ">40")
h1n1$year <- as.factor(h1n1$year)

# duration
dur_h1 <- lmer(tot_duration ~ log(flua_h1n1pdm_gm, base = 10) + age_cat + sex + 
                 cleanhiv + bmicat + site + (1|year), 
              data = h1n1) 
summary(dur_h1)

# with hh cluster
dur_h1_hh <- lmer(tot_duration ~ log(flua_h1n1pdm_gm, base = 10) + age_cat + sex + 
                 cleanhiv + bmicat + site + (1|year) + (1|hh_id), 
               data = h1n1) 
summary(dur_h1_hh)

h1durp <- plot_model(dur_h1, show.values = TRUE, value.offset = .35, colors = c("black"),
           order.terms = c(1, 2, 5, 3, 4, 6, 9, 7, 8, 12, 11, 10, 13),
            axis.labels = c("Site Klerkdorp",  "BMI Obese", "BMI Overweight", 
                "BMI Underweight",
                "PLWH CD4 >=200", "PLWH CD4 <200", "HIV Status Unknown",
                "Sex Male",
                "Age 19-40", "Age 12-18", "Age 5-11", "Age <5",
                "H1N1 Pre-Season Titer - Log Scale")) + 
  theme(panel.background = element_rect(fill = '#e4f0f6')) + labs(title = "A(H1N1)pdm09")
export_plot(h1durp, "dur plot h1", 4, 4.5)
h1durpsmall <- plot_model(dur_h1, show.values = TRUE, value.offset = .35, colors = c("black"),
                          rm.terms = c("sex [Male]", "cleanhiv [Positive CD4 <200]", 
                                       "cleanhiv [Positive CD4 >=200]", "cleanhiv [Unknown]",
                                       "site [Klerksdorp]", "bmicat [Obese]", "bmicat [Overweight]",
                                       "bmicat [Underweight]", "hh_income [R1601-3200]", 
                                       "hh_income [R3201-6400]", "hh_income [R6401-12800]",
                                       "hh_income [R801-1600]"),
                          order.terms = c(1, 2, 5, 3, 4),
                          axis.labels = c("Age 19-40", "Age 12-18", "Age 5-11", "Age <5",
                                          "H1N1 Titer (Log Scale) ")) + 
  theme(panel.background = element_rect(fill = '#e4f0f6')) + labs(title = "A(H1N1)pdm09")
export_plot(h1durpsmall, "small dur plot h1", 2.5, 2.5)
h1durphh <- plot_model(dur_h1_hh, show.values = TRUE, value.offset = .35, colors = c("black"),
                     order.terms = c(1, 2, 5, 3, 4, 6, 9, 7, 8, 12, 11, 10, 13),
                     axis.labels = c("Site Klerkdorp",  "BMI Obese", "BMI Overweight", 
                                     "BMI Underweight",
                                     "PLWH CD4 >=200", "PLWH CD4 <200", "HIV Status Unknown",
                                     "Sex Male",
                                     "Age 19-40", "Age 12-18", "Age 5-11", "Age <5",
                                     "H1N1 Pre-Season Titer - Log Scale")) + 
  theme(panel.background = element_rect(fill = '#e4f0f6')) + labs(title = "A(H1N1)pdm09")
export_plot(h1durphh, "dur plot h1 hh cluster", 4, 4.5)

# minimum Ct

#ct_h1 <- lm(meanct ~ log(flua_h1n1pdm_gm, base = 10) + age_cat + sex + 
#              cleanhiv + bmicat + site, 
#             data = h1n1)
#bptest(ct_h1)

# generalized least squares in the presence of heteroskedasticity - year as controlled effect! 
gls.h1 <- gls(meanct ~ log(flua_h1n1pdm_gm, base = 10) + age_cat + sex + cleanhiv + 
                bmicat + site + year,
              data = h1n1)
summary(gls.h1)

h1ctp <- plot_model(gls.h1,show.values = TRUE, value.offset = .35, colors = c("black"),
                    order.terms = c(1, 2, 5, 3, 4, 6, 9, 7, 8, 12, 11, 10, 13, 14, 15),
                    axis.labels = c("Year 2018", "Year 2017", "Site Klerkdorp",  "BMI Obese", "BMI Overweight", 
                                    "BMI Underweight",
                                    "PLWH CD4 >=200", "PLWH CD4 <200", "HIV Status Unknown",
                                    "Sex Male",
                                    "Age 19-40", "Age 12-18", "Age 5-11", "Age <5",
                                    "H1N1 Pre-Season Titer - Log Scale")) + 
  theme(panel.background = element_rect(fill = '#e4f0f6')) + labs(title = "A(H1N1)pdm09")
export_plot(h1ctp, "ct plot h1", 4, 4.5)


# h3n2  ###############################
h3n2 <- all_dat %>% filter(TYPE == "AH3")
h3n2$age_cat <- relevel(as.factor(h3n2$age_cat), ref = ">40")
h3n2$year <- as.factor(h3n2$year)

# total duration
dur_h3 <- lmer(tot_duration ~ log(flua_h3n2_gm, base = 10) + age_cat + sex + 
               cleanhiv + bmicat + site + (1|year), 
              data = h3n2) 
summary(dur_h3)
cbind("estimates" = coef(dur_h3), confint.default(dur_h3, level = 0.95))


dur_h3_hh <- lmer(tot_duration ~ log(flua_h3n2_gm, base = 10) + age_cat + sex + 
                    cleanhiv + bmicat + site + (1|year) + (1|hh_id), 
                  data = h3n2) 
summary(dur_h3_hh)


h3durp <- plot_model(dur_h3, show.values = TRUE, value.offset = .35, colors = c("black"),
                     order.terms = c(1, 2, 5, 3, 4, 6, 9, 7, 8, 12, 11, 10, 13),
                     axis.labels = c("Site Klerkdorp",  "BMI Obese", "BMI Overweight", 
                                     "BMI Underweight",
                                     "PLWH CD4 >=200", "PLWH CD4 <200", "HIV Status Unknown",
                                     "Sex Male",
                                     "Age 19-40", "Age 12-18", "Age 5-11", "Age <5",
                                     "H3N2 Pre-Season Titer - Log Scale")) + 
  theme(panel.background = element_rect(fill = '#bacfde')) + labs(title = "A(H3N2)")
export_plot(h3durp, "dur plot h3", 4, 4.5)
h3durp <- plot_model(dur_h3, show.values = TRUE, value.offset = .35, colors = c("black"),
                     rm.terms = c("sex [Male]", "cleanhiv [Positive CD4 <200]", 
                                  "cleanhiv [Positive CD4 >=200]", "cleanhiv [Unknown]",
                                  "site [Klerksdorp]", "bmicat [Obese]", "bmicat [Overweight]",
                                  "bmicat [Underweight]", "hh_income [R1601-3200]", 
                                  "hh_income [R3201-6400]", "hh_income [R6401-12800]",
                                  "hh_income [R801-1600]"),
                     order.terms = c(1, 2, 5, 3, 4),
                     axis.labels = c("Age 19-40", "Age 12-18", "Age 5-11", "Age <5",
                                     "H3N2 Titer")) + 
  theme(panel.background = element_rect(fill = '#bacfde')) + labs(title = "A(H3N2)")
export_plot(h3durp, "small dur plot h3", 2.5, 2.5)
h3durphh <- plot_model(dur_h3_hh, show.values = TRUE, value.offset = .35, colors = c("black"),
                     order.terms = c(1, 2, 5, 3, 4, 6, 9, 7, 8, 12, 11, 10, 13),
                     axis.labels = c("Site Klerkdorp",  "BMI Obese", "BMI Overweight", 
                                     "BMI Underweight",
                                     "PLWH CD4 >=200", "PLWH CD4 <200", "HIV Status Unknown",
                                     "Sex Male",
                                     "Age 19-40", "Age 12-18", "Age 5-11", "Age <5",
                                     "H3N2 Pre-Season Titer - Log Scale")) + 
  theme(panel.background = element_rect(fill = '#bacfde')) + labs(title = "A(H3N2)")
export_plot(h3durphh, "dur plot h3 hh cluster", 4, 4.5)


# minimum Ct
#ct_h3 <- lm(meanct ~ log(flua_h3n2_gm, base = 10) + age_cat + sex + cleanhiv + 
#              bmicat + site, 
#            data = h3n2) 
#bptest(ct_h3)

gls.h3 <- gls(meanct ~ log(flua_h3n2_gm, base = 10) + age_cat + sex + cleanhiv + 
                bmicat + site + year,
              data = h3n2)
summary(gls.h3)
cbind("estimates" = coef(gls.h3), confint.default(gls.h3, level = 0.95))


h3ctp <- plot_model(gls.h3, show.values = TRUE, value.offset = .35, colors = c("black"),
                    order.terms = c(1, 2, 5, 3, 4, 6, 9, 7, 8, 12, 11, 10, 13:15),
                    axis.labels = c("Year 2018", "Year 2017", "Site Klerkdorp",  "BMI Obese", 
                                    "BMI Overweight", "BMI Underweight",
                                    "PLWH CD4 >=200", "PLWH CD4 <200", "HIV Status Unknown",
                                    "Sex Male",
                                    "Age 19-40", "Age 12-18", "Age 5-11", "Age <5",
                                    "H3N2 Pre-Season Titer - Log Scale")) + 
  theme(panel.background = element_rect(fill = '#bacfde')) + labs(title = "A(H3N2)")
export_plot(h3ctp, "ct plot h3", 4, 4.5)


#  VICTORIA ####################################
vic <- all_dat %>% filter(TYPE == "BVic")
vic$age_cat <- relevel(as.factor(vic$age_cat), ref = ">40")
vic$year <- as.factor(vic$year)

dur_vic <- lmer(tot_duration ~ log(flub_victoria_gm, base = 10) + age_cat + sex + 
                  cleanhiv + bmicat + site + (1|year), 
               data = vic) 
summary(dur_vic)

dur_vic_hh <- lmer(tot_duration ~ log(flub_victoria_gm, base = 10) + age_cat + sex + 
                  cleanhiv + bmicat + site + (1|year) + (1|hh_id), 
                data = vic) 
summary(dur_vic_hh)

vicdurp <- plot_model(dur_vic, show.values = TRUE, value.offset = .35, colors = c("black"),
                      order.terms = c(1, 2, 5, 3, 4, 6, 9, 7, 8, 12, 11, 10, 13),
                      axis.labels = c("Site Klerkdorp",  "BMI Obese", "BMI Overweight", 
                                      "BMI Underweight",
                                      "PLWH CD4 >=200", "PLWH CD4 <200", "HIV Status Unknown",
                                      "Sex Male",
                                      "Age 19-40", "Age 12-18", "Age 5-11", "Age <5",
                                      "Victoria Pre-Season Titer - Log Scale")) + 
  theme(panel.background = element_rect(fill = '#a7c6a5')) + labs(title = "B/Victoria")
export_plot(vicdurp, "dur plot vic", 4, 4.5)
vicdurp <- plot_model(dur_vic, show.values = TRUE, value.offset = .35, colors = c("black"),
                     rm.terms = c("sex [Male]", "cleanhiv [Positive CD4 <200]", 
                                  "cleanhiv [Positive CD4 >=200]", "cleanhiv [Unknown]",
                                  "site [Klerksdorp]", "bmicat [Obese]", "bmicat [Overweight]",
                                  "bmicat [Underweight]", "hh_income [R1601-3200]", 
                                  "hh_income [R3201-6400]", "hh_income [R6401-12800]",
                                  "hh_income [R801-1600]"),
                     order.terms = c(1, 2, 5, 3, 4),
                     axis.labels = c("Age 19-40", "Age 12-18", "Age 5-11", "Age <5",
                                     "Victoria Titer")) + 
  theme(panel.background = element_rect(fill = '#a7c6a5')) + labs(title = "B/Victoria")
export_plot(vicdurp, "small dur plot vic", 2.5, 2.5)
vicdurphh <- plot_model(dur_vic_hh, show.values = TRUE, value.offset = .35, colors = c("black"),
                      order.terms = c(1, 2, 5, 3, 4, 6, 9, 7, 8, 12, 11, 10, 13),
                      axis.labels = c("Site Klerkdorp",  "BMI Obese", "BMI Overweight", 
                                      "BMI Underweight",
                                      "PLWH CD4 >=200", "PLWH CD4 <200", "HIV Status Unknown",
                                      "Sex Male",
                                      "Age 19-40", "Age 12-18", "Age 5-11", "Age <5",
                                      "Victoria Pre-Season Titer - Log Scale")) + 
  theme(panel.background = element_rect(fill = '#a7c6a5')) + labs(title = "B/Victoria")
export_plot(vicdurphh, "dur plot vic hh cluster", 4, 4.5)



# minimum Ct 
gls.v <- gls(meanct ~ log(flub_victoria_gm, base = 10) + age_cat + sex + cleanhiv + 
                bmicat + site + year,
              data = vic)
summary(gls.v)


vicctp <- plot_model(gls.v, show.values = TRUE, value.offset = .35, colors = c("black"),
                     order.terms = c(1, 2, 5, 3, 4, 6, 9, 7, 8, 12, 11, 10, 13, 14),
                     axis.labels = c("Year 2018", "Site Klerkdorp",  "BMI Obese", "BMI Overweight", 
                                     "BMI Underweight",
                                     "PLWH CD4 >=200", "PLWH CD4 <200", "HIV Status Unknown",
                                     "Sex Male",
                                     "Age 19-40", "Age 12-18", "Age 5-11", "Age <5",
                                     "Victoria Pre-Season Titer - Log Scale")) + 
  theme(panel.background = element_rect(fill = '#a7c6a5')) + labs(title = "B/Victoria")
export_plot(vicctp, "ct plot vic", 4, 4.5)


# YAMAGATA ####################
yam <- all_dat %>% filter(TYPE == "BYam")
yam$age_cat <- relevel(as.factor(yam$age_cat), ref = ">40")
yam$year <- as.factor(yam$year)

dur_yam <- lmer(tot_duration ~ log(flub_yamagata_gm, base = 10) + age_cat + sex + 
                cleanhiv + bmicat + site + (1|year), 
                data = yam) 
summary(dur_yam)

dur_yamhh <- lmer(tot_duration ~ log(flub_yamagata_gm, base = 10) + age_cat + sex + 
                  cleanhiv + bmicat + site + (1|year) + (1|hh_id), 
                data = yam) 
summary(dur_yamhh)

yamdurp <- plot_model(dur_yam, show.values = TRUE, value.offset = .35, colors = c("black"),
                      order.terms = c(1, 2, 5, 3, 4, 6, 9, 7, 8, 12, 11, 10, 13),
                      axis.labels = c("Site Klerkdorp",  "BMI Obese", "BMI Overweight", 
                                      "BMI Underweight",
                                      "PLWH CD4 >=200", "PLWH CD4 <200", "HIV Status Unknown",
                                      "Sex Male",
                                      "Age 19-40", "Age 12-18", "Age 5-11", "Age <5",
                                      "Yamgata Pre-Season Titer - Log Scale")) + 
  theme(panel.background = element_rect(fill = '#def1ce')) + labs(title = "B/Yamagata")
export_plot(yamdurp, "dur plot yam", 4, 4.5)
yamdurp <- plot_model(dur_yam, show.values = TRUE, value.offset = .35, colors = c("black"),
                      rm.terms = c("sex [Male]", "cleanhiv [Positive CD4 <200]", 
                                   "cleanhiv [Positive CD4 >=200]", "cleanhiv [Unknown]",
                                   "site [Klerksdorp]", "bmicat [Obese]", "bmicat [Overweight]",
                                   "bmicat [Underweight]", "hh_income [R1601-3200]", 
                                   "hh_income [R3201-6400]", "hh_income [R6401-12800]",
                                   "hh_income [R801-1600]"),
                      order.terms = c(1, 2, 5, 3, 4),
                      axis.labels = c("Age 19-40", "Age 12-18", "Age 5-11", "Age <5",
                                      "Yamagata Titer")) + 
  theme(panel.background = element_rect(fill = '#def1ce')) + labs(title = "B/Yamagata")
export_plot(yamdurp, "small dur plot yam", 2.5, 2.5)
yamdurphh <- plot_model(dur_yamhh, show.values = TRUE, value.offset = .35, colors = c("black"),
                      order.terms = c(1, 2, 5, 3, 4, 6, 9, 7, 8, 12, 11, 10, 13),
                      axis.labels = c("Site Klerkdorp",  "BMI Obese", "BMI Overweight", 
                                      "BMI Underweight",
                                      "PLWH CD4 >=200", "PLWH CD4 <200", "HIV Status Unknown",
                                      "Sex Male",
                                      "Age 19-40", "Age 12-18", "Age 5-11", "Age <5",
                                      "Yamgata Pre-Season Titer - Log Scale")) + 
  theme(panel.background = element_rect(fill = '#def1ce')) + labs(title = "B/Yamagata")
export_plot(yamdurphh, "dur plot yam hh cluster", 4, 4.5)


# minimum ct

gls.y <- gls(meanct ~ log(flub_yamagata_gm, base = 10) + age_cat + sex + cleanhiv + 
               bmicat + site,
             data = yam)
summary(gls.y)

yamctp <- plot_model(gls.y, show.values = TRUE, value.offset = .35, colors = c("black"),
                     order.terms = c(1, 2, 5, 3, 4, 6, 9, 7, 8, 12, 11, 10, 13),
                     axis.labels = c("Year 2017", "Site Klerkdorp",  "BMI Obese", "BMI Overweight", 
                                     "BMI Underweight",
                                     "PLWH CD4 >=200", "PLWH CD4 <200", "HIV Status Unknown",
                                     "Sex Male",
                                     "Age 19-40", "Age 12-18", "Age 5-11", "Age <5",
                                     "Yamgata Pre-Season Titer - Log Scale")) + 
  theme(panel.background = element_rect(fill = '#def1ce')) + labs(title = "B/Yamagata")
export_plot(yamctp, "ct plot yam", 4, 4.5)



## SAMPLE SIZE ANALYSIS ######################################
# if it is an issue of sample size - H3N2 should lose significance when sampled 
# down to 80 for H1N1 or 43 for yamagata

c <- data.frame()
for (i in 1:100) {
  rand_h3_80 <- h3n2[sample(nrow(h3n2), size=80), ]
  sampledown80_h3 <- lmer(tot_duration ~ log(flua_h3n2_gm, base = 10) + age_cat + sex + 
                            cleanhiv + bmicat + site + (1|year), 
                          data = rand_h3_80) 
 c <- rbind(c, confint(sampledown80_h3, "log(flua_h3n2_gm, base = 10)"))
}
c <- c %>% filter(!(is.na(`2.5 %`)))
mean(c$`2.5 %`)
mean(c$`97.5 %`)

c2 <- data.frame()
for (i in 1:100) {
  rand_h3_43 <- h3n2[sample(nrow(h3n2), size=43), ]
  sampledown43_h3 <- lm(tot_duration ~ log(flua_h3n2_gm, base = 10) + age_cat + sex + 
                            cleanhiv + bmicat + site, 
                          data = rand_h3_43) 
  c2 <- rbind(c2, confint(sampledown43_h3, "log(flua_h3n2_gm, base = 10)"))
}
c2 <- c2 %>% filter(!(is.na(`2.5 %`)))
mean(c2$`2.5 %`)
mean(c2$`97.5 %`)
