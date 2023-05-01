
# CREATING MOSAIC PLOT/HEAT MAP FOR RAW DATA ###################

library(tidyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(stringr)
library(quantreg)
library(car)
library(lmtest)
source('export.R')

rm(list = ls())

shed_dat <- read.csv('transmission/PHIRST flu shedding and profile 2016-2018.csv') %>% 
  dplyr::select(-X) %>%
  mutate(tot_duration = meanwp + meanwc)

# 2016
d2016 <- shed_dat %>% filter(year == 2016)
mosaic16 <- ggplot(d2016, aes(x = funum, y = indid, fill = TYPE, alpha = CT)) +
  geom_tile() +
  theme_minimal() +
  labs(x = "Collection Number", y = "Individual", 
       title = "2016") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_alpha(range = c(1, .2),limits=c(11, 37)) +
  scale_x_continuous(expand=c(0,0), n.breaks = 52, labels = c("", "", "2016-05-02", rep("", 49), "2016-10-30")) +
  scale_fill_manual(name = element_blank(), labels = c("A(H1N1)pdm09", "A(H3N2)", 
                                                       "B/Yamagata", "B/Victoria", "Negative", "Missing"),
                    values=c(AH1 = "#79b5d5", AH3 = "#196293",
                             BYam = "#94d25c", BVic = "#24721f", Negative = "#F0F0F0", Missing = "#FFFFFF")) 
#print(mosaic16)
export_plot(mosaic16, "mosaic16", 7, 5)

# mosaic plot for Ct values - 2017
d2017 <- shed_dat %>% filter(year == 2017)
mosaic17 <- ggplot(d2017, aes(x = funum, y = indid, fill = TYPE, alpha = CT)) +
  geom_tile() +
  theme_minimal() +
  labs(x = "Collection Number", y = "Individual", 
       title = "2017") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_alpha(range = c(1, .2),limits=c(11, 37)) +
  scale_x_continuous(expand=c(0,0), n.breaks = 82, labels = c("", "", "2017-01-16", rep("", 79), "2017-10-28")) +
  scale_fill_manual(name = element_blank(), labels = c("A(H1N1)pdm09", "A(H3N2)", 
                                                       "B/Yamagata", "B/Victoria", "Negative", "Missing"),
                    values=c(AH1 = "#79b5d5", AH3 = "#196293",
                             BYam = "#94d25c", BVic = "#24721f", Negative = "#F0F0F0", Missing = "#FFFFFF")) 
#print(mosaic17)
export_plot(mosaic17, "mosaic17tall", 5, 6)

# 2018 
d2018 <- shed_dat %>% filter(year == 2018)
mosaic18 <- ggplot(d2018, aes(x = funum, y = indid, fill = TYPE, alpha = CT)) +
  geom_tile() +
  theme_minimal() +
  labs(x = "Collection Number", y = "Individual", 
       title = "2018") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_alpha(range = c(1, .2),limits=c(11, 37)) +
  scale_x_continuous(expand=c(0,0), n.breaks = 83, labels = c("", "", "2018-01-15", rep("", 80), "2018-10-31")) +
  scale_fill_manual(name = element_blank(), labels = c("A(H1N1)pdm09", "A(H3N2)", 
                                                       "B/Yamagata", "B/Victoria", "Negative", "Missing"),
                    values=c(AH1 = "#79b5d5", AH3 = "#196293",
                             BYam = "#94d25c", BVic = "#24721f", Negative = "#F0F0F0", Missing = "#FFFFFF")) 
#print(mosaic18)
export_plot(mosaic18, "mosaic18", 7, 5)

shed_dat <- read.csv('transmission/PHIRST flu shedding and profile 2016-2018.csv') %>% 
  dplyr::select(-X) %>%
  mutate(tot_duration = meanwp + meanwc)

# ATTACK RATE FIGURE AND INDIVIDUAL CHARACTERISTICS ########################

# AGE AND SUBTYPE/LINEAGE
indid_inf <- unique(shed_dat$indid_inf)
indid_inf <- indid_inf[which(!is.na(indid_inf))]
types <- c()
ages <- c()

for (n in indid_inf) {
  sub <- shed_dat %>%
    filter(indid_inf == n)
  types <- c(types, sub$TYPE[1])
  ages <- c(ages, sub$age_cat[1])
}

type_age <- data.frame(ages, types)
df <- data.frame(xtabs(~types + ages, data = type_age))
df$Freq <- ifelse(df$ages == "<5", df$Freq / length(unique(shed_dat$indid[which(shed_dat$age_at_consent < 5)])), df$Freq)
df$Freq <- ifelse(df$ages == "5-11", df$Freq / length(unique(shed_dat$indid[which(shed_dat$age_at_consent >= 5 & shed_dat$age_at_consent <=11)])), df$Freq)
df$Freq <- ifelse(df$ages == "12-18", df$Freq / length(unique(shed_dat$indid[which(shed_dat$age_at_consent >= 12 & shed_dat$age_at_consent <=18)])), df$Freq)
df$Freq <- ifelse(df$ages == "19-40", df$Freq / length(unique(shed_dat$indid[which(shed_dat$age_at_consent >= 19 & shed_dat$age_at_consent <=40)])), df$Freq)
df$Freq <- ifelse(df$ages == ">40", df$Freq / length(unique(shed_dat$indid[which(shed_dat$age_at_consent > 40)])), df$Freq)

df$types <- factor(df$types, levels = c("AH1", "AH3", "BYam", "BVic"), order = T)

arate <- ggplot(data = df, aes(x = ages, y = Freq, fill = types)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(labels = c("A(H1N1)pdm09", "A(H3N2)", "B/Yamagata", "B/Victoria"), 
                    values = c(AH1 = "#79b5d5", AH3 = "#196293",
                               BYam = "#94d25c", BVic = "#24721f"),
                    name = "Subtype/Lineage") +
  scale_x_discrete(limits = c("<5", "5-11", "12-18", "19-40", ">40")) +
  theme_bw() +
  labs(x = "Age", y = "Attack Rate 
       (# of Episodes/Age Group Size)")
export_plot(arate, "attackrate", 5, 3)

# individual characteristics 
# included are those that fulfilled enrollment criteria and had all relevant data included
# this excludes 2018 without titer
hai_dat <- read.csv('data/PHIRST qry_flu_hai_2016-2018_2023-02-13.csv') 
hai_dat <- hai_dat %>% 
  filter(draw == 1) %>% 
  dplyr::select(2,8:19) %>%
  rename(indid = ind_id)
shed_dat <- read.csv('transmission/PHIRST flu shedding and profile 2016-2018.csv') %>%  
  left_join(hai_dat, by = "indid")
shed_dat <- distinct(shed_dat, indid, .keep_all = TRUE)

shed_dat <- shed_dat %>% filter(!is.na(flua_h1n1pdm_gm))
shed_dat <- shed_dat %>% filter(bmicat != "")

hhdf <- distinct(shed_dat, hh_id, .keep_all = TRUE) 


### TITER DISTRIBUTION PLOT  #####################
shed_dat <- read.csv('transmission/PHIRST flu shedding and profile 2016-2018.csv') %>% 
  dplyr::select(-X)

hai_dat <- read.csv('data/qry_flu_hai_2016-2018_23Mar2023.csv') %>% 
  filter(draw == 1) %>% 
  dplyr::select(2,8:19) %>%
  rename(indid = ind_id)

all_dat <- shed_dat %>% 
  left_join(hai_dat, by = "indid") 
all_dat <- distinct(all_dat, indid, .keep_all = TRUE)

all_dat <- all_dat %>% 
  filter(!(is.na(flua_h1n1pdm_gm)) & !(is.na(age_at_consent)))

h1_plot <- ggplot(data = all_dat, aes(x = age_at_consent, y = log(flua_h1n1pdm_gm))) +
  geom_point() + 
  geom_smooth(color = "white", fill = "black") + 
  theme(panel.background = element_rect(fill = '#e4f0f6')) + 
  labs(x = "Age at Enrollment", y = "Pre-season A(H1N1)pdm09 HAI Titer", title = "A(H1N1)pdm09")
export_plot(h1_plot, "titer age h1", 4.5, 3)

h3_plot <- ggplot(data = all_dat, aes(x = age_at_consent, y = log(flua_h3n2_gm))) +
  geom_point() + 
  geom_smooth(color = "white", fill = "black") + 
  theme(panel.background = element_rect(fill = '#bacfde')) + 
  labs(x = "Age at Enrollment", y = "Pre-season A(H3N2) HAI Titer", title = "A(H3N2)")
export_plot(h3_plot, "titer age h3", 4.5, 3)

yam_plot <- ggplot(data = all_dat, aes(x = age_at_consent, y = log(flub_yamagata_gm))) +
  geom_point() + 
  geom_smooth(color = "white", fill = "black") + 
  theme(panel.background = element_rect(fill = '#def1ce')) + 
  labs(x = "Age at Enrollment", y = "Pre-season B/Yamagata HAI Titer", title = "B/Yamagata")
export_plot(yam_plot, "titer age yam", 4.5, 3)

vic_plot <- ggplot(data = all_dat, aes(x = age_at_consent, y = log(flub_victoria_gm))) +
  geom_point() + 
  geom_smooth(color = "white", fill = "black") + 
  theme(panel.background = element_rect(fill = '#a7c6a5')) + 
  labs(x = "Age at Enrollment", y = "Pre-season B/Victoria HAI Titer", title = "B/Victoria")
export_plot(vic_plot, "titer age vic", 4.5, 3)


# SHEDDING FIGURE BOXPLOTS AND BETWEEN STRAIN ANALYSIS


all <- read.csv('transmission/PHIRST flu shedding and profile 2016-2018.csv')

individual <- distinct(all, indid_inf, .keep_all = TRUE)
individual <- individual %>% filter(!is.na(meanwc))

# take one point at 37 with the -(meanwp) and then the other point at 37 with (meanwc) 
# and a point at 0 for meanct

start <- individual %>% 
  mutate(date = -meanwp) %>%
  dplyr::select(date, TYPE)
start$ct <- 37 
mid <- individual %>% 
  mutate(ct = meanct) %>%
  dplyr::select(ct, TYPE)
mid$date <- 0 
end <- individual %>% 
  mutate(date = meanwc) %>%
  dplyr::select(date, TYPE)
end$ct <- 37 

fig_df <- rbind(start, mid, end)

shedfig <- ggplot(data = fig_df, aes(x = date, y = ct, color = TYPE, fill = TYPE)) +
  stat_smooth(data = subset(fig_df, date <= 0), method = "lm") +
  stat_smooth(data = subset(fig_df, date >= 0), method = "lm") +
  scale_y_reverse() +
  coord_cartesian(ylim = c(37,20), xlim = c(-5,12)) +
  theme_minimal() +
  scale_color_manual(labels = c("A(H1N1)pdm09", "A(H3N2)", "B/Yamagata", "B/Victoria"), 
                     values = c(AH1 = "#79b5d5", AH3 = "#196293",
                                BYam = "#94d25c", BVic = "#24721f")) +
  scale_fill_manual(labels = c("A(H1N1)pdm09", "A(H3N2)", "B/Yamagata", "B/Victoria"), 
                    values = c(AH1 = "#A6CEE3", AH3 = "#1F78B4",
                               BYam = "#B2DF8A", BVic = "#33A02C")) +
  xlab("Time since minimum Ct (days)") + ylab("Cycle threshold (Ct)")
export_plot(shedfig, "shed char", 5, 3) 


# boxplots ##########


pbox <- ggplot(data = individual, aes(x = TYPE, y = meanwp, fill = TYPE)) + 
  geom_boxplot() + 
  scale_fill_manual(values =c(AH1 = "#79b5d5", AH3 = "#196293",
                              BYam = "#94d25c", BVic = "#24721f")) +
  scale_x_discrete(limits = c("AH1", "AH3", "BYam", "BVic")) + 
  theme_bw() + 
  labs(title = "Proliferation Duration", y = "Proliferation Duration (Days)", 
       x = NULL) +
  guides(fill="none")
export_plot(pbox, "box prolif", 3.5, 3.5)

cbox <- ggplot(data = individual, aes(x = TYPE, y = meanwc, fill = TYPE)) + 
  geom_boxplot() + 
  scale_fill_manual(values =c(AH1 = "#79b5d5", AH3 = "#196293",
                              BYam = "#94d25c", BVic = "#24721f")) +
  scale_x_discrete(limits = c("AH1", "AH3", "BYam", "BVic")) + 
  theme_bw() + 
  labs(title = "Clearance Duration", y = "Clearance Duration (Days)", 
       x = NULL)+
  guides(fill="none")
export_plot(cbox, "box clear", 3.5, 3.5)

tbox <- ggplot(data = individual, aes(x = TYPE, y = tot_duration, fill = TYPE)) + 
  geom_boxplot() + 
  scale_fill_manual(values =c(AH1 = "#79b5d5", AH3 = "#196293",
                              BYam = "#94d25c", BVic = "#24721f")) +
  scale_x_discrete(limits = c("AH1", "AH3", "BYam", "BVic")) + 
  theme_bw() + 
  labs(title = "Total Duration", y = "Total Duration (Days)", 
       x = NULL)+
  guides(fill="none")
export_plot(tbox, "box total", 3.5, 3.5)

ctbox <- ggplot(data = individual, aes(x = TYPE, y = meanct, fill = TYPE)) + 
  geom_boxplot() + 
  scale_fill_manual(values =c(AH1 = "#79b5d5", AH3 = "#196293",
                              BYam = "#94d25c", BVic = "#24721f")) +
  scale_x_discrete(limits = c("AH1", "AH3", "BYam", "BVic")) + 
  theme_bw() + scale_y_reverse() +
  labs(title = "Minimum Ct (Peak Viral Load)", y = "Minimum Ct", 
       x = NULL)+
  guides(fill="none")
export_plot(ctbox, "box ct", 3.5, 3.5)



# regression to compare proliferation/clearance/total durations and then Ct 

individual$age_cat <- factor(individual$age_cat)
individual$sex <- factor(individual$sex)
individual$cleanhiv <- factor(individual$cleanhiv)
individual$bmicat <- factor(individual$bmicat)

# prolif 
summary(individual$meanwp)
summary(individual$meanwp[which(individual$TYPE == "AH1")])
summary(individual$meanwp[which(individual$TYPE == "AH3")])
summary(individual$meanwp[which(individual$TYPE == "BYam")])
summary(individual$meanwp[which(individual$TYPE == "BVic")])
#prolif <- lm(meanwp ~ TYPE + age_cat + sex + cleanhiv + bmicat, data = individual)
#summary(prolif)
#bptest(prolif)
#shapiro.test(rstandard(prolif))

# correction for weight of age cat and TYPE
library(nlme)
weights = varIdent(form= ~1|age_cat * TYPE)
gls.wp <- gls(meanwp ~ TYPE + age_cat + sex + cleanhiv + bmicat, weight = weights, data = individual)
summary(gls.wp)

prolifplot <- plot_model(gls.wp, show.values = TRUE, value.offset = .35, colors = c("black"),
                         vline.color = "red", 
                         rm.terms = c("sex [Male]", "cleanhiv [Positive CD4 <200]", 
                                      "cleanhiv [Positive CD4 >=200]", "cleanhiv [Unknown]",
                                      "site [Klerksdorp]", "bmicat [Obese]", "bmicat [Overweight]",
                                      "bmicat [Underweight]", "age_cat [5-11]", "age_cat [19-40]",
                                      "age_cat [12-18]", "age_cat [>40]"),
                         axis.labels = c("B/Yamagata", "B/Victoria", "A(H3N2)")) + 
  theme_bw() +
  labs(title = "Proliferation Duration")
#export_plot(prolifplot, "strain compare prolif", 3, 3)



# clear 
summary(individual$meanwc)
summary(individual$meanwc[which(individual$TYPE == "AH1")])
summary(individual$meanwc[which(individual$TYPE == "AH3")])
summary(individual$meanwc[which(individual$TYPE == "BYam")])
summary(individual$meanwc[which(individual$TYPE == "BVic")])
#clear <- lm(meanwc ~ TYPE + age_cat + sex + cleanhiv + bmicat, data = individual) 
#bptest(clear)
#shapiro.test(rstandard(clear))
#summary(clear)

weights = varIdent(form= ~1|age_cat * TYPE)
gls.wc <- gls(meanwc ~ TYPE + age_cat + sex + cleanhiv + bmicat,weights = weights, data = individual)
summary(gls.wc)

clearplot <- plot_model(gls.wc, show.values = TRUE, value.offset = .35, colors = c("black"),
                        vline.color = "red", 
                        rm.terms = c("sex [Male]", "cleanhiv [Positive CD4 <200]", 
                                     "cleanhiv [Positive CD4 >=200]", "cleanhiv [Unknown]",
                                     "site [Klerksdorp]", "bmicat [Obese]", "bmicat [Overweight]",
                                     "bmicat [Underweight]", "age_cat [5-11]", "age_cat [19-40]",
                                     "age_cat [12-18]", "age_cat [>40]"),
                        axis.labels = c("B/Yamagata", "B/Victoria", "A(H3N2)")) + 
  theme_bw() +
  labs(title = "Clearance Duration")
#export_plot(clearplot, "strain compare clear", 3, 3)

# total 
summary(individual$tot_duration)
summary(individual$tot_duration[which(individual$TYPE == "AH1")])
summary(individual$tot_duration[which(individual$TYPE == "AH3")])
summary(individual$tot_duration[which(individual$TYPE == "BYam")])
summary(individual$tot_duration[which(individual$TYPE == "BVic")])
mean_se(individual$tot_duration[which(individual$TYPE == "BVic")])
#tot <- lm(tot_duration ~ TYPE + age_cat + sex + cleanhiv + bmicat, data = individual) 
#bptest(tot)
#shapiro.test(rstandard(tot))
#summary(tot)

weights = varIdent(form= ~1|age_cat * TYPE)
gls.tot <- gls(tot_duration ~ TYPE + age_cat + sex + cleanhiv + bmicat, weights = weights, data = individual)
summary(gls.tot)

totplot <- plot_model(gls.tot, show.values = TRUE, value.offset = .35, colors = c("black"),
                      vline.color = "red", 
                      rm.terms = c("sex [Male]", "cleanhiv [Positive CD4 <200]", 
                                   "cleanhiv [Positive CD4 >=200]", "cleanhiv [Unknown]",
                                   "site [Klerksdorp]", "bmicat [Obese]", "bmicat [Overweight]",
                                   "bmicat [Underweight]", "age_cat [5-11]", "age_cat [19-40]",
                                   "age_cat [12-18]", "age_cat [>40]"),
                      axis.labels = c("B/Yamagata", "B/Victoria", "A(H3N2)")) + 
  theme_bw() +
  labs(title = "Total Duration")
#export_plot(totplot, "strain compare total", 3, 3)



# ct 
summary(individual$meanct)
summary(individual$meanct[which(individual$TYPE == "AH1")])
summary(individual$meanct[which(individual$TYPE == "AH3")])
summary(individual$meanct[which(individual$TYPE == "BYam")])
summary(individual$meanct[which(individual$TYPE == "BVic")])

weights = varIdent(form= ~1|age_cat * TYPE)
gls.ct <- gls(meanct ~ TYPE + age_cat + sex + cleanhiv + bmicat, weights = weights, data = individual)
summary(gls.ct)

ctplot <- plot_model(gls.ct, show.values = TRUE, value.offset = .35, colors = c("black"),
                     vline.color = "red", 
                     rm.terms = c("sex [Male]", "cleanhiv [Positive CD4 <200]", 
                                  "cleanhiv [Positive CD4 >=200]", "cleanhiv [Unknown]",
                                  "site [Klerksdorp]", "bmicat [Obese]", "bmicat [Overweight]",
                                  "bmicat [Underweight]", "age_cat [5-11]", "age_cat [19-40]",
                                  "age_cat [12-18]", "age_cat [>40]"),
                     axis.labels = c("B/Yamagata", "B/Victoria", "A(H3N2)")) + 
  theme_bw() +
  labs(title = "Minimum Ct (Peak Viral Load)")
#export_plot(totplot, "strain compare ct", 3, 3)

