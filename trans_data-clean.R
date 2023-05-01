#
# Cleaning data for the transmission model
#

rm(list = ls())

library(tidyverse)
library(ggplot2)


raw_flu <- read.csv('data/PHIRST Nasopharyngeal specimens 2016-2018 Flu 2022-05-03.csv')

# CLEAN VARIABLE FOR FLU TYPE/SUBTYPE ##########################################
raw_flu$TYPE <- NA

# Types
raw_flu$TYPE <- ifelse(raw_flu$npsinfa == "1" & raw_flu$npsh3 == "1" 
                       & raw_flu$npsinfb == "0" & raw_flu$npsh1 == "0", "AH3", raw_flu$TYPE)
raw_flu$TYPE <- ifelse(raw_flu$npsinfa == "1" & raw_flu$npsh1 == "1" 
                       & raw_flu$npsinfb == "0" & raw_flu$npsh3 == "0", "AH1", raw_flu$TYPE)
raw_flu$TYPE <- ifelse(raw_flu$npsinfb == "1" & raw_flu$npsyam == "1" 
                       & raw_flu$npsinfa == "0" & raw_flu$npsvic == "0", "BYam", raw_flu$TYPE)
raw_flu$TYPE <- ifelse(raw_flu$npsinfb == "1" & raw_flu$npsvic == "1" 
                       & raw_flu$npsinfa == "0" & raw_flu$npsyam == "0", "BVic", raw_flu$TYPE)
raw_flu$TYPE <- ifelse((raw_flu$npsinfa == "1" & raw_flu$npsh1 == "0" &  raw_flu$npsh3 == "0") &
                         raw_flu$npsinfb == "0" 
                       | (raw_flu$npsinfa == "1" & is.na(raw_flu$npsh1) & is.na(raw_flu$npsh3)), 
                       "ANoSub", raw_flu$TYPE)
raw_flu$TYPE <- ifelse((raw_flu$npsinfb == "1" & raw_flu$npsinfa == "0" & raw_flu$npsyam == "0" &  raw_flu$npsvic == "0") 
                       | (raw_flu$npsinfb == "1" & raw_flu$npsinfa == "0" & is.na(raw_flu$npsyam) & is.na(raw_flu$npsvic)), 
                       "BNoSub", raw_flu$TYPE)

# Variable to denote unsubtyped sample
raw_flu[ ,'unsub_sam'] <- NA 
raw_flu$unsub_sam <- ifelse(raw_flu$TYPE == "ANoSub", 1, raw_flu$unsub_sam)
raw_flu$unsub_sam <- ifelse(raw_flu$TYPE == "BNoSub", 1, raw_flu$unsub_sam)

# Check ANoSub
raw_flu$TYPE <- ifelse(raw_flu$TYPE == "ANoSub" & lag(raw_flu$TYPE, n=1) != "ANoSub" 
                       & (lag(raw_flu$npsinfb, n=1) != 1) & !is.na(lag(raw_flu$TYPE, n=1)), 
                       lag(raw_flu$TYPE, n=1), raw_flu$TYPE) 
raw_flu$TYPE <- ifelse(raw_flu$TYPE == "ANoSub" & lag(raw_flu$TYPE, n=2) != "ANoSub" 
                       & (lag(raw_flu$npsinfb, n=2) != 1) & !is.na(lag(raw_flu$TYPE, n=2)), 
                       lag(raw_flu$TYPE, n=2), raw_flu$TYPE) 
raw_flu$TYPE <- ifelse(raw_flu$TYPE == "ANoSub" & lead(raw_flu$TYPE, n=1) != "ANoSub" 
                       & (lead(raw_flu$npsinfb, n=1) != 1) & !is.na(lead(raw_flu$TYPE, n=1)), 
                       lead(raw_flu$TYPE, n=1), raw_flu$TYPE) 
raw_flu$TYPE <- ifelse(raw_flu$TYPE == "ANoSub" & lead(raw_flu$TYPE, n=2) != "ANoSub" 
                       & (lead(raw_flu$npsinfb, n=2) != 1) & !is.na(lead(raw_flu$TYPE, n=2)), 
                       lead(raw_flu$TYPE, n=2), raw_flu$TYPE) 

# Check BNoSub
raw_flu$TYPE <- ifelse(raw_flu$TYPE == "BNoSub" & lag(raw_flu$TYPE, n=1) != "BNoSub" 
                       & (lag(raw_flu$npsinfa, n=1) != 1) & !is.na(lag(raw_flu$TYPE, n=1)), 
                       lag(raw_flu$TYPE, n=1), raw_flu$TYPE) 
raw_flu$TYPE <- ifelse(raw_flu$TYPE == "BNoSub" & lag(raw_flu$TYPE, n=2) != "BNoSub" 
                       & (lag(raw_flu$npsinfa, n=2) != 1) & !is.na(lag(raw_flu$TYPE, n=2)), 
                       lag(raw_flu$TYPE, n=2), raw_flu$TYPE) 
raw_flu$TYPE <- ifelse(raw_flu$TYPE == "BNoSub" & lead(raw_flu$TYPE, n=1) != "BNoSub" 
                       & (lead(raw_flu$npsinfa, n=1) != 1) & !is.na(lead(raw_flu$TYPE, n=1)), 
                       lead(raw_flu$TYPE, n=1), raw_flu$TYPE) 
raw_flu$TYPE <- ifelse(raw_flu$TYPE == "BNoSub" & lead(raw_flu$TYPE, n=2) != "BNoSub" 
                       & (lead(raw_flu$npsinfa, n=2) != 1) & !is.na(lead(raw_flu$TYPE, n=2)), 
                       lead(raw_flu$TYPE, n=2), raw_flu$TYPE) 

# Coinfection with AH3BVic -- only 1, it occurs in Agincourt in 2016
raw_flu$TYPE <- ifelse(raw_flu$npsinfa == "1" & raw_flu$npsh3 == "1" & raw_flu$npsinfb == "1" 
                       & (raw_flu$npsyam == "0" | is.na(raw_flu$npsyam)), "AH3BVic", raw_flu$TYPE)
# Coinfection with AH3BYam - 2, occur in Klerksdorp, in 2017
raw_flu$TYPE <- ifelse(raw_flu$npsinfa == "1" & raw_flu$npsh3 == "1" & raw_flu$npsinfb == "1" 
                       & raw_flu$npsyam == "1", "AH3BYam", raw_flu$TYPE)
# Coinfection with AH3H1 - 22, occurred in both 2016 and 2017, both sites
raw_flu$TYPE <- ifelse(raw_flu$npsinfa == "1" & raw_flu$npsh3 == "1" & raw_flu$npsh1 == "1" 
                       & !is.na(raw_flu$npsh3) & !is.na(raw_flu$npsh1), 
                       "AH3H1", raw_flu$TYPE)
# Coinfection with BYamVic - 2, occurred in 2016, in both sites
raw_flu$TYPE <- ifelse(raw_flu$npsinfb == "1" & raw_flu$npsyam == "1" & raw_flu$npsvic == "1"
                       & !is.na(raw_flu$npsyam) & !is.na(raw_flu$npsvic), 
                       "BYamVic", raw_flu$TYPE)

raw_flu$TYPE <- ifelse(raw_flu$flu == "Neg", "Negative", raw_flu$TYPE)
raw_flu$TYPE <- ifelse(raw_flu$flu == "", "Missing", raw_flu$TYPE)


# ASSIGNING SUBTYPES BASED ON DOMINANT CIRCULATING #############################

# Here assign further based on dominate circulating strain
raw_flu$npsdatecol <- as.Date(raw_flu$npsdatecol, format = "%d/%m/%Y")
raw_flu$year <- format(as.Date(raw_flu$npsdatecol, format = "%d/%m/%Y"), "%Y")


raw_flu$TYPE <- ifelse(raw_flu$TYPE == "ANoSub" & raw_flu$year == 2017, "AH3", raw_flu$TYPE)
raw_flu$TYPE <- ifelse(raw_flu$TYPE == "ANoSub" & raw_flu$year == 2018, "AH1", raw_flu$TYPE)
raw_flu$TYPE <- ifelse(raw_flu$TYPE == "BNoSub" & raw_flu$year == 2016, "BVic", raw_flu$TYPE)
raw_flu$TYPE <- ifelse(raw_flu$TYPE == "BNoSub" & raw_flu$year == 2017, "BYam", raw_flu$TYPE)
raw_flu$TYPE <- ifelse(raw_flu$TYPE == "BNoSub" & raw_flu$year == 2018, "BVic", raw_flu$TYPE)

# removing the rest of ANoSub for now (will need to exclude these households later on in 2016)
raw_flu <- raw_flu %>%
  filter(TYPE != "ANoSub")

# table(raw_flu$TYPE)
# how many exclusions: 14 samples (14 episodes)

# COINFECTIONS  ###########

co_dat <- raw_flu %>% 
  filter(TYPE == "AH3H1" | TYPE == "AH3BVic" | TYPE == "BYamVic" | TYPE == "AH3BYam")
for (i in 1:nrow(co_dat)) {
  co_dat <- co_dat %>%
    rbind(co_dat[i,])
}
co_dat <- co_dat %>%
  arrange(fuid) %>%
  # also going to add identifier to the new fuid
  mutate(fuid_co = fuid)

# go by each type of coinfection
h3h1 <- co_dat %>%
  filter(TYPE == "AH3H1")
h3h1$TYPE <- ifelse(h3h1$fuid_co == lead(h3h1$fuid_co), "AH3", h3h1$TYPE) 
h3h1$TYPE <- ifelse(h3h1$fuid_co == lag(h3h1$fuid_co) & !is.na(lag(h3h1$fuid_co)), "AH1", h3h1$TYPE)
h3h1$fuid_co <- ifelse(h3h1$fuid_co == lag(h3h1$fuid_co) & !is.na(lag(h3h1$fuid_co)), paste(h3h1$fuid_co, "-1"), h3h1$fuid_co)

byamvic <- co_dat %>%
  filter(TYPE == "BYamVic")
byamvic$TYPE <- ifelse(byamvic$fuid_co == lead(byamvic$fuid_co), "BYam", byamvic$TYPE)
byamvic$TYPE <- ifelse(byamvic$fuid_co == lag(byamvic$fuid_co) & !is.na(lag(byamvic$fuid_co)), "BVic", byamvic$TYPE)
byamvic$fuid_co <- ifelse(byamvic$fuid_co == lag(byamvic$fuid_co) & !is.na(lag(byamvic$fuid_co)), paste(byamvic$fuid_co, "-1"), byamvic$fuid_co)

h3bns <- co_dat %>%
  filter(TYPE == "AH3BVic")
h3bns$TYPE <- ifelse(h3bns$fuid_co == lead(h3bns$fuid_co), "AH3", h3bns$TYPE)
h3bns$TYPE <- ifelse(h3bns$fuid_co == lag(h3bns$fuid_co) & !is.na(lag(h3bns$fuid_co)), "BVic", h3bns$TYPE)
h3bns$fuid_co <- ifelse(h3bns$fuid_co == lag(h3bns$fuid_co) & !is.na(lag(h3bns$fuid_co)), paste(h3bns$fuid_co, "-1"), h3bns$fuid_co)

h3byam <- co_dat %>%
  filter(TYPE == "AH3BYam")
h3byam$TYPE <- ifelse(h3byam$fuid_co == lead(h3byam$fuid_co), "AH3", h3byam$TYPE)
h3byam$TYPE <- ifelse(h3byam$fuid_co == lag(h3byam$fuid_co) & !is.na(lag(h3byam$fuid_co)), "BYam", h3byam$TYPE)
h3byam$fuid_co <- ifelse(h3byam$fuid_co == lag(h3byam$fuid_co) & !is.na(lag(h3byam$fuid_co)), paste(h3byam$fuid_co, "-1"), h3byam$fuid_co)

coinf_refine <- h3h1 %>%
  rbind(byamvic, h3bns, h3byam) %>%
  # create variable to note the co-infection is occurring
  mutate(coinfect = 1)

# combine with other data to make refined data 
non_coinf <- raw_flu %>%
  filter(TYPE == "Missing" | TYPE == "Negative" | TYPE == "AH1" | TYPE == "AH3" |TYPE == "BYam" | TYPE == "BVic" | TYPE == "ANoSub" | TYPE == "BNoSub") %>%
  mutate(coinfect = 0,fuid_co = fuid) 

raw_flu <- non_coinf %>%
  rbind(coinf_refine) %>%
  arrange(fuid)

rm(byamvic, co_dat, coinf_refine, h3bns, h3byam, h3h1, non_coinf, i)

# CLUSTER SAMPLES FOR INFECTION ################################################

raw_flu$date <- NA
raw_flu$date <- ifelse(!is.na(raw_flu$npsdatecol) & raw_flu$year == 2016, 
                       as.integer(raw_flu$npsdatecol - as.Date("2016-05-02")), raw_flu$date)
raw_flu$date <- ifelse(!is.na(raw_flu$npsdatecol) & raw_flu$year == 2017, 
                       as.integer(raw_flu$npsdatecol - as.Date("2017-01-16")), raw_flu$date)
raw_flu$date <- ifelse(!is.na(raw_flu$npsdatecol) & raw_flu$year == 2018, 
                       as.integer(raw_flu$npsdatecol - as.Date("2018-01-15")), raw_flu$date)

clus <- raw_flu %>% filter(flu == "Pos")
clus$type_num <- NA
clus$type_num <- ifelse(clus$TYPE == "AH3", 0, clus$type_num)
clus$type_num <- ifelse(clus$TYPE == "AH1", 15, clus$type_num)
clus$type_num <- ifelse(clus$TYPE == "BYam", 30, clus$type_num)
clus$type_num <- ifelse(clus$TYPE == "BVic", 45, clus$type_num)

indid_l <- unique(c(clus$indid))
clusters <- c()

for (n in indid_l) {
  #subset data for that individual
  sub <- clus %>%
    filter(indid == n) 
  if (nrow(sub) == 1) {
    clusters <- c(clusters, 1)
  }
  if (nrow(sub) > 1) {
    sub <- sub %>%
      arrange(fuid) %>%
      dplyr::select(date, type_num)
    sub_dist <- dist(sub, method = "euclidean")
    hc_sub <- hclust(sub_dist, method = "single")
    
    clusters <- c(clusters, cutree(hc_sub, h = 14))
  }
}

clus$infcluster <- clusters
clus <- clus %>% 
  # making new variable to combine indid and the cluster/infection number (easier to work with later)
  mutate(indid_inf = paste(indid, "_I", infcluster, sep = "")) %>%
  dplyr::select(-type_num)
neg <- raw_flu %>%
  filter(flu != "Pos") %>%
  mutate(infcluster = NA) %>% 
  mutate(indid_inf = NA)
raw_flu <- clus %>%
  rbind(neg) %>%
  arrange(fuid)

rm(clus, hc_sub, neg, sub, clusters, indid_l, n, sub_dist)

# CT VALUE CLEAN ###############################################################

raw_flu$CT <- NA
raw_flu$CT <- ifelse(raw_flu$TYPE == "AH3" | raw_flu$TYPE == "AH1", raw_flu$npsinfact, raw_flu$CT)
raw_flu$CT <- ifelse(raw_flu$TYPE == "BYam" | raw_flu$TYPE == "BVic", raw_flu$npsinfbct, raw_flu$CT)
raw_flu$CT <- ifelse(raw_flu$TYPE == "Negative", 37, raw_flu$CT)
raw_flu$CT <- 37 - raw_flu$CT

# ADD IN PROFILE DATA - MISSING HAI TITER ########################

raw_flu <- raw_flu %>% 
  dplyr::select(1:5, 21, 24, 26:29)

indid_flu <- c(unique(raw_flu$indid))
profile <- read.csv('data/PHIRST 2016-2018 Metadata 2023-01-26.csv') %>% 
  rename(indid = ï..indid) %>%
  filter(indid %in% indid_flu)
flu_prof <- raw_flu %>% 
  left_join(profile, by = "indid") %>% 
  arrange(indid)

# Clean age categories
flu_prof$age_cat <- NA 
flu_prof$age_cat <- ifelse(flu_prof$age_at_consent <= 4, "<5", flu_prof$age_cat)
flu_prof$age_cat <- ifelse(flu_prof$age_at_consent >= 5 & flu_prof$age_at_consent <= 11, "5-11", flu_prof$age_cat)
flu_prof$age_cat <- ifelse(flu_prof$age_at_consent >= 12 & flu_prof$age_at_consent <= 18, "12-18", flu_prof$age_cat)
flu_prof$age_cat <- ifelse(flu_prof$age_at_consent >= 19 & flu_prof$age_at_consent <= 40, "19-40", flu_prof$age_cat)
flu_prof$age_cat <- ifelse(flu_prof$age_at_consent >= 41, ">40", flu_prof$age_cat)

# clean HIV variable
flu_prof$HIV <- NA
flu_prof$HIV <- ifelse(flu_prof$hiv_status == "Negative", "Negative", flu_prof$HIV)
flu_prof$HIV <- ifelse(flu_prof$hiv_status == "Unknown", "Unknown", flu_prof$HIV)
flu_prof$HIV <- ifelse(flu_prof$hiv_status == "Positive" & flu_prof$cd4_cat == "200-500", 
                       "Positive; 200-500", flu_prof$HIV)
flu_prof$HIV <- ifelse(flu_prof$hiv_status == "Positive" & flu_prof$cd4_cat == ">500", 
                       "Positive; >500", flu_prof$HIV)
flu_prof$HIV <- ifelse(flu_prof$hiv_status == "Positive" & flu_prof$cd4_cat == "<200", 
                       "Positive; <200", flu_prof$HIV)
flu_prof$HIV <- ifelse(flu_prof$hiv_status == "Positive" & flu_prof$cd4_cat == "", 
                       "Positive; <200", flu_prof$HIV)

# add a cleaner HIV variable 
flu_prof$cleanhiv <- NA
flu_prof$cleanhiv <- ifelse(flu_prof$HIV == "Negative", "Negative", flu_prof$cleanhiv)
flu_prof$cleanhiv <- ifelse(flu_prof$HIV == "Unknown", "Unknown", flu_prof$cleanhiv)
flu_prof$cleanhiv <- ifelse(flu_prof$HIV == "Positive; 200-500" | flu_prof$HIV ==  "Positive; >500", "Positive CD4 >=200", flu_prof$cleanhiv)
flu_prof$cleanhiv <- ifelse(flu_prof$HIV == "Positive; <200", "Positive CD4 <200", flu_prof$cleanhiv)


# continuous hh size
profile2 <- read.csv('data/PHIRST 2016-2018 Metadata 2023-02-02.csv') %>%
  dplyr::select(hh_id, true_hh_size) %>% 
  distinct(.keep_all = TRUE)

flu_prof <- flu_prof %>%
  left_join(profile2, by = "hh_id") %>%
  arrange(indid)

# ADDING IN SHEDDING DATA ######################################################

data_list = c("h1_5", "h1_5-11", "h1_12-18", "h1_19-40", "h1_40",
              "h3_5", "h3_5-11", "h3_12-18", "h3_19-40", "h3_40",
              "byam_5", "byam_5-11", "byam_12-18", "byam_19-40", "byam_40",
              "bvic_5", "bvic_5-11", "bvic_12-18", "bvic_19-40", "bvic_40")
shedding <- data.frame()
for(i in data_list){
  data <- read.csv(paste0("model/alldata_fits/data/", i, ".csv", sep = "")) %>%
    dplyr::select(iteration, indid, indid_inf, tp, x, wp, wc)
  shedding <- rbind(shedding, data)
}
rm(data, i, data_list)

indid_l <- c(unique(shedding$indid_inf))
means <- data.frame()
for(i in indid_l) {
  subset <- shedding %>% 
    filter(indid_inf == i) %>%
    mutate(meantp = mean(tp), 
           meanct = mean(37 - x), 
           meanwp = mean(wp),
           meanwc = mean(wc)) %>%
    filter(iteration == 1)
  means <- rbind(means, subset)
}
rm(subset, i, indid_l) 

id_list <- c(unique(flu_prof$indid_inf))
id_list <- id_list[!is.na(id_list)]
date_dat <- data.frame()
for(i in id_list) {
  sub <- flu_prof %>%
    filter(indid_inf == i) %>% 
    dplyr::select(indid, indid_inf, date, npsdatecol, CT, TYPE)
  tp <- c(sub$date[which(sub$CT == max(sub$CT))])
  sub$center_date <- sub$date[which(sub$date == tp[1])]
  date_dat <- rbind(date_dat, sub[1,])
}
rm(id_list, sub, i, tp)

date_dat <- date_dat %>% 
  dplyr::select(indid_inf, center_date)
means <- means %>%
  left_join(date_dat, by = "indid_inf")

means <- means %>% 
  mutate(center_point = center_date + meantp) %>%
  mutate(start_point = center_point - meanwp) %>%
  mutate(end_point = center_point + meanwc) %>% 
  mutate(round_start = round(start_point)) %>% 
  mutate(tot_duration = meanwp + meanwc) %>%
  arrange(indid)

# calculating proliferation/clearance slope to use later to get day by day estimate of viral load
means <- means %>%
  mutate(prolif_slope = meanct / meanwp) %>%
  mutate(clear_slope = -meanct / meanwc) %>% 
  dplyr::select(indid_inf, center_date, center_point, start_point, end_point, round_start, 
                meanct, prolif_slope, clear_slope, meanwp, meanwc, tot_duration)

# combine with master data
flu_prof <- flu_prof %>% 
  left_join(means, by = "indid_inf")

write.csv(flu_prof, 'transmission/PHIRST flu shedding and profile 2016-2018.csv')

















