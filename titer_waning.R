

library(tidyr)
library(dplyr)
library(ggplot2)


hai_dat <- read.csv('data/PHIRST qry_flu_hai_2016-2018_2023-02-13.csv') %>% 
  rename(year = ï..year)

all <- read.csv('transmission/PHIRST flu shedding and profile 2016-2018.csv')

ind_w_infect <- c(unique(all$indid[which(all$flu == "Pos")]))

hai_no_inf <- hai_dat %>% filter(!(ind_id %in% ind_w_infect))

hai <- hai_no_inf %>% 
  select(year, ind_id, draw, flua_h1n1pdm_gm, flua_h3n2_gm, flub_victoria_gm, flub_yamagata_gm)
d2016 <- subset(hai, year == 2016) %>%
  filter(draw == 1 | draw == 2)
d20178 <- subset(hai, year == 2017 | year == 2018) %>%
  filter(draw == 1 | draw == 2 | draw == 3)
hai <- rbind(d2016, d20178)

ages <- all %>% 
  select(indid, age_cat) %>% 
  rename(ind_id = indid) %>% 
  distinct(ind_id, .keep_all = TRUE)

hai <- hai %>% left_join(ages, by = "ind_id")



# 2016 - only two points
h2016 <- hai %>% 
  filter(year == 2016 & !is.na(flua_h1n1pdm_gm))
incl1 <- unique(h2016$ind_id[which(h2016$draw == 2)])
incl2 <- unique(h2016$ind_id[which(h2016$draw == 1)])
incl <- intersect(incl1, incl2)
h12016 <- h2016 %>% filter(ind_id %in% incl)
 
h2016$changeh1 <- NA
h2016$changeh3 <- NA
h2016$changeyam <- NA
h2016$changevic <- NA
for(i in incl) {
  h2016$changeh1[which(h2016$ind_id == i)] <- h2016$flua_h1n1pdm_gm[which(h2016$ind_id == i & h2016$draw == 2)] /
    h2016$flua_h1n1pdm_gm[which(h2016$ind_id == i & h2016$draw == 1)]
  h2016$changeh3[which(h2016$ind_id == i)] <- h2016$flua_h3n2_gm[which(h2016$ind_id == i & h2016$draw == 2)] /
    h2016$flua_h3n2_gm[which(h2016$ind_id == i & h2016$draw == 1)]
  h2016$changeyam[which(h2016$ind_id == i)] <- h2016$flub_yamagata_gm[which(h2016$ind_id == i & h2016$draw == 2)] /
    h2016$flub_yamagata_gm[which(h2016$ind_id == i & h2016$draw == 1)]
  h2016$changevic[which(h2016$ind_id == i)] <- h2016$flub_victoria_gm[which(h2016$ind_id == i & h2016$draw == 2)] /
    h2016$flub_victoria_gm[which(h2016$ind_id == i & h2016$draw == 1)]
}

h1718 <- hai %>% 
  filter((year == 2017 | year == 2018) & !is.na(flua_h1n1pdm_gm))
incl1 <- unique(h1718$ind_id[which(h1718$draw == 3)])
incl2 <- unique(h1718$ind_id[which(h1718$draw == 1)])
incl <- intersect(incl1, incl2)
h1718 <- h1718 %>% filter(ind_id %in% incl)

h1718$changeh1 <- NA
h1718$changeh3 <- NA
h1718$changeyam <- NA
h1718$changevic <- NA
for(i in incl) {
  h1718$changeh1[which(h1718$ind_id == i)] <- h1718$flua_h1n1pdm_gm[which(h1718$ind_id == i & h1718$draw == 3)] /
    h1718$flua_h1n1pdm_gm[which(h1718$ind_id == i & h1718$draw == 1)]
  h1718$changeh3[which(h1718$ind_id == i)] <- h1718$flua_h3n2_gm[which(h1718$ind_id == i & h1718$draw == 3)] /
    h1718$flua_h3n2_gm[which(h1718$ind_id == i & h1718$draw == 1)]
  h1718$changeyam[which(h1718$ind_id == i)] <- h1718$flub_yamagata_gm[which(h1718$ind_id == i & h1718$draw == 3)] /
    h1718$flub_yamagata_gm[which(h1718$ind_id == i & h1718$draw == 1)]
  h1718$changevic[which(h1718$ind_id == i)] <- h1718$flub_victoria_gm[which(h1718$ind_id == i & h1718$draw == 3)] /
    h1718$flub_victoria_gm[which(h1718$ind_id == i & h1718$draw == 1)]
}

hai_clean <- rbind(h2016, h1718)

# exclude over 4 because that suggests infection
hai_clean$age_cat <- factor(hai_clean$age_cat, levels = c("<5", "5-11", "12-18", "19-40", ">40"), ordered = TRUE)

h1 <- hai_clean %>% filter(changeh1 <= 1)
h1p2 <- ggplot(h1, aes(x = age_cat, y = changeh1)) +
  geom_violin(aes(fill = age_cat)) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., digits = 3)), vjust = -0.7) +
  scale_fill_brewer(palette="Blues") +
  theme_minimal() + 
  labs(title = "H1N1 Pre vs. Post Season (all increase excluded)", xlab = "Age Group", 
       ylab = "Fractional Change in Titer (Post-Season / Pre-Seaon Titer)")

h3 <- hai_clean %>% filter(changeh3 <= 1)
h3p2 <- ggplot(h3, aes(x = age_cat, y = changeh3)) +
  geom_violin(aes(fill = age_cat)) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., digits = 3)), vjust = -0.7) +
  scale_fill_brewer(palette="Blues") +
  theme_minimal() + 
  labs(title = "H3N2 Pre vs. Post Season (all increase excluded)", xlab = "Age Group", 
       ylab = "Fractional Change in Titer (Post-Season / Pre-Seaon Titer)")

yam <- hai_clean %>% filter(changeyam <= 1)
y1p2 <- ggplot(yam, aes(x = age_cat, y = changeyam)) +
  geom_violin(aes(fill = age_cat)) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., digits = 3)), vjust = -0.7) +
  scale_fill_brewer(palette="Blues") +
  theme_minimal() + 
  labs(title = "Yamagata Pre vs. Post Season (all increase excluded)", xlab = "Age Group", 
       ylab = "Fractional Change in Titer (Post-Season / Pre-Seaon Titer)")

vic <- hai_clean %>% filter(changevic <= 1)
v1p2 <- ggplot(vic, aes(x = age_cat, y = changevic)) +
  geom_violin(aes(fill = age_cat)) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., digits = 3)), vjust = -0.7) +
  scale_fill_brewer(palette="Blues") +
  theme_minimal() + 
  labs(title = "Victoria Pre vs. Post Season (all increase excluded)", xlab = "Age Group", 
       ylab = "Fractional Change in Titer (Post-Season / Pre-Seaon Titer)")

pdf(file = 'figures/titer_waning.pdf')

print(h1p1)
print(h3p1)
print(y1p1)
print(v1p1)
print(h1p2)
print(h3p2)
print(y1p2)
print(v1p2)

dev.off()





