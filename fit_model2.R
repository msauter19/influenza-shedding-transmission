
# FITTING THE MODEL WITH INSTANCES OF REINFECTION, COINFECTION, AND ASSIGNED UN-SUBTYPED


rm(list = ls())

.rs.restartR()

library(rstan)
library(tidyverse)
library(dplyr)
library(purrr)
library(shinystan)
library(bayesplot)
source('export.R')
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

mod_dat <- read.csv('data/model_data.csv') 
one_dat <- read.csv('data/one_data.csv')

# prep data for model ##########################################################

# csv messed up age categories 
mod_dat$age_cat <- NA 
mod_dat$age_cat <- ifelse(mod_dat$age <= 4, "<5", mod_dat$age_cat)
mod_dat$age_cat <- ifelse(mod_dat$age >= 5 & mod_dat$age <= 11, "5-11", mod_dat$age_cat)
mod_dat$age_cat <- ifelse(mod_dat$age >= 12 & mod_dat$age <= 18, "12-18", mod_dat$age_cat)
mod_dat$age_cat <- ifelse(mod_dat$age >= 19 & mod_dat$age <= 40, "19-40", mod_dat$age_cat)
mod_dat$age_cat <- ifelse(mod_dat$age >= 41, ">40", mod_dat$age_cat)

one_dat$age_cat <- NA 
one_dat$age_cat <- ifelse(one_dat$age <= 4, "<5", one_dat$age_cat)
one_dat$age_cat <- ifelse(one_dat$age >= 5 & one_dat$age <= 11, "5-11", one_dat$age_cat)
one_dat$age_cat <- ifelse(one_dat$age >= 12 & one_dat$age <= 18, "12-18", one_dat$age_cat)
one_dat$age_cat <- ifelse(one_dat$age >= 19 & one_dat$age <= 40, "19-40", one_dat$age_cat)
one_dat$age_cat <- ifelse(one_dat$age >= 41, ">40", one_dat$age_cat)

# subset
type <- "BVic"
age <- "12-18"
indid_l <- c(unique(mod_dat$indid_inf[which(mod_dat$TYPE == type & mod_dat$age_cat == age)]))
subset <- mod_dat %>%
  filter(indid_inf %in% indid_l)
person_dat <- one_dat %>% 
  filter(indid_inf %in% indid_l)

# must change clean_id to consecutive list of individuals
id_list <- sort(unique(subset$indid_inf))
id_df <- data.frame(indid_inf = id_list, id_clean = 1:length(id_list), stringsAsFactors=FALSE)
subset <- left_join(subset, id_df, by="indid_inf")
person_dat <- left_join(person_dat, id_df, by = "indid_inf") %>%
  dplyr::select(id_clean, everything())

rm(id_df, id_list, indid_l)

# parameters ###################################################################


if (type == "AH3" | type == "AH1") {
  # for influenza a
  subset$pars <- list(N = nrow(subset), 
             n_id = length(unique(subset$id_clean)),
             id = subset$id_clean,
             t = subset$date, 
             ct = subset$CT,
             lod = as.numeric(37),
             tpsd_p = as.numeric(2),
             dctpmean_p = as.numeric(37/2),
             dctpsd_p = as.numeric(37/6),
             wpmax = as.numeric(15),
             wpmean_p = as.numeric(2),
             wpsd_p = as.numeric(1), 
             wcmax = as.numeric(30),
             wcmean_p = as.numeric(4),
             wcsd_p = as.numeric(3),
             sigma = as.numeric(1)
             )
}
if (type == "BVic" | type == "BYam") {
  # for influenza b
  pars <- list(N = nrow(subset), 
               n_id = length(unique(subset$id_clean)),
               id = subset$id_clean,
               t = subset$date, 
               ct = subset$CT,
               lod = as.numeric(37),
               tpsd_p = as.numeric(2),
               dctpmean_p = as.numeric(37/2),
               dctpsd_p = as.numeric(37/6),
               wpmax = as.numeric(15),
               wpmean_p = as.numeric(2),
               wpsd_p = as.numeric(1), 
               wcmax = as.numeric(40),
               wcmean_p = as.numeric(5),
               wcsd_p = as.numeric(3),
               sigma = as.numeric(1)
  )
}

# running model ################################################################

c_mod <- stanc("model/flu_ct_model.stan")
ct_model <- stan_model(stanc_ret = c_mod)

fit_startq <- Sys.time()
ct_fit <- sampling(ct_model, 
                   data = list(
                     N = as.list(pars)$N, 
                     n_id = as.list(pars)$n_id,
                     id = as.list(pars)$id,
                     t = as.list(pars)$t, 
                     ct = as.list(pars)$ct, 
                     lod = as.list(pars)$lod,
                     tpsd_p = as.list(pars)$tpsd_p,
                     dctpmean_p = as.list(pars)$dctpmean_p,
                     dctpsd_p = as.list(pars)$dctpsd_p,
                     wpmax = as.list(pars)$wpmax,
                     wpmean_p = as.list(pars)$wpmean_p,
                     wpsd_p = as.list(pars)$wpsd_p,
                     wcmax = as.list(pars)$wcmax,
                     wcmean_p = as.list(pars)$wcmean_p,
                     wcsd_p = as.list(pars)$wcsd_p,
                     sigma = as.list(pars)$sigma), 
                   iter=1000, chains=4) #, control = list(adapt_delta = 0.9))
fit_endq <- Sys.time()
print(paste0("Fit time: ",difftime(fit_endq, fit_startq, units="min")," mins"))
rm(fit_endq, fit_startq)

# saving data ##################################################################

load(file = "model/alldata_fits/ct_fits/bvic_12-18")
c <- get_sampler_params(ct_fit)
summary(do.call(rbind, c), digits = 2)

print(ct_fit, pars = c("dctpmean", "dctpsd", "wpmean", "wpsd", "wcmean", "wcsd"))
print(ct_fit, pars = c("tp"))
print(ct_fit, pars = c("x"))
print(ct_fit, pars = c("wp"))
print(ct_fit, pars = c("wc"))
print(ct_fit, pars = c("mu"))

# functions adapted from code from kissler ####
make_indiv_params_df <- function(extracted_params, parnames, n_indiv){
  out <- reduce(lapply(parnames, function(x) parseparam(extracted_params,x,n_indiv)), cbind) %>% 
    as_tibble %>% 
    mutate(iteration=1:n()) %>% 
    pivot_longer(-iteration) %>% 
    separate(name, c("param","id"), sep="_") %>%
    pivot_wider(c("id","iteration"), names_from="param", values_from="value")
}
parseparam <- function(extracted_params, parname, n_indiv){
  as_tibble(setNames(as.data.frame(extracted_params[[parname]]), makenames(parname,n_indiv)))
}
makenames <- function(parname, n_indiv){
  unlist(lapply(1:n_indiv, function(x) paste0(parname, "_", x)))
}
make_shared_params_df <- function(extracted_params, parnames){
   out <- reduce(lapply(parnames, function(x) 
    as_tibble(setNames(as.data.frame(extracted_params[[x]]),x))
  ), cbind) %>%
    as_tibble() %>%
    mutate(iteration=1:n())
}


params <- rstan::extract(ct_fit)
indiv_params_df <- make_indiv_params_df(params, c("tp","x","wp","wc"), pars$n_id) %>%
  rename(id_clean = id)
shared_params_df <- make_shared_params_df(params, c("dctpmean","wpmean","wcmean","dctpsd","wpsd","wcsd"))
params_df <- indiv_params_df %>% 
  left_join(shared_params_df, by="iteration")

person_dat$id_clean <- as.character(person_dat$id_clean)
params_df <- params_df %>%
  left_join(person_dat, by="id_clean") 
#params_df$center_date = as.Date(params_df$center_date, format = "%m/%d/%Y")
#params_df$tp_trans = params_df$center_date + round(params_df$tp, digits = 0)

save(ct_fit, file = paste("model/alldata_fits/ct_fits/", type, age, ".csv", sep = ""))
write.csv(params_df, file = paste("model/alldata_fits/data/", type, age, ".csv", sep = ""))

# plot data ####################################################################

traceplot(ct_fit, pars = c("dctpmean", "dctpsd", "wpmean", "wpsd", "wcmean", "wcsd"), 
          inc_warmup = TRUE, nrow = 2)
tpsampling <- traceplot(ct_fit, pars = c("tp"), inc_warmup = FALSE)
export_plot(tpsampling, "fits tp sampling", 15, 10)



# plot for fit estimates compared to data
params_df2 <- indiv_params_df 
params_df2$tpmean <- NA
params_df2$xmean <- NA
params_df2$wpmean <- NA 
params_df2$wcmean <- NA

for(i in 1:pars$n_id) {
  params_df2$tpmean[which(params_df2$id_clean == i)] <- mean(params_df2$tp[which(params_df2$id_clean == i)])
  params_df2$xmean[which(params_df2$id_clean == i)] <- mean(params_df2$x[which(params_df2$id_clean == i)])
  params_df2$wpmean[which(params_df2$id_clean == i)] <- mean(params_df2$wp[which(params_df2$id_clean == i)])
  params_df2$wcmean[which(params_df2$id_clean == i)] <- mean(params_df2$wc[which(params_df2$id_clean == i)])
}

plot_ct_fit <- function(params_df, lod, indiv_data, ntraces){
  params_df %>% 
    sample_n(ntraces) %>% 
    ggplot() + 
    # Plot traces:
    geom_segment(aes(x=-Inf, y=lod, xend=tp-wp, yend=lod), alpha=0.3, col="#73D055FF") + 
    geom_segment(aes(x=tp-wp, y=lod, xend=tp, yend=lod-x), alpha=0.3, col="#73D055FF") + 
    geom_segment(aes(x=tp, y=lod-x, xend=tp+wc, yend=lod), alpha=0.3, col="#73D055FF") + 
    geom_segment(aes(x=tp+wc, y=lod, xend=Inf, yend=lod), alpha=0.3, col="#73D055FF") + 
    geom_segment(aes(x=-Inf, y=lod, xend=tpmean-wpmean, yend=lod), size = 1.5, col="#2D708EFF") + 
    geom_segment(aes(x=tpmean-wpmean, y=lod, xend=tpmean, yend=lod-xmean), size = 1.5, col="#2D708EFF") + 
    geom_segment(aes(x=tpmean, y=lod-xmean, xend=tpmean+wcmean, yend=lod), size = 1.5, col="#2D708EFF") + 
    geom_segment(aes(x=tpmean+wcmean, y=lod, xend=Inf, yend=lod), size = 1.5, col="#2D708EFF") + 
    # Plot data:
    geom_point(data=indiv_data, aes(x=date, y=CT), size=3) + 
    theme_classic() + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
    labs(x="Time (days)", y="Ct") + 
    scale_y_reverse() + 
    facet_wrap(~id_clean)
}
params_df2$id_clean <- as.numeric(params_df2$id_clean)
fitsplot <- plot_ct_fit(params_df2, 37, subset, ntraces = 500)
export_plot(fitsplot, "fits plot", 15, 10)


rm(ct_fit, subset, indiv_params_df, shared_params_df, params, params_df, person_dat, pars)





