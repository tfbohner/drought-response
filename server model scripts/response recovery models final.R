# models for ch2 draft results

library(tidyverse)
library(brms)
library(outliers)


sitechar <- read.csv("Raw Data/site characteristic table.csv")
l_leg <- read.csv("Processed Data/legacy_1post.csv") %>% 
  group_by(Species) %>%
  mutate(log_hegyi=log(hegyi),
         dfx_score=scores(drought_fx, type="z"),
         recov_score=scores(recovery, type="z"),
         drought_fx=ifelse(dfx_score<=-3, NA, drought_fx),
         simple_comp=ifelse(simple_comp>900, NA, simple_comp),
         comp=BA,
         mean_spei12=mean(spei12, na.rm=T),
         sd_spei12=sd(spei12, na.rm=T),
         mean_dfx=mean(drought_fx, na.rm=T),
         sd_dfx=sd(drought_fx, na.rm=T),
         mean_dbh=mean(DBH, na.rm=T),
         sd_dbh=sd(DBH, na.rm=T),
         mean_comp=mean(comp, na.rm=T),
         sd_comp=sd(comp, na.rm=T)) %>% 
  left_join(sitechar)

write.csv(l_leg, "Processed Data/legacy_modeling_data.csv")

newdat <- l_leg %>% 
  mutate_at(vars(spei12, drought_fx, DBH, comp), scale)



summary <- group_by(l_leg, tree.uniqueID) %>% 
  summarize(n_events=length(drought_fx))

# for(i in 1:9){
#   assign(paste("resp", i, sep=""), readRDS(paste("Brms models/final rr models/resp", i, ".rds", sep="")), envir = .GlobalEnv)
# }
# 
# for(i in 1:10){
#   assign(paste("recov", i, sep=""), readRDS(paste("Brms models/final rr models/recov", i, ".rds", sep="")), envir = .GlobalEnv)
# }

#### Q1: Species level differences----
## response
resp0 <- brm(bf(drought_fx~Species + (1|tree.uniqueID), sigma~Species), family=gaussian(), data=l_leg, cores=4)
hyp <- "exp(sigma_Intercept + sigma_SpeciesPJ) > exp(sigma_Intercept)"
(hyp <- hypothesis(resp0, hyp))
plot(hyp)

hyp <- "exp(sigma_Intercept + sigma_SpeciesPL) > exp(sigma_Intercept)"
(hyp <- hypothesis(resp0, hyp))
plot(hyp)

hyp <- "exp(sigma_SpeciesPL) > exp(sigma_SpeciesPJ)"
(hyp <- hypothesis(resp0, hyp))
plot(hyp)

plot(marginal_effects(resp0), points = TRUE)

# resp1 <- brm(drought_fx~Species + (1|tree.uniqueID), data=l_leg, cores=4)
resp2 <- update(resp0, formula. = ~ . - Species + spei12, newdata=l_leg, cores=4)
resp3 <- update(resp0, formula. = ~ . + spei12, newdata=l_leg, cores=4)
resp4 <- update(resp0, formula. = ~ . + spei12*Species, newdata=l_leg, cores=4)


for(i in 2:4){
  mod <- get(paste("resp", i, sep="")) 
  saveRDS(mod, paste("Brms models/final rr models/resp", i, ".rds", sep=""))
  
  mod <- get(paste("resp", i, sep="")) 
  lx <- loo(mod)
  assign(paste0("l", i), lx)
}

loo_model_weights(list(l1, l2, l3, l4), method="pseudobma")

loo_model_weights(list(l3, l4), method="pseudobma")

model_weights(resp1, resp2, resp3, resp4, weights = "waic")
# b1 <- add_criterion(resp1, "waic")
# b2 <- add_criterion(resp2, "waic")
# b3 <- add_criterion(resp3, "waic")
# b4 <- add_criterion(resp4, "waic")
# 
# w <- loo_compare(b1, b2, b3, b4, criterion = "waic")
# print(w, simplify = F)
# cbind(waic_diff = w[, 1] * -2,
#       se        = w[, 2] *  2)

## recovery
recov0 <- brm(bf(recovery~Species + (1|tree.uniqueID), sigma~Species), family=gaussian(), data=l_leg, cores=4)
hyp <- "exp(sigma_Intercept + sigma_SpeciesPJ) > exp(sigma_Intercept)"
(hyp <- hypothesis(recov0, hyp))
plot(hyp)

hyp <- "exp(sigma_Intercept + sigma_SpeciesPL) > exp(sigma_Intercept)"
(hyp <- hypothesis(recov0, hyp))
plot(hyp)

hyp <- "exp(sigma_SpeciesPL) > exp(sigma_SpeciesPJ)"
(hyp <- hypothesis(recov, hyp))
plot(hyp)



# recov1 <- brm(recovery~Species + (1|tree.uniqueID), data=l_leg, cores=4)
recov2 <- update(recov0, formula. = ~ . + drought_fx - Species, newdata=l_leg, cores=4)
recov3 <- update(recov0, formula. = ~ . + drought_fx, newdata=l_leg, cores=4)
recov4 <- update(recov0, formula. = ~ . + drought_fx + drought_fx:Species, newdata=l_leg, cores=4)

plot(marginal_effects(recov4), points = TRUE)

for(i in 2:4){
  mod <- get(paste("recov", i, sep="")) 
  saveRDS(mod, paste("Brms models/final rr models/recov", i, ".rds", sep=""))
  
  mod <- get(paste("recov", i, sep="")) 
  lx <- loo(mod)
  assign(paste0("l", i, "a"), lx)
}

loo_model_weights(list(l1a, l2a, l3a, l4a))
loo_model_weights(list(l3a, l4a))
model_weights(recov1, recov2, recov3, recov4, weights = "waic")
# b1a <- add_criterion(recov1, "waic")
# b2a <- add_criterion(recov2, "waic")
# b3a <- add_criterion(recov3, "waic")
# b4a <- add_criterion(recov4, "waic")
# 
# w <- loo_compare(b1a, b2a, b3a, b4a, criterion = "waic")
# print(w, simplify = F)
# cbind(waic_diff = w[, 1] * -2,
#       se        = w[, 2] *  2)

#### Q2 Regional differences----
## response

# resp5 <- brm(drought_fx~spei12 + (1|tree.uniqueID) + (1|Region/hilo), data=l_leg, cores=4, control=list(adapt_delta=0.99))
# resp6 <- brm(drought_fx~spei12 + (1|tree.uniqueID) + (spei12|Region/hilo), data=l_leg, cores=4, control=list(adapt_delta=0.99))
# resp7 <- brm(drought_fx~spei12 + (1|tree.uniqueID) + (1|Region), data=l_leg, cores=4, control=list(adapt_delta=0.99))
# resp8 <- brm(drought_fx~spei12 + (1|tree.uniqueID) + (spei12|Region), data=l_leg, cores=4, control=list(adapt_delta=0.99))
# 
# for(i in 5:8){
#   mod <- get(paste("resp", i, sep="")) 
#   saveRDS(mod, paste("Brms models/final rr models/resp", i, ".rds", sep=""))
#   
#   mod <- get(paste("resp", i, sep="")) 
#   lx <- loo(mod)
#   assign(paste0("l", i), lx)
# }
# 
# loo_model_weights(list(l5, l6, l7, l8))
# model_weights(resp5, resp6, resp7, resp8, weights = "waic")

spp_modeler <- function(data, spp, form) {
  sub <- filter(data, Species==spp)
  brm(form, sub, cores=4, control=list(adapt_delta=0.99))
}

spp_modeler2 <- function(data, spp, model) {
  sub <- filter(data, Species==spp)
  update(model, newdata=sub, cores=4, control=list(adapt_delta=0.99))
}

resp4 <- bf(drought_fx~spei12 + (1|tree.uniqueID))
resp4_ac <- spp_modeler(l_leg, "AC", resp4)
resp4_pj <- spp_modeler2(l_leg, "PJ", resp4_ac)
resp4_pl <- spp_modeler2(l_leg, "PL", resp4_ac)

resp5 <- bf(drought_fx~spei12 + (1|tree.uniqueID) + (1|Region/hilo))
resp5_ac <- spp_modeler(l_leg, "AC", resp5)
resp5_pj <- spp_modeler2(l_leg, "PJ", resp5_ac)
resp5_pl <- spp_modeler2(l_leg, "PL", resp5_ac)

resp6 <- bf(drought_fx~spei12 + (1|tree.uniqueID) + (spei12|Region/hilo))
resp6_ac <- spp_modeler(l_leg, "AC", resp6)
resp6_pj <- spp_modeler2(l_leg, "PJ", resp6_ac)
resp6_pl <- spp_modeler2(l_leg, "PL", resp6_ac)

resp7 <- bf(drought_fx~spei12 + (1|tree.uniqueID) + (1|Region))
resp7_ac <- spp_modeler(l_leg, "AC", resp7)
resp7_pj <- spp_modeler2(l_leg, "PJ", resp7_ac)
resp7_pl <- spp_modeler2(l_leg, "PL", resp7_ac)

resp8 <- bf(drought_fx~spei12 + (1|tree.uniqueID) + (spei12|Region))
resp8_ac <- spp_modeler(l_leg, "AC", resp8)
resp8_pj <- spp_modeler2(l_leg, "PJ", resp8_ac)
resp8_pl <- spp_modeler2(l_leg, "PL", resp8_ac)

mod_list <- c("resp4_ac", "resp4_pj", "resp4_pl","resp5_ac", "resp5_pj", "resp5_pl", "resp6_ac", "resp6_pj", "resp6_pl", "resp7_ac", "resp7_pj", "resp7_pl", "resp8_ac", "resp8_pj", "resp8_pl")

for(i in mod_list) {
  mod <- get(mod_list[which(mod_list==i)])
  saveRDS(mod, paste("Brms models/final rr models/", i, ".rds", sep=""))
  
  lx <- loo(mod)
  assign(paste0("l_", i), lx)
  
}

## recovery

recov4 <- bf(recovery~drought_fx + (1|tree.uniqueID))
recov4_ac <- spp_modeler(l_leg, "AC", recov4)
recov4_pj <- spp_modeler2(l_leg, "PJ", recov4_ac)
recov4_pl <- spp_modeler2(l_leg, "PL", recov4_ac)

recov5 <- bf(recovery~drought_fx + (1|tree.uniqueID) + (1|Region/hilo))
recov5_ac <- spp_modeler(l_leg, "AC", recov5)
recov5_pj <- spp_modeler2(l_leg, "PJ", recov5_ac)
recov5_pl <- spp_modeler2(l_leg, "PL", recov5_ac)

recov6 <- bf(recovery~drought_fx + (1|tree.uniqueID) + (drought_fx|Region/hilo))
recov6_ac <- spp_modeler(l_leg, "AC", recov6)
recov6_pj <- spp_modeler2(l_leg, "PJ", recov6_ac)
recov6_pl <- spp_modeler2(l_leg, "PL", recov6_ac)

recov7 <- bf(recovery~drought_fx + (1|tree.uniqueID) + (1|Region))
recov7_ac <- spp_modeler(l_leg, "AC", recov7)
recov7_pj <- spp_modeler2(l_leg, "PJ", recov7_ac)
recov7_pl <- spp_modeler2(l_leg, "PL", recov7_ac)

recov8 <- bf(recovery~drought_fx + (1|tree.uniqueID) + (drought_fx|Region))
recov8_ac <- spp_modeler(l_leg, "AC", recov8)
recov8_pj <- spp_modeler2(l_leg, "PJ", recov8_ac)
recov8_pl <- spp_modeler2(l_leg, "PL", recov8_ac)

mod_list <- c("recov4_ac", "recov4_pj", "recov4_pl","recov5_ac", "recov5_pj", "recov5_pl", "recov6_ac", "recov6_pj", "recov6_pl", "recov7_ac", "recov7_pj", "recov7_pl", "recov8_ac", "recov8_pj", "recov8_pl")

for(i in mod_list) {
  mod <- get(mod_list[which(mod_list==i)])
  saveRDS(mod, paste("Brms models/final rr models/", i, ".rds", sep=""))

  lx <- loo(mod)
  assign(paste0("l_", i), lx)

}

# #### Q3 Competition and diameter----
spp_modeler <- function(data, spp, form) {
  sub <- filter(data, Species==spp)
  brm(form, sub, cores=4, #prior=set_prior(horseshoe(df = 4)),
      iter=2000, control=list(adapt_delta=0.999, max_treedepth=12)) 
}

spp_modeler2 <- function(data, spp, model) {
  sub <- filter(data, Species==spp)
  update(model, newdata=sub, cores=4, #prior=set_prior(horseshoe(df = 4)), 
         iter=2000, control=list(adapt_delta=0.999, max_treedepth=12))
}


### Full models ----
            
# resp_full <- brm(drought_fx~spei12 + DBH + log_hegyi + spei12:DBH + spei12:log_hegyi +
#                    (1|tree.uniqueID) + (spei12+DBH+log_hegyi|Region + Species),
#                  data=l_leg, cores=4, prior=set_prior(horseshoe(df = 3, par_ratio = 0.1)), iter=4000, control=list(adapt_delta=0.999, max_treedepth=14))
# saveRDS(resp_full, "Brms models/final rr models/resp_full.rds")
# 
# recov_full <- brm(recovery~drought_fx + DBH + log_hegyi + drought_fx:DBH + drought_fx:log_hegyi +
#                    (1|tree.uniqueID) + (drought_fx+DBH+log_hegyi|Region + Species),
#                  data=l_leg, cores=4, prior=set_prior(horseshoe(df = 3, par_ratio = 0.1)),  iter=4000, control=list(adapt_delta=0.999, max_treedepth=14))
# saveRDS(recov_full, "Brms models/final rr models/recov_full.rds")



### Species Full models ----

# resp_mod <- bf(drought_fx~spei12 + DBH + log_hegyi + spei12:DBH + spei12:log_hegyi + 
#                   (1|tree.uniqueID) + 
#                   (spei12+DBH+log_hegyi + spei12:DBH + spei12:log_hegyi|Region))
resp_mod <- bf(drought_fx~spei12 + DBH + comp + spei12:DBH + spei12:comp + 
                 (1|tree.uniqueID) + 
                 (spei12|Region))

resp_full_ac <- spp_modeler(l_leg, "AC", resp_mod)
saveRDS(resp_full_ac, "Brms models/final rr models/resp_full_ac.rds")
resp_full_pj <- spp_modeler2(l_leg, "PJ", resp_full_ac)
saveRDS(resp_full_pj, "Brms models/final rr models/resp_full_pj.rds")
resp_full_pl <- spp_modeler2(l_leg, "PL", resp_full_ac)
saveRDS(resp_full_pl, "Brms models/final rr models/resp_full_pl.rds")


# recov_mod <- bf(recovery~drought_fx + DBH + log_hegyi + drought_fx:DBH + drought_fx:log_hegyi + 
#                   (1|tree.uniqueID) + 
#                   (drought_fx+DBH+log_hegyi + drought_fx:DBH + drought_fx:log_hegyi|Region))
# recov_mod <- bf(recovery~drought_fx + DBH + log_hegyi + drought_fx:DBH + drought_fx:log_hegyi + 
#                   (1|tree.uniqueID) + 
#                   (drought_fx|Region))

recov_mod <- bf(recovery~drought_fx + DBH + comp + drought_fx:DBH + drought_fx:comp + 
                  (1|tree.uniqueID) + 
                  (drought_fx|Region))

recov_full_ac <- spp_modeler(l_leg, "AC", recov_mod)
saveRDS(recov_full_ac, "Brms models/final rr models/recov_full_ac.rds")
recov_full_pj <- spp_modeler2(l_leg, "PJ", recov_full_ac)
saveRDS(recov_full_pj, "Brms models/final rr models/recov_full_pj.rds")
recov_full_pl <- spp_modeler2(l_leg, "PL", recov_full_ac)
saveRDS(recov_full_pl, "Brms models/final rr models/recov_full_pl.rds")


#### Scaled data ----
resp_mod <- bf(drought_fx~0 + Intercept + spei12 + DBH + comp + spei12:DBH + spei12:comp + 
                 (1|tree.uniqueID) + 
                 (spei12|Region))

resp_full_ac <- spp_modeler(newdat, "AC", resp_mod)
resp_full_pj <- spp_modeler2(newdat, "PJ", resp_full_ac)
resp_full_pl <- spp_modeler2(newdat, "PL", resp_full_ac)

recov_mod <- bf(recovery~0 + Intercept + drought_fx + DBH + comp + drought_fx:DBH + drought_fx:comp + 
                  (1|tree.uniqueID) + 
                  (drought_fx|Region))

recov_full_ac <- spp_modeler(newdat, "AC", recov_mod)
recov_full_pj <- spp_modeler2(newdat, "PJ", recov_full_ac)
recov_full_pl <- spp_modeler2(newdat, "PL", recov_full_ac)
