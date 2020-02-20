### predict growth without post-drought years

library(tidyverse)
library(shiny)
library(brms)
library(bayesplot)

tree_rings <- read.csv("Processed Data/matched_rings_field.csv")
climate <- read.csv("Processed Data/precip_temp_spei.csv") %>% 
  mutate(Site=str_sub(Name, 1,2),
         Neighborhood=as.numeric(str_sub(Name, 3,3))) %>% 
  select(-c(Name, hydroyear, siteno, Region, hilo)) #

alldata <- left_join(tree_rings, climate, by=c("Site", "Neighborhood", "Year")) %>% 
  mutate(unique.nb = paste0(Site, Neighborhood, sep=""))

datasub <- filter(alldata, zeroRing=="no")
datasub <- filter(alldata, zeroRing=="no", Species!="AM") 
datasub$Species[str_which(datasub$Species, "PP")] <- "PJ" 
datasub$Species <- droplevels(datasub$Species)
newdata <- datasub

datasub <- filter(newdata, Year>2008)


mod0 <- brm(spline_growth ~ (spei12|Species), data=datasub, family=gaussian(), cores=4, control=list(adapt_delta=0.99))
saveRDS(mod0, "Brms models/spline_growth models/spline_growth0.rds")
mod1 <- brm(spline_growth ~ Species + (spei12|Species), data=datasub, family=gaussian(), cores=4, control=list(adapt_delta=0.99))
saveRDS(mod1, "Brms models/spline_growth models/spline_growth1.rds")
mod2 <- brm(spline_growth ~ spei12 * Species, data=datasub, family=gaussian(), cores=4, control=list(adapt_delta=0.99))
saveRDS(mod2, "Brms models/spline_growth models/spline_growth2.rds")


mod3 <- brm(spline_growth ~ Species + (spei12|Species/Region), data=datasub, family=gaussian(), cores=4, control=list(adapt_delta=0.99))
saveRDS(mod3, "Brms models/spline_growth models/spline_growth3.rds")
mod4 <- brm(spline_growth ~ (spei12|Species/Region), data=datasub, family=gaussian(), cores=4, control=list(adapt_delta=0.99))
saveRDS(mod4, "Brms models/spline_growth models/spline_growth4.rds")
mod5 <- brm(spline_growth ~ (spei12*Species|Region), data=datasub, family=gaussian(), cores=4, control=list(adapt_delta=0.99))
saveRDS(mod5, "Brms models/spline_growth models/spline_growth5.rds")
mod6 <- brm(spline_growth ~ spei12 * Species + (spei12 + Species|Region), data=datasub, family=gaussian(), cores=4, control=list(adapt_delta=0.99))
saveRDS(mod6, "Brms models/spline_growth models/spline_growth6.rds")


mod7 <- brm(spline_growth ~ (spei12*Species|Region/hilo), data=datasub, family=gaussian(), cores=4, control=list(adapt_delta=0.99))
saveRDS(mod7, "Brms models/spline_growth models/spline_growth7.rds")

mod8 <- brm(spline_growth ~ (spei12|Species/Region/hilo), data=datasub, family=gaussian(), cores=4, control=list(adapt_delta=0.99))
saveRDS(mod8, "Brms models/spline_growth models/spline_growth8.rds")

mod9 <- brm(spline_growth ~ (spei12*Species|Region/hilo) + (1|tree.uniqueID), data=datasub, family=gaussian(), cores=4, control=list(adapt_delta=0.99))
saveRDS(mod9, "Brms models/spline_growth models/spline_growth9.rds")

mod10 <- brm(spline_growth ~ (spei12*Species|Site) + (1|tree.uniqueID), data=datasub, family=gaussian(), cores=4, control=list(adapt_delta=0.99))
saveRDS(mod10, "Brms models/spline_growth models/spline_growth10.rds")

mod11 <- brm(spline_growth ~ (spei12*Species|unique.nb) + (1|tree.uniqueID), data=datasub, family=gaussian(), cores=4, control=list(adapt_delta=0.99))
saveRDS(mod11, "Brms models/spline_growth models/spline_growth11.rds")

mod12 <- brm(spline_growth ~ (spei12*Species|unique.nb) + (1|tree.uniqueID), autocor = cor_ar(formula = ~Year|tree.uniqueID, p = 1, cov = FALSE), data=datasub, family=gaussian(), cores=4, control=list(adapt_delta=0.99))
saveRDS(mod12, "Brms models/spline_growth models/spline_growth12.rds")


model_weights(mod11, mod12, weights = "waic")
bayes_R2(mod12)
bayes_R2(mod11)
bayes_R2(mod10)
sbayes_R2(mod9)
bayes_R2(mod7)


for(i in 0:10){
  assign(paste("mod", i, sep=""), readRDS(paste("Brms models/spline_growth models/spline_growth", i, ".rds", sep="")), envir = .GlobalEnv)
}

## alldata run----
mod11a <- brm(spline_growth ~ (spei12*Species|unique.nb) + (1|tree.uniqueID), data=newdata, family=gaussian(), cores=4, control=list(adapt_delta=0.99))
saveRDS(mod11a, "Brms models/spline_growth models/spline_growth11_alldata.rds")

mod12a <- brm(spline_growth ~ (spei12*Species|unique.nb) + (1|tree.uniqueID), autocor = cor_ar(formula = ~Year|tree.uniqueID, p = 1, cov = FALSE), data=newdata, family=gaussian(), cores=4, control=list(adapt_delta=0.99))
saveRDS(mod12a, "Brms models/spline_growth models/spline_growth12_alldata.rds")

## Predicted data
dfx_data_raw <- read.csv("Processed Data/dfx_yearly.csv")

dfx_data <- filter(dfx_data_raw, zeroRing=="no", Species!="AM", Site!="Mammoth") 
dfx_data$Species[str_which(dfx_data$Species, "PP")] <- "PJ" 
dfx_data$Species <- droplevels(dfx_data$Species)

for(i in unique(dfx_data$tree.uniqueID)){ # unique(dfx_data$tree.uniqueID)
  for(j in unique(dfx_data$Year[which(dfx_data$tree.uniqueID==i)])) {
    newdata <- filter(dfx_data, tree.uniqueID==i&Year==j) %>% 
      select(c(spei_post1, Species, Region, hilo, tree.uniqueID, post1)) %>% 
      rename(spei12=spei_post1)
    
    pred_grow <- fitted(mod8, newdata=newdata, summary=FALSE, nsamples = 1000)
    newy_df <- data.frame(pred_grow) %>% 
      mutate(tree.uniqueID=i,
             Year=j,
             post1=newdata$post1,
             diff=pred_grow-post1)
    if(i==first(unique(datasub$tree.uniqueID))&j==min(dfx_data$Year[which(dfx_data$tree.uniqueID==i)])) {
      big_df <- newy_df
    } else {big_df <- bind_rows(big_df, newy_df)}
  }
}

for(i in unique(dfx_data$tree.uniqueID)){ # unique(dfx_data$tree.uniqueID)
  for(j in unique(dfx_data$Year[which(dfx_data$tree.uniqueID==i)])) {
    newdata <- filter(dfx_data, tree.uniqueID==i&Year==j) %>% 
      select(c(spei_post2, Species, Region, hilo, tree.uniqueID, post2)) %>% 
      rename(spei12=spei_post2)
    
    pred_grow <- fitted(mod8, newdata=newdata, summary=FALSE, nsamples = 1000)
    newy_df <- data.frame(pred_grow) %>% 
      mutate(tree.uniqueID=i,
             Year=j,
             post2=newdata$post2,
             diff=pred_grow-post2)
    if(i==first(unique(datasub$tree.uniqueID))&j==min(dfx_data$Year[which(dfx_data$tree.uniqueID==i)])) {
      big_df2 <- newy_df
    } else {big_df2 <- bind_rows(big_df2, newy_df)}
  }
}

legacy <- big_df %>% 
  group_by(tree.uniqueID, Year) %>% 
  summarise(pred_mean = mean(pred_grow, na.rm=T),
            pred_50 = quantile(pred_grow, .5, na.rm=T),
            pred_025 = quantile(pred_grow, .025, na.rm=T),
            pred_975 = quantile(pred_grow, .975, na.rm=T),
            leg1_mean = mean(diff, na.rm=T),
            leg1_50 = quantile(diff, .5, na.rm=T),
            leg1_025 = quantile(diff, .025, na.rm=T),
            leg1_975 = quantile(diff, .975, na.rm=T)) %>% 
  left_join(select(dfx_data, c(Year, Site, Species, Neighborhood, Region, hilo, tree.uniqueID, hegyi, DBH, reduction,
                               spei12, post1, post2, spei_post1, spei_post2, Year_post1, Year_post2, drought_fx)))

legacy2 <- big_df2 %>% 
  group_by(tree.uniqueID, Year) %>% 
  summarise(pred_mean = mean(pred_grow, na.rm=T),
            pred_50 = quantile(pred_grow, .5, na.rm=T),
            pred_025 = quantile(pred_grow, .025, na.rm=T),
            pred_975 = quantile(pred_grow, .975, na.rm=T),
            leg2_mean = mean(diff, na.rm=T),
            leg2_50 = quantile(diff, .5, na.rm=T),
            leg2_025 = quantile(diff, .025, na.rm=T),
            leg2_975 = quantile(diff, .975, na.rm=T)) %>% 
  select(-c(pred_mean, pred_50, pred_025, pred_975)) %>% 
  left_join(legacy)
# write.csv(legacy2, "Processed Data/legacy.csv", row.names = F)
