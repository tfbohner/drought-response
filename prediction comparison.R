mod <- readRDS("Brms models/spline_growth models/spline_growth11_alldata.rds")

tree <- read.csv("Processed Data/dfx_allyears.csv") %>% 
  filter(tree.uniqueID=="9pj1_1") %>% 
  mutate(unique.nb = paste0(Site, Neighborhood, sep=""),
         Species=droplevels(Species))

ggplot(tree, aes(Year, spline_growth)) +
  geom_line()
  

for(i in unique(tree$Year)){
  newdata <- filter(tree, Year==i) %>% 
    dplyr::select(c(spei12, Species, unique.nb, tree.uniqueID, Year))
  pred_grow <- predict(mod, newdata=newdata, summary=FALSE, nsamples = 1000)
  newy_df <- data.frame(pred_grow) %>%
    mutate(tree.uniqueID="9pj1_1",
           Year=i)
  if(i==min(tree$Year)) {
    big_df <- newy_df 
  } else {big_df <- bind_rows(big_df, newy_df)}
  
}


concept_dat <- big_df %>%
  group_by(Year) %>%
  summarise(pred_mean = mean(pred_grow, na.rm=T),
            pred_50 = quantile(pred_grow, .5, na.rm=T),
            pred_025 = quantile(pred_grow, .025, na.rm=T),
            pred_975 = quantile(pred_grow, .975, na.rm=T)) %>%
  left_join(dplyr::select(tree, c(Year, Site, Species, Neighborhood, Region, hilo, tree.uniqueID, hegyi, simple_comp, BA, DBH, avg_growth, drought_fx, droughtyear, spline_growth, spei12)))

mod <- readRDS("Brms models/spline_growth models/spline_growth12_alldata.rds")
for(i in unique(tree$Year)){
  newdata <- filter(tree, Year==i) %>% 
    dplyr::select(c(spei12, Species, unique.nb, tree.uniqueID, Year))
  fit_grow <- predict(mod, newdata=newdata, summary=FALSE, nsamples = 1000)
  newy_df <- data.frame(fit_grow) %>%
    mutate(tree.uniqueID="9pj1_1",
           Year=i)
  if(i==min(tree$Year)) {
    big_df <- newy_df 
  } else {big_df <- bind_rows(big_df, newy_df)}
  
}


concept_dat2 <- big_df %>%
  group_by(Year) %>%
  summarise(fit_mean = mean(fit_grow, na.rm=T),
            fit_50 = quantile(fit_grow, .5, na.rm=T),
            fit_025 = quantile(fit_grow, .025, na.rm=T),
            fit_975 = quantile(fit_grow, .975, na.rm=T)) %>%
  left_join(dplyr::select(tree, c(Year, Site, Species, Neighborhood, Region, hilo, tree.uniqueID, hegyi, simple_comp, BA, DBH, avg_growth, drought_fx, droughtyear, spline_growth, spei12)))

concept_dat3 <- left_join(concept_dat, concept_dat2)
ggplot(concept_dat3, aes(Year, spline_growth)) +
  geom_line() +
  geom_ribbon(aes(Year, ymin=pred_025, ymax=pred_975), alpha=0.3) +
  geom_line(aes(Year, pred_50), color="red") +
  geom_ribbon(aes(Year, ymin=fit_025, ymax=fit_975), alpha=0.3) +
  geom_line(aes(Year, fit_50), color="blue") 

ggplot(concept_dat3, aes(spei12, spline_growth)) +
  geom_point() +
  geom_smooth(method="lm") +
  geom_point(aes(spei12, pred_50), color="red") +
  geom_point(aes(spei12, fit_50), color="blue")


ggplot(concept_dat, aes(Year, spline_growth)) +
  geom_line() +
  geom_ribbon(aes(Year, ymin=pred_025, ymax=pred_975), alpha=0.3) +
  geom_line(aes(Year, pred_50), color="red") 

ggplot(concept_dat, aes(spei12, spline_growth)) +
  geom_point() +
  geom_smooth(method="lm") +
  geom_point(aes(spei12, pred_50), color="red")

ggplot(concept_dat2, aes(Year, spline_growth)) +
  geom_line() +
  geom_ribbon(aes(Year, ymin=fit_025, ymax=fit_975), alpha=0.3) +
  geom_line(aes(Year, fit_50), color="red") 

ggplot(concept_dat2, aes(spei12, spline_growth)) +
  geom_point() +
  geom_smooth(method="lm") +
  geom_point(aes(spei12, fit_50), color="red")

plot(concept_dat2$fit_50, concept_dat$pred_50)
abline(0,1)

plot(concept_dat$pred_mean, concept_dat$pred_50)
abline(0,1)


plot(concept_dat2$fit_mean, concept_dat2$fit_50)
abline(0,1)

ggplot(concept_dat, aes(spei12, fit_mean)) +
  geom_point() +
  geom_smooth(method="lm") +
  geom_point(aes(spei12, pred_mean), color="red") +
  geom_point(aes(spei12, spline_growth), color="gray")

ggplot(concept_dat, aes(pred_mean, spline_growth)) +
  geom_point() +
  geom_smooth(method='lm')

ggplot(concept_dat, aes(fit_mean, spline_growth)) +
  geom_point() +
  geom_smooth(method='lm') +
  xlim(c(0.98, 0.99))

