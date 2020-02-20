---
title: "brms Model taxonomy"
author: "Teresa Bohner"
date: "11/22/2019"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Brms models

This document details the difference in parameter estimates, model convergence, and speed when parameters are included in the population or group level portions of the model (inside or out of the parentheses).

### Start with spei and species only.
Compare varying slope model with simple interaction

```{r interaction, eval=F}
mod1 <- brm(spline_growth ~ Species + (spei12|Species), data=datasub, family=gaussian(), cores=4)
mod2 <- brm(spline_growth ~ spei12 * Species, data=datasub, family=gaussian(), cores=4)
```

```{r interaction plot, echo=F}
mod1 <- readRDS("Brms models/test models/spline_growth1")
mod2 <- readRDS("Brms models/test models/spline_growth2")
param1 <- coef(mod1, summary=F)
int1 <- data.frame(param1$Species[,,"Intercept"]) %>% 
  mutate(AM = AM + param1$Species[,"AM","SpeciesAM"],
         PJ = PJ + param1$Species[,"PJ","SpeciesPJ"],
         PL = PL + param1$Species[,"PL","SpeciesPL"],
         PP = PP + param1$Species[,"PP","SpeciesPP"]) %>% 
  pivot_longer(everything(), names_to="Species", values_to="Sample") %>% 
  group_by(Species) %>% 
  summarise(mean=mean(Sample),
            med=median(Sample),
            upper=quantile(Sample, 0.975),
            lower=quantile(Sample, 0.025),
            mod="mod1")

param2 <- fixef(mod2, summary=F)
int2 <- data.frame(param2[,1]) %>% 
  rename(AC = param2...1.) %>% 
  mutate(AM = AC + param2[,3],
         PJ = AC + param2[,4],
         PL = AC + param2[,5],
         PP = AC + param2[,6]) %>% 
  pivot_longer(everything(), names_to="Species", values_to="Sample") %>% 
  group_by(Species) %>% 
  summarise(mean=mean(Sample),
            med=median(Sample),
            upper=quantile(Sample, 0.975),
            lower=quantile(Sample, 0.025),
            mod="mod2")

int <- bind_rows(int1, int2)
ggplot(int, aes(Species, mean, ymin=lower, ymax=upper, group=mod, color=mod)) +
  geom_pointrange(position = position_dodge(width=0.5)) +
  ggtitile("Intercept comparison")

slope1 <- data.frame(param1$Species[,,"spei12"]) %>% 
  pivot_longer(everything(), names_to="Species", values_to="Sample") %>% 
  group_by(Species) %>% 
  summarise(mean=mean(Sample),
            med=median(Sample),
            upper=quantile(Sample, 0.975),
            lower=quantile(Sample, 0.025),
            mod="mod1")

slope2 <- data.frame(param2[,2]) %>% 
  rename(AC = param2...2.) %>% 
  mutate(AM = AC + param2[,7],
         PJ = AC + param2[,1],
         PL = AC + param2[,9],
         PP = AC + param2[,10]) %>% 
  pivot_longer(everything(), names_to="Species", values_to="Sample") %>% 
  group_by(Species) %>% 
  summarise(mean=mean(Sample),
            med=median(Sample),
            upper=quantile(Sample, 0.975),
            lower=quantile(Sample, 0.025),
            mod="mod2")

slope <- bind_rows(slope1, slope2)
ggplot(slope, aes(Species, mean, ymin=lower, ymax=upper, group=mod, color=mod)) +
  geom_pointrange(position = position_dodge(width=0.5)) +
  ggtitle("Slope comparison")
```


