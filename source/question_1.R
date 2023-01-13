rm(list = ls(all = TRUE)) 

library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(car)
library(lme4)
library(lmerTest)


# Question #1 -------------------------------------------------------------

### (1) Do key IAA pathway gene insertion mutant lines vary in their response to exogenous
### IAA treatments? We predicted that lines with inserts in IAA pathway genes will vary in 
### their response to exogenous IAA treatments showing differences in below ground and above
### ground phenotypes paired with differences in plant fitness and fruit production.

df <- read.csv("./data/combined_data_clean.csv") 
df$mutant <- ifelse(df$insert.location == "phytometer", "N", "Y")
df$ttl.biomass <- df$above.ground.mass.mg + df$below.ground.mass.mg

### create df_SALK df with just mutants

df_SALK <- df %>% filter(mutant == "Y")%>%
  select(Stock_number,treatment, tray:column, collection.time, rosette.diameter, germination.day, Length.cm., primary.length:lower.LR.length,
         LR.count:above.ground.mass.mg, day.to.bolt:day.to.flower, day.to.mature, basal.fruit.length:flw.buds, ttl.biomass)

# Analysis of EARLY stage -------------------------------------------------

early <- df_SALK %>% filter(treatment != "promix") %>% filter(collection.time == "14d") %>%
  select(Stock_number, treatment, tray:column, rosette.diameter:above.ground.mass.mg,ttl.biomass) %>% 
  select(-LR.count) %>% select(-LR.density) %>% select(-ttl.LR.count)

early$ttl.LR <- early$lower.LR.count + early$mid.LR.count + early$upper.LR.count
early$LR.density <- early$ttl.LR / early$primary.length

# complete transformations
 
early$upper.LR.count <- sqrt(early$upper.LR.count) + .5
early$lower.LR.count <- sqrt(early$lower.LR.count) + .5
early$lower.LR.length <- sqrt(early$lower.LR.length) + .5
early$upper.LR.length <- sqrt(early$upper.LR.length) + .5
early$LR.density <- sqrt(early$LR.density) + .5
early$upper.LR.density <- sqrt(early$upper.LR.density) + .5
early$mid.LR.density <- sqrt(early$mid.LR.density) + .5
early$lower.LR.density <- sqrt(early$lower.LR.density) + .5
early$below.ground.mass.mg <- sqrt(early$below.ground.mass.mg) + 0.5
early$above.ground.mass.mg <- log10(early$above.ground.mass.mg) + 1


# run early models

pdf(file = './figs/question_1_early_models.pdf')
for(i in 6:ncol(early)) {
  column <- names(early[i])
  mod <- lm(early[,i] ~ treatment + Stock_number*treatment, data = early)
  print(column)
  result_fixed <- anova(mod)
  print(result_fixed)
  plot(fitted(mod),resid(mod))
  mtext(names(early[i]))
  mtext('lmer(early[,i]~treatment+Stock_number*treatment',
        side = 1)
  abline(h = 0)
  qqnorm(resid(mod))
  qqline(resid(mod))
  mtext(names(early[i]))
  mtext('lmer(early[,i]~treatment+Stock_number*treatment',
        side = 1)
}
dev.off()

# Analysis of LATE stage --------------------------------------------------

late <- df_SALK %>% filter(treatment != "promix") %>% filter(collection.time == "21d") %>% 
  select(Stock_number, treatment, tray:column, rosette.diameter:above.ground.mass.mg, ttl.biomass) %>% 
  select(-LR.count, -upper.LR.length,-upper.LR.count, -ttl.LR.count, -upper.LR.density)

### add LR.count column

late$LR.count <- late$mid.LR.count + late$lower.LR.count

# complete transformations

late$rosette.diameter <- log10(late$rosette.diameter) + 1
late$primary.length <- log10(late$primary.length) + 1
late$mid.LR.count <- sqrt(late$mid.LR.count) + .5
late$lower.LR.count <- sqrt(late$lower.LR.count) + .5
late$lower.LR.length <- sqrt(late$lower.LR.length) + .5
late$mid.LR.length <- sqrt(late$mid.LR.length) + .5
late$LR.density <- sqrt(late$LR.density) + 0.5
late$mid.LR.density <- sqrt(late$mid.LR.density) + 0.5
late$LR.count <- sqrt(late$LR.count) + .5
late$below.ground.mass.mg <- sqrt(late$below.ground.mass.mg) + .5
late$above.ground.mass.mg <- sqrt(late$above.ground.mass.mg) + .5
late$ttl.biomass <- log10(late$ttl.biomass) + 1

# run late models

pdf(file = './figs/question_1_late_models.pdf')
for(i in 6:ncol(late)) {
  column <- names(late[i])
  mod <- lm(late[,i] ~ treatment + Stock_number*treatment, data = late)
  print(column)
  result_fixed <- anova(mod)
  print(result_fixed)
  plot(fitted(mod),resid(mod))
  mtext(names(late[i]))
  mtext('lm(late[,i]~treatment+Stock_number*treatment',
        side = 1)
  abline(h = 0)
  qqnorm(resid(mod))
  qqline(resid(mod))
  mtext(names(late[i]))
  mtext('lm(late[,i]~treatment+Stock_number*treatment',
        side = 1)
}
dev.off()

# Analysis of MATURE stage ------------------------------------------------

mature <- df_SALK %>% filter(treatment != "promix") %>% filter(collection.time == "fruit") %>%
  select(Stock_number, treatment, tray:column, germination.day:Length.cm., below.ground.mass.mg:ttl.biomass)

mature$avg.fruit.length = (mature$basal.fruit.length + mature$mid.fruit.length + mature$upper.fruit.length) / 3
mature$fitness = mature$avg.fruit.length * mature$fruit.num

# complete transformations

mature$day.to.bolt <- sqrt(mature$day.to.bolt) + .5
mature$day.to.flower <- sqrt(mature$day.to.flower) + .5
mature$basal.fruit.length <- log10(mature$basal.fruit.length) + 1
mature$mid.fruit.length <- log10(mature$mid.fruit.length) + 1
mature$upper.fruit.length <- log10(mature$upper.fruit.length) + 1
mature$inflorescence.height <- log10(mature$inflorescence.height) + 1
mature$basal.branch <- sqrt(mature$basal.branch) + .5
mature$ttl.maininfl.branch <- sqrt(mature$ttl.maininfl.branch) + .5
mature$branch.basalbranch <- sqrt(mature$branch.basalbranch) + .5
mature$days.to.bolt <- sqrt(mature$days.to.bolt) + .5
mature$days.to.flower <- sqrt(mature$days.to.flower) + .5
mature$ttl.maininfl.branch <- sqrt(mature$ttl.maininfl.branch) + 0.5
mature$fruit.num <- sqrt(mature$fruit.num) + .5
mature$flw.num <- sqrt(mature$flw.num) +.5
mature$aborted.fruits <- sqrt(mature$aborted.fruits) + .5
mature$flw.buds <- sqrt(mature$flw.buds) + .5

# run mature models

pdf(file = './figs/question_1_mature_models.pdf')
for(i in 6:ncol(mature)) {
  column <- names(mature[i])
  mod <- lmer(mature[,i] ~ treatment + Stock_number*treatment + (1|tray), data = mature)
  print(column)
  result_fixed <- anova(mod)
  result_rand <- ranova(mod)
  print(result_fixed)
  print(result_rand)
  plot(fitted(mod),resid(mod))
  mtext(names(mature[i]))
  mtext('lmer(mature[,i]~treatment+Stock_number*treatment+(1|tray))',
        side = 1)
  abline(h = 0)
  qqnorm(resid(mod))
  qqline(resid(mod))
  mtext(names(mature[i]))
  mtext('lmer(mature[,i]~treatment+Stock_number*treatment+(1|tray))',
        side = 1)
}
dev.off()

