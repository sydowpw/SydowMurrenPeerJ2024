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

# Define additional traits

df$mutant <- ifelse(df$insert.location == "phytometer", "N", "Y")
df$LR.density <- ifelse(df$collection.time == "14d", (df$upper.LR.count + df$mid.LR.count + df$lower.LR.count) / df$primary.length, 
                        ifelse(df$collection.time =="21d", (df$mid.LR.count + df$lower.LR.count) / ((2/3)*df$primary.length), NA))
df$avg.fruit.length <- ifelse(df$collection.time == "fruit", (df$basal.fruit.length + df$mid.fruit.length + df$upper.fruit.length) / 3, NA)
df$fitness <- ifelse(df$collection.time == 'fruit', df$avg.fruit.length * df$fruit.num, NA)

### create df_SALK df with just mutants

df_SALK <- df %>% filter(mutant == "Y") %>%
  select(Stock_number,treatment, tray:column, collection.time, rosette.diameter, LR.density, diameter.at.bolt, 
         Length.cm., below.ground.mass.mg, above.ground.mass.mg, inflorescence.height, ttl.branch, 
         fruit.num, avg.fruit.length, fitness, germination.day, day.to.bolt, day.to.flower, days.to.collect, ttl.branch)

# Analysis of EARLY stage -------------------------------------------------

early <- df_SALK %>% filter(treatment != "promix") %>% filter(collection.time == "14d") %>%
  select(Stock_number, treatment, tray:column, rosette.diameter, Length.cm., below.ground.mass.mg, 
         above.ground.mass.mg, LR.density)

# complete transformations
 
early$LR.density <- sqrt(early$LR.density) + .5
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
  select(Stock_number, treatment, tray:column, rosette.diameter, Length.cm., below.ground.mass.mg, 
         above.ground.mass.mg, LR.density)

# complete transformations

late$rosette.diameter <- log10(late$rosette.diameter) + 1
late$LR.density <- sqrt(late$LR.density) + 0.5
late$below.ground.mass.mg <- sqrt(late$below.ground.mass.mg) + .5
late$above.ground.mass.mg <- sqrt(late$above.ground.mass.mg) + .5

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
  select(Stock_number, treatment, tray:column, diameter.at.bolt, Length.cm., 
         below.ground.mass.mg, above.ground.mass.mg, inflorescence.height,
         ttl.branch, fruit.num, avg.fruit.length, fitness, germination.day, day.to.bolt, 
         day.to.flower, days.to.collect, ttl.branch)

# complete transformations

mature$day.to.bolt <- sqrt(mature$day.to.bolt) + .5
mature$day.to.flower <- sqrt(mature$day.to.flower) + .5
mature$inflorescence.height <- log10(mature$inflorescence.height) + 1
mature$fruit.num <- sqrt(mature$fruit.num) + .5

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

