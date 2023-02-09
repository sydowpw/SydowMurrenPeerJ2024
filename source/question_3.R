rm(list = ls(all = TRUE)) 

library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(car)
library(lme4)
library(lmerTest)


# Question #2 -------------------------------------------------------------

### Does the root phenotype response of mutant lines of IAA genes correspond with 
### the number of alternative polyadenylation sites, and is this pattern the same as 
### seedlings and reproductive adult plants? We predicted that genes with fewer poly(A) 
### sites would function similarly across plant development with the root trait phenotype 
### of their respective insert mutants varying to a greater extent by exhibiting both fewer 
### lateral roots and less overall root length when compared to mutants with inserts in genes 
### with many poly(A) sites that may confer functional variation across the lifecycle of the plant

df <- read.csv("./data/combined_data_clean.csv") 
df$mutant <- ifelse(df$insert.location == "phytometer", "N", "Y")

# Define additional traits

df$mutant <- ifelse(df$insert.location == "phytometer", "N", "Y")
df$LR.density <- ifelse(df$collection.time == "14d", (df$mid.LR.count + df$lower.LR.count) / ((2/3)*df$primary.length),
                        ifelse(df$collection.time =="21d", (df$mid.LR.count + df$lower.LR.count) / ((2/3)*df$primary.length), NA))

### create df_SALK df with just mutants

df_SALK <- df %>% filter(mutant == "Y") %>%
  select(PA.group, Stock_number,treatment, tray:column, collection.time, rosette.diameter, LR.density, 
         Length.cm., below.ground.mass.mg, above.ground.mass.mg)

### filter out two lines with multiple inserts without PA.group

df_SALK <- df_SALK %>% filter(!is.na(PA.group))

# Scaled Development Models 14d - fruit -----------------------------------------------

### filter df to traits across all three dev. stages and variables for model

all_SALK <- df_SALK %>% filter(treatment != "promix") %>%
  select(PA.group, Stock_number, treatment, tray, collection.time,
         Length.cm., below.ground.mass.mg, above.ground.mass.mg)

### create function to scale all columns that are numeric or integer class and 
### create a new column in the df scaled_"trait"

scale_traits <- function(x){
  for(i in 6:ncol(x)) {
    if(class(x[,i]) == "numeric" | class(x[,i]) == "integer") {
      x[,ncol(x)+ 1] <- scale(x[,i])
      names(x)[ncol(x)] <- paste("scaled_", names(x)[i], sep ='')
    }
  }
  return(x)
}

### break up df and scaled traits by their collection time and treatment

early_scaled_IAA <- all_SALK %>% filter(collection.time == '14d' & treatment == "IAA") %>% scale_traits()
late_scaled_IAA <- all_SALK %>% filter(collection.time == '21d' & treatment == "IAA") %>% scale_traits()
fruit_scaled_IAA <- all_SALK %>% filter(collection.time == 'fruit' & treatment == "IAA") %>% scale_traits()

early_scaled_control <- all_SALK %>% filter(collection.time == '14d' & treatment == "control") %>% scale_traits()
late_scaled_control <- all_SALK %>% filter(collection.time == '21d' & treatment == "control") %>% scale_traits()
fruit_scaled_control <- all_SALK %>% filter(collection.time == 'fruit' & treatment == "control") %>% scale_traits()

### update tray column so that tray 1 of seedling is not considered tray 1 of mature and vice versa

fruit_scaled_IAA$tray = fruit_scaled_IAA$tray + 3
fruit_scaled_control$tray = fruit_scaled_control$tray + 3

### throw all time points back together

all_SALK_scaled_control <- early_scaled_control %>% bind_rows(late_scaled_control) %>% bind_rows(fruit_scaled_control)
all_SALK_scaled_IAA <- early_scaled_IAA %>% bind_rows(late_scaled_IAA) %>% bind_rows(fruit_scaled_IAA)


### run models on control data and save plot figs in _models.pdf

pdf(file = './figs/question_3_all_scaled_control_models.pdf')
for(i in 1:ncol(all_SALK_scaled_control)) {
  #if there is an error, the whole loop will stop
  #tryCatch attempts to prevent that
  tryCatch({
    if(grepl("scaled",names(all_SALK_scaled_control)[i], fixed=TRUE)){
      column <- names(all_SALK_scaled_control[i])
      mod <- lmer(all_SALK_scaled_control[,i] ~ PA.group + collection.time*PA.group + (1|tray), data = all_SALK_scaled_control)
      print(column)
      result_fixed <- anova(mod)
      result_rand <- ranova(mod)
      print(result_fixed)
      print(result_rand)
      plot(fitted(mod),resid(mod))
      mtext(names(all_SALK_scaled_control[i]))
      mtext('lmer(all_SALK_scaled_control[,i]~PA.group+collection.time*PA.group+(1|tray))',
            side = 1)
      abline(h = 0)
      qqnorm(resid(mod))
      qqline(resid(mod))
      mtext(names(all_SALK_scaled_control[i]))
      mtext('lmer(all_SALK_scaled_control[,i]~PA.group+collection.time*PA.group+(1|tray))',
            side = 1)
    }
  }, error=function(e){cat("error in",i,":",conditionMessage(e),"\n")})
}
dev.off()

### run models on IAA data and save plot figs in _models.pdf

pdf(file = './figs/question_3_all_scaled_IAA_models.pdf')
for(i in 6:ncol(all_SALK_scaled_IAA)) {
  #if there is an error, the whole loop will stop
  #tryCatch attempts to prevent that
  tryCatch({
    if(grepl("scaled",names(all_SALK_scaled_IAA)[i], fixed=TRUE)){
      column <- names(all_SALK_scaled_IAA[i])
      mod <- lmer(all_SALK_scaled_IAA[,i] ~ PA.group + collection.time*PA.group + (1|tray), data = all_SALK_scaled_IAA)
      print(column)
      result_fixed <- anova(mod)
      result_rand <- ranova(mod)
      print(result_fixed)
      print(result_rand)
      plot(fitted(mod),resid(mod))
      mtext(names(all_SALK_scaled_IAA[i]))
      mtext('lmer(all_SALK_scaled_IAA[,i]~PA.group+collection.time*PA.group+(1|tray))',
            side = 1)
      abline(h = 0)
      qqnorm(resid(mod))
      qqline(resid(mod))
      mtext(names(all_SALK_scaled_IAA[i]))
      mtext('lmer(all_SALK_scaled_IAA[,i]~PA.group+collection.time*PA.group+(1|tray))',
            side = 1)
    }
  }, error=function(e){cat("error in",i,":",conditionMessage(e),"\n")})
}
dev.off()

# emmeans example plot ----------------------------------------------------
library(emmeans)
mod <- lm(Length.cm. ~ collection.time + Stock_number*collection.time, data = all_SALK_NAT_scaled_control)
means <- emmeans(mod, list(pairwise ~ Stock_number|collection.time), adjust = 'tukey')
plot(means)

mod <- lm(Length.cm. ~ collection.time + Stock_number*collection.time , data = seedling_SALK_NAT_control)
aov(mod)
anova(mod)

# Emmeans work to determine which lines vary ------------------------------
options(max.print=1000000)

# All emmeans --------------------------------------------------------

scaled_IAA_all_above.ground.mass.mg  <- lmer(scaled_above.ground.mass.mg ~ collection.time + PA.group*collection.time + (1|tray), data = all_SALK_scaled_IAA)
scaled_IAA_all_above.ground.mass.mg_means <- emmeans(scaled_IAA_all_above.ground.mass.mg, list(pairwise ~ PA.group|collection.time), adjust = 'tukey')
plot(scaled_IAA_all_above.ground.mass.mg_means)
scaled_IAA_all_above.ground.mass.mg_means
scaled_IAA_all_above.ground.mass.mg_test<-as.data.frame(contrast(scaled_IAA_all_above.ground.mass.mg_means))
scaled_IAA_all_above.ground.mass.mg_test <- scaled_IAA_all_above.ground.mass.mg_test[,c(8,9,14)] %>% filter(scaled_IAA_all_above.ground.mass.mg_test[,14] <= 0.05)
scaled_IAA_all_above.ground.mass.mg_test

scaled_IAA_all_below.ground.mass.mg  <- lmer(scaled_below.ground.mass.mg ~ collection.time + PA.group*collection.time + (1|tray), data = all_SALK_scaled_IAA)
scaled_IAA_all_below.ground.mass.mg_means <- emmeans(scaled_IAA_all_below.ground.mass.mg, list(pairwise ~ PA.group|collection.time), adjust = 'tukey')
plot(scaled_IAA_all_below.ground.mass.mg_means)
scaled_IAA_all_below.ground.mass.mg_means
scaled_IAA_all_below.ground.mass.mg_test<-as.data.frame(contrast(scaled_IAA_all_below.ground.mass.mg_means))
scaled_IAA_all_below.ground.mass.mg_test <- scaled_IAA_all_below.ground.mass.mg_test[,c(8,9,14)] %>% filter(scaled_IAA_all_below.ground.mass.mg_test[,14] <= 0.05)
scaled_IAA_all_below.ground.mass.mg_test

scaled_IAA_all_Length.cm.  <- lmer(scaled_Length.cm. ~ collection.time + PA.group*collection.time + (1|tray), data = all_SALK_scaled_IAA)
scaled_IAA_all_Length.cm._means <- emmeans(scaled_IAA_all_Length.cm., list(pairwise ~ PA.group|collection.time), adjust = 'tukey')
plot(scaled_IAA_all_Length.cm._means)
scaled_IAA_all_Length.cm._means
scaled_IAA_all_Length.cm._test<-as.data.frame(contrast(scaled_IAA_all_Length.cm._means))
scaled_IAA_all_Length.cm._test <- scaled_IAA_all_Length.cm._test[,c(8,9,14)] %>% filter(scaled_IAA_all_Length.cm._test[,14] <= 0.05)
scaled_IAA_all_Length.cm._test

scaled_control_all_above.ground.mass.mg  <- lmer(scaled_above.ground.mass.mg ~ collection.time + PA.group*collection.time + (1|tray), data = all_SALK_scaled_control)
scaled_control_all_above.ground.mass.mg_means <- emmeans(scaled_control_all_above.ground.mass.mg, list(pairwise ~ PA.group|collection.time), adjust = 'tukey')
plot(scaled_control_all_above.ground.mass.mg_means)
scaled_control_all_above.ground.mass.mg_means
scaled_control_all_above.ground.mass.mg_test<-as.data.frame(contrast(scaled_control_all_above.ground.mass.mg_means))
scaled_control_all_above.ground.mass.mg_test <- scaled_control_all_above.ground.mass.mg_test[,c(8,9,14)] %>% filter(scaled_control_all_above.ground.mass.mg_test[,14] <= 0.05)
scaled_control_all_above.ground.mass.mg_test

scaled_control_all_below.ground.mass.mg  <- lmer(scaled_below.ground.mass.mg ~ collection.time + PA.group*collection.time + (1|tray), data = all_SALK_scaled_control)
scaled_control_all_below.ground.mass.mg_means <- emmeans(scaled_control_all_below.ground.mass.mg, list(pairwise ~ PA.group|collection.time), adjust = 'tukey')
plot(scaled_control_all_below.ground.mass.mg_means)
scaled_control_all_below.ground.mass.mg_means
scaled_control_all_below.ground.mass.mg_test<-as.data.frame(contrast(scaled_control_all_below.ground.mass.mg_means))
scaled_control_all_below.ground.mass.mg_test <- scaled_control_all_below.ground.mass.mg_test[,c(8,9,14)] %>% filter(scaled_control_all_below.ground.mass.mg_test[,14] <= 0.05)
scaled_control_all_below.ground.mass.mg_test

scaled_control_all_Length.cm.  <- lmer(scaled_Length.cm. ~ collection.time + PA.group*collection.time + (1|tray), data = all_SALK_scaled_control)
scaled_control_all_Length.cm._means <- emmeans(scaled_control_all_Length.cm., list(pairwise ~ PA.group|collection.time), adjust = 'tukey')
plot(scaled_control_all_Length.cm._means)
scaled_control_all_Length.cm._means
scaled_control_all_Length.cm._test<-as.data.frame(contrast(scaled_control_all_Length.cm._means))
scaled_control_all_Length.cm._test <- scaled_control_all_Length.cm._test[,c(8,9,14)] %>% filter(scaled_control_all_Length.cm._test[,14] <= 0.05)
scaled_control_all_Length.cm._test