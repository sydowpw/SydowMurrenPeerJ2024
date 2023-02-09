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

### (1) Do differences in root traits of mutant lines and natural accessions
### change across developmental stages? 

df <- read.csv("./data/combined_data_clean.csv") 
df$mutant <- ifelse(df$insert.location == "phytometer", "N", "Y")

# Define additional traits

df$mutant <- ifelse(df$insert.location == "phytometer", "N", "Y")
df$LR.density <- ifelse(df$collection.time == "14d", (df$mid.LR.count + df$lower.LR.count) / ((2/3)*df$primary.length),
                        ifelse(df$collection.time =="21d", (df$mid.LR.count + df$lower.LR.count) / ((2/3)*df$primary.length), NA))

### create df_SALK_NAT df with all genotypes

df_SALK_NAT <- df %>%
  select(Stock_number,treatment, tray:column, collection.time, rosette.diameter, LR.density, 
         Length.cm., below.ground.mass.mg, above.ground.mass.mg)

# Scaled Development Models 14d - fruit -----------------------------------------------

### filter df to traits across all three dev. stages and variables for model

all_SALK_NAT <- df_SALK_NAT %>% filter(treatment != "promix") %>%
  select(Stock_number, treatment, tray, Length.cm., below.ground.mass.mg, above.ground.mass.mg, collection.time)

### create function to scale all columns that are numeric or integer class and 
### create a new column in the df scaled_"trait"

scale_traits <- function(x){
  for(i in 4:ncol(x)) {
    if(class(x[,i]) == "numeric" | class(x[,i]) == "integer") {
      x[,ncol(x)+ 1] <- scale(x[,i])
      names(x)[ncol(x)] <- paste("scaled_", names(x)[i], sep ='')
    }
  }
  return(x)
}

### break up df and scaled traits by their collection time and treatment

early_scaled_IAA <- all_SALK_NAT %>% filter(collection.time == '14d' & treatment == "IAA") %>% scale_traits()
late_scaled_IAA <- all_SALK_NAT %>% filter(collection.time == '21d' & treatment == "IAA") %>% scale_traits()
fruit_scaled_IAA <- all_SALK_NAT %>% filter(collection.time == 'fruit' & treatment == "IAA") %>% scale_traits()

early_scaled_control <- all_SALK_NAT %>% filter(collection.time == '14d' & treatment == "control") %>% scale_traits()
late_scaled_control <- all_SALK_NAT %>% filter(collection.time == '21d' & treatment == "control") %>% scale_traits()
fruit_scaled_control <- all_SALK_NAT %>% filter(collection.time == 'fruit' & treatment == "control") %>% scale_traits()

### update tray column so that tray 1 of seedling is not considered tray 1 of mature and vice versa

fruit_scaled_IAA$tray = fruit_scaled_IAA$tray + 3
fruit_scaled_control$tray = fruit_scaled_control$tray + 3

### throw all time points back together

all_SALK_NAT_scaled_control <- early_scaled_control %>% bind_rows(late_scaled_control) %>% bind_rows(fruit_scaled_control)
all_SALK_NAT_scaled_IAA <- early_scaled_IAA %>% bind_rows(late_scaled_IAA) %>% bind_rows(fruit_scaled_IAA)

### run models on control data and save plot figs in _models.pdf

pdf(file = './figs/question_2_all_scaled_control_models.pdf')
for(i in 1:ncol(all_SALK_NAT_scaled_control)) {
  #if there is an error, the whole loop will stop
  #tryCatch attempts to prevent that
  tryCatch({
    if(grepl("scaled",names(all_SALK_NAT_scaled_control)[i], fixed=TRUE)){
      column <- names(all_SALK_NAT_scaled_control[i])
      mod <- lmer(all_SALK_NAT_scaled_control[,i] ~ collection.time + Stock_number*collection.time + (1|tray), data = all_SALK_NAT_scaled_control)
      print(column)
      result_fixed <- anova(mod)
      result_rand <- ranova(mod)
      print(result_fixed)
      print(result_rand)
      plot(fitted(mod),resid(mod))
      mtext(names(all_SALK_NAT_scaled_control[i]))
      mtext('lmer(all_SALK_NAT_scaled_control[,i]~collection.time+Stock_number*collection.time+(1|tray))',
            side = 1)
      abline(h = 0)
      qqnorm(resid(mod))
      qqline(resid(mod))
      mtext(names(all_SALK_NAT_scaled_control[i]))
      mtext('lmer(all_SALK_NAT_scaled_control[,i]~collection.time+Stock_number*collection.time+(1|tray))',
            side = 1)
  }
 }, error=function(e){cat("error in",i,":",conditionMessage(e),"\n")})
}
dev.off()

### run models on IAA data and save plot figs in _models.pdf

pdf(file = './figs/question_2_all_scaled_IAA_models.pdf')
for(i in 4:ncol(all_SALK_NAT_scaled_IAA)) {
  #if there is an error, the whole loop will stop
  #tryCatch attempts to prevent that
  tryCatch({
    if(grepl("scaled",names(all_SALK_NAT_scaled_IAA)[i], fixed=TRUE)){
      column <- names(all_SALK_NAT_scaled_IAA[i])
      mod <- lmer(all_SALK_NAT_scaled_IAA[,i] ~ collection.time + Stock_number*collection.time + (1|tray), data = all_SALK_NAT_scaled_IAA)
      print(column)
      result_fixed <- anova(mod)
      result_rand <- ranova(mod)
      print(result_fixed)
      print(result_rand)
      plot(fitted(mod),resid(mod))
      mtext(names(all_SALK_NAT_scaled_IAA[i]))
      mtext('lmer(all_SALK_NAT_scaled_IAA[,i]~collection.time+Stock_number*collection.time+(1|tray))',
            side = 1)
      abline(h = 0)
      qqnorm(resid(mod))
      qqline(resid(mod))
      mtext(names(all_SALK_NAT_scaled_IAA[i]))
      mtext('lmer(all_SALK_NAT_scaled_IAA[,i]~collection.time+Stock_number*collection.time+(1|tray))',
            side = 1)
    }
  }, error=function(e){cat("error in",i,":",conditionMessage(e),"\n")})
}
dev.off()

# Scaled Development Models 14d - 21d -------------------------------------

# select traits shared across 14d and 21d plants

seedling_SALK_NAT <- df_SALK_NAT %>% filter(treatment != "promix") %>%
  select(Stock_number:tray, collection.time, rosette.diameter, Length.cm., below.ground.mass.mg,
         above.ground.mass.mg, LR.density) %>% filter(collection.time == "14d" | collection.time == "21d")

# break up df and scaled traits by their collection time

early_seedling_SALK_NAT_IAA <- seedling_SALK_NAT %>% filter(collection.time == "14d" & treatment == "IAA") %>% scale_traits()
early_seedling_SALK_NAT_control <- seedling_SALK_NAT %>% filter(collection.time == "14d" & treatment == "control") %>% scale_traits()

late_seedling_SALK_NAT_IAA <- seedling_SALK_NAT %>% filter(collection.time == "21d" & treatment == "IAA") %>% scale_traits()
late_seedling_SALK_NAT_control <- seedling_SALK_NAT %>% filter(collection.time == "21d" & treatment == "control") %>% scale_traits()

# throw timepoints together 

seedling_SALK_NAT_scaled_IAA <- late_seedling_SALK_NAT_IAA %>% bind_rows(early_seedling_SALK_NAT_IAA)
seedling_SALK_NAT_scaled_control <- late_seedling_SALK_NAT_control %>% bind_rows(early_seedling_SALK_NAT_control)

### run models on control data and save plot figs in _models.pdf

pdf(file = './figs/question_2_seedling_scaled_control_models.pdf')
for(i in 1:ncol(seedling_SALK_NAT_scaled_control)) {
  #if there is an error, the whole loop will stop
  #tryCatch attempts to prevent that
  tryCatch({
    if(grepl("scaled",names(seedling_SALK_NAT_scaled_control)[i], fixed=TRUE)){
      column <- names(seedling_SALK_NAT_scaled_control[i])
      mod <- lm(seedling_SALK_NAT_scaled_control[,i] ~ collection.time + Stock_number*collection.time, data = seedling_SALK_NAT_scaled_control)
      print(column)
      result_fixed <- anova(mod)
      print(result_fixed)
      plot(fitted(mod),resid(mod))
      mtext(names(seedling_SALK_NAT_scaled_control[i]))
      mtext('lm(seedling_SALK_NAT_scaled_control[,i]~collection.time+Stock_number*collection.time)',
            side = 1)
      abline(h = 0)
      qqnorm(resid(mod))
      qqline(resid(mod))
      mtext(names(seedling_SALK_NAT_scaled_control[i]))
      mtext('lm(seedling_SALK_NAT_scaled_control[,i]~collection.time+Stock_number*collection.time',
            side = 1)
    }
  }, error=function(e){cat("error in",i,":",conditionMessage(e),"\n")})
}
dev.off()

### run models on control data and save plot figs in _models.pdf

pdf(file = './figs/question_2_seedling_scaled_control_models.pdf')
for(i in 1:ncol(seedling_SALK_NAT_scaled_control)) {
  #if there is an error, the whole loop will stop
  #tryCatch attempts to prevent that
  tryCatch({
    if(grepl("scaled",names(seedling_SALK_NAT_scaled_control)[i], fixed=TRUE)){
      column <- names(seedling_SALK_NAT_scaled_control[i])
      mod <- lm(seedling_SALK_NAT_scaled_control[,i] ~ collection.time + Stock_number*collection.time, data = seedling_SALK_NAT_scaled_control)
      print(column)
      result_fixed <- anova(mod)
      print(result_fixed)
      plot(fitted(mod),resid(mod))
      mtext(names(seedling_SALK_NAT_scaled_control[i]))
      mtext('lm(seedling_SALK_NAT_scaled_control[,i]~collection.time+Stock_number*collection.time',
            side = 1)
      abline(h = 0)
      qqnorm(resid(mod))
      qqline(resid(mod))
      mtext(names(seedling_SALK_NAT_scaled_control[i]))
      mtext('lm(seedling_SALK_NAT_scaled_control[,i]~collection.time+Stock_number*collection.time',
            side = 1)
    }
  }, error=function(e){cat("error in",i,":",conditionMessage(e),"\n")})
}
dev.off()

### run models on IAA data and save plot figs in _models.pdf

pdf(file = './figs/question_2_seedling_scaled_IAA_models.pdf')
for(i in 1:ncol(seedling_SALK_NAT_scaled_IAA)) {
  #if there is an error, the whole loop will stop
  #tryCatch attempts to prevent that
  tryCatch({
    if(grepl("scaled",names(seedling_SALK_NAT_scaled_IAA)[i], fixed=TRUE)){
      column <- names(seedling_SALK_NAT_scaled_IAA[i])
      mod <- lm(seedling_SALK_NAT_scaled_IAA[,i] ~ collection.time + Stock_number*collection.time, data = seedling_SALK_NAT_scaled_IAA)
      print(column)
      result_fixed <- anova(mod)
      print(result_fixed)
      plot(fitted(mod),resid(mod))
      mtext(names(seedling_SALK_NAT_scaled_IAA[i]))
      mtext('lmer(seedling_SALK_NAT_scaled_IAA[,i]~collection.time+Stock_number*collection.time',
            side = 1)
      abline(h = 0)
      qqnorm(resid(mod))
      qqline(resid(mod))
      mtext(names(seedling_SALK_NAT_scaled_IAA[i]))
      mtext('lmer(seedling_SALK_NAT_scaled_IAA[,i]~collection.time+Stock_number*collection.time',
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

# Seedling emmeans --------------------------------------------------------

#seedling unscaled
unscaled_IAA_LR.density  <- lm(LR.density ~ collection.time + Stock_number*collection.time, data = seedling_SALK_NAT_IAA)
unscaled_IAA_LR.density_means <- emmeans(unscaled_IAA_LR.density, list(pairwise ~ Stock_number|collection.time), adjust = 'tukey')
plot(unscaled_IAA_LR.density_means)
unscaled_IAA_LR.density_means
unscaled_IAA_LR.density_test<-as.data.frame(contrast(unscaled_IAA_LR.density_means))
unscaled_IAA_LR.density_test <- unscaled_IAA_LR.density_test[,c(8,9,14)] %>% filter(unscaled_IAA_LR.density_test[,14] <= 0.05)
unscaled_IAA_LR.density_test

#seedling scaled
scaled_IAA_LR.density  <- lm(scaled_LR.density ~ collection.time + Stock_number*collection.time, data = seedling_SALK_NAT_scaled_IAA)
scaled_IAA_LR.density_means <- emmeans(scaled_IAA_LR.density, list(pairwise ~ Stock_number|collection.time), adjust = 'tukey')
plot(scaled_IAA_LR.density_means)
scaled_IAA_LR.density_means
scaled_IAA_LR.density_test<-as.data.frame(contrast(scaled_IAA_LR.density_means))
scaled_IAA_LR.density_test <- scaled_IAA_LR.density_test[,c(8,9,14)] %>% filter(scaled_IAA_LR.density_test[,14] <= 0.05)
scaled_IAA_LR.density_test

scaled_control_LR.density  <- lm(scaled_LR.density ~ collection.time + Stock_number*collection.time, data = seedling_SALK_NAT_scaled_control)
scaled_control_LR.density_means <- emmeans(scaled_control_LR.density, list(pairwise ~ Stock_number|collection.time), adjust = 'tukey')
plot(scaled_control_LR.density_means)
scaled_control_LR.density_means
scaled_control_LR.density_test<-as.data.frame(contrast(scaled_control_LR.density_means))
scaled_control_LR.density_test <- scaled_control_LR.density_test[,c(8,9,14)] %>% filter(scaled_control_LR.density_test[,14] <= 0.05)
scaled_control_LR.density_test

scaled_IAA_Length.cm.  <- lm(scaled_Length.cm. ~ collection.time + Stock_number*collection.time, data = seedling_SALK_NAT_scaled_IAA)
scaled_IAA_Length.cm._means <- emmeans(scaled_IAA_Length.cm., list(pairwise ~ Stock_number|collection.time), adjust = 'tukey')
plot(scaled_IAA_Length.cm._means)
scaled_IAA_Length.cm._means
scaled_IAA_Length.cm._test<-as.data.frame(contrast(scaled_IAA_Length.cm._means))
scaled_IAA_Length.cm._test <- scaled_IAA_Length.cm._test[,c(8,9,14)] %>% filter(scaled_IAA_Length.cm._test[,14] <= 0.05)
scaled_IAA_Length.cm._test

scaled_control_Length.cm.  <- lm(scaled_Length.cm. ~ collection.time + Stock_number*collection.time, data = seedling_SALK_NAT_scaled_control)
scaled_control_Length.cm._means <- emmeans(scaled_control_Length.cm., list(pairwise ~ Stock_number|collection.time), adjust = 'tukey')
plot(scaled_control_Length.cm._means)
scaled_control_Length.cm._means
scaled_control_Length.cm._test<-as.data.frame(contrast(scaled_control_Length.cm._means))
scaled_control_Length.cm._test <- scaled_control_Length.cm._test[,c(8,9,14)] %>% filter(scaled_control_Length.cm._test[,14] <= 0.05)
scaled_control_Length.cm._test

scaled_IAA_below.ground.mass.mg  <- lm(scaled_below.ground.mass.mg ~ collection.time + Stock_number*collection.time, data = seedling_SALK_NAT_scaled_IAA)
scaled_IAA_below.ground.mass.mg_means <- emmeans(scaled_IAA_below.ground.mass.mg, list(pairwise ~ Stock_number|collection.time), adjust = 'tukey')
plot(scaled_IAA_below.ground.mass.mg_means)
scaled_IAA_below.ground.mass.mg_means
scaled_IAA_below.ground.mass.mg_test<-as.data.frame(contrast(scaled_IAA_below.ground.mass.mg_means))
scaled_IAA_below.ground.mass.mg_test <- scaled_IAA_below.ground.mass.mg_test[,c(8,9,14)] %>% filter(scaled_IAA_below.ground.mass.mg_test[,14] <= 0.05)
scaled_IAA_below.ground.mass.mg_test

scaled_control_below.ground.mass.mg  <- lm(scaled_below.ground.mass.mg ~ collection.time + Stock_number*collection.time, data = seedling_SALK_NAT_scaled_control)
scaled_control_below.ground.mass.mg_means <- emmeans(scaled_control_below.ground.mass.mg, list(pairwise ~ Stock_number|collection.time), adjust = 'tukey')
plot(scaled_control_below.ground.mass.mg_means)
scaled_control_below.ground.mass.mg_means
scaled_control_below.ground.mass.mg_test<-as.data.frame(contrast(scaled_control_below.ground.mass.mg_means))
scaled_control_below.ground.mass.mg_test <- scaled_control_below.ground.mass.mg_test[,c(8,9,14)] %>% filter(scaled_control_below.ground.mass.mg_test[,14] <= 0.05)
scaled_control_below.ground.mass.mg_test

scaled_IAA_above.ground.mass.mg  <- lm(scaled_above.ground.mass.mg ~ collection.time + Stock_number*collection.time, data = seedling_SALK_NAT_scaled_IAA)
scaled_IAA_above.ground.mass.mg_means <- emmeans(scaled_IAA_above.ground.mass.mg, list(pairwise ~ Stock_number|collection.time), adjust = 'tukey')
plot(scaled_IAA_above.ground.mass.mg_means)
scaled_IAA_above.ground.mass.mg_means
scaled_IAA_above.ground.mass.mg_test<-as.data.frame(contrast(scaled_IAA_above.ground.mass.mg_means))
scaled_IAA_above.ground.mass.mg_test <- scaled_IAA_above.ground.mass.mg_test[,c(8,9,14)] %>% filter(scaled_IAA_above.ground.mass.mg_test[,14] <= 0.05)
scaled_IAA_above.ground.mass.mg_test

scaled_control_above.ground.mass.mg  <- lm(scaled_above.ground.mass.mg ~ collection.time + Stock_number*collection.time, data = seedling_SALK_NAT_scaled_control)
scaled_control_above.ground.mass.mg_means <- emmeans(scaled_control_above.ground.mass.mg, list(pairwise ~ Stock_number|collection.time), adjust = 'tukey')
plot(scaled_control_above.ground.mass.mg_means)
scaled_control_above.ground.mass.mg_means
scaled_control_above.ground.mass.mg_test<-as.data.frame(contrast(scaled_control_above.ground.mass.mg_means))
scaled_control_above.ground.mass.mg_test <- scaled_control_above.ground.mass.mg_test[,c(8,9,14)] %>% filter(scaled_control_above.ground.mass.mg_test[,14] <= 0.05)
scaled_control_above.ground.mass.mg_test

scaled_IAA_rosette.diameter  <- lm(scaled_rosette.diameter ~ collection.time + Stock_number*collection.time, data = seedling_SALK_NAT_scaled_IAA)
scaled_IAA_rosette.diameter_means <- emmeans(scaled_IAA_rosette.diameter, list(pairwise ~ Stock_number|collection.time), adjust = 'tukey')
plot(scaled_IAA_rosette.diameter_means)
scaled_IAA_rosette.diameter_means
scaled_IAA_rosette.diameter_test<-as.data.frame(contrast(scaled_IAA_rosette.diameter_means))
scaled_IAA_rosette.diameter_test <- scaled_IAA_rosette.diameter_test[,c(8,9,14)] %>% filter(scaled_IAA_rosette.diameter_test[,14] <= 0.05)
scaled_IAA_rosette.diameter_test

scaled_control_rosette.diameter  <- lm(scaled_rosette.diameter ~ collection.time + Stock_number*collection.time, data = seedling_SALK_NAT_scaled_control)
scaled_control_rosette.diameter_means <- emmeans(scaled_control_rosette.diameter, list(pairwise ~ Stock_number|collection.time), adjust = 'tukey')
plot(scaled_control_rosette.diameter_means)
scaled_control_rosette.diameter_means
scaled_control_rosette.diameter_test<-as.data.frame(contrast(scaled_control_rosette.diameter_means))
scaled_control_rosette.diameter_test <- scaled_control_rosette.diameter_test[,c(8,9,14)] %>% filter(scaled_control_rosette.diameter_test[,14] <= 0.05)
scaled_control_rosette.diameter_test

# All emmeans --------------------------------------------------------

scaled_IAA_all_above.ground.mass.mg  <- lmer(scaled_above.ground.mass.mg ~ collection.time + Stock_number*collection.time + (1|tray), data = all_SALK_NAT_scaled_IAA)
scaled_IAA_all_above.ground.mass.mg_means <- emmeans(scaled_IAA_all_above.ground.mass.mg, list(pairwise ~ Stock_number|collection.time), adjust = 'tukey')
plot(scaled_IAA_all_above.ground.mass.mg_means)
scaled_IAA_all_above.ground.mass.mg_means
scaled_IAA_all_above.ground.mass.mg_test<-as.data.frame(contrast(scaled_IAA_all_above.ground.mass.mg_means))
scaled_IAA_all_above.ground.mass.mg_test <- scaled_IAA_all_above.ground.mass.mg_test[,c(8,9,14)] %>% filter(scaled_IAA_all_above.ground.mass.mg_test[,14] <= 0.05)
scaled_IAA_all_above.ground.mass.mg_test

scaled_IAA_all_below.ground.mass.mg  <- lmer(scaled_below.ground.mass.mg ~ collection.time + Stock_number*collection.time + (1|tray), data = all_SALK_NAT_scaled_IAA)
scaled_IAA_all_below.ground.mass.mg_means <- emmeans(scaled_IAA_all_below.ground.mass.mg, list(pairwise ~ Stock_number|collection.time), adjust = 'tukey')
plot(scaled_IAA_all_below.ground.mass.mg_means)
scaled_IAA_all_below.ground.mass.mg_means
scaled_IAA_all_below.ground.mass.mg_test<-as.data.frame(contrast(scaled_IAA_all_below.ground.mass.mg_means))
scaled_IAA_all_below.ground.mass.mg_test <- scaled_IAA_all_below.ground.mass.mg_test[,c(8,9,14)] %>% filter(scaled_IAA_all_below.ground.mass.mg_test[,14] <= 0.05)
scaled_IAA_all_below.ground.mass.mg_test

scaled_IAA_all_Length.cm.  <- lmer(scaled_Length.cm. ~ collection.time + Stock_number*collection.time + (1|tray), data = all_SALK_NAT_scaled_IAA)
scaled_IAA_all_Length.cm._means <- emmeans(scaled_IAA_all_Length.cm., list(pairwise ~ Stock_number|collection.time), adjust = 'tukey')
plot(scaled_IAA_all_Length.cm._means)
scaled_IAA_all_Length.cm._means
scaled_IAA_all_Length.cm._test<-as.data.frame(contrast(scaled_IAA_all_Length.cm._means))
scaled_IAA_all_Length.cm._test <- scaled_IAA_all_Length.cm._test[,c(8,9,14)] %>% filter(scaled_IAA_all_Length.cm._test[,14] <= 0.05)
scaled_IAA_all_Length.cm._test

scaled_control_all_above.ground.mass.mg  <- lmer(scaled_above.ground.mass.mg ~ collection.time + Stock_number*collection.time + (1|tray), data = all_SALK_NAT_scaled_control)
scaled_control_all_above.ground.mass.mg_means <- emmeans(scaled_control_all_above.ground.mass.mg, list(pairwise ~ Stock_number|collection.time), adjust = 'tukey')
plot(scaled_control_all_above.ground.mass.mg_means)
scaled_control_all_above.ground.mass.mg_means
scaled_control_all_above.ground.mass.mg_test<-as.data.frame(contrast(scaled_control_all_above.ground.mass.mg_means))
scaled_control_all_above.ground.mass.mg_test <- scaled_control_all_above.ground.mass.mg_test[,c(8,9,14)] %>% filter(scaled_control_all_above.ground.mass.mg_test[,14] <= 0.05)
scaled_control_all_above.ground.mass.mg_test

scaled_control_all_below.ground.mass.mg  <- lmer(scaled_below.ground.mass.mg ~ collection.time + Stock_number*collection.time + (1|tray), data = all_SALK_NAT_scaled_control)
scaled_control_all_below.ground.mass.mg_means <- emmeans(scaled_control_all_below.ground.mass.mg, list(pairwise ~ Stock_number|collection.time), adjust = 'tukey')
plot(scaled_control_all_below.ground.mass.mg_means)
scaled_control_all_below.ground.mass.mg_means
scaled_control_all_below.ground.mass.mg_test<-as.data.frame(contrast(scaled_control_all_below.ground.mass.mg_means))
scaled_control_all_below.ground.mass.mg_test <- scaled_control_all_below.ground.mass.mg_test[,c(8,9,14)] %>% filter(scaled_control_all_below.ground.mass.mg_test[,14] <= 0.05)
scaled_control_all_below.ground.mass.mg_test

scaled_control_all_Length.cm.  <- lmer(scaled_Length.cm. ~ collection.time + Stock_number*collection.time + (1|tray), data = all_SALK_NAT_scaled_control)
scaled_control_all_Length.cm._means <- emmeans(scaled_control_all_Length.cm., list(pairwise ~ Stock_number|collection.time), adjust = 'tukey')
plot(scaled_control_all_Length.cm._means)
scaled_control_all_Length.cm._means
scaled_control_all_Length.cm._test<-as.data.frame(contrast(scaled_control_all_Length.cm._means))
scaled_control_all_Length.cm._test <- scaled_control_all_Length.cm._test[,c(8,9,14)] %>% filter(scaled_control_all_Length.cm._test[,14] <= 0.05)
scaled_control_all_Length.cm._test