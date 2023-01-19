rm(list = ls(all = TRUE)) 

library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(car)
library(lme4)
library(lmerTest)
library(gridExtra)
library(vistime)

# Load data
df <- read.csv("./data/combined_data_clean.csv") 
df$mutant <- ifelse(df$insert.location == "phytometer", "N", "Y")

df$avg.fruit.length = (df$basal.fruit.length + df$mid.fruit.length + df$upper.fruit.length) / 3
df$fitness = df$avg.fruit.length * df$fruit.num

# Change collection time names
df$plant.age <- ifelse(df$collection.time == '14d', 14, ifelse(
  df$collection.time == '21d', 21, 44))

# Order PA.groups

df$PA.group <- as.factor(df$PA.group)
df$PA.group <- ordered(df$PA.group, levels = c('0', '1+2', '>2'))


# PA.group by dev by treatment plot ----------------------------------------------------

sumstats_dev_treat__group <- df %>% filter(PA.group == '0' | PA.group == '1+2' | PA.group == '>2') %>% 
  filter(treatment != "promix") %>%
  group_by(plant.age, treatment, PA.group) %>%
  summarise(mean_Length.cm. = mean(Length.cm., na.rm=T),
            Length.cm._sd = sd(Length.cm.,na.rm=T), Length.cm._samp = n(),
            mean_above.ground.mass.mg = mean(above.ground.mass.mg, na.rm=T),
            above.ground.mass.mg_sd = sd(above.ground.mass.mg,na.rm=T), above.ground.mass.mg_samp = n(),
            mean_below.ground.mass.mg = mean(below.ground.mass.mg, na.rm=T),
            below.ground.mass.mg_sd = sd(below.ground.mass.mg,na.rm=T), below.ground.mass.mg_samp = n(),
            mean_LR.density = mean(LR.density, na.rm=T),
            LR.density_sd = sd(LR.density,na.rm=T), LR.density_samp = n()) %>%
  mutate(Length.cm..stderr = Length.cm._sd/sqrt(Length.cm._samp)) %>%
  mutate(above.ground.mass.mg.stderr = above.ground.mass.mg_sd/sqrt(above.ground.mass.mg_samp)) %>%
  mutate(below.ground.mass.mg.stderr = below.ground.mass.mg_sd/sqrt(below.ground.mass.mg_samp)) %>%
  mutate(LR.density.stderr = LR.density_sd/sqrt(LR.density_samp)) 

# need to change plant.age means to match treatment (control and IAA differ)

# determine means

day.to.mature_means <- df %>% filter(mutant == 'Y') %>% filter(treatment != "promix") %>%
  group_by(treatment) %>%
  summarise(mean_day.to.mature = mean(day.to.mature, na.rm=T),
            day.to.mature_sd = sd(day.to.mature,na.rm=T), day.to.mature_samp = n()) %>%
  mutate(day.to.mature.stderr = day.to.mature_sd/day.to.mature_samp) 

# replace values in sumstats table\

sumstats_dev_treat__group$plant.age <- ifelse(sumstats_dev_treat__group$plant.age == 44 &
                                                sumstats_dev_treat__group$treatment == "IAA",
                                              47, ifelse(
                                                sumstats_dev_treat__group$plant.age == 44 &
                                                  sumstats_dev_treat__group$treatment == "control",
                                                44.5,sumstats_dev_treat__group$plant.age  
                                              )
)

# plot

# Create new treatment labels

sumstats_dev_treat__group$treatment <- factor(sumstats_dev_treat__group$treatment, levels = c('control', 'IAA'),
                                              labels = c("Control", "Auxin"))

p1 <- sumstats_dev_treat__group %>% 
  ggplot(aes(x = plant.age, y = mean_Length.cm., color = PA.group)) + geom_point(size = 4) +
  facet_wrap(~ treatment) +
  scale_color_discrete("APA Group:") +
  scale_color_brewer(palette = 'Dark2') +
  geom_errorbar((aes(ymin = mean_Length.cm. - Length.cm..stderr,
                     ymax = mean_Length.cm. + Length.cm..stderr)), width = 1) +
  stat_summary(aes(group = PA.group), geom = "line", fun.y = mean, size = 2) +
  ylab("Mean Root Length (cm)") +
  theme_classic() +
  theme(legend.position = "top",
        text=element_text(size=20, face='bold'),
        axis.line = element_line(size=2)) + xlab("Plant Age (Days)")


# plot boxplot to be flipped and pasted on the side

p2 <- df %>% filter(treatment != 'promix') %>% filter(PA.group == '0' | PA.group == '1+2' | PA.group == '>2') %>%
  filter(plant.age == 44) %>%
  ggplot(aes(y = Length.cm., x = treatment)) + geom_boxplot(lwd=1.25) + 
  coord_cartesian(ylim = c(30,450)) +
  stat_summary(fun.y=mean, geom='point', shape=19, size=5, color="black", fill='black')+
  scale_x_discrete(labels = c('Control', "Auxin")) +
  theme_classic() +
  xlab('Treatment') + ylab('') +
  ggtitle("Mature") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        text=element_text(size=20, face='bold'),
        axis.line = element_line(size=2))

# throw together

grid.arrange(p1, p2, nrow = 1, widths=c(5,1))



# phenology by treatment plot ----------------------------------------------------


pheno <- df %>% filter(treatment != 'promix') %>% filter(plant.age == 44) %>%
  select(treatment, germination.day, day.to.bolt, day.to.flower, day.to.mature)



pheno_long <- gather(pheno, phase, day, germination.day:day.to.mature, factor_key = TRUE)
pheno_long$phase <- ifelse(pheno_long$phase == "germination.day", 'germination',
                           ifelse(pheno_long$phase == "day.to.bolt", 'bolting',
                                  ifelse(pheno_long$phase == 'day.to.flower', 'flowering', 'fruiting')))

sumstats_pheno <- pheno_long  %>%
  group_by(treatment, phase) %>%
  summarise(mean_day = mean(day, na.rm=T),
            day_sd = sd(day,na.rm=T), day_samp = n()) %>%
  mutate(day.stderr = day_sd/sqrt(day_samp))

sumstats_pheno$phase <- factor(sumstats_pheno$phase, levels = c('germination', 'bolting', 'flowering', 'fruiting'))

# plot time line

pheno_1 <- sumstats_pheno %>% ggplot(aes(x = mean_day, y = treatment, color = phase)) + geom_point(size = 6) +
  scale_y_discrete(labels = c('Control', "Auxin")) +
  scale_color_discrete("Event:") +
  theme_classic() +
  theme(legend.position = "top",
        text=element_text(size=20, face='bold'),
        axis.line = element_line(size=2))+
  xlab('Day') + ylab('Treatment')

# plot fruit.num boxplot

pheno_2 <- df %>% filter(treatment != 'promix') %>% filter(plant.age == 44) %>%
  ggplot(aes(x = treatment, y = fitness)) + geom_boxplot(lwd = 1.25) +
  scale_x_discrete(labels = c('Control', "Auxin")) +
  theme_classic() +
  ggtitle("") +
  xlab('Treatment') + ylab('Fitness') +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        text=element_text(size=20, face='bold'),
        axis.line = element_line(size=2)) +
  coord_flip()

# add together

grid.arrange(pheno_1, pheno_2, nrow = 1, widths=c(2.75,1))

test <- df %>% filter(treatment != 'promix') %>% filter(plant.age == 44) %>% 
  filter(mutant == 'Y')

summary(test$Length.cm.)


# reaction norm for natural accessions and SALK lines ---------------------

react <- df %>% filter(plant.age == 21) %>%
  select(Stock_number, LR.density, treatment, mutant, plant.age)

react$line = ifelse(react$mutant == "N", "Natural Accession", "Insert Mutant")

sumstats_LR.density_react <- react %>% 
  filter(treatment != "promix") %>%
  group_by(treatment, line, Stock_number) %>%
  summarise(mean_LR.density = mean(LR.density, na.rm=T),
            LR.density_sd = sd(LR.density,na.rm=T), LR.density_samp = n()) %>%
  mutate(LR.density.stderr = LR.density_sd/sqrt(LR.density_samp))

sumstats_LR.density_react %>% 
  ggplot(aes(x = treatment, y = mean_LR.density, color = Stock_number)) + geom_point(size = 4) +
  facet_wrap(~ line) +
  scale_x_discrete(labels = c('Control', "Auxin")) +
  geom_errorbar((aes(ymin = mean_LR.density - LR.density.stderr,
                     ymax = mean_LR.density + LR.density.stderr)), width = .15) +
  stat_summary(aes(group = Stock_number), geom = "line", fun.y = mean, size = 2) +
  ylab("Mean Lateral Root density") +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size=20, face = "bold"),
        axis.line = element_line(size=2)) + xlab("Treatment")

# reaction norm for natural accessions and SALK lines GXEXDEV ---------------------

react <- df %>% filter(plant.age == 21 | plant.age == 14) %>%
  select(Stock_number, LR.density, treatment, mutant, plant.age)

react$line = ifelse(react$mutant == "N", "Natural Accession", "Insert Mutant")
react$plant.age = ifelse(react$plant.age == 14, "Early Seedling", "Late Seedling")
sumstats_LR.density_react <- react %>% 
  filter(treatment != "promix") %>%
  group_by(treatment, line, Stock_number, plant.age) %>%
  summarise(mean_LR.density = mean(LR.density, na.rm=T),
            LR.density_sd = sd(LR.density,na.rm=T), LR.density_samp = n()) %>%
  mutate(LR.density.stderr = LR.density_sd/sqrt(LR.density_samp))

sumstats_LR.density_react %>% 
  ggplot(aes(x = treatment, y = mean_LR.density, color = Stock_number)) + geom_point(size = 3) +
  facet_wrap(line ~ plant.age) +
  scale_x_discrete(labels = c('Control', "Auxin")) +
  geom_errorbar((aes(ymin = mean_LR.density - LR.density.stderr,
                     ymax = mean_LR.density + LR.density.stderr)), width = .15) +
  stat_summary(aes(group = Stock_number), geom = "line", fun.y = mean, size = 1) +
  ylab("Mean Lateral Root density") +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size=10, face = "bold"),
        axis.line = element_line(size=1)) + xlab("Treatment")


# Histograms to show NAT vs SALK mutant variation -------------------------
library(scales)

df$Line = ifelse(is.na(df$PA.group), "Natural Accession", "Insert Mutant")

R1 <- df %>% filter(treatment == 'IAA') %>% filter(collection.time == "14d") %>% 
  ggplot(aes(x=Length.cm., fill=Line, color=Line)) +
  geom_histogram(alpha=.6) +
  theme_classic() +
  ggtitle("") +
  xlab('') + ylab('Count') +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.position = 'none',
        text=element_text(size=10, face='bold'),
        axis.line = element_line(size=1)) 

R2 <- df %>% filter(treatment == 'IAA') %>% filter(collection.time == "21d") %>% 
  ggplot(aes(x=Length.cm., fill=Line, color=Line)) +
  geom_histogram(alpha=.6) +
  theme_classic() +
  ggtitle("") +
  xlab('Root Length (cm)') + ylab('Count') +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.position = 'none',
        text=element_text(size=10, face='bold'),
        axis.title.y = element_blank(),
        axis.line = element_line(size=1))

R3 <- df %>% filter(treatment == 'IAA') %>% filter(collection.time == "fruit") %>% 
  ggplot(aes(x=Length.cm., fill=Line, color=Line)) +
  geom_histogram(alpha=.6) +
  theme_classic() +
  ggtitle("") +
  xlab('') + ylab('Count') +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.position = 'none',
        text=element_text(size=10, face='bold'),
        axis.title.y = element_blank(),
        axis.line = element_line(size=1)) 

B1 <- df %>% filter(treatment == 'IAA') %>% filter(collection.time == "14d") %>% 
  ggplot(aes(x=below.ground.mass.mg, fill=Line, color=Line)) +
  geom_histogram(alpha=.6) +
  theme_classic() +
  ggtitle("") +
  xlab('') + ylab('Count') +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.position = 'none',
        text=element_text(size=10, face='bold'),
        axis.line = element_line(size=1)) 

B2 <- df %>% filter(treatment == 'IAA') %>% filter(collection.time == "21d") %>% 
  ggplot(aes(x=below.ground.mass.mg, fill=Line, color=Line)) +
  geom_histogram(alpha=.6) +
  theme_classic() +
  ggtitle("") +
  xlab('Below Mass (mg)') + ylab('Count') +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.position = 'none',
        text=element_text(size=10, face='bold'),
        axis.title.y = element_blank(),
        axis.line = element_line(size=1))


B3 <- df %>% filter(treatment == 'IAA') %>% filter(collection.time == "fruit") %>% 
  ggplot(aes(x=below.ground.mass.mg, fill=Line, color=Line)) +
  geom_histogram(alpha=.6) +
  theme_classic() +
  ggtitle("") +
  xlab('') + ylab('Count') +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.position = 'none',
        text=element_text(size=10, face='bold'),
        axis.title.y = element_blank(),
        axis.line = element_line(size=1)) 

A1 <- df %>% filter(treatment == 'IAA') %>% filter(collection.time == "14d") %>% 
  ggplot(aes(x=above.ground.mass.mg, fill=Line, color=Line)) +
  geom_histogram(alpha=.6) +
  theme_classic() +
  ggtitle("") +
  xlab('') + ylab('Count') +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.position = 'none',
        text=element_text(size=10, face='bold'),
        axis.line = element_line(size=1)) 

A2 <- df %>% filter(treatment == 'IAA') %>% filter(collection.time == "21d") %>% 
  ggplot(aes(x=above.ground.mass.mg, fill=Line, color=Line)) +
  geom_histogram(alpha=.6) +
  theme_classic() +
  ggtitle("") +
  xlab('Above Mass (mg)') + ylab('Count') +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.position = 'none',
        text=element_text(size=10, face='bold'),
        axis.title.y = element_blank(),
        axis.line = element_line(size=1))


A3 <- df %>% filter(treatment == 'IAA') %>% filter(collection.time == "fruit") %>% 
  ggplot(aes(x=above.ground.mass.mg, fill=Line, color=Line)) +
  geom_histogram(alpha=.6) +
  theme_classic() +
  ggtitle("") +
  xlab('') + ylab('Count') +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.position = 'none',
        text=element_text(size=10, face='bold'),
        axis.title.y = element_blank(),
        axis.line = element_line(size=1)) 

CR1 <- df %>% filter(treatment == 'control') %>% filter(collection.time == "14d") %>% 
  ggplot(aes(x=Length.cm., fill=Line, color=Line)) +
  geom_histogram(alpha=.6) +
  theme_classic() +
  ggtitle("") +
  xlab('') + ylab('Count') +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.position = 'none',
        text=element_text(size=10, face='bold'),
        axis.line = element_line(size=1)) 

CR2 <- df %>% filter(treatment == 'control') %>% filter(collection.time == "21d") %>% 
  ggplot(aes(x=Length.cm., fill=Line, color=Line)) +
  geom_histogram(alpha=.6) +
  theme_classic() +
  ggtitle("") +
  xlab('Root Length (cm)') + ylab('Count') +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.position = 'none',
        text=element_text(size=10, face='bold'),
        axis.title.y = element_blank(),
        axis.line = element_line(size=1))

CR3 <- df %>% filter(treatment == 'control') %>% filter(collection.time == "fruit") %>% 
  ggplot(aes(x=Length.cm., fill=Line, color=Line)) +
  geom_histogram(alpha=.6) +
  theme_classic() +
  ggtitle("") +
  xlab('') + ylab('Count') +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.position = 'none',
        text=element_text(size=10, face='bold'),
        axis.title.y = element_blank(),
        axis.line = element_line(size=1)) 

CB1 <- df %>% filter(treatment == 'control') %>% filter(collection.time == "14d") %>% 
  ggplot(aes(x=below.ground.mass.mg, fill=Line, color=Line)) +
  geom_histogram(alpha=.6) +
  theme_classic() +
  ggtitle("") +
  xlab('') + ylab('Count') +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.position = 'none',
        text=element_text(size=10, face='bold'),
        axis.line = element_line(size=1)) 

CB2 <- df %>% filter(treatment == 'control') %>% filter(collection.time == "21d") %>% 
  ggplot(aes(x=below.ground.mass.mg, fill=Line, color=Line)) +
  geom_histogram(alpha=.6) +
  theme_classic() +
  ggtitle("") +
  xlab('Below Mass (mg)') + ylab('Count') +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.position = 'none',
        text=element_text(size=10, face='bold'),
        axis.title.y = element_blank(),
        axis.line = element_line(size=1))


CB3 <- df %>% filter(treatment == 'control') %>% filter(collection.time == "fruit") %>% 
  ggplot(aes(x=below.ground.mass.mg, fill=Line, color=Line)) +
  geom_histogram(alpha=.6) +
  theme_classic() +
  ggtitle("") +
  xlab('') + ylab('Count') +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.position = 'none',
        text=element_text(size=10, face='bold'),
        axis.title.y = element_blank(),
        axis.line = element_line(size=1)) 

CA1 <- df %>% filter(treatment == 'control') %>% filter(collection.time == "14d") %>% 
  ggplot(aes(x=above.ground.mass.mg, fill=Line, color=Line)) +
  geom_histogram(alpha=.6) +
  theme_classic() +
  ggtitle("") +
  xlab('') + ylab('Count') +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.position = 'none',
        text=element_text(size=10, face='bold'),
        axis.line = element_line(size=1)) 

CA2 <- df %>% filter(treatment == 'control') %>% filter(collection.time == "21d") %>% 
  ggplot(aes(x=above.ground.mass.mg, fill=Line, color=Line)) +
  geom_histogram(alpha=.6) +
  theme_classic() +
  ggtitle("") +
  xlab('Above Mass (mg)') + ylab('Count') +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.position = 'none',
        text=element_text(size=10, face='bold'),
        axis.title.y = element_blank(),
        axis.line = element_line(size=1))


CA3 <- df %>% filter(treatment == 'control') %>% filter(collection.time == "fruit") %>% 
  ggplot(aes(x=above.ground.mass.mg, fill=Line, color=Line)) +
  geom_histogram(alpha=.6) +
  theme_classic() +
  ggtitle("") +
  xlab('') + ylab('Count') +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.position = 'none',
        text=element_text(size=10, face='bold'),
        axis.title.y = element_blank(),
        axis.line = element_line(size=1)) 

# IAA treatment

grid.arrange(R1, R2, R3, B1, B2, B3, A1, A2, A3, ncol = 3, widths=c(1.08,1,1))

# Control treatment

grid.arrange(CR1, CR2, CR3, CB1, CB2, CB3, CA1, CA2, CA3, ncol = 3, widths=c(1.08,1,1))



#SD values

sumstats_dev_line_group <- df %>%
  filter(treatment != "promix") %>%
  group_by(collection.time, treatment, Line) %>%
  summarise(mean_Length.cm. = mean(Length.cm., na.rm=T),
            Length.cm._sd = sd(Length.cm.,na.rm=T), Length.cm._samp = n(),
            mean_above.ground.mass.mg = mean(above.ground.mass.mg, na.rm=T),
            above.ground.mass.mg_sd = sd(above.ground.mass.mg,na.rm=T), above.ground.mass.mg_samp = n(),
            mean_below.ground.mass.mg = mean(below.ground.mass.mg, na.rm=T),
            below.ground.mass.mg_sd = sd(below.ground.mass.mg,na.rm=T), below.ground.mass.mg_samp = n()) %>%
  mutate(Length.cm..stderr = Length.cm._sd/sqrt(Length.cm._samp)) %>%
  mutate(above.ground.mass.mg.stderr = above.ground.mass.mg_sd/sqrt(above.ground.mass.mg_samp)) %>%
  mutate(below.ground.mass.mg.stderr = below.ground.mass.mg_sd/sqrt(below.ground.mass.mg_samp))



#levene's  testing

early.IAA <- df %>% filter(treatment == 'IAA' & plant.age == 14)
early.control <- df %>% filter(treatment == 'control' & plant.age == 14)

early.IAA.Length.cm. = leveneTest(Length.cm. ~ Line, early.IAA)
early.IAA.Length.cm.
early.IAA.below.ground.mass.mg = leveneTest(below.ground.mass.mg ~ Line, early.IAA)
early.IAA.below.ground.mass.mg
early.IAA.above.ground.mass.mg = leveneTest(above.ground.mass.mg ~ Line, early.IAA)
early.IAA.above.ground.mass.mg

early.control.Length.cm. = leveneTest(Length.cm. ~ Line, early.control)
early.control.Length.cm.
early.control.below.ground.mass.mg = leveneTest(below.ground.mass.mg ~ Line, early.control)
early.control.below.ground.mass.mg
early.control.above.ground.mass.mg = leveneTest(above.ground.mass.mg ~ Line, early.control)
early.control.above.ground.mass.mg

late.IAA <- df %>% filter(treatment == 'IAA' & plant.age == 21)
late.control <- df %>% filter(treatment == 'control' & plant.age == 21)

late.IAA.Length.cm. = leveneTest(Length.cm. ~ Line, late.IAA)
late.IAA.Length.cm.
late.IAA.below.ground.mass.mg = leveneTest(below.ground.mass.mg ~ Line, late.IAA)
late.IAA.below.ground.mass.mg
late.IAA.above.ground.mass.mg = leveneTest(above.ground.mass.mg ~ Line, late.IAA)
late.IAA.above.ground.mass.mg

late.control.Length.cm. = leveneTest(Length.cm. ~ Line, late.control)
late.control.Length.cm.
late.control.below.ground.mass.mg = leveneTest(below.ground.mass.mg ~ Line, late.control)
late.control.below.ground.mass.mg
late.control.above.ground.mass.mg = leveneTest(above.ground.mass.mg ~ Line, late.control)
late.control.above.ground.mass.mg

mature.IAA <- df %>% filter(treatment == 'IAA' & collection.time == 'fruit')
mature.control <- df %>% filter(treatment == 'control' & collection.time == 'fruit')

mature.IAA.Length.cm. = leveneTest(Length.cm. ~ Line, mature.IAA)
mature.IAA.Length.cm.
mature.IAA.below.ground.mass.mg = leveneTest(below.ground.mass.mg ~ Line, mature.IAA)
mature.IAA.below.ground.mass.mg
mature.IAA.above.ground.mass.mg = leveneTest(above.ground.mass.mg ~ Line, mature.IAA)
mature.IAA.above.ground.mass.mg

mature.control.Length.cm. = leveneTest(Length.cm. ~ Line, mature.control)
mature.control.Length.cm.
mature.control.below.ground.mass.mg = leveneTest(below.ground.mass.mg ~ Line, mature.control)
mature.control.below.ground.mass.mg
mature.control.above.ground.mass.mg = leveneTest(above.ground.mass.mg ~ Line, mature.control)
mature.control.above.ground.mass.mg


# "Dev figure" ------------------------------------------------------------

#get gene names on to mature plants

temp <- df %>% filter(collection.time == "14d" & gene != "phytometer") %>%
  count(gene, locus_capt)

temp$mutant_gene <- temp$gene

df <- temp %>% select(mutant_gene, locus_capt) %>% right_join(df)


df_scaled <- df %>% group_by(collection.time, treatment)%>%
  mutate(scaled_above.ground=scale(above.ground.mass.mg))

scaled_dev_plot <- df_scaled %>%
  filter(treatment != "promix") %>%
  group_by(plant.age, treatment, Line, mutant_gene) %>%
  summarise(mean_scaled_above.ground = mean(scaled_above.ground, na.rm=T),
            scaled_above.ground_sd = sd(scaled_above.ground,na.rm=T), scaled_above.ground_samp = n()) %>%
  mutate(scaled_above.ground.stderr = scaled_above.ground_sd/sqrt(scaled_above.ground_samp))

scaled_dev_plot %>% filter(Line == "Insert Mutant") %>%
ggplot(aes(x = plant.age, y = mean_scaled_above.ground, color = mutant_gene)) + geom_point(size = 4) +
  facet_wrap(~ treatment) +
  geom_errorbar((aes(ymin = mean_scaled_above.ground - scaled_above.ground.stderr,
                     ymax = mean_scaled_above.ground + scaled_above.ground.stderr)), width = 1) +
  stat_summary(aes(group = mutant_gene), geom = "line", fun.y = mean, size = 2) +
  scale_color_discrete(name = "Genotype") +
  ylab("Mean Scaled Aboveground Mass") +
  theme_classic() +
  theme(legend.position = "top",
        text=element_text(size=17, face='bold'),
        legend.text=element_text(face="italic"),
        legend.title=element_blank(),
        axis.line = element_line(size=2)) + xlab("Plant Age (Days)")

# Q1 14d ------------------------------------------------------------------

df$gene <- str_to_lower(df$gene)

Q1A_sum <- df %>%
  filter(treatment != "promix") %>% filter(collection.time == "14d") %>%
  filter(gene != "phytometer") %>%
  group_by(treatment, gene) %>%
  summarise(mean_Length.cm. = mean(Length.cm., na.rm=T),
            Length.cm._sd = sd(Length.cm.,na.rm=T), Length.cm._samp = n()) %>%
  mutate(Length.cm..stderr = Length.cm._sd/sqrt(Length.cm._samp))

Q1A <- Q1A_sum %>% 
  ggplot(aes(x = treatment, y = mean_Length.cm., color = gene)) + geom_point(size = 4) +
  geom_errorbar((aes(ymin = mean_Length.cm. - Length.cm..stderr,
                     ymax = mean_Length.cm. + Length.cm..stderr)), width = .2) +
  stat_summary(aes(group = gene), geom = "line", fun.y = mean, size = 1.2) +
  scale_color_discrete(name = "Genotype") +
  ylab("Mean Root Length (cm)") +
  scale_x_discrete(labels = c('Control', "Auxin")) +
  theme_classic() +
  ggtitle("A)") +
  theme(legend.position = "none",
        text=element_text(size=17, face='bold'),
        legend.text=element_text(face="italic"),
        legend.title=element_blank(),
        axis.line = element_line(size=2)) + xlab("Treatment")

Q1A

Q1B_sum <- df %>%
  filter(treatment != "promix") %>% filter(collection.time == "14d") %>%
  filter(gene != "phytometer") %>%
  group_by(treatment, gene) %>%
  summarise(mean_LR.density = mean(LR.density, na.rm=T),
            LR.density_sd = sd(LR.density,na.rm=T), LR.density_samp = n()) %>%
  mutate(LR.density.stderr = LR.density_sd/sqrt(LR.density_samp))

Q1B <- Q1B_sum %>% 
  ggplot(aes(x = treatment, y = mean_LR.density, color = gene)) + geom_point(size = 4) +
  geom_errorbar((aes(ymin = mean_LR.density - LR.density.stderr,
                     ymax = mean_LR.density + LR.density.stderr)), width = .2) +
  stat_summary(aes(group = gene), geom = "line", fun.y = mean, size = 1.2) +
  scale_color_discrete(name = "Genotype") +
  ylab("Mean LR Density") +
  scale_x_discrete(labels = c('Control', "Auxin")) +
  theme_classic() +
  ggtitle("B)") +
  theme(legend.position = "right",
        text=element_text(size=17, face='bold'),
        legend.text=element_text(face="italic"),
        legend.title=element_blank(),
        axis.line = element_line(size=2)) + xlab("Treatment")

Q1B

grid.arrange(Q1A, Q1B, nrow = 1, widths=c(1,1.25))

# Q1 21d ------------------------------------------------------------------

df$gene <- str_to_lower(df$gene)

Q1C_sum <- df %>%
  filter(treatment != "promix") %>% filter(collection.time == "21d") %>%
  filter(gene != "phytometer") %>%
  group_by(treatment, gene) %>%
  summarise(mean_above.ground.mass.mg = mean(above.ground.mass.mg, na.rm=T),
            above.ground.mass.mg_sd = sd(above.ground.mass.mg,na.rm=T), above.ground.mass.mg_samp = n()) %>%
  mutate(above.ground.mass.mg.stderr = above.ground.mass.mg_sd/sqrt(above.ground.mass.mg_samp))

Q1C <- Q1C_sum %>% 
  ggplot(aes(x = treatment, y = mean_above.ground.mass.mg, color = gene)) + geom_point(size = 4) +
  geom_errorbar((aes(ymin = mean_above.ground.mass.mg - above.ground.mass.mg.stderr,
                     ymax = mean_above.ground.mass.mg + above.ground.mass.mg.stderr)), width = .2) +
  stat_summary(aes(group = gene), geom = "line", fun.y = mean, size = 1.2) +
  scale_color_discrete(name = "Genotype") +
  ylab("Mean Aboveground Mass (mg)") +
  scale_x_discrete(labels = c('Control', "Auxin")) +
  theme_classic() +
  ggtitle("A)") +
  theme(legend.position = "none",
        text=element_text(size=17, face='bold'),
        legend.text=element_text(face="italic"),
        legend.title=element_blank(),
        axis.line = element_line(size=2)) + xlab("Treatment")

Q1D_sum <- df %>%
  filter(treatment != "promix") %>% filter(collection.time == "21d") %>%
  filter(gene != "phytometer") %>%
  group_by(treatment, gene) %>%
  summarise(mean_Length.cm. = mean(Length.cm., na.rm=T),
            Length.cm._sd = sd(Length.cm.,na.rm=T), Length.cm._samp = n()) %>%
  mutate(Length.cm..stderr = Length.cm._sd/sqrt(Length.cm._samp))

Q1D <- Q1D_sum %>% 
  ggplot(aes(x = treatment, y = mean_Length.cm., color = gene)) + geom_point(size = 4) +
  geom_errorbar((aes(ymin = mean_Length.cm. - Length.cm..stderr,
                     ymax = mean_Length.cm. + Length.cm..stderr)), width = .2) +
  stat_summary(aes(group = gene), geom = "line", fun.y = mean, size = 1.2) +
  scale_color_discrete(name = "Genotype") +
  ylab("Mean Root Length (cm)") +
  scale_x_discrete(labels = c('Control', "Auxin")) +
  theme_classic() +
  ggtitle("B)") +
  theme(legend.position = "right",
        text=element_text(size=17, face='bold'),
        legend.text=element_text(face="italic"),
        legend.title=element_blank(),
        axis.line = element_line(size=2)) + xlab("Treatment")

Q1D

grid.arrange(Q1C, Q1D, nrow = 1, widths=c(1,1.25))

# Q1 Mature ------------------------------------------------------------------

# Get gene names

temp <- df %>% filter(collection.time == "14d" & gene != "phytometer") %>%
  count(gene, locus_capt)
temp$mutant_gene <- temp$gene
df <- temp %>% select(mutant_gene, locus_capt) %>% right_join(df)
df$gene <- df$mutant_gene

Q1E_sum <- df %>%
  filter(treatment != "promix") %>% filter(collection.time == "fruit") %>%
  filter(gene != "phytometer") %>%
  group_by(treatment, gene) %>%
  summarise(mean_day.to.flower = mean(day.to.flower, na.rm=T),
            day.to.flower_sd = sd(day.to.flower,na.rm=T), day.to.flower_samp = n()) %>%
  mutate(day.to.flower.stderr = day.to.flower_sd/sqrt(day.to.flower_samp))
Q1E_sum$gene <- str_to_lower(Q1E_sum$gene)


Q1E <- Q1E_sum %>% 
  ggplot(aes(x = treatment, y = mean_day.to.flower, color = gene)) + geom_point(size = 4) +
  geom_errorbar((aes(ymin = mean_day.to.flower - day.to.flower.stderr,
                     ymax = mean_day.to.flower + day.to.flower.stderr)), width = .2) +
  stat_summary(aes(group = gene), geom = "line", fun.y = mean, size = 1.2) +
  scale_color_discrete(name = "Genotype") +
  ylab("Mean Days to Flower") +
  scale_x_discrete(labels = c('Control', "Auxin")) +
  theme_classic() +
  ggtitle("A)") +
  theme(legend.position = "none",
        text=element_text(size=17, face='bold'),
        legend.text=element_text(face="italic"),
        legend.title=element_blank(),
        axis.line = element_line(size=2)) + xlab("Treatment")

Q1F_sum <- df %>%
  filter(treatment != "promix") %>% filter(collection.time == "fruit") %>%
  filter(gene != "phytometer") %>%
  group_by(treatment, gene) %>%
  summarise(mean_Length.cm. = mean(Length.cm., na.rm=T),
            Length.cm._sd = sd(Length.cm.,na.rm=T), Length.cm._samp = n()) %>%
  mutate(Length.cm..stderr = Length.cm._sd/sqrt(Length.cm._samp))
Q1F_sum$gene <- str_to_lower(Q1F_sum$gene)

Q1F <- Q1F_sum %>% 
  ggplot(aes(x = treatment, y = mean_Length.cm., color = gene)) + geom_point(size = 4) +
  geom_errorbar((aes(ymin = mean_Length.cm. - Length.cm..stderr,
                     ymax = mean_Length.cm. + Length.cm..stderr)), width = .2) +
  stat_summary(aes(group = gene), geom = "line", fun.y = mean, size = 1.2) +
  scale_color_discrete(name = "Genotype") +
  ylab("Mean Root Length (cm)") +
  scale_x_discrete(labels = c('Control', "Auxin")) +
  theme_classic() +
  ggtitle("B)") +
  theme(legend.position = "right",
        text=element_text(size=17, face='bold'),
        legend.text=element_text(face="italic"),
        legend.title=element_blank(),
        axis.line = element_line(size=2)) + xlab("Treatment")

grid.arrange(Q1E, Q1F, nrow = 1, widths=c(1,1.25))


# Q2 Seedling ----------------------------------------------------------------------

df$mutant <- ifelse(df$insert.location == "phytometer", "N", "Y")

df$gene <- str_to_lower(df$gene)

Q2A_scaled <- df %>% filter(collection.time != "fruit") %>%
  group_by(collection.time, treatment)%>%
  mutate(scaled_LR.density=scale(LR.density))

Q2A <- Q2A_scaled %>%
  filter(treatment != "promix") %>%
  group_by(mutant,plant.age, treatment, Stock_number,gene) %>%
  summarise(mean_scaled_LR.density = mean(scaled_LR.density, na.rm=T),
            scaled_LR.density_sd = sd(scaled_LR.density,na.rm=T), scaled_LR.density_samp = n()) %>%
  mutate(scaled_LR.density.stderr = scaled_LR.density_sd/sqrt(scaled_LR.density_samp))

Q2A$gene <- str_to_lower(Q2A$gene)

Q2A$gene <- ifelse(Q2A$gene == "phytometer", Q2A$Stock_number, Q2A$gene)

Q2A$treatment <- factor(Q2A$treatment, levels = c('control', 'IAA'),
                                              labels = c("Control", "Auxin"))

Q2A %>%
  ggplot(aes(x = plant.age, y = mean_scaled_LR.density, color = gene)) + geom_point(size = 4) +
  facet_wrap(~ treatment) +
  geom_errorbar((aes(ymin = mean_scaled_LR.density - scaled_LR.density.stderr,
                     ymax = mean_scaled_LR.density + scaled_LR.density.stderr)), width = 1) +
  stat_summary(aes(group = Stock_number, linetype = mutant), geom = "line", fun.y = mean, size = 1.2) +
  scale_linetype_manual(values=c(1, 5)) +
  scale_x_continuous(breaks=c(14,21)) +
  ylab("Mean Scaled LR Density") +
  theme_classic() +
  theme(legend.position = "right",
        text=element_text(size=17, face='bold'),
        legend.text=element_text(face="italic"),
        legend.title=element_blank(),
        axis.line = element_line(size=2)) + xlab("Plant Age (Days)") +
        guides(linetype = FALSE)

# Q2 Adult ----------------------------------------------------------------------

# root length

df$mutant <- ifelse(is.na(df$PA.group), "N", "Y")

Q2B_scaled <- df %>%
  group_by(collection.time, treatment)%>%
  mutate(scaled_Length.cm.=scale(Length.cm.))

Q2B <- Q2B_scaled %>%
  filter(treatment != "promix") %>%
  group_by(mutant,plant.age, treatment, Stock_number, gene) %>%
  summarise(mean_scaled_Length.cm. = mean(scaled_Length.cm., na.rm=T),
            scaled_Length.cm._sd = sd(scaled_Length.cm.,na.rm=T), scaled_Length.cm._samp = n()) %>%
  mutate(scaled_Length.cm..stderr = scaled_Length.cm._sd/sqrt(scaled_Length.cm._samp))

Q2B$gene <- ifelse(is.na(Q2B$gene), Q2B$Stock_number, Q2B$gene)

# replace values in sumstats table

Q2B$plant.age <- ifelse(Q2B$plant.age == 44 &
                          Q2B$treatment == "IAA",
                                              47, ifelse(
                                                Q2B$plant.age == 44 &
                                                  Q2B$treatment == "control",
                                                44.5,Q2B$plant.age  
                                              )
)
                        
Q2B$treatment <- factor(Q2B$treatment, levels = c('control', 'IAA'),
                        labels = c("Control", "Auxin"))

Q2B %>%
  ggplot(aes(x = plant.age, y = mean_scaled_Length.cm., color = gene)) +
  geom_point(size = 4) +
  facet_wrap(~ treatment) +
  geom_errorbar((aes(ymin = mean_scaled_Length.cm. - scaled_Length.cm..stderr,
                     ymax = mean_scaled_Length.cm. + scaled_Length.cm..stderr)), width = 1) +
  stat_summary(aes(group = Stock_number, linetype = mutant), geom = "line", fun.y = mean, size = 1.2) +
  scale_linetype_manual(values=c(1,6)) +
  scale_x_continuous(breaks=c(14,21,44, 47)) +
  ylab("Mean Scaled Root Length") +
  theme_classic() +
  theme(legend.position = "right",
        text=element_text(size=17, face='bold'),
        legend.text=element_text(face="italic"),
        legend.title=element_blank(),
        axis.line = element_line(size=2)) + xlab("Plant Age (Days)") +
  guides(linetype = FALSE)

# Below mass

Q2C_scaled <- df %>%
  group_by(collection.time, treatment)%>%
  mutate(scaled_below.ground.mass.mg=scale(below.ground.mass.mg))

Q2C <- Q2C_scaled %>%
  filter(treatment != "promix") %>%
  group_by(mutant,plant.age, treatment, Stock_number,gene) %>%
  summarise(mean_scaled_below.ground.mass.mg = mean(scaled_below.ground.mass.mg, na.rm=T),
            scaled_below.ground.mass.mg_sd = sd(scaled_below.ground.mass.mg,na.rm=T), scaled_below.ground.mass.mg_samp = n()) %>%
  mutate(scaled_below.ground.mass.mg.stderr = scaled_below.ground.mass.mg_sd/sqrt(scaled_below.ground.mass.mg_samp))

Q2C$gene <- ifelse(is.na(Q2C$gene), Q2C$Stock_number, Q2C$gene)


# replace values in sumstats table for plant age to be accurate for both treat

Q2C$plant.age <- ifelse(Q2C$plant.age == 44 &
                          Q2C$treatment == "IAA",
                        47, ifelse(
                          Q2C$plant.age == 44 &
                            Q2C$treatment == "control",
                          44.5,Q2C$plant.age  
                        )
)

Q2C$treatment <- factor(Q2C$treatment, levels = c('control', 'IAA'),
                        labels = c("Control", "Auxin"))

Q2C %>%
  ggplot(aes(x = plant.age, y = mean_scaled_below.ground.mass.mg, color = gene)) +
  geom_point(size = 4) +
  facet_wrap(~ treatment) +
  geom_errorbar((aes(ymin = mean_scaled_below.ground.mass.mg - scaled_below.ground.mass.mg.stderr,
                     ymax = mean_scaled_below.ground.mass.mg + scaled_below.ground.mass.mg.stderr)), width = 1) +
  stat_summary(aes(group = Stock_number, linetype = mutant), geom = "line", fun.y = mean, size = 1.2) +
  scale_linetype_manual(values=c(1,6)) +
  scale_x_continuous(breaks=c(14,21,44, 47)) +
  ylab("Mean Scaled Belowground Mass") +
  theme_classic() +
  theme(legend.position = "right",
        text=element_text(size=17, face='bold'),
        legend.text=element_text(face="italic"),
        legend.title=element_blank(),
        axis.line = element_line(size=2)) + xlab("Plant Age (Days)")+
  guides(linetype = FALSE)