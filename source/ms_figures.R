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
library(ggpubr)

# Load data
df <- read.csv("./data/combined_data_clean.csv") 
df$mutant <- ifelse(is.na(df$PA.group), "N", "Y")

# Create new variables
df$LR.density <- ifelse(df$collection.time == "14d", (df$upper.LR.count + df$mid.LR.count + df$lower.LR.count) / df$primary.length, 
                        ifelse(df$collection.time =="21d", (df$mid.LR.count + df$lower.LR.count) / ((2/3)*df$primary.length), NA))
df$avg.fruit.length = (df$basal.fruit.length + df$mid.fruit.length + df$upper.fruit.length) / 3
df$fitness = df$avg.fruit.length * df$fruit.num

# Change collection time names
df$plant.age <- ifelse(df$collection.time == '14d', 14, ifelse(
  df$collection.time == '21d', 21, 44))

# Order PA.groups
df$PA.group <- as.factor(df$PA.group)
df$PA.group <- ordered(df$PA.group, levels = c('0', '1+2', '>2'))

# Make gene names lowercase
df$gene <- str_to_lower(df$gene)

# Define color scale for plots with just mutants
scale.mutants <- c('#CC79A7','#D55E00','#0072B2','#0072B2','#0072B2','#0072B2',
                   '#0072B2','#F0E442','#F0E442','#F0E442','#F0E442','#009E73',
                   '#56B4E9','#56B4E9','#E69F00')

# Define color scale for plots with mutants and WT
#scale.all <- c('#E53935','#9C27B0','#0D47A1','#1565C0','#1E88E5','#42A5F5',
 #              '#90CAF9','#2E7D32','#43A047','#66BB6A','#A5D6A7','#00796B',
  #             '#FBC02D','#FFEB3B','#FF6F00',"#78909C","#546E7A",'#455A64',
   #            '#37474F','#263238','black')

#From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

scale.all <- c('#CC79A7','#D55E00','#0072B2','#0072B2','#0072B2','#0072B2',
               '#0072B2','#F0E442','#F0E442','#F0E442','#F0E442','#009E73',
               '#56B4E9','#56B4E9','#E69F00',"#78909C","#546E7A",'#455A64',
               '#37474F','#263238','black')

# Define shape scale for plots with mutants
scale.shape.mutants <- c(19,19,19,15,17,18,8,19,15,17,18,19,19,15,19)

# Define shape scale for plots with mutants and WT
scale.shape.all <- c(19,19,19,15,17,18,8,19,15,17,18,19,19,15,19,19,19,19,
                     19,19,19)

# Figure 1 Mutant Phenotypes by Treatment and Developmental Stage -----------------

Fig1A_sum <- df %>%
  filter(treatment != "promix") %>% filter(collection.time == "14d") %>%
  filter(gene != "phytometer") %>%
  group_by(treatment, gene) %>%
  summarise(mean_Length.cm. = mean(Length.cm., na.rm=T),
            Length.cm._sd = sd(Length.cm.,na.rm=T), Length.cm._samp = n()) %>%
  mutate(Length.cm..stderr = Length.cm._sd/sqrt(Length.cm._samp))

# Order gene names
Fig1A_sum$gene <- as.factor(Fig1A_sum$gene)
Fig1A_sum$gene <- ordered(Fig1A_sum$gene, levels = c("abp1", "afb3", "arf1", "arf4", "arf11",
                                                     "arf14", "arf21", "iaa5", "iaa8", "iaa15",
                                                     "iaa33", "lax2", "pin7", "pin8", "sur1"))

Fig1A <- Fig1A_sum %>% 
  ggplot(aes(x = treatment, y = mean_Length.cm., color = gene, shape = gene)) + geom_point(size = 4) +
  #geom_errorbar((aes(ymin = mean_Length.cm. - Length.cm..stderr,
  #                   ymax = mean_Length.cm. + Length.cm..stderr)), width = .2) +
  stat_summary(aes(group = gene), geom = "line", fun.y = mean, size = 1.2) +
  scale_color_manual(values = scale.mutants) +
  scale_shape_manual(values = scale.shape.mutants) +
  ylab("Mean Root Length (cm)") +
  scale_x_discrete(labels = c('Control', "Auxin")) +
  theme_classic() +
  ggtitle("A)        Early Seedling") +
  theme(legend.position = "none",
        text=element_text(size=17, face='bold'),
        legend.text=element_text(face="italic"),
        legend.title=element_blank(),
        axis.line = element_line(size=2)) + xlab("Treatment")

Fig1B_sum <- df %>%
  filter(treatment != "promix") %>% filter(collection.time == "21d") %>%
  filter(gene != "phytometer") %>%
  group_by(treatment, gene) %>%
  summarise(mean_Length.cm. = mean(Length.cm., na.rm=T),
            Length.cm._sd = sd(Length.cm.,na.rm=T), Length.cm._samp = n()) %>%
  mutate(Length.cm..stderr = Length.cm._sd/sqrt(Length.cm._samp))

# Order gene names
Fig1B_sum$gene <- as.factor(Fig1B_sum$gene)
Fig1B_sum$gene <- ordered(Fig1B_sum$gene, levels = c("abp1", "afb3", "arf1", "arf4", "arf11",
                                                     "arf14", "arf21", "iaa5", "iaa8", "iaa15",
                                                     "iaa33", "lax2", "pin7", "pin8", "sur1"))

Fig1B <- Fig1B_sum %>% 
  ggplot(aes(x = treatment, y = mean_Length.cm., color = gene, shape = gene)) + geom_point(size = 4) +
  #geom_errorbar((aes(ymin = mean_Length.cm. - Length.cm..stderr,
  #                   ymax = mean_Length.cm. + Length.cm..stderr)), width = .2) +
  stat_summary(aes(group = gene), geom = "line", fun.y = mean, size = 1.2) +
  scale_color_manual(values = scale.mutants) +
  scale_shape_manual(values = scale.shape.mutants) +
  ylab("Mean Root Length (cm)") +
  scale_x_discrete(labels = c('Control', "Auxin")) +
  theme_classic() +
  ggtitle("B)        Late Seedling") +
  theme(legend.position = "none",
        text=element_text(size=17, face='bold'),
        legend.text=element_text(face="italic"),
        legend.title=element_blank(),
        axis.line = element_line(size=2)) + xlab("Treatment")


Fig1D_sum <- df %>%
  filter(treatment != "promix") %>% filter(collection.time == "14d") %>%
  filter(gene != "phytometer") %>%
  group_by(treatment, gene) %>%
  summarise(mean_LR.density = mean(LR.density, na.rm=T),
            LR.density_sd = sd(LR.density,na.rm=T), LR.density_samp = n()) %>%
  mutate(LR.density.stderr = LR.density_sd/sqrt(LR.density_samp))

# Order gene names
Fig1D_sum$gene <- as.factor(Fig1D_sum$gene)
Fig1D_sum$gene <- ordered(Fig1D_sum$gene, levels = c("abp1", "afb3", "arf1", "arf4", "arf11",
                                                     "arf14", "arf21", "iaa5", "iaa8", "iaa15",
                                                     "iaa33", "lax2", "pin7", "pin8", "sur1"))

Fig1D <- Fig1D_sum %>% 
  ggplot(aes(x = treatment, y = mean_LR.density, color = gene, shape = gene)) + geom_point(size = 4) +
  #geom_errorbar((aes(ymin = mean_LR.density - LR.density.stderr,
  #                   ymax = mean_LR.density + LR.density.stderr)), width = .2) +
  stat_summary(aes(group = gene), geom = "line", fun.y = mean, size = 1.2) +
  scale_color_manual(values = scale.mutants) +
  scale_shape_manual(values = scale.shape.mutants) +
  ylab("Mean LR Density") +
  scale_x_discrete(labels = c('Control', "Auxin")) +
  theme_classic() +
  ggtitle("D)        Early Seedling") +
  theme(legend.position = "right",
        text=element_text(size=17, face='bold'),
        legend.text=element_text(face="italic"),
        legend.title=element_blank(),
        axis.line = element_line(size=2)) + xlab("Treatment")

Fig1E_sum <- df %>%
  filter(treatment != "promix") %>% filter(collection.time == "21d") %>%
  filter(gene != "phytometer") %>%
  group_by(treatment, gene) %>%
  summarise(mean_above.ground.mass.mg = mean(above.ground.mass.mg, na.rm=T),
            above.ground.mass.mg_sd = sd(above.ground.mass.mg,na.rm=T), above.ground.mass.mg_samp = n()) %>%
  mutate(above.ground.mass.mg.stderr = above.ground.mass.mg_sd/sqrt(above.ground.mass.mg_samp))

# Order gene names
Fig1E_sum$gene <- as.factor(Fig1E_sum$gene)
Fig1E_sum$gene <- ordered(Fig1E_sum$gene, levels = c("abp1", "afb3", "arf1", "arf4", "arf11",
                                                     "arf14", "arf21", "iaa5", "iaa8", "iaa15",
                                                     "iaa33", "lax2", "pin7", "pin8", "sur1"))

Fig1E <- Fig1E_sum %>% 
  ggplot(aes(x = treatment, y = mean_above.ground.mass.mg, color = gene, shape = gene)) + geom_point(size = 4) +
  #geom_errorbar((aes(ymin = mean_above.ground.mass.mg - above.ground.mass.mg.stderr,
  #                   ymax = mean_above.ground.mass.mg + above.ground.mass.mg.stderr)), width = .2) +
  stat_summary(aes(group = gene), geom = "line", fun.y = mean, size = 1.2) +
  scale_color_manual(values = scale.mutants) +
  scale_shape_manual(values = scale.shape.mutants) +
  ylab("Mean Abovegroundmass (mg)") +
  scale_x_discrete(labels = c('Control', "Auxin")) +
  theme_classic() +
  ggtitle("E)       Late Seedling") +
  theme(legend.position = "none",
        text=element_text(size=17, face='bold'),
        legend.text=element_text(face="italic"),
        legend.title=element_blank(),
        axis.line = element_line(size=2)) + xlab("Treatment")

# Get gene names for mature plants
temp <- df %>% filter(collection.time == "14d" & gene != "phytometer") %>%
  count(gene, locus_capt)
temp$mutant_gene <- temp$gene
df <- temp %>% select(mutant_gene, locus_capt) %>% right_join(df)
df$gene <- df$mutant_gene

Fig1F_sum <- df %>%
  filter(treatment != "promix") %>% filter(collection.time == "fruit") %>%
  filter(gene != "phytometer") %>%
  group_by(treatment, gene) %>%
  summarise(mean_day.to.flower = mean(day.to.flower, na.rm=T),
            day.to.flower_sd = sd(day.to.flower,na.rm=T), day.to.flower_samp = n()) %>%
  mutate(day.to.flower.stderr = day.to.flower_sd/sqrt(day.to.flower_samp))

# Order gene names
Fig1F_sum$gene <- as.factor(Fig1F_sum$gene)
Fig1F_sum$gene <- ordered(Fig1F_sum$gene, levels = c("abp1", "afb3", "arf1", "arf4", "arf11",
                                                     "arf14", "arf21", "iaa5", "iaa8", "iaa15",
                                                     "iaa33", "lax2", "pin7", "pin8", "sur1"))

Fig1F <- Fig1F_sum %>% 
  ggplot(aes(x = treatment, y = mean_day.to.flower, color = gene, shape = gene)) + geom_point(size = 4) +
  #geom_errorbar((aes(ymin = mean_day.to.flower - day.to.flower.stderr,
  #                   ymax = mean_day.to.flower + day.to.flower.stderr)), width = .2) +
  stat_summary(aes(group = gene), geom = "line", fun.y = mean, size = 1.2) +
  scale_color_manual(values = scale.mutants) +
  scale_shape_manual(values = scale.shape.mutants) +
  ylab("Mean Days to Flower") +
  scale_x_discrete(labels = c('Control', "Auxin")) +
  theme_classic() +
  ggtitle("F)             Mature") +
  theme(legend.position = "none",
        text=element_text(size=17, face='bold'),
        legend.text=element_text(face="italic"),
        legend.title=element_blank(),
        axis.line = element_line(size=2)) + xlab("Treatment")

Fig1C_sum <- df %>%
  filter(treatment != "promix") %>% filter(collection.time == "fruit") %>%
  filter(gene != "phytometer") %>%
  group_by(treatment, gene) %>%
  summarise(mean_Length.cm. = mean(Length.cm., na.rm=T),
            Length.cm._sd = sd(Length.cm.,na.rm=T), Length.cm._samp = n()) %>%
  mutate(Length.cm..stderr = Length.cm._sd/sqrt(Length.cm._samp))

# Order gene names
Fig1C_sum$gene <- as.factor(Fig1C_sum$gene)
Fig1C_sum$gene <- ordered(Fig1C_sum$gene, levels = c("abp1", "afb3", "arf1", "arf4", "arf11",
                                                     "arf14", "arf21", "iaa5", "iaa8", "iaa15",
                                                     "iaa33", "lax2", "pin7", "pin8", "sur1"))

Fig1C <- Fig1C_sum %>% 
  ggplot(aes(x = treatment, y = mean_Length.cm., color = gene, shape = gene)) + geom_point(size = 4) +
  #geom_errorbar((aes(ymin = mean_Length.cm. - Length.cm..stderr,
  #                   ymax = mean_Length.cm. + Length.cm..stderr)), width = .2) +
  stat_summary(aes(group = gene), geom = "line", fun.y = mean, size = 1.2) +
  scale_color_manual(values = scale.mutants) +
  scale_shape_manual(values = scale.shape.mutants) +
  ylab("Mean Root Length (cm)") +
  scale_x_discrete(labels = c('Control', "Auxin")) +
  theme_classic() +
  ggtitle("C)             Mature") +
  theme(legend.position = "none",
        text=element_text(size=17, face='bold'),
        legend.text=element_text(face="italic"),
        legend.title=element_blank(),
        axis.line = element_line(size=2)) + xlab("Treatment")

pdf(file = './figs/Figure1.pdf', height = 10, width = 14)
ggarrange(Fig1A, Fig1B, Fig1C, Fig1D, Fig1E, Fig1F, ncol = 3, nrow = 2, common.legend = TRUE, legend = "right")
dev.off()

# Figure 2 Plant Phenology by Treatment ------------------------------------------

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
        text=element_text(size=18, face='bold'),
        axis.line = element_line(size=2))+
  xlab('Day') + ylab('Treatment')

# plot fruit.num boxplot

pheno_2 <- df %>% filter(treatment != 'promix') %>% filter(plant.age == 44) %>%
  ggplot(aes(x = treatment, y = fitness)) + geom_boxplot(lwd = 1.25) +
  scale_x_discrete(labels = c('Control', "Auxin")) +
  scale_y_continuous(breaks = c(0,100,300,500))+
  theme_classic() +
  ggtitle("") +
  xlab('Treatment') + ylab('Fitness') +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        text=element_text(size=18, face='bold'),
        axis.line = element_line(size=2)) +
  coord_flip()

# add together

pdf(file = './figs/Figure2.pdf', height = 7, width = 10)
grid.arrange(pheno_1, pheno_2, nrow = 1, widths=c(2.75,1))
dev.off()

# Figure 3 Flowering time of ARF and IAA gene families --------------------

Fig3_iaa <- df %>% filter(treatment != "promix") %>% filter(collection.time == "fruit") %>%
  filter(grepl('iaa', gene)) %>% mutate(family="IAA")
Fig3_arf <- df %>% filter(treatment != "promix") %>% filter(collection.time == "fruit") %>%
  filter(grepl('arf', gene)) %>% mutate(family="ARF")
Fig3_df <- bind_rows(Fig3_arf, Fig3_iaa)

Fig3_sum <- Fig3_df %>%
  group_by(treatment, family) %>%
  summarise(mean_day.to.flower = mean(day.to.flower, na.rm=T),
            day.to.flower_sd = sd(day.to.flower,na.rm=T), day.to.flower_samp = n()) %>%
  mutate(day.to.flower.stderr = day.to.flower_sd/sqrt(day.to.flower_samp))

Fig3 <- Fig3_sum %>% 
  ggplot(aes(x = treatment, y = mean_day.to.flower, color = family)) +
  geom_errorbar((aes(ymin = mean_day.to.flower - day.to.flower.stderr,
                     ymax = mean_day.to.flower + day.to.flower.stderr)), width = .2) +
  stat_summary(aes(group = family), geom = "line", fun.y = mean, size = 1.2) +
  geom_point(size = 4) +
  scale_color_discrete(name = "Genotype") +
  ylab("Mean Days to Flower") +
  scale_x_discrete(labels = c('Control', "Auxin")) +
  scale_color_manual(values = c("black", "grey"))+
  theme_classic() +
  theme(legend.position = "right",
        text=element_text(size=17, face='bold'),
        legend.text=element_text(face="italic"),
        legend.title=element_blank(),
        axis.line = element_line(size=2)) + xlab("Treatment")

pdf(file = './figs/Figure3.pdf', height = 5, width = 5)
Fig3
dev.off()

# Figure 4 LR Density Across Seedlings by Treatment -----------------------

# Redefine LR.Density to compare across stages
df$LR.density <- ifelse(df$collection.time == "14d", (df$mid.LR.count + df$lower.LR.count) / ((2/3)*df$primary.length),
                        ifelse(df$collection.time =="21d", (df$mid.LR.count + df$lower.LR.count) / ((2/3)*df$primary.length), NA))

df$mutant <- ifelse(df$insert.location == "phytometer", "N", "Y")
df$gene <- str_to_lower(df$gene)

Fig4_scaled <- df %>% filter(collection.time != "fruit") %>%
  group_by(collection.time, treatment)%>%
  mutate(scaled_LR.density=scale(LR.density))

Fig4 <- Fig4_scaled %>%
  filter(treatment != "promix") %>%
  group_by(mutant,plant.age, treatment, Stock_number,gene) %>%
  summarise(mean_scaled_LR.density = mean(scaled_LR.density, na.rm=T),
            scaled_LR.density_sd = sd(scaled_LR.density,na.rm=T), scaled_LR.density_samp = n()) %>%
  mutate(scaled_LR.density.stderr = scaled_LR.density_sd/sqrt(scaled_LR.density_samp))

Fig4$gene <- str_to_lower(Fig4$gene)
Fig4$gene <- ifelse(is.na(Fig4$gene), Fig4$Stock_number, Fig4$gene)

Fig4$treatment <- factor(Fig4$treatment, levels = c('control', 'IAA'),
                         labels = c("Control", "Auxin"))

# Order gene / accession names
Fig4$gene <- as.factor(Fig4$gene)
Fig4$gene <- ordered(Fig4$gene, levels = c("abp1", "afb3", "arf1", "arf4", "arf11",
                                           "arf14", "arf21", "iaa5", "iaa8", "iaa15",
                                           "iaa33", "lax2", "pin7", "pin8", "sur1",
                                           "CS22596", "CS22617", "CS22636", "CS22647",
                                           "CS22658", "CS70000"))

# New facet label names for mutant variable
mutant.labs <- c("Mutant", "Natural Accession")
names(mutant.labs) <- c("Y", "N")

Fig4 <- Fig4 %>%
  ggplot(aes(x = plant.age, y = mean_scaled_LR.density, color = gene, shape = gene)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size = 0.8) +
  facet_grid(mutant ~ treatment,  labeller = labeller(mutant = mutant.labs)) +
  stat_summary(aes(group = Stock_number), geom = "line", fun = mean, size = 1.2) +
  geom_point(size = 4) +
  scale_x_continuous(breaks=c(14,21,44), label = c("14", "21", "Mature")) +
  scale_y_continuous(breaks=c(-1,0,1), limits = c(-1.67, 1.67)) +
  scale_color_manual(values = scale.all) +
  scale_shape_manual(values = scale.shape.all) +
  ylab("Mean Scaled LR Density") +
  theme_classic() +
  theme(legend.position = "bottom",
        text=element_text(size=17, face='bold'),
        legend.text=element_text(face="italic"),
        legend.title=element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_line(size=2)) + xlab("Plant Age (Days)") 

pdf(file = './figs/Figure4.pdf', height = 7, width = 6)
Fig4
dev.off()

# Figure 5 Belowground Mass Across Stages by Treatment ----------------------
Fig5A_scaled <- df %>%
  group_by(collection.time, treatment)%>%
  mutate(scaled_below.ground.mass.mg=scale(below.ground.mass.mg))

Fig5A <- Fig5A_scaled %>%
  filter(treatment != "promix") %>%
  group_by(mutant,plant.age, treatment, Stock_number, gene) %>%
  summarise(mean_scaled_below.ground.mass.mg = mean(scaled_below.ground.mass.mg, na.rm=T),
            scaled_below.ground.mass.mg_sd = sd(scaled_below.ground.mass.mg,na.rm=T), scaled_below.ground.mass.mg_samp = n()) %>%
  mutate(scaled_below.ground.mass.mg.stderr = scaled_below.ground.mass.mg_sd/sqrt(scaled_below.ground.mass.mg_samp))

# replace values in sumstats table

Fig5A$plant.age <- ifelse(Fig5A$plant.age == 44 &
                           Fig5A$treatment == "IAA",
                         47, ifelse(
                           Fig5A$plant.age == 44 &
                             Fig5A$treatment == "control",
                           44.5,Fig5A$plant.age  
                         )
)

Fig5A$treatment <- factor(Fig5A$treatment, levels = c('control', 'IAA'),
                         labels = c("Control", "Auxin"))

Fig5A$gene <- ifelse(is.na(Fig5A$gene), Fig5A$Stock_number, Fig5A$gene)
Fig5A$mutant <- ifelse(is.na(Fig5A$mutant), "N", Fig5A$mutant)

# Order gene / accession names
Fig5A$gene <- as.factor(Fig5A$gene)
Fig5A$gene <- ordered(Fig5A$gene, levels = c("abp1", "afb3", "arf1", "arf4", "arf11",
                                           "arf14", "arf21", "iaa5", "iaa8", "iaa15",
                                           "iaa33", "lax2", "pin7", "pin8", "sur1",
                                           "CS22596", "CS22617", "CS22636", "CS22647",
                                           "CS22658", "CS70000"))

# New facet label names for mutant variable
mutant.labs <- c("Mutant", "Natural Accession")
names(mutant.labs) <- c("Y", "N")

Fig5A <- Fig5A %>% select(-7,-8, -9) %>%
  ggplot(aes(x = plant.age, y = mean_scaled_below.ground.mass.mg, color = gene, shape = gene)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size = 0.8) +
  facet_grid(mutant ~ treatment,  labeller = labeller(mutant = mutant.labs)) +
  stat_summary(aes(group = Stock_number), geom = "line", fun = mean, size = 1.2) +
  geom_point(size = 4) +
  scale_x_continuous(breaks=c(14,21,44), label = c("14", "21", "Mature")) +
  scale_y_continuous(breaks=c(-1,0,1), limits = c(-1.6, 1.6)) +
  scale_color_manual(values = scale.all) +
  scale_shape_manual(values = scale.shape.all) +
  ylab("Mean Scaled Belowground Mass") +
  theme_classic() +
  ggtitle("A)") +
  theme(legend.position = "right",
        text=element_text(size=17, face='bold'),
        legend.text=element_text(face="italic"),
        legend.title=element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_line(size=2)) + xlab("Plant Age (Days)") 

Fig5B_scaled <- df %>%
  group_by(collection.time, treatment)%>%
  mutate(scaled_Length.cm.=scale(Length.cm.))

Fig5B <- Fig5B_scaled %>%
  filter(treatment != "promix") %>%
  group_by(mutant,plant.age, treatment, Stock_number, gene) %>%
  summarise(mean_scaled_Length.cm. = mean(scaled_Length.cm., na.rm=T),
            scaled_Length.cm._sd = sd(scaled_Length.cm.,na.rm=T), scaled_Length.cm._samp = n()) %>%
  mutate(scaled_Length.cm..stderr = scaled_Length.cm._sd/sqrt(scaled_Length.cm._samp))

# replace values in sumstats table

Fig5B$plant.age <- ifelse(Fig5B$plant.age == 44 &
                            Fig5B$treatment == "IAA",
                          47, ifelse(
                            Fig5B$plant.age == 44 &
                              Fig5B$treatment == "control",
                            44.5,Fig5B$plant.age  
                          )
)

Fig5B$treatment <- factor(Fig5B$treatment, levels = c('control', 'IAA'),
                          labels = c("Control", "Auxin"))

Fig5B$gene <- ifelse(is.na(Fig5B$gene), Fig5B$Stock_number, Fig5B$gene)
Fig5B$mutant <- ifelse(is.na(Fig5B$mutant), "N", Fig5B$mutant)

# Order gene / accession names
Fig5B$gene <- as.factor(Fig5B$gene)
Fig5B$gene <- ordered(Fig5B$gene, levels = c("abp1", "afb3", "arf1", "arf4", "arf11",
                                             "arf14", "arf21", "iaa5", "iaa8", "iaa15",
                                             "iaa33", "lax2", "pin7", "pin8", "sur1",
                                             "CS22596", "CS22617", "CS22636", "CS22647",
                                             "CS22658", "CS70000"))

# New facet label names for mutant variable
mutant.labs <- c("Mutant", "Natural Accession")
names(mutant.labs) <- c("Y", "N")

Fig5B <- Fig5B %>%
  ggplot(aes(x = plant.age, y = mean_scaled_Length.cm., color = gene, shape = gene)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size = 0.8) +
  facet_grid(mutant ~ treatment,  labeller = labeller(mutant = mutant.labs)) +
  stat_summary(aes(group = Stock_number), geom = "line", fun = mean, size = 1.2) +
  geom_point(size = 4) +
  scale_x_continuous(breaks=c(14,21,44), label = c("14", "21", "Mature")) +
  scale_y_continuous(breaks=c(-1,0,1), limits = c(-1.6, 1.6)) +
  scale_color_manual(values = scale.all) +
  scale_shape_manual(values = scale.shape.all) +
  ylab("Mean Scaled Root Length") +
  theme_classic() +
  ggtitle("B)") +
  theme(legend.position = "right",
        text=element_text(size=17, face='bold'),
        legend.text=element_text(face="italic"),
        legend.title=element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_line(size=2)) + xlab("Plant Age (Days)") 

pdf(file = './figs/Figure5.pdf', height = 7, width = 14)
ggarrange(Fig5A, Fig5B, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

# Figure 6 Root Length by PA.group by Stage by Treatment -----------------------

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

# replace values in sumstats table

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

# Define color scale PA groups
scale.PA <- c('#fb6a4a','#de2d26','#a50f15')

p1 <- sumstats_dev_treat__group %>% 
  ggplot(aes(x = plant.age, y = mean_Length.cm., color = PA.group)) + geom_point(size = 4) +
  facet_wrap(~ treatment) +
  scale_color_manual(values = scale.PA) +
  geom_errorbar((aes(ymin = mean_Length.cm. - Length.cm..stderr,
                     ymax = mean_Length.cm. + Length.cm..stderr)), width = 1) +
  stat_summary(aes(group = PA.group), geom = "line", fun.y = mean, size = 2) +
  ylab("Mean Root Length (cm)") +
  theme_classic() +
  theme(legend.position = "top",
        text=element_text(size=20, face='bold'),
        axis.line = element_line(size=2)) + xlab("Plant Age (Days)") +
  guides(color=guide_legend(title="APA Group:")) 

# plot boxplot to be flipped and pasted on the side

p2 <- df %>% filter(treatment != 'promix') %>% filter(PA.group == '0' | PA.group == '1+2' | PA.group == '>2') %>%
  filter(plant.age == 44) %>%
  ggplot(aes(y = Length.cm., x = treatment)) + geom_boxplot(lwd=1.25) + 
  coord_cartesian(ylim = c(30,450)) +
  stat_summary(fun.y=mean, geom='point', shape=19, size=5, color="black", fill='black')+
  scale_x_discrete(labels = c('Control', "Auxin")) +
  theme_classic() +
  xlab('Treatment') + ylab('') +
  ggtitle("") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        text=element_text(size=20, face='bold'),
        axis.line = element_line(size=2))

pdf(file = './figs/Figure6.pdf', height = 8, width = 12)
grid.arrange(p1, p2, nrow = 1, widths=c(5,1))
dev.off()

# Figure S5 Histograms to show NAT vs SALK mutant variation -------------------------
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

Fig7a <- grid.arrange(R1, R2, R3, B1, B2, B3, A1, A2, A3, ncol = 3, widths=c(1.08,1,1))

# Control treatment

Fig7b <- grid.arrange(CR1, CR2, CR3, CB1, CB2, CB3, CA1, CA2, CA3, ncol = 3, widths=c(1.08,1,1))

pdf(file = './figs/FigureS5.pdf', height = 6.5, width = 10)
grid.arrange(Fig7b, Fig7a, ncol = 2)
dev.off()

# Levene's testing --------------------------------------------------------

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