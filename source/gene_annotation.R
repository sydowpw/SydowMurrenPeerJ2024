rm(list = ls(all = TRUE)) 

library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

### load up datasets
exp_df <- read.csv("./data/combined_data_clean.csv")
exp_df_salk <- exp_df %>% filter(!is.na(PA.site.number))
salk_loci <- unique(exp_df_salk$locus_capt)
salk_lines <- unique(exp_df_salk$Stock_number)

salk_lines_edit <-substr(salk_lines,1,11)

### load up insert location data and gene annotation data

salk_insert <- read.csv("./data/APA_files/SALK_insert_loc.csv")

AT1G34410 <- read.csv("./data/APA_files/AT1G34410.csv")

AT1G34410$height <- ifelse(AT1G34410$feature == "exon", 1, ifelse
                           (AT1G34410$feature == "CDS", 1, 0.05))

### RAN 1/25/22 PAC_df <- read.csv("plantapadb_arabidopsis_PAC_downloaded_12_15_21.csv")
### RAN 1/25/22 filtered_PAC_df <- PAC_df %>% filter(gene_id == salk_loci)

### RAN 1/25/22 write.csv(filtered_PAC_df, "filtered_PAC_df.csv") 

PAC_df <- read.csv("./data/filtered_PAC_df.csv")
PAC_df$locus <- PAC_df$gene_id
PAC_df$PAC_begin <- PAC_df$start
PAC_df$PAC_end <- PAC_df$end
PAC_df$feature = "PAC"
PAC_df <- select(PAC_df, PAC_begin, PAC_end, locus, feature)

### all loci with insert mutants from the experiment have PAC data!
unique(PAC_df$gene_id)
salk_loci

### test some graphical stuffs

### set target for plotting order
target = c("mRNA", "CDS", "exon", "five_prime_UTR", "three_prime_UTR")

###AT1G34410

AT1G34410 <- read.csv("./data/APA_files/AT1G34410.csv")
AT1G34410$height <- ifelse(AT1G34410$feature == "exon", 0.75, ifelse
                           (AT1G34410$feature == "CDS", 0.77, 0.05))
AT1G34410$feature <- factor(AT1G34410$feature, levels = c("mRNA", "CDS", "exon", "five_prime_UTR", "three_prime_UTR"))

###order features
AT1G34410 <- AT1G34410 %>% arrange(factor(feature, levels = target))

AT1G34410_salk_insert <- salk_insert %>% filter(locus == "AT1G34410")
AT1G34410_PAC_df <- PAC_df %>% filter(locus == "AT1G34410")

ggplot()+ geom_rect(data = AT1G34410, aes(xmin = begin, xmax = end, ymax = height, ymin = -1 *(height),
                                         fill = feature)) + ggtitle("AT1G34410") + xlab("coordinate") + ylab("") +
  scale_fill_manual(breaks = c("CDS", "exon", "five_prime_UTR", "mRNA", "three_prime_UTR"),
                    values = c("#073763", "#6FA8DC", "#93C47D", "#434343", "#E06666")) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y =  element_blank(),
        axis.ticks.y = element_blank()) +
  xlim(min(AT1G34410$begin - ((AT1G34410$end - AT1G34410$begin)/4)), max(AT1G34410$end +((AT1G34410$end - AT1G34410$begin)/4))) +
  ylim(-1,1) +
  geom_point(data = AT1G34410_PAC_df, aes(x = PAC_begin, y = 0.7, size = 1.5, color = "")) +
  geom_segment(data = AT1G34410_salk_insert, aes(x = begin, xend = end, y = 0.7, yend = 0.7, size = 1, color = "")) +
  facet_wrap(~transcript)

###AT2G46530

AT2G46530 <- read.csv("./data/APA_files/AT2G46530.csv")
AT2G46530$height <- ifelse(AT2G46530$feature == "exon", 0.75, ifelse
                           (AT2G46530$feature == "CDS", 0.77, 0.05))
AT2G46530$feature <- factor(AT2G46530$feature, levels = c("mRNA", "CDS", "exon", "five_prime_UTR", "three_prime_UTR"))

###order features
AT2G46530 <- AT2G46530 %>% arrange(factor(feature, levels = target))

AT2G46530_salk_insert <- salk_insert %>% filter(locus == "AT2G46530")
AT2G46530_PAC_df <- PAC_df %>% filter(locus == "AT2G46530")

ggplot()+ geom_rect(data = AT2G46530, aes(xmin = begin, xmax = end, ymax = height, ymin = -1 *(height),
                                          fill = feature)) + ggtitle("AT2G46530") + xlab("coordinate") + ylab("") +
  scale_fill_manual(breaks = c("CDS", "exon", "five_prime_UTR", "mRNA", "three_prime_UTR"),
                    values = c("#073763", "#6FA8DC", "#93C47D", "#434343", "#E06666")) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y =  element_blank(),
        axis.ticks.y = element_blank()) +
  xlim(min(AT2G46530$begin - ((AT2G46530$end - AT2G46530$begin)/4)), max(AT2G46530$end +((AT2G46530$end - AT2G46530$begin)/4))) +
  ylim(-1,1) +
  geom_point(data = AT2G46530_PAC_df, aes(x = PAC_begin, y = 0.7, size = 1.5, color = "")) +
  geom_segment(data = AT2G46530_salk_insert, aes(x = begin, xend = end, y = 0.7, yend = 0.7, size = 1, color = "")) +
  facet_wrap(~transcript)

###AT1G12820

AT1G12820 <- read.csv("./data/APA_files/AT1G12820.csv")
AT1G12820$height <- ifelse(AT1G12820$feature == "exon", 0.75, ifelse
                           (AT1G12820$feature == "CDS", 0.77, 0.05))
AT1G12820$feature <- factor(AT1G12820$feature, levels = c("mRNA", "CDS", "exon", "five_prime_UTR", "three_prime_UTR"))

###order features
AT1G12820 <- AT1G12820 %>% arrange(factor(feature, levels = target))

AT1G12820_salk_insert <- salk_insert %>% filter(locus == "AT1G12820")
AT1G12820_PAC_df <- PAC_df %>% filter(locus == "AT1G12820")

ggplot()+ geom_rect(data = AT1G12820, aes(xmin = begin, xmax = end, ymax = height, ymin = -1 *(height),
                                          fill = feature)) + ggtitle("AT1G12820") + xlab("coordinate") + ylab("") +
  scale_fill_manual(breaks = c("CDS", "exon", "five_prime_UTR", "mRNA", "three_prime_UTR"),
                    values = c("#073763", "#6FA8DC", "#93C47D", "#434343", "#E06666")) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y =  element_blank(),
        axis.ticks.y = element_blank()) +
  xlim(min(AT1G12820$begin - ((AT1G12820$end - AT1G12820$begin)/4)), max(AT1G12820$end +((AT1G12820$end - AT1G12820$begin)/4))) +
  ylim(-1,1) +
  geom_point(data = AT1G12820_PAC_df, aes(x = PAC_begin, y = 0.7, size = 1.5, color = "")) +
  geom_segment(data = AT1G12820_salk_insert, aes(x = begin, xend = end, y = 0.7, yend = 0.7, size = 1, color = "")) +
  facet_wrap(~transcript)

###AT2G21050

AT2G21050 <- read.csv("./data/APA_files/AT2G21050.csv")
AT2G21050$height <- ifelse(AT2G21050$feature == "exon", 0.75, ifelse
                           (AT2G21050$feature == "CDS", 0.77, 0.05))
AT2G21050$feature <- factor(AT2G21050$feature, levels = c("mRNA", "CDS", "exon", "five_prime_UTR", "three_prime_UTR"))

###order features
AT2G21050 <- AT2G21050 %>% arrange(factor(feature, levels = target))

AT2G21050_salk_insert <- salk_insert %>% filter(locus == "AT2G21050")
AT2G21050_PAC_df <- PAC_df %>% filter(locus == "AT2G21050")

ggplot()+ geom_rect(data = AT2G21050, aes(xmin = begin, xmax = end, ymax = height, ymin = -1 *(height),
                                          fill = feature)) + ggtitle("AT2G21050") + xlab("coordinate") + ylab("") +
  scale_fill_manual(breaks = c("CDS", "exon", "five_prime_UTR", "mRNA", "three_prime_UTR"),
                    values = c("#073763", "#6FA8DC", "#93C47D", "#434343", "#E06666")) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y =  element_blank(),
        axis.ticks.y = element_blank()) +
  xlim(min(AT2G21050$begin - ((AT2G21050$end - AT2G21050$begin)/4)), max(AT2G21050$end +((AT2G21050$end - AT2G21050$begin)/4))) +
  ylim(-1,1) +
  geom_point(data = AT2G21050_PAC_df, aes(x = PAC_begin, y = 0.7, size = 1.5, color = "")) +
  geom_segment(data = AT2G21050_salk_insert, aes(x = begin, xend = end, y = 0.7, yend = 0.7, size = 1, color = "")) +
  facet_wrap(~transcript)

###AT4G02980

AT4G02980 <- read.csv("./data/APA_files/AT4G02980.csv")
AT4G02980$height <- ifelse(AT4G02980$feature == "exon", 0.75, ifelse
                           (AT4G02980$feature == "CDS", 0.77, 0.05))
AT4G02980$feature <- factor(AT4G02980$feature, levels = c("mRNA", "CDS", "exon", "five_prime_UTR", "three_prime_UTR"))

###order features
AT4G02980 <- AT4G02980 %>% arrange(factor(feature, levels = target))

AT4G02980_salk_insert <- salk_insert %>% filter(locus == "AT4G02980")
AT4G02980_PAC_df <- PAC_df %>% filter(locus == "AT4G02980")

ggplot()+ geom_rect(data = AT4G02980, aes(xmin = begin, xmax = end, ymax = height, ymin = -1 *(height),
                                          fill = feature)) + ggtitle("AT4G02980") + xlab("coordinate") + ylab("") +
  scale_fill_manual(breaks = c("CDS", "exon", "five_prime_UTR", "mRNA", "three_prime_UTR"),
                    values = c("#073763", "#6FA8DC", "#93C47D", "#434343", "#E06666")) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y =  element_blank(),
        axis.ticks.y = element_blank()) +
  xlim(min(AT4G02980$begin - ((AT4G02980$end - AT4G02980$begin)/4)), max(AT4G02980$end +((AT4G02980$end - AT4G02980$begin)/4))) +
  ylim(-1,1) +
  geom_point(data = AT4G02980_PAC_df, aes(x = PAC_begin, y = 0.7, size = 1.5, color = "")) +
  geom_segment(data = AT4G02980_salk_insert, aes(x = begin, xend = end, y = 0.7, yend = 0.7, size = 1, color = "")) +
  facet_wrap(~transcript)

###AT1G15580

AT1G15580 <- read.csv("./data/APA_files/AT1G15580.csv")
AT1G15580$height <- ifelse(AT1G15580$feature == "exon", 0.75, ifelse
                           (AT1G15580$feature == "CDS", 0.77, 0.05))
AT1G15580$feature <- factor(AT1G15580$feature, levels = c("mRNA", "CDS", "exon", "five_prime_UTR", "three_prime_UTR"))

###order features
AT1G15580 <- AT1G15580 %>% arrange(factor(feature, levels = target))

AT1G15580_salk_insert <- salk_insert %>% filter(locus == "AT1G15580")
AT1G15580_PAC_df <- PAC_df %>% filter(locus == "AT1G15580")

ggplot()+ geom_rect(data = AT1G15580, aes(xmin = begin, xmax = end, ymax = height, ymin = -1 *(height),
                                          fill = feature)) + ggtitle("AT1G15580") + xlab("coordinate") + ylab("") +
  scale_fill_manual(breaks = c("CDS", "exon", "five_prime_UTR", "mRNA", "three_prime_UTR"),
                    values = c("#073763", "#6FA8DC", "#93C47D", "#434343", "#E06666")) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y =  element_blank(),
        axis.ticks.y = element_blank()) +
  xlim(min(AT1G15580$begin - ((AT1G15580$end - AT1G15580$begin)/4)), max(AT1G15580$end +((AT1G15580$end - AT1G15580$begin)/4))) +
  ylim(-1,1) +
  geom_point(data = AT1G15580_PAC_df, aes(x = PAC_begin, y = 0.7, size = 1.5, color = "")) +
  geom_segment(data = AT1G15580_salk_insert, aes(x = begin, xend = end, y = 0.7, yend = 0.7, size = 1, color = "")) +
  facet_wrap(~transcript)

###AT5G15100

AT5G15100 <- read.csv("./data/APA_files/AT5G15100.csv")
AT5G15100$height <- ifelse(AT5G15100$feature == "exon", 0.75, ifelse
                           (AT5G15100$feature == "CDS", 0.77, 0.05))
AT5G15100$feature <- factor(AT5G15100$feature, levels = c("mRNA", "CDS", "exon", "five_prime_UTR", "three_prime_UTR"))

###order features
AT5G15100 <- AT5G15100 %>% arrange(factor(feature, levels = target))

AT5G15100_salk_insert <- salk_insert %>% filter(locus == "AT5G15100")
AT5G15100_PAC_df <- PAC_df %>% filter(locus == "AT5G15100")

ggplot()+ geom_rect(data = AT5G15100, aes(xmin = begin, xmax = end, ymax = height, ymin = -1 *(height),
                                          fill = feature)) + ggtitle("AT5G15100") + xlab("coordinate") + ylab("") +
  scale_fill_manual(breaks = c("CDS", "exon", "five_prime_UTR", "mRNA", "three_prime_UTR"),
                    values = c("#073763", "#6FA8DC", "#93C47D", "#434343", "#E06666")) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y =  element_blank(),
        axis.ticks.y = element_blank()) +
  xlim(min(AT5G15100$begin - ((AT5G15100$end - AT5G15100$begin)/4)), max(AT5G15100$end +((AT5G15100$end - AT5G15100$begin)/4))) +
  ylim(-1,1) +
  geom_point(data = AT5G15100_PAC_df, aes(x = PAC_begin, y = 0.7, size = 1.5, color = "")) +
  geom_segment(data = AT5G15100_salk_insert, aes(x = begin, xend = end, y = 0.7, yend = 0.7, size = 1, color = "")) +
  facet_wrap(~transcript)

###AT2G20610

AT2G20610 <- read.csv("./data/APA_files/AT2G20610.csv")
AT2G20610$height <- ifelse(AT2G20610$feature == "exon", 0.75, ifelse
                           (AT2G20610$feature == "CDS", 0.77, 0.05))
AT2G20610$feature <- factor(AT2G20610$feature, levels = c("mRNA", "CDS", "exon", "five_prime_UTR", "three_prime_UTR"))

###order features
AT2G20610 <- AT2G20610 %>% arrange(factor(feature, levels = target))

AT2G20610_salk_insert <- salk_insert %>% filter(locus == "AT2G20610")
AT2G20610_PAC_df <- PAC_df %>% filter(locus == "AT2G20610")

ggplot()+ geom_rect(data = AT2G20610, aes(xmin = begin, xmax = end, ymax = height, ymin = -1 *(height),
                                          fill = feature)) + ggtitle("AT2G20610") + xlab("coordinate") + ylab("") +
  scale_fill_manual(breaks = c("CDS", "exon", "five_prime_UTR", "mRNA", "three_prime_UTR"),
                    values = c("#073763", "#6FA8DC", "#93C47D", "#434343", "#E06666")) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y =  element_blank(),
        axis.ticks.y = element_blank()) +
  xlim(min(AT2G20610$begin - ((AT2G20610$end - AT2G20610$begin)/4)), max(AT2G20610$end +((AT2G20610$end - AT2G20610$begin)/4))) +
  ylim(-1,1) +
  geom_point(data = AT2G20610_PAC_df, aes(x = PAC_begin, y = 0.7, size = 1.5, color = "")) +
  geom_segment(data = AT2G20610_salk_insert, aes(x = begin, xend = end, y = 0.7, yend = 0.7, size = 1, color = "")) +
  facet_wrap(~transcript)

###AT1G80390

AT1G80390 <- read.csv("./data/APA_files/AT1G80390.csv")
AT1G80390$height <- ifelse(AT1G80390$feature == "exon", 0.75, ifelse
                           (AT1G80390$feature == "CDS", 0.77, 0.05))
AT1G80390$feature <- factor(AT1G80390$feature, levels = c("mRNA", "CDS", "exon", "five_prime_UTR", "three_prime_UTR"))

###order features
AT1G80390 <- AT1G80390 %>% arrange(factor(feature, levels = target))

AT1G80390_salk_insert <- salk_insert %>% filter(locus == "AT1G80390")
AT1G80390_PAC_df <- PAC_df %>% filter(locus == "AT1G80390")

ggplot()+ geom_rect(data = AT1G80390, aes(xmin = begin, xmax = end, ymax = height, ymin = -1 *(height),
                                          fill = feature)) + ggtitle("AT1G80390") + xlab("coordinate") + ylab("") +
  scale_fill_manual(breaks = c("CDS", "exon", "five_prime_UTR", "mRNA", "three_prime_UTR"),
                    values = c("#073763", "#6FA8DC", "#93C47D", "#434343", "#E06666")) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y =  element_blank(),
        axis.ticks.y = element_blank()) +
  xlim(min(AT1G80390$begin - ((AT1G80390$end - AT1G80390$begin)/4)), max(AT1G80390$end +((AT1G80390$end - AT1G80390$begin)/4))) +
  ylim(-1,1) +
  geom_point(data = AT1G80390_PAC_df, aes(x = PAC_begin, y = 0.7, size = 1.5, color = "")) +
  geom_segment(data = AT1G80390_salk_insert, aes(x = begin, xend = end, y = 0.7, yend = 0.7, size = 1, color = "")) +
  facet_wrap(~transcript)

###AT2G22670

AT2G22670 <- read.csv("./data/APA_files/AT2G22670.csv")
AT2G22670$height <- ifelse(AT2G22670$feature == "exon", 0.75, ifelse
                           (AT2G22670$feature == "CDS", 0.77, 0.05))
AT2G22670$feature <- factor(AT2G22670$feature, levels = c("mRNA", "CDS", "exon", "five_prime_UTR", "three_prime_UTR"))

###order features
AT2G22670 <- AT2G22670 %>% arrange(factor(feature, levels = target))

AT2G22670_salk_insert <- salk_insert %>% filter(locus == "AT2G22670")
AT2G22670_PAC_df <- PAC_df %>% filter(locus == "AT2G22670")

ggplot()+ geom_rect(data = AT2G22670, aes(xmin = begin, xmax = end, ymax = height, ymin = -1 *(height),
                                          fill = feature)) + ggtitle("AT2G22670") + xlab("coordinate") + ylab("") +
  scale_fill_manual(breaks = c("CDS", "exon", "five_prime_UTR", "mRNA", "three_prime_UTR"),
                    values = c("#073763", "#6FA8DC", "#93C47D", "#434343", "#E06666")) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y =  element_blank(),
        axis.ticks.y = element_blank()) +
  xlim(min(AT2G22670$begin - ((AT2G22670$end - AT2G22670$begin)/4)), max(AT2G22670$end +((AT2G22670$end - AT2G22670$begin)/4))) +
  ylim(-1,1) +
  geom_point(data = AT2G22670_PAC_df, aes(x = PAC_begin, y = 0.7, size = 1.5, color = "")) +
  geom_segment(data = AT2G22670_salk_insert, aes(x = begin, xend = end, y = 0.7, yend = 0.7, size = 1, color = "")) +
  facet_wrap(~transcript)

###AT1G35540

AT1G35540 <- read.csv("./data/APA_files/AT1G35540.csv")
AT1G35540$height <- ifelse(AT1G35540$feature == "exon", 0.75, ifelse
                           (AT1G35540$feature == "CDS", 0.77, 0.05))
AT1G35540$feature <- factor(AT1G35540$feature, levels = c("mRNA", "CDS", "exon", "five_prime_UTR", "three_prime_UTR"))

###order features
AT1G35540 <- AT1G35540 %>% arrange(factor(feature, levels = target))

AT1G35540_salk_insert <- salk_insert %>% filter(locus == "AT1G35540")
AT1G35540_PAC_df <- PAC_df %>% filter(locus == "AT1G35540")

ggplot()+ geom_rect(data = AT1G35540, aes(xmin = begin, xmax = end, ymax = height, ymin = -1 *(height),
                                          fill = feature)) + ggtitle("AT1G35540") + xlab("coordinate") + ylab("") +
  scale_fill_manual(breaks = c("CDS", "exon", "five_prime_UTR", "mRNA", "three_prime_UTR"),
                    values = c("#073763", "#6FA8DC", "#93C47D", "#434343", "#E06666")) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y =  element_blank(),
        axis.ticks.y = element_blank()) +
  xlim(min(AT1G35540$begin - ((AT1G35540$end - AT1G35540$begin)/4)), max(AT1G35540$end +((AT1G35540$end - AT1G35540$begin)/4))) +
  ylim(-1,1) +
  geom_point(data = AT1G35540_PAC_df, aes(x = PAC_begin, y = 0.7, size = 1.5, color = "")) +
  geom_segment(data = AT1G35540_salk_insert, aes(x = begin, xend = end, y = 0.7, yend = 0.7, size = 1, color = "")) +
  facet_wrap(~transcript)

###AT1G59750

AT1G59750 <- read.csv("./data/APA_files/AT1G59750.csv")
AT1G59750$height <- ifelse(AT1G59750$feature == "exon", 0.75, ifelse
                           (AT1G59750$feature == "CDS", 0.77, 0.05))
AT1G59750$feature <- factor(AT1G59750$feature, levels = c("mRNA", "CDS", "exon", "five_prime_UTR", "three_prime_UTR"))

###order features
AT1G59750 <- AT1G59750 %>% arrange(factor(feature, levels = target))

AT1G59750_salk_insert <- salk_insert %>% filter(locus == "AT1G59750")
AT1G59750_PAC_df <- PAC_df %>% filter(locus == "AT1G59750")

ggplot()+ geom_rect(data = AT1G59750, aes(xmin = begin, xmax = end, ymax = height, ymin = -1 *(height),
                                          fill = feature)) + ggtitle("AT1G59750") + xlab("coordinate") + ylab("") +
  scale_fill_manual(breaks = c("CDS", "exon", "five_prime_UTR", "mRNA", "three_prime_UTR"),
                    values = c("#073763", "#6FA8DC", "#93C47D", "#434343", "#E06666")) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y =  element_blank(),
        axis.ticks.y = element_blank()) +
  xlim(min(AT1G59750$begin - ((AT1G59750$end - AT1G59750$begin)/4)), max(AT1G59750$end +((AT1G59750$end - AT1G59750$begin)/4))) +
  ylim(-1,1) +
  geom_point(data = AT1G59750_PAC_df, aes(x = PAC_begin, y = 0.7, size = 1.5, color = "")) +
  geom_segment(data = AT1G59750_salk_insert, aes(x = begin, xend = end, y = 0.7, yend = 0.7, size = 1, color = "")) +
  facet_wrap(~transcript)

###AT1G23080

AT1G23080 <- read.csv("./data/APA_files/AT1G23080.csv")
AT1G23080$height <- ifelse(AT1G23080$feature == "exon", 0.75, ifelse
                           (AT1G23080$feature == "CDS", 0.77, 0.05))
AT1G23080$feature <- factor(AT1G23080$feature, levels = c("mRNA", "CDS", "exon", "five_prime_UTR", "three_prime_UTR"))

###order features
AT1G23080 <- AT1G23080 %>% arrange(factor(feature, levels = target))

AT1G23080_salk_insert <- salk_insert %>% filter(locus == "AT1G23080")
AT1G23080_PAC_df <- PAC_df %>% filter(locus == "AT1G23080")

ggplot()+ geom_rect(data = AT1G23080, aes(xmin = begin, xmax = end, ymax = height, ymin = -1 *(height),
                                          fill = feature)) + ggtitle("AT1G23080") + xlab("coordinate") + ylab("") +
  scale_fill_manual(breaks = c("CDS", "exon", "five_prime_UTR", "mRNA", "three_prime_UTR"),
                    values = c("#073763", "#6FA8DC", "#93C47D", "#434343", "#E06666")) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y =  element_blank(),
        axis.ticks.y = element_blank()) +
  xlim(min(AT1G23080$begin - ((AT1G23080$end - AT1G23080$begin)/4)), max(AT1G23080$end +((AT1G23080$end - AT1G23080$begin)/4))) +
  ylim(-1,1) +
  geom_point(data = AT1G23080_PAC_df, aes(x = PAC_begin, y = 0.7, size = 1.5, color = "")) +
  geom_segment(data = AT1G23080_salk_insert, aes(x = begin, xend = end, y = 0.7, yend = 0.7, size = 1, color = "")) +
  facet_wrap(~transcript)

###AT5G60450

AT5G60450 <- read.csv("./data/APA_files/AT5G60450.csv")
AT5G60450$height <- ifelse(AT5G60450$feature == "exon", 0.75, ifelse
                           (AT5G60450$feature == "CDS", 0.77, 0.05))
AT5G60450$feature <- factor(AT5G60450$feature, levels = c("mRNA", "CDS", "exon", "five_prime_UTR", "three_prime_UTR"))

###order features
AT5G60450 <- AT5G60450 %>% arrange(factor(feature, levels = target))

AT5G60450_salk_insert <- salk_insert %>% filter(locus == "AT5G60450")
AT5G60450_PAC_df <- PAC_df %>% filter(locus == "AT5G60450")

ggplot()+ geom_rect(data = AT5G60450, aes(xmin = begin, xmax = end, ymax = height, ymin = -1 *(height),
                                          fill = feature)) + ggtitle("AT5G60450") + xlab("coordinate") + ylab("") +
  scale_fill_manual(breaks = c("CDS", "exon", "five_prime_UTR", "mRNA", "three_prime_UTR"),
                    values = c("#073763", "#6FA8DC", "#93C47D", "#434343", "#E06666")) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y =  element_blank(),
        axis.ticks.y = element_blank()) +
  xlim(min(AT5G60450$begin - ((AT5G60450$end - AT5G60450$begin)/4)), max(AT5G60450$end +((AT5G60450$end - AT5G60450$begin)/4))) +
  ylim(-1,1) +
  geom_point(data = AT5G60450_PAC_df, aes(x = PAC_begin, y = 0.7, size = 1.5, color = "")) +
  geom_segment(data = AT5G60450_salk_insert, aes(x = begin, xend = end, y = 0.7, yend = 0.7, size = 1, color = "")) +
  facet_wrap(~transcript)

###AT5G57420

AT5G57420 <- read.csv("./data/APA_files/AT5G57420.csv")
AT5G57420$height <- ifelse(AT5G57420$feature == "exon", 0.75, ifelse
                           (AT5G57420$feature == "CDS", 0.77, 0.05))
AT5G57420$feature <- factor(AT5G57420$feature, levels = c("mRNA", "CDS", "exon", "five_prime_UTR", "three_prime_UTR"))

###order features
AT5G57420 <- AT5G57420 %>% arrange(factor(feature, levels = target))

AT5G57420_salk_insert <- salk_insert %>% filter(locus == "AT5G57420")
AT5G57420_PAC_df <- PAC_df %>% filter(locus == "AT5G57420")

ggplot()+ geom_rect(data = AT5G57420, aes(xmin = begin, xmax = end, ymax = height, ymin = -1 *(height),
                                          fill = feature)) + ggtitle("AT5G57420") + xlab("coordinate") + ylab("") +
  scale_fill_manual(breaks = c("CDS", "exon", "five_prime_UTR", "mRNA", "three_prime_UTR"),
                    values = c("#073763", "#6FA8DC", "#93C47D", "#434343", "#E06666")) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y =  element_blank(),
        axis.ticks.y = element_blank()) +
  xlim(min(23269349), max(AT5G57420$end +((AT5G57420$end - AT5G57420$begin)/4))) +
  ylim(-1,1) +
  geom_point(data = AT5G57420_PAC_df, aes(x = PAC_begin, y = 0.7, size = 1.5, color = "")) +
  geom_segment(data = AT5G57420_salk_insert, aes(x = begin, xend = end, y = 0.7, yend = 0.7, size = 1, color = "")) +
  facet_wrap(~transcript)





