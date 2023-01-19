rm(list = ls(all = TRUE)) 

library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# load up insert location and PAC data

salk_insert <- read.csv("APA_files/SALK_insert_loc.csv")

PAC_df <- read.csv("filtered_PAC_df.csv")
PAC_df$locus <- PAC_df$gene_id
PAC_df$PAC_begin <- PAC_df$start
PAC_df$PAC_end <- PAC_df$end
PAC_df$feature = "PAC"
PAC_df <- select(PAC_df, PAC_begin, PAC_end, locus, feature, strand)

# create midpoint in insert and PAC

salk_insert$insert_mid = salk_insert$begin + (salk_insert$end - salk_insert$begin)

PAC_df$PAC_mid = PAC_df$PAC_begin + (PAC_df$PAC_end - PAC_df$PAC_begin)

# identify and remove SALK lines with mutiple inserts notes in SIGNAL df.
# unfortunately will not be able to determine which is the single insert present
# within our set of confirmed single insert lines

count(salk_insert, locus)

salk_insert <- salk_insert %>% filter(locus != "AT1G12820") %>% 
  filter(locus != "AT1G80390")

# PAC_df <- PAC_df %>% filter(locus != "AT1G12820") %>% 
  # filter(locus != "AT1G80390")

PAC_df <- PAC_df 
  
# select mid and locus columns and join data sets

PAC_insert_dist <- PAC_df %>% left_join(salk_insert, by = 'locus') %>% 
  select(locus, insert_mid, PAC_mid, strand)

# plot location of insert on pos vs negative strand

PAC_insert_dist %>% ggplot(aes(x = PAC_mid, y = strand, color = strand)) + 
  geom_point() + facet_wrap(~locus, scales = "free_x")

# create new column for PAC no.

PAC_count <- count(PAC_df, locus) 
PAC_count$PAC_count = PAC_count$n
PAC_count <- PAC_count %>% select(locus, PAC_count)

# calc distance
PAC_insert_dist$PAC_insert_dist <- abs(PAC_insert_dist$PAC_mid - PAC_insert_dist$insert_mid)

# calc means
means_all <- PAC_insert_dist %>% group_by(locus) %>%
  summarise(mean_PAC_insert_dist = mean(PAC_insert_dist, na.rm = T))

# do for pos and neg strand as well
PAC_count_pos <- PAC_df %>% filter(strand == "+") %>% count(locus)
PAC_count_pos$PAC_count_pos = PAC_count_pos$n
PAC_count_pos <- PAC_count_pos %>% select(locus, PAC_count_pos)

PAC_count_neg <- PAC_df %>% filter(strand == "-") %>% count(locus)
PAC_count_neg$PAC_count_neg = PAC_count_neg$n
PAC_count_neg <- PAC_count_neg %>% select(locus, PAC_count_neg)

# join means to PAC_df

PAC_df <- PAC_df %>% left_join(PAC_insert_dist)

means_pos <- PAC_df %>% filter(strand == "+") %>% group_by(locus) %>%
  summarise(mean_PAC_insert_dist_pos = mean(PAC_insert_dist, na.rm = T))
means_neg <- PAC_df %>% filter(strand == "-") %>% group_by(locus) %>%
  summarise(mean_PAC_insert_dist_neg = mean(PAC_insert_dist, na.rm = T))

# join
means <- means_all %>% left_join(means_pos) %>% left_join(means_neg) %>%
  left_join(PAC_count) %>% left_join(PAC_count_pos) %>%
  left_join(PAC_count_neg)

means = means %>% select(PAC_count, PAC_count_neg, PAC_count_pos, mean_PAC_insert_dist,
                         mean_PAC_insert_dist_pos, mean_PAC_insert_dist_neg, locus)

#replace PAC_count NAs with 0
means$PAC_count_neg = ifelse(is.na(means$PAC_count_neg), 0,means$PAC_count_neg)
means$PAC_count_pos = ifelse(is.na(means$PAC_count_pos), 0,means$PAC_count_pos)

#create PAC groups

means$PAC_strand_group <- ifelse(means$PAC_count_neg == 0, "+", "-")

# put in 

#distribution of PAC_count is relatively even

ggplot(data = mean, aes(x = PAC_count)) + geom_histogram()

# join with larger dataset

df <- read.csv("../combined_data_1_25_22.csv")

df$locus = df$locus_capt

df <- left_join(df, means)

# make sure AT5G15100 rows have 0 and not NA for PAC count columns

df$PAC_count <- ifelse(df$locus == "AT5G15100", 0, df$PAC_count)
df$PAC_count_pos <- ifelse(df$locus == "AT5G15100", 0, df$PAC_count_pos)
df$PAC_count_neg <- ifelse(df$locus == "AT5G15100", 0, df$PAC_count_neg)

# mmmm
count(df, locus, PA.site.number)

test <- df %>% select(locus, PA.site.number, PAC_count, PAC_count_pos, PAC_count_neg)
test$PAC_diff = test$PAC_count - test$PA.site.number
test$PAC_diff_pos = test$PAC_count_pos - test$PA.site.number
test$PAC_diff_neg = test$PAC_count_neg - test$PA.site.number

test <- count(test, locus, PAC_diff, PAC_diff_pos, PAC_diff_neg) %>%
  select(-n)
test


### ran 3/3/2022 write.csv(df, "combined_data_clean.csv")


