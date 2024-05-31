
####### partial coverage ratio of 5'end 300 nt to 3'end 300 nt
####### genes filtered by both ends coverage, total >= 1500 (avg. 5 coverage per nt)
####### further remove genes with length <= 300 nt
####### separately for each condition and time point

library(tidyverse)

####### CDS length 
CDS_length <- read.table('../index/msmeg_CDS_length.txt', header = TRUE) %>% 
        rownames_to_column('gene')

####### get filtered coverage ratios
###### con_Atc
con_Atc_0 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/con_Atc_0_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
con_Atc_1 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/con_Atc_1_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
con_Atc_2 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/con_Atc_2_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
con_Atc_4 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/con_Atc_4_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
con_Atc_8 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/con_Atc_8_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
con_Atc_16 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/con_Atc_16_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))

###### con_noAtc
con_noAtc_0 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/con_noAtc_0_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
con_noAtc_1 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/con_noAtc_1_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
con_noAtc_2 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/con_noAtc_2_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
con_noAtc_4 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/con_noAtc_4_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
con_noAtc_8 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/con_noAtc_8_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
con_noAtc_16 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/con_noAtc_16_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
con_noAtc_32 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/con_noAtc_32_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))

###### EKD_Atc
EKD_Atc_0 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/EKD_Atc_0_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
EKD_Atc_2 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/EKD_Atc_2_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
EKD_Atc_4 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/EKD_Atc_4_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
EKD_Atc_8 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/EKD_Atc_8_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
EKD_Atc_16 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/EKD_Atc_16_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
EKD_Atc_32 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/EKD_Atc_32_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))

###### EKD_noAtc
EKD_noAtc_0 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/EKD_noAtc_0_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
EKD_noAtc_1 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/EKD_noAtc_1_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
EKD_noAtc_2 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/EKD_noAtc_2_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
EKD_noAtc_4 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/EKD_noAtc_4_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
EKD_noAtc_8 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/EKD_noAtc_8_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
EKD_noAtc_16 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/EKD_noAtc_16_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))

###### hypoxia
hyp_0 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/hyp_0_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
hyp_3 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/hyp_3_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
hyp_6 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/hyp_6_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
hyp_9 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/hyp_9_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
hyp_15 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/hyp_15_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
hyp_30 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/hyp_30_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
hyp_60 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/hyp_60_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))

####### get ratio density distributions
###### con_Atc
con_Atc_0_filterByLength <- con_Atc_0 %>% 
        mutate(ratio_5pTo3p_log2 = log2(con_Atc_0$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(con_Atc_0_filterByLength$ratio_5pTo3p_log2)
con_Atc_0_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

con_Atc_1_filterByLength <- con_Atc_1 %>% 
        mutate(ratio_5pTo3p_log2 = log2(con_Atc_1$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(con_Atc_1_filterByLength$ratio_5pTo3p_log2)
con_Atc_1_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

con_Atc_2_filterByLength <- con_Atc_2 %>% 
        mutate(ratio_5pTo3p_log2 = log2(con_Atc_2$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(con_Atc_2_filterByLength$ratio_5pTo3p_log2)
con_Atc_2_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

con_Atc_4_filterByLength <- con_Atc_4 %>% 
        mutate(ratio_5pTo3p_log2 = log2(con_Atc_4$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(con_Atc_4_filterByLength$ratio_5pTo3p_log2)
con_Atc_4_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

con_Atc_8_filterByLength <- con_Atc_8 %>% 
        mutate(ratio_5pTo3p_log2 = log2(con_Atc_8$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(con_Atc_8_filterByLength$ratio_5pTo3p_log2)
con_Atc_8_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

con_Atc_16_filterByLength <- con_Atc_16 %>% 
        mutate(ratio_5pTo3p_log2 = log2(con_Atc_16$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(con_Atc_16_filterByLength$ratio_5pTo3p_log2)
con_Atc_16_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

###### con_noAtc
con_noAtc_0_filterByLength <- con_noAtc_0 %>% 
        mutate(ratio_5pTo3p_log2 = log2(con_noAtc_0$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(con_noAtc_0_filterByLength$ratio_5pTo3p_log2)
con_noAtc_0_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

con_noAtc_1_filterByLength <- con_noAtc_1 %>% 
        mutate(ratio_5pTo3p_log2 = log2(con_noAtc_1$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(con_noAtc_1_filterByLength$ratio_5pTo3p_log2)
con_noAtc_1_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

con_noAtc_2_filterByLength <- con_noAtc_2 %>% 
        mutate(ratio_5pTo3p_log2 = log2(con_noAtc_2$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(con_noAtc_2_filterByLength$ratio_5pTo3p_log2)
con_noAtc_2_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

con_noAtc_4_filterByLength <- con_noAtc_4 %>% 
        mutate(ratio_5pTo3p_log2 = log2(con_noAtc_4$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(con_noAtc_4_filterByLength$ratio_5pTo3p_log2)
con_noAtc_4_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

con_noAtc_8_filterByLength <- con_noAtc_8 %>% 
        mutate(ratio_5pTo3p_log2 = log2(con_noAtc_8$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(con_noAtc_8_filterByLength$ratio_5pTo3p_log2)
con_noAtc_8_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

con_noAtc_16_filterByLength <- con_noAtc_16 %>% 
        mutate(ratio_5pTo3p_log2 = log2(con_noAtc_16$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(con_noAtc_16_filterByLength$ratio_5pTo3p_log2)
con_noAtc_16_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

con_noAtc_32_filterByLength <- con_noAtc_32 %>% 
        mutate(ratio_5pTo3p_log2 = log2(con_noAtc_32$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(con_noAtc_32_filterByLength$ratio_5pTo3p_log2)
con_noAtc_32_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

###### EKD_Atc
EKD_Atc_0_filterByLength <- EKD_Atc_0 %>% 
        mutate(ratio_5pTo3p_log2 = log2(EKD_Atc_0$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(EKD_Atc_0_filterByLength$ratio_5pTo3p_log2)
EKD_Atc_0_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

EKD_Atc_2_filterByLength <- EKD_Atc_2 %>% 
        mutate(ratio_5pTo3p_log2 = log2(EKD_Atc_2$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(EKD_Atc_2_filterByLength$ratio_5pTo3p_log2)
EKD_Atc_2_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

EKD_Atc_4_filterByLength <- EKD_Atc_4 %>% 
        mutate(ratio_5pTo3p_log2 = log2(EKD_Atc_4$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(EKD_Atc_4_filterByLength$ratio_5pTo3p_log2)
EKD_Atc_4_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

EKD_Atc_8_filterByLength <- EKD_Atc_8 %>% 
        mutate(ratio_5pTo3p_log2 = log2(EKD_Atc_8$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(EKD_Atc_8_filterByLength$ratio_5pTo3p_log2)
EKD_Atc_8_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

EKD_Atc_16_filterByLength <- EKD_Atc_16 %>% 
        mutate(ratio_5pTo3p_log2 = log2(EKD_Atc_16$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(EKD_Atc_16_filterByLength$ratio_5pTo3p_log2)
EKD_Atc_16_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

EKD_Atc_32_filterByLength <- EKD_Atc_32 %>% 
        mutate(ratio_5pTo3p_log2 = log2(EKD_Atc_32$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(EKD_Atc_32_filterByLength$ratio_5pTo3p_log2)
EKD_Atc_32_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

###### EKD_noAtc
EKD_noAtc_0_filterByLength <- EKD_noAtc_0 %>% 
        mutate(ratio_5pTo3p_log2 = log2(EKD_noAtc_0$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(EKD_noAtc_0_filterByLength$ratio_5pTo3p_log2)
EKD_noAtc_0_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

EKD_noAtc_1_filterByLength <- EKD_noAtc_1 %>% 
        mutate(ratio_5pTo3p_log2 = log2(EKD_noAtc_1$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(EKD_noAtc_1_filterByLength$ratio_5pTo3p_log2)
EKD_noAtc_1_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

EKD_noAtc_2_filterByLength <- EKD_noAtc_2 %>% 
        mutate(ratio_5pTo3p_log2 = log2(EKD_noAtc_2$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(EKD_noAtc_2_filterByLength$ratio_5pTo3p_log2)
EKD_noAtc_2_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

EKD_noAtc_4_filterByLength <- EKD_noAtc_4 %>% 
        mutate(ratio_5pTo3p_log2 = log2(EKD_noAtc_4$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(EKD_noAtc_4_filterByLength$ratio_5pTo3p_log2)
EKD_noAtc_4_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

EKD_noAtc_8_filterByLength <- EKD_noAtc_8 %>% 
        mutate(ratio_5pTo3p_log2 = log2(EKD_noAtc_8$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(EKD_noAtc_8_filterByLength$ratio_5pTo3p_log2)
EKD_noAtc_8_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

EKD_noAtc_16_filterByLength <- EKD_noAtc_16 %>% 
        mutate(ratio_5pTo3p_log2 = log2(EKD_noAtc_16$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(EKD_noAtc_16_filterByLength$ratio_5pTo3p_log2)
EKD_noAtc_16_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

###### hypoxia
hyp_0_filterByLength <- hyp_0 %>% 
        mutate(ratio_5pTo3p_log2 = log2(hyp_0$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(hyp_0_filterByLength$ratio_5pTo3p_log2)
hyp_0_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

hyp_3_filterByLength <- hyp_3 %>% 
        mutate(ratio_5pTo3p_log2 = log2(hyp_3$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(hyp_3_filterByLength$ratio_5pTo3p_log2)
hyp_3_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

hyp_6_filterByLength <- hyp_6 %>% 
        mutate(ratio_5pTo3p_log2 = log2(hyp_6$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(hyp_6_filterByLength$ratio_5pTo3p_log2)
hyp_6_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

hyp_9_filterByLength <- hyp_9 %>% 
        mutate(ratio_5pTo3p_log2 = log2(hyp_9$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(hyp_9_filterByLength$ratio_5pTo3p_log2)
hyp_9_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

hyp_15_filterByLength <- hyp_15 %>% 
        mutate(ratio_5pTo3p_log2 = log2(hyp_15$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(hyp_15_filterByLength$ratio_5pTo3p_log2)
hyp_15_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

hyp_30_filterByLength <- hyp_30 %>% 
        mutate(ratio_5pTo3p_log2 = log2(hyp_30$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(hyp_30_filterByLength$ratio_5pTo3p_log2)
hyp_30_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

hyp_60_filterByLength <- hyp_60 %>% 
        mutate(ratio_5pTo3p_log2 = log2(hyp_60$ratio_5pTo3p)) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)
summary(hyp_60_filterByLength$ratio_5pTo3p_log2)
hyp_60_filterByLength %>%
        ggplot(aes(ratio_5pTo3p_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())
