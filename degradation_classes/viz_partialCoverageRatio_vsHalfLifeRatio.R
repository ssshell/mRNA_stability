
####### evaluate the correlation between partial coverage ratio, half-life and half-life ratio
####### partial coverage ratios of 5'end 300 nt to 3'end 300 nt
####### half-life ratios of 5'end 300 nt to entire gene
####### further remove genes with length <= 300 nt

library(tidyverse)

####### CDS length 
CDS_length <- read.table('../index/msmeg_CDS_length.txt', header = TRUE) %>% 
        rownames_to_column('gene')

####### get half-life and half-life ratios
###### entire gene with CDS length > 300 nt
hl_CDS_complete_class <- read.table('./half_life/CDS_halfLife_completeClass.txt', header = TRUE)

hl_CDS_logPhase <- hl_CDS_complete_class %>% 
        select(gene, HalfLife_vals_logPhase) %>% 
        drop_na() %>%  
        rename(HalfLife_entireGene_logPhase = HalfLife_vals_logPhase) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)

hl_CDS_hypoxia <- hl_CDS_complete_class %>% 
        select(gene, HalfLife_vals_hypoxia) %>% 
        drop_na() %>% 
        rename(HalfLife_entireGene_hypoxia = HalfLife_vals_hypoxia) %>% 
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)

hl_rne_p1 <- read.table('./half_life/RNaseE_halfLife.txt', header = TRUE, sep = "\t", col.names = c("gene", "HalfLife_vals_rneAtc", "HalfLife_vals_rneNoAtc"))
hl_rne <- hl_rne_p1 %>% 
        drop_na() %>%
        inner_join(CDS_length, by = "gene") %>% 
        filter(CDS_length > 300)

###### 5'end 300 nt
hl_CDS5p300nt_logPhase <- read.table('../normalization/halfLife_5p300nt_con_noAtc.txt', col.names = c('gene', 'HalfLife_5p300nt_logPhase'))
hl_CDS5p300nt_hypoxia <- read.table('../normalization/halfLife_5p300nt_hypoxia.txt', col.names = c('gene', 'HalfLife_5p300nt_hypoxia'))

###### half-life ratios of 5'end vs. entire gene
##### log phase
hl_ratio_5pVScds_logPhase <- hl_CDS5p300nt_logPhase %>% 
        inner_join(hl_CDS_logPhase, by = 'gene') %>% 
        mutate(halfLifeRatio_5pVScds = HalfLife_5p300nt_logPhase / HalfLife_entireGene_logPhase) %>% 
        mutate(halfLifeRatio_5pVScds_log2 = log2(halfLifeRatio_5pVScds))
summary(hl_ratio_5pVScds_logPhase$halfLifeRatio_5pVScds_log2)

hl_ratio_5pVScds_logPhase %>%
        ggplot(aes(halfLifeRatio_5pVScds_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

##### hypoxia
hl_ratio_5pVScds_hypoxia <- hl_CDS5p300nt_hypoxia %>% 
        inner_join(hl_CDS_hypoxia, by = 'gene') %>% 
        mutate(halfLifeRatio_5pVScds = HalfLife_5p300nt_hypoxia / HalfLife_entireGene_hypoxia) %>% 
        mutate(halfLifeRatio_5pVScds_log2 = log2(halfLifeRatio_5pVScds))
summary(hl_ratio_5pVScds_hypoxia$halfLifeRatio_5pVScds_log2)

hl_ratio_5pVScds_hypoxia %>%
        ggplot(aes(halfLifeRatio_5pVScds_log2)) +
        geom_density() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

####### evaluate correlations:
####### between coverage ratios and half-life ratios
####### between entire gene half-life and coverage ratios of 5'end 300 nt to 3'end 300 nt
###### get coverage ratios of 5'end vs. 3'end
con_noAtc_0 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/con_noAtc_0_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
hyp_0 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/hyp_0_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
EKD_Atc_0 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/EKD_Atc_0_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
EKD_noAtc_0 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/EKD_noAtc_0_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))

###### get correlations between coverage ratios and half-life ratios
##### log phase
con_noAtc_0$coverageRatio_5pTo3p_log2 <- log2(con_noAtc_0$ratio_5pTo3p)

corr_con_noAtc_0 <- con_noAtc_0 %>% 
        inner_join(hl_ratio_5pVScds_logPhase, by = 'gene') %>% 
        drop_na()
cor.test(corr_con_noAtc_0$coverageRatio_5pTo3p_log2, corr_con_noAtc_0$halfLifeRatio_5pVScds_log2)

corr_con_noAtc_0 %>% 
        ggplot(aes(x = coverageRatio_5pTo3p_log2, y = halfLifeRatio_5pVScds_log2)) + 
        geom_point(size = 4, stroke = 0, alpha = 0.3) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

###### hypoxia
hyp_0$coverageRatio_5pTo3p_log2 <- log2(hyp_0$ratio_5pTo3p)

corr_hyp_0 <- hyp_0 %>% 
        inner_join(hl_ratio_5pVScds_hypoxia, by = 'gene') %>% 
        drop_na()
cor.test(corr_hyp_0$coverageRatio_5pTo3p_log2, corr_hyp_0$halfLifeRatio_5pVScds_log2)

corr_hyp_0 %>% 
        ggplot(aes(x = coverageRatio_5pTo3p_log2, y = halfLifeRatio_5pVScds_log2)) + 
        geom_point(size = 4, stroke = 0, alpha = 0.3) +
        ylim(-5, 2) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

###### get correlations entire gene half-life and coverage ratios of 5'end 300 nt to 3'end 300 nt
##### log phase
corr_con_noAtc_0$HalfLife_entireGene_logPhase_log2 <- log2(corr_con_noAtc_0$HalfLife_entireGene_logPhase)
cor.test(corr_con_noAtc_0$coverageRatio_5pTo3p_log2, corr_con_noAtc_0$HalfLife_entireGene_logPhase_log2)

corr_con_noAtc_0 %>% 
        ggplot(aes(x = coverageRatio_5pTo3p_log2, y = HalfLife_entireGene_logPhase_log2)) + 
        geom_point(size = 4, stroke = 0, alpha = 0.3) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

###### hypoxia
corr_hyp_0$HalfLife_entireGene_hypoxia_log2 <- log2(corr_hyp_0$HalfLife_entireGene_hypoxia)
cor.test(corr_hyp_0$coverageRatio_5pTo3p_log2, corr_hyp_0$HalfLife_entireGene_hypoxia_log2)

corr_hyp_0 %>% 
        ggplot(aes(x = coverageRatio_5pTo3p_log2, y = HalfLife_entireGene_hypoxia_log2)) + 
        geom_point(size = 4, stroke = 0, alpha = 0.3) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

###### EKD_Atc
corr_EKDatc_0 <- EKD_Atc_0 %>% 
        inner_join(hl_rne, by = "gene") %>% 
        drop_na() %>% 
        filter(HalfLife_vals_rneAtc != "stable") %>% 
        mutate(ratio_5pTo3p = as.numeric(ratio_5pTo3p)) %>% 
        mutate(HalfLife_vals_rneAtc = as.numeric(HalfLife_vals_rneAtc)) %>% 
        mutate(coverageRatio_5pTo3p_log2 = log2(ratio_5pTo3p)) %>% 
        mutate(HalfLife_entireGene_EKDatc_log2 = log2(HalfLife_vals_rneAtc))
cor.test(corr_EKDatc_0$coverageRatio_5pTo3p_log2, corr_EKDatc_0$HalfLife_entireGene_EKDatc_log2)

corr_EKDatc_0 %>% 
        ggplot(aes(x = coverageRatio_5pTo3p_log2, y = HalfLife_entireGene_EKDatc_log2)) + 
        geom_point(size = 4, stroke = 0, alpha = 0.3) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

###### EKD_noAtc
corr_EKDnoAtc_0 <- EKD_noAtc_0 %>% 
        inner_join(hl_rne, by = "gene") %>% 
        drop_na() %>% 
        filter(HalfLife_vals_rneAtc != "stable") %>% 
        mutate(ratio_5pTo3p = as.numeric(ratio_5pTo3p)) %>% 
        mutate(HalfLife_vals_rneNoAtc = as.numeric(HalfLife_vals_rneNoAtc)) %>% 
        mutate(coverageRatio_5pTo3p_log2 = log2(ratio_5pTo3p)) %>% 
        mutate(HalfLife_entireGene_EKDnoAtc_log2 = log2(HalfLife_vals_rneNoAtc))
cor.test(corr_EKDnoAtc_0$coverageRatio_5pTo3p_log2, corr_EKDnoAtc_0$HalfLife_entireGene_EKDnoAtc_log2)

corr_EKDnoAtc_0 %>% 
        ggplot(aes(x = coverageRatio_5pTo3p_log2, y = HalfLife_entireGene_EKDnoAtc_log2)) + 
        geom_point(size = 4, stroke = 0, alpha = 0.3) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

####### get distribution of half-life ratios among coverage ratio subsets
###### split into subsets by similar coverage ratios
##### log phase
summary(corr_con_noAtc_0$coverageRatio_5pTo3p_log2)

corr_con_noAtc_0_coverageSubset <- corr_con_noAtc_0 %>% 
        mutate(coverageRatio_5pTo3p_log2_subset = case_when(coverageRatio_5pTo3p_log2 <= -1 ~ "logPhase_G1",
                                                            coverageRatio_5pTo3p_log2 > -1 & coverageRatio_5pTo3p_log2 <= -0.5 ~ "logPhase_G2",
                                                            coverageRatio_5pTo3p_log2 > -0.5 & coverageRatio_5pTo3p_log2 <= 0.5 ~ "logPhase_G3",
                                                            coverageRatio_5pTo3p_log2 > 0.5 & coverageRatio_5pTo3p_log2 <= 1 ~ "logPhase_G4",
                                                            coverageRatio_5pTo3p_log2 > 1 & coverageRatio_5pTo3p_log2 <= 2 ~ "logPhase_G5",
                                                            coverageRatio_5pTo3p_log2 > 2 ~ "logPhase_G6"))

corr_con_noAtc_0_coverageSubset_G3 <- corr_con_noAtc_0_coverageSubset %>% filter(coverageRatio_5pTo3p_log2_subset == "logPhase_G3")
write.table(corr_con_noAtc_0_coverageSubset_G3, "./half_life/conNoAtc_t0_endCoverageG3.txt", sep = "\t", quote = FALSE, row.names = FALSE)

corr_con_noAtc_0_coverageSubset %>% 
        group_by(coverageRatio_5pTo3p_log2_subset) %>% 
        summarize(halfLifeRatio_median = median(halfLifeRatio_5pVScds_log2),
                  gene_number = n())

corr_con_noAtc_0_coverageSubset %>% 
        ggplot(aes(x = coverageRatio_5pTo3p_log2_subset, y = halfLifeRatio_5pVScds_log2)) +
        geom_violin() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text.y = element_text(size = 20, face = "bold"),
              axis.text.x = element_text(size = 20, face = "bold", angle = 30, hjust = 1),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())

##### hypoxia
summary(corr_hyp_0$coverageRatio_5pTo3p_log2)

corr_hyp_0_coverageSubset <- corr_hyp_0 %>% 
        mutate(coverageRatio_5pTo3p_log2_subset = case_when(coverageRatio_5pTo3p_log2 <= -1 ~ "hypoxia_G1",
                                                            coverageRatio_5pTo3p_log2 > -1 & coverageRatio_5pTo3p_log2 <= -0.5 ~ "hypoxia_G2",
                                                            coverageRatio_5pTo3p_log2 > -0.5 & coverageRatio_5pTo3p_log2 <= 0.5 ~ "hypoxia_G3",
                                                            coverageRatio_5pTo3p_log2 > 0.5 & coverageRatio_5pTo3p_log2 <= 1 ~ "hypoxia_G4",
                                                            coverageRatio_5pTo3p_log2 > 1 ~ "hypoxia_G5"))

corr_hyp_0_coverageSubset_G3 <- corr_hyp_0_coverageSubset %>% filter(coverageRatio_5pTo3p_log2_subset == "hypoxia_G3")
write.table(corr_hyp_0_coverageSubset_G3, "./half_life/hypoxia_t0_endCoverageG3.txt", sep = "\t", quote = FALSE, row.names = FALSE)

corr_hyp_0_coverageSubset %>% 
        group_by(coverageRatio_5pTo3p_log2_subset) %>% 
        summarize(halfLifeRatio_median = median(halfLifeRatio_5pVScds_log2),
                  gene_number = n())

corr_hyp_0_coverageSubset %>% 
        ggplot(aes(x = coverageRatio_5pTo3p_log2_subset, y = halfLifeRatio_5pVScds_log2)) +
        geom_violin() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text.y = element_text(size = 20, face = "bold"),
              axis.text.x = element_text(size = 20, face = "bold", angle = 30, hjust = 1),
              axis.title =  element_text(size = 25, face = "bold"),
              strip.background = element_blank())
