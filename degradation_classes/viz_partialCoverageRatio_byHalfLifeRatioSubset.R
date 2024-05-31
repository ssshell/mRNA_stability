
####### get distribution of coverage ratios among half-life ratio subsets
####### split into subsets by similar half-life ratio
####### partial coverage ratios of 5'end 300 nt to 3'end 300 nt
####### half-life ratios of 5'end 300 nt to entire gene
####### only use genes with length > 300 nt

library(tidyverse)

###### CDS length 
CDS_length <- read.table('../index/msmeg_CDS_length.txt', header = TRUE) %>% 
        rownames_to_column('gene')

####### get genes with CDS length > 300 nt
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

####### get half-life ratios of 5'end vs. entire gene
###### 5'end 300 nt half-life
hl_CDS5p300nt_logPhase <- read.table('../normalization/halfLife_5p300nt_con_noAtc.txt', col.names = c('gene', 'HalfLife_5p300nt_logPhase'))
hl_CDS5p300nt_hypoxia <- read.table('../normalization/halfLife_5p300nt_hypoxia.txt', col.names = c('gene', 'HalfLife_5p300nt_hypoxia'))

##### log phase
hl_ratio_5pVScds_logPhase <- hl_CDS5p300nt_logPhase %>% 
        inner_join(hl_CDS_logPhase, by = 'gene') %>% 
        mutate(halfLifeRatio_5pVScds = HalfLife_5p300nt_logPhase / HalfLife_entireGene_logPhase) %>% 
        mutate(halfLifeRatio_5pVScds_log2 = log2(halfLifeRatio_5pVScds))

##### hypoxia
hl_ratio_5pVScds_hypoxia <- hl_CDS5p300nt_hypoxia %>% 
        inner_join(hl_CDS_hypoxia, by = 'gene') %>% 
        mutate(halfLifeRatio_5pVScds = HalfLife_5p300nt_hypoxia / HalfLife_entireGene_hypoxia) %>% 
        mutate(halfLifeRatio_5pVScds_log2 = log2(halfLifeRatio_5pVScds))

####### get coverage ratios of 5'end vs. 3'end
con_noAtc_0 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/con_noAtc_0_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))
hyp_0 <- read.table('../normalization/coverage_ratio_splitByConditionAndTimePoint_afterFilter/hyp_0_filtered.txt', col.names = c('gene', 'ratio_5pTo3p'))

con_noAtc_0$coverageRatio_5pTo3p_log2 <- log2(con_noAtc_0$ratio_5pTo3p)
corr_con_noAtc_0 <- con_noAtc_0 %>% 
        inner_join(hl_ratio_5pVScds_logPhase, by = 'gene') %>% 
        drop_na()

hyp_0$coverageRatio_5pTo3p_log2 <- log2(hyp_0$ratio_5pTo3p)
corr_hyp_0 <- hyp_0 %>% 
        inner_join(hl_ratio_5pVScds_hypoxia, by = 'gene') %>% 
        drop_na()

####### get distribution of coverage ratios among half-life ratio subsets
###### log phase
summary(corr_con_noAtc_0$halfLifeRatio_5pVScds_log2)

##### split into subsets by similar half-life ratios
corr_con_noAtc_0_halfLifeSubset <- corr_con_noAtc_0 %>% 
        mutate(halfLifeRatio_5pVScds_log2_subset = case_when(halfLifeRatio_5pVScds_log2 <= -0.181224 ~ "logPhase_halfLifeRatio_G1",
                                                             halfLifeRatio_5pVScds_log2 > -0.181224 & halfLifeRatio_5pVScds_log2 <= -0.05 ~ "logPhase_halfLifeRatio_G2",
                                                             halfLifeRatio_5pVScds_log2 > -0.05 & halfLifeRatio_5pVScds_log2 <= 0.05 ~ "logPhase_halfLifeRatio_G3",
                                                             halfLifeRatio_5pVScds_log2 > 0.05 & halfLifeRatio_5pVScds_log2 <= 0.120642 ~ "logPhase_halfLifeRatio_G4",
                                                             halfLifeRatio_5pVScds_log2 > 0.120642 ~ "logPhase_halfLifeRatio_G5"))

##### get group3 with most similar half-life (log2 ratio around 0)
table(corr_con_noAtc_0_halfLifeSubset$halfLifeRatio_5pVScds_log2_subset)
corr_con_noAtc_0_halfLifeSubset_G3 <- corr_con_noAtc_0_halfLifeSubset %>% filter(halfLifeRatio_5pVScds_log2_subset == "logPhase_halfLifeRatio_G3")
write.table(corr_con_noAtc_0_halfLifeSubset_G3, "./half_life/conNoAtc_t0_halfLifeRatioG3.txt", sep = "\t", quote = FALSE, row.names = FALSE)

##### compare gene composition with group3 from splitting by coverage ratio
#### get subset of genes with similar 5' to 3' coverage ratio (G3)
coverageSubset_G3_logPhase <- read.table("./half_life/conNoAtc_t0_endCoverageG3.txt", header = TRUE)

corr_con_noAtc_0_halfLifeSubset_G3_intersect <- inner_join(corr_con_noAtc_0_halfLifeSubset_G3, coverageSubset_G3_logPhase, by = "gene")
coverageSubset_G3_logPhase_intersect <- inner_join(coverageSubset_G3_logPhase, corr_con_noAtc_0_halfLifeSubset_G3, by = "gene")

##### get summary and visualization of coverage ratio in half-life subsets
corr_con_noAtc_0_halfLifeSubset %>% 
        group_by(halfLifeRatio_5pVScds_log2_subset) %>% 
        summarize(coverageRatio_median = median(coverageRatio_5pTo3p_log2),
                  gene_number = n())

corr_con_noAtc_0_halfLifeSubset %>% 
        ggplot(aes(x = halfLifeRatio_5pVScds_log2_subset, y = coverageRatio_5pTo3p_log2,
                   fill = halfLifeRatio_5pVScds_log2_subset)) +
        geom_violin() +
        scale_fill_manual(values = c("white", "white", "#EB636A", "white", "white")) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text.y = element_text(size = 17, face = "bold"),
              axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1),
              axis.title =  element_text(size = 17, face = "bold"),
              strip.background = element_blank(),
              legend.position = 'None'
              )
ggsave("./viz_CDSlengthCorr_byHalfLifeRatio/corr_con_noAtc_0_coverageWithHalfLifeSubset.png", width = 25, height = 15, units = "cm", dpi = 600)

###### hypoxia
summary(corr_hyp_0$halfLifeRatio_5pVScds_log2)

##### split into subsets by similar half-life ratios
corr_hyp_0_halfLifeSubset <- corr_hyp_0 %>% 
        mutate(halfLifeRatio_5pVScds_log2_subset = case_when(halfLifeRatio_5pVScds_log2 <= -0.344476 ~ "hypoxia_halfLifeRatio_G1",
                                                             halfLifeRatio_5pVScds_log2 > -0.344476 & halfLifeRatio_5pVScds_log2 <= -0.128058 ~ "hypoxia_halfLifeRatio_G2",
                                                             halfLifeRatio_5pVScds_log2 > -0.128058 & halfLifeRatio_5pVScds_log2 <= -0.05 ~ "hypoxia_halfLifeRatio_G3",
                                                             halfLifeRatio_5pVScds_log2 > -0.05 & halfLifeRatio_5pVScds_log2 <= 0.05 ~ "hypoxia_halfLifeRatio_G4",
                                                             halfLifeRatio_5pVScds_log2 > 0.05 ~ "hypoxia_halfLifeRatio_G5"))

##### get group3 with most similar half-life (log2 ratio around 0)
table(corr_hyp_0_halfLifeSubset$halfLifeRatio_5pVScds_log2_subset)
corr_hyp_0_halfLifeSubset_G4 <- corr_hyp_0_halfLifeSubset %>% filter(halfLifeRatio_5pVScds_log2_subset == "hypoxia_halfLifeRatio_G4")
write.table(corr_hyp_0_halfLifeSubset_G4, "./half_life/hypoxia_t0_halfLifeRatioG4.txt", sep = "\t", quote = FALSE, row.names = FALSE)

##### compare gene composition with group3 from splitting by coverage ratio
#### get subset of genes with similar 5' to 3' coverage ratio (G3)
coverageSubset_G3_hypoxia <- read.table("./half_life/hypoxia_t0_endCoverageG3.txt", header = TRUE)

corr_hyp_0_halfLifeSubset_G4_intersect <- inner_join(corr_hyp_0_halfLifeSubset_G4, coverageSubset_G3_hypoxia, by = "gene")
coverageSubset_G4_hypoxia_intersect <- inner_join(coverageSubset_G3_hypoxia, corr_hyp_0_halfLifeSubset_G4, by = "gene")

##### get summary and visualization of coverage ratio in half-life subsets
corr_hyp_0_halfLifeSubset %>% 
        group_by(halfLifeRatio_5pVScds_log2_subset) %>% 
        summarize(coverageRatio_median = median(coverageRatio_5pTo3p_log2),
                  gene_number = n())

corr_hyp_0_halfLifeSubset %>% 
        ggplot(aes(x = halfLifeRatio_5pVScds_log2_subset, y = coverageRatio_5pTo3p_log2,
                   fill = halfLifeRatio_5pVScds_log2_subset)) +
        geom_violin() +
        scale_fill_manual(values = c("white", "white", "white", "#EB636A", "white")) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text.y = element_text(size = 17, face = "bold"),
              axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1),
              axis.title =  element_text(size = 17, face = "bold"),
              strip.background = element_blank(),
              legend.position = 'None'
              )
ggsave("./viz_CDSlengthCorr_byHalfLifeRatio/corr_hyp_0_coverageWithHalfLifeSubset.png", width = 25, height = 15, units = "cm", dpi = 600)
