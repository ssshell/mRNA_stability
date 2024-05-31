
####### check correlation between CDS half-life and length for different gene sets
####### gene sets are defined by the combinations of gene scope, gene type, condition and half-life restriction
####### gene scope: genes with entire gene half-life; genes with 5'end half-life
####### condition: log phase; hypoxia
####### gene type: all leadered; all leaderless (all leadered in ml models; all leaderless in ml models maybe later)
####### half-life restriction: try with only similar half-life coverage ratio (coverage subset G3 or G4; CDS length > 300 nt)
####### visualize correlation by distribution over half-life class

library(tidyverse)
library(RColorBrewer)
library(PupillometryR)
library(scales)

####### get color scheme for feature distributions
col_hl_cls <- c("#7DC462", "#0D95D0", "#774FA0", "#E72F52")

####### CDS length 
CDS_length <- read.table('../index/msmeg_CDS_length.txt', header = TRUE) %>% 
        rownames_to_column('gene')

####### get 5'end 300 nt half-life
hl_CDS5p300nt_logPhase <- read.table('../normalization/halfLife_5p300nt_con_noAtc.txt', col.names = c('gene', 'HalfLife_5p300nt_logPhase'))
hl_CDS5p300nt_hypoxia <- read.table('../normalization/halfLife_5p300nt_hypoxia.txt', col.names = c('gene', 'HalfLife_5p300nt_hypoxia'))

####### get subset of genes with similar half-life ratio (G3 in log phase, G4 in hypoxia)
halfLifeSubset_G3_logPhase <- read.table("./half_life/conNoAtc_t0_halfLifeRatioG3.txt", header = TRUE)
halfLifeSubset_G4_hypoxia <- read.table("./half_life/hypoxia_t0_halfLifeRatioG4.txt", header = TRUE)

####### get leadered and leaderless genes
gene_leadered <- read.table("../index/msmeg_CombinedAnnotation_CDS_leadered.bed", col.names = c("gene", "start", "end", "strand"))
gene_leaderless <- read.table("../index/msmeg_CombinedAnnotation_CDS_leaderless.bed", col.names = c("gene", "start", "end", "strand"))

####### get half-life of entire gene
hl_CDS_complete_class <- read.table('./half_life/CDS_halfLife_completeClass.txt', header = TRUE)

###### log phase
hl_CDS_logPhase <- hl_CDS_complete_class %>% 
        select(gene, HalfLife_vals_logPhase, HalfLife_cls_logPhase) %>% 
        drop_na()

###### hypoxia 
hl_CDS_hypoxia <- hl_CDS_complete_class %>% 
        select(gene, HalfLife_vals_hypoxia, HalfLife_cls_hypoxia) %>% 
        drop_na()

####### all genes with similar half-life ratio (G3 or G4)
###### log phase
halfLifeSubset_G3_logPhase_CDSlength <- halfLifeSubset_G3_logPhase %>%
        select(gene, CDS_length) %>%
        inner_join(hl_CDS_logPhase, by = "gene") %>% 
        mutate(feature_toViz = CDS_length)

halfLifeSubset_G3_logPhase_CDSlength_den <- halfLifeSubset_G3_logPhase_CDSlength %>% 
        select(CDS_length, HalfLife_cls_logPhase) %>%
        mutate_at("HalfLife_cls_logPhase", as.factor) %>%
        group_by(HalfLife_cls_logPhase) %>%
        rename(feature_toViz = CDS_length)

halfLifeSubset_G3_logPhase_CDSlength %>%
        mutate_at("HalfLife_cls_logPhase", as.factor) %>%
        ggplot(aes(colour = HalfLife_cls_logPhase)) +
        scale_colour_manual(values = col_hl_cls) +
        scale_fill_manual(values = col_hl_cls) +
        stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
                     fun = "median", geom = "line", color = "#424242", linewidth = 6) +
        stat_summary(
                mapping = aes(HalfLife_cls_logPhase, feature_toViz),
                fun.min = function(z) {quantile(z, 0.25)},
                fun.max = function(z) {quantile(z, 0.75)},
                fun = median, size = 5,
                linewidth = 5) +
        geom_flat_violin(data = halfLifeSubset_G3_logPhase_CDSlength_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                         alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
        theme_bw() +
        theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.text.x = element_text(size = 50, face = "bold"),
              axis.text.y = element_text(size = 20, face = "bold"),
              axis.ticks = element_blank(), axis.title = element_blank(),
              legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
        scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
        ylim(-500, 2500) +
        coord_flip()
ggsave("./viz_CDSlengthCorr_byHalfLifeRatio/HalfLifeRatio_logPhaseG3_allGenes_distributionOverHalfLifeClass.png", width = 35, height = 18, units = "cm", dpi = 600)

###### hypoxia
halfLifeSubset_G4_hypoxia_CDSlength <- halfLifeSubset_G4_hypoxia %>%
        select(gene, CDS_length) %>%
        inner_join(hl_CDS_hypoxia, by = "gene") %>% 
        mutate(feature_toViz = CDS_length)

halfLifeSubset_G4_hypoxia_CDSlength_den <- halfLifeSubset_G4_hypoxia_CDSlength %>% 
        select(CDS_length, HalfLife_cls_hypoxia) %>%
        mutate_at("HalfLife_cls_hypoxia", as.factor) %>%
        group_by(HalfLife_cls_hypoxia) %>%
        rename(feature_toViz = CDS_length)

halfLifeSubset_G4_hypoxia_CDSlength %>%
        mutate_at("HalfLife_cls_hypoxia", as.factor) %>%
        ggplot(aes(colour = HalfLife_cls_hypoxia)) +
        scale_colour_manual(values = col_hl_cls) +
        scale_fill_manual(values = col_hl_cls) +
        stat_summary(mapping = aes(HalfLife_cls_hypoxia, feature_toViz, group = 1),
                     fun = "median", geom = "line", color = "#424242", linewidth = 6) +
        stat_summary(
                mapping = aes(HalfLife_cls_hypoxia, feature_toViz),
                fun.min = function(z) {quantile(z, 0.25)},
                fun.max = function(z) {quantile(z, 0.75)},
                fun = median, size = 5,
                linewidth = 5) +
        geom_flat_violin(data = halfLifeSubset_G4_hypoxia_CDSlength_den, aes(HalfLife_cls_hypoxia, feature_toViz, fill = HalfLife_cls_hypoxia),
                         alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
        theme_bw() +
        theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.text.x = element_text(size = 50, face = "bold"),
              axis.text.y = element_text(size = 20, face = "bold"),
              axis.ticks = element_blank(), axis.title = element_blank(),
              legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
        scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
        ylim(-500, 2500) +
        coord_flip()
ggsave("./viz_CDSlengthCorr_byHalfLifeRatio/HalfLifeRatio_hypoxiaG4_allGenes_distributionOverHalfLifeClass.png", width = 35, height = 18, units = "cm", dpi = 600)

####### all leadered genes with similar half-life ratio (G3 or G4)
###### log phase
halfLifeSubset_G3_logPhase_CDSlength_leadered <- halfLifeSubset_G3_logPhase %>%
        select(gene, CDS_length) %>%
        inner_join(gene_leadered, by = "gene") %>%
        inner_join(hl_CDS_logPhase, by = "gene") %>% 
        mutate(feature_toViz = CDS_length)

halfLifeSubset_G3_logPhase_CDSlength_leadered_den <- halfLifeSubset_G3_logPhase_CDSlength_leadered %>% 
        select(CDS_length, HalfLife_cls_logPhase) %>%
        mutate_at("HalfLife_cls_logPhase", as.factor) %>%
        group_by(HalfLife_cls_logPhase) %>%
        rename(feature_toViz = CDS_length)

halfLifeSubset_G3_logPhase_CDSlength_leadered %>%
        mutate_at("HalfLife_cls_logPhase", as.factor) %>%
        ggplot(aes(colour = HalfLife_cls_logPhase)) +
        scale_colour_manual(values = col_hl_cls) +
        scale_fill_manual(values = col_hl_cls) +
        stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
                     fun = "median", geom = "line", color = "#424242", linewidth = 6) +
        stat_summary(
                mapping = aes(HalfLife_cls_logPhase, feature_toViz),
                fun.min = function(z) {quantile(z, 0.25)},
                fun.max = function(z) {quantile(z, 0.75)},
                fun = median, size = 5,
                linewidth = 5) +
        geom_flat_violin(data = halfLifeSubset_G3_logPhase_CDSlength_leadered_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                         alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
        theme_bw() +
        theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.text.x = element_text(size = 50, face = "bold"),
              axis.text.y = element_text(size = 20, face = "bold"),
              axis.ticks = element_blank(), axis.title = element_blank(),
              legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
        scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
        ylim(-500, 2500) +
        coord_flip()
ggsave("./viz_CDSlengthCorr_byHalfLifeRatio/HalfLifeRatio_logPhaseG3_leadered_distributionOverHalfLifeClass.png", width = 35, height = 18, units = "cm", dpi = 600)

###### hypoxia
halfLifeSubset_G4_hypoxia_CDSlength_leadered <- halfLifeSubset_G4_hypoxia %>%
        select(gene, CDS_length) %>%
        inner_join(gene_leadered, by = "gene") %>%
        inner_join(hl_CDS_hypoxia, by = "gene") %>% 
        mutate(feature_toViz = CDS_length)

halfLifeSubset_G4_hypoxia_CDSlength_leadered_den <- halfLifeSubset_G4_hypoxia_CDSlength_leadered %>% 
        select(CDS_length, HalfLife_cls_hypoxia) %>%
        mutate_at("HalfLife_cls_hypoxia", as.factor) %>%
        group_by(HalfLife_cls_hypoxia) %>%
        rename(feature_toViz = CDS_length)

halfLifeSubset_G4_hypoxia_CDSlength_leadered %>%
        mutate_at("HalfLife_cls_hypoxia", as.factor) %>%
        ggplot(aes(colour = HalfLife_cls_hypoxia)) +
        scale_colour_manual(values = col_hl_cls) +
        scale_fill_manual(values = col_hl_cls) +
        stat_summary(mapping = aes(HalfLife_cls_hypoxia, feature_toViz, group = 1),
                     fun = "median", geom = "line", color = "#424242", linewidth = 6) +
        stat_summary(
                mapping = aes(HalfLife_cls_hypoxia, feature_toViz),
                fun.min = function(z) {quantile(z, 0.25)},
                fun.max = function(z) {quantile(z, 0.75)},
                fun = median, size = 5,
                linewidth = 5) +
        geom_flat_violin(data = halfLifeSubset_G4_hypoxia_CDSlength_leadered_den, aes(HalfLife_cls_hypoxia, feature_toViz, fill = HalfLife_cls_hypoxia),
                         alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
        theme_bw() +
        theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.text.x = element_text(size = 50, face = "bold"),
              axis.text.y = element_text(size = 20, face = "bold"),
              axis.ticks = element_blank(), axis.title = element_blank(),
              legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
        scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
        ylim(-500, 2500) +
        coord_flip()
ggsave("./viz_CDSlengthCorr_byHalfLifeRatio/HalfLifeRatio_hypoxiaG4_leadered_distributionOverHalfLifeClass.png", width = 35, height = 18, units = "cm", dpi = 600)

####### all leaderless genes with similar half-life ratio (G3 or G4)
###### log phase
halfLifeSubset_G3_logPhase_CDSlength_leaderless <- halfLifeSubset_G3_logPhase %>%
        select(gene, CDS_length) %>%
        inner_join(gene_leaderless, by = "gene") %>%
        inner_join(hl_CDS_logPhase, by = "gene") %>% 
        mutate(feature_toViz = CDS_length)

halfLifeSubset_G3_logPhase_CDSlength_leaderless_den <- halfLifeSubset_G3_logPhase_CDSlength_leaderless %>% 
        select(CDS_length, HalfLife_cls_logPhase) %>%
        mutate_at("HalfLife_cls_logPhase", as.factor) %>%
        group_by(HalfLife_cls_logPhase) %>%
        rename(feature_toViz = CDS_length)

halfLifeSubset_G3_logPhase_CDSlength_leaderless %>%
        mutate_at("HalfLife_cls_logPhase", as.factor) %>%
        ggplot(aes(colour = HalfLife_cls_logPhase)) +
        scale_colour_manual(values = col_hl_cls) +
        scale_fill_manual(values = col_hl_cls) +
        stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
                     fun = "median", geom = "line", color = "#424242", linewidth = 6) +
        stat_summary(
                mapping = aes(HalfLife_cls_logPhase, feature_toViz),
                fun.min = function(z) {quantile(z, 0.25)},
                fun.max = function(z) {quantile(z, 0.75)},
                fun = median, size = 5,
                linewidth = 5) +
        geom_flat_violin(data = halfLifeSubset_G3_logPhase_CDSlength_leaderless_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                         alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
        theme_bw() +
        theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.text.x = element_text(size = 50, face = "bold"),
              axis.text.y = element_text(size = 20, face = "bold"),
              axis.ticks = element_blank(), axis.title = element_blank(),
              legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
        scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
        ylim(-500, 2500) +
        coord_flip()
ggsave("./viz_CDSlengthCorr_byHalfLifeRatio/HalfLifeRatio_logPhaseG3_leaderless_distributionOverHalfLifeClass.png", width = 35, height = 18, units = "cm", dpi = 600)

###### hypoxia
halfLifeSubset_G4_hypoxia_CDSlength_leaderless <- halfLifeSubset_G4_hypoxia %>%
        select(gene, CDS_length) %>%
        inner_join(gene_leaderless, by = "gene") %>%
        inner_join(hl_CDS_hypoxia, by = "gene") %>% 
        mutate(feature_toViz = CDS_length)

halfLifeSubset_G4_hypoxia_CDSlength_leaderless_den <- halfLifeSubset_G4_hypoxia_CDSlength_leaderless %>% 
        select(CDS_length, HalfLife_cls_hypoxia) %>%
        mutate_at("HalfLife_cls_hypoxia", as.factor) %>%
        group_by(HalfLife_cls_hypoxia) %>%
        rename(feature_toViz = CDS_length)

halfLifeSubset_G4_hypoxia_CDSlength_leaderless %>%
        mutate_at("HalfLife_cls_hypoxia", as.factor) %>%
        ggplot(aes(colour = HalfLife_cls_hypoxia)) +
        scale_colour_manual(values = col_hl_cls) +
        scale_fill_manual(values = col_hl_cls) +
        stat_summary(mapping = aes(HalfLife_cls_hypoxia, feature_toViz, group = 1),
                     fun = "median", geom = "line", color = "#424242", linewidth = 6) +
        stat_summary(
                mapping = aes(HalfLife_cls_hypoxia, feature_toViz),
                fun.min = function(z) {quantile(z, 0.25)},
                fun.max = function(z) {quantile(z, 0.75)},
                fun = median, size = 5,
                linewidth = 5) +
        geom_flat_violin(data = halfLifeSubset_G4_hypoxia_CDSlength_leaderless_den, aes(HalfLife_cls_hypoxia, feature_toViz, fill = HalfLife_cls_hypoxia),
                         alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
        theme_bw() +
        theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.text.x = element_text(size = 50, face = "bold"),
              axis.text.y = element_text(size = 20, face = "bold"),
              axis.ticks = element_blank(), axis.title = element_blank(),
              legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
        scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
        ylim(-500, 2500) +
        coord_flip()
ggsave("./viz_CDSlengthCorr_byHalfLifeRatio/HalfLifeRatio_hypoxiaG4_leaderless_distributionOverHalfLifeClass.png", width = 35, height = 18, units = "cm", dpi = 600)
