
####### check essentiality classes using CRISPRi overlap with half-life classes
####### Bosch, Barbara et al. Cell. 2021

library(tidyverse)
library(RColorBrewer)
library(scales)

#######
# display.brewer.all(type = 'qual')
# display.brewer.pal(n = 12, name = 'Paired')
# brewer.pal(n = 12, name = 'Paired')
# show_col("#1F78B4")

####### get color scheme for half-life classes
col_hl_cls <- c("#7DC462", "#0D95D0", "#774FA0", "#E72F52")

####### get half-life values
hl_CDS_complete_class <- read.table('./half_life/CDS_halfLife_completeClass.txt', header = TRUE)

####### get essentiality class
ess_class_p1 <- read.csv('./essentiality_class/msmeg_ess_class.csv')
ess_class <- ess_class_p1 %>% 
        select(locus_tag, crispr_ess) %>% 
        rename(gene = locus_tag) %>% 
        mutate(gene = str_replace(gene, str_sub(gene, 1, 5), "MSMEG_"))

####### combine essentiality class with half-life classes
hl_ess_CDS <- hl_CDS_complete_class %>% 
        inner_join(ess_class, by = "gene") %>% 
        select(gene, HalfLife_cls_logPhase, HalfLife_cls_hypoxia, 
               HalfLife_FCcls_hypoxiaToLogPhase, crispr_ess) %>% 
        mutate(across(c(HalfLife_cls_logPhase, HalfLife_cls_hypoxia, 
                        HalfLife_FCcls_hypoxiaToLogPhase, crispr_ess), factor))

####### essentiality distributions of frequency count for half-life class of log phase, hypoxia, fold change
###### with complete available genes in each half-life condition class
###### Vetoed
##### log phase
# hl_ess_CDS_logPhase <- hl_ess_CDS %>%
#         select(gene, HalfLife_cls_logPhase, crispr_ess) %>%
#         drop_na()
# hl_ess_CDS_logPhase %>%
#         group_by(HalfLife_cls_logPhase, crispr_ess) %>%
#         summarise(freq = n()) %>%
#         ggplot(aes(x = HalfLife_cls_logPhase, y = freq, fill = HalfLife_cls_logPhase)) +
#         geom_bar(stat="identity", width = 0.7) +
#         scale_fill_manual(values = col_hl_cls) +
#         geom_text(aes(label = freq)) +
#         facet_wrap(~ crispr_ess, ncol = 1, scales = "free_y") +
#         theme_bw() +
#         theme(
#                 strip.background = element_blank(),
#                 panel.grid.major = element_blank(),
#                 panel.grid.minor = element_blank(),
#                 panel.border = element_blank(),
#                 axis.ticks = element_blank(),
#                 strip.text.x = element_blank(),
#                 axis.text = element_blank(),
#                 axis.title = element_blank(),
#                 legend.position = "None"
#         )
# ggsave("./viz_essentialGeneInHalfLifeClass/CDS_halfLifeEss_logPhase_complete.png", width = 3, height = 4, units = "cm", dpi = 600)

##### hypoxia
# hl_ess_CDS_hypoxia <- hl_ess_CDS %>% 
#         select(gene, HalfLife_cls_hypoxia, crispr_ess) %>% 
#         drop_na()
# hl_ess_CDS_hypoxia %>% 
#         group_by(HalfLife_cls_hypoxia, crispr_ess) %>% 
#         summarise(freq = n()) %>% 
#         ggplot(aes(x = HalfLife_cls_hypoxia, y = freq, fill = HalfLife_cls_hypoxia)) +
#         geom_bar(stat="identity", width = 0.7) +
#         scale_fill_manual(values = col_hl_cls) +
#         geom_text(aes(label = freq)) +
#         facet_wrap(~ crispr_ess, ncol = 1, scales = "free_y") +
#         theme_bw() +
#         theme(
#                 strip.background = element_blank(),
#                 panel.grid.major = element_blank(),
#                 panel.grid.minor = element_blank(),
#                 panel.border = element_blank(),
#                 axis.ticks = element_blank(),
#                 strip.text.x = element_blank(),
#                 axis.text = element_blank(),
#                 axis.title = element_blank(),
#                 legend.position = "None"
#         )
# ggsave("./viz_essentialGeneInHalfLifeClass/CDS_halfLifeEss_hypoxia_complete.png", width = 3, height = 4, units = "cm", dpi = 600)

##### half-life change from log phase to hypoxia 
# hl_ess_CDS_FC <- hl_ess_CDS %>% 
#         select(gene, HalfLife_FCcls_hypoxiaToLogPhase, crispr_ess) %>% 
#         drop_na() %>% 
#         mutate(HalfLife_FCcls_hypoxiaToLogPhase = factor(HalfLife_FCcls_hypoxiaToLogPhase, 
#                                                          levels = c("Small", "Med-small", "Med-large", "Large")))
# hl_ess_CDS_FC %>% 
#         group_by(HalfLife_FCcls_hypoxiaToLogPhase, crispr_ess) %>% 
#         summarise(freq = n()) %>% 
#         ggplot(aes(x = HalfLife_FCcls_hypoxiaToLogPhase, y = freq, fill = HalfLife_FCcls_hypoxiaToLogPhase)) +
#         geom_bar(stat="identity", width = 0.7) +
#         scale_fill_manual(values = col_hl_cls) +
#         geom_text(aes(label = freq)) +
#         facet_wrap(~ crispr_ess, ncol = 1, scales = "free_y") +
#         theme_bw() +
#         theme(
#                 strip.background = element_blank(),
#                 panel.grid.major = element_blank(),
#                 panel.grid.minor = element_blank(),
#                 panel.border = element_blank(),
#                 axis.ticks = element_blank(),
#                 strip.text.x = element_blank(),
#                 axis.text = element_blank(),
#                 axis.title = element_blank(),
#                 legend.position = "None"
#         )
# ggsave("./viz_essentialGeneInHalfLifeClass/CDS_halfLifeFcEss_complete.png", width = 3, height = 4, units = "cm", dpi = 600)

###### with only available genes in all class
###### Vetoed
# hl_ess_CDS_common <- hl_ess_CDS %>% 
#         drop_na()

##### log phase
# hl_ess_CDS_common %>% 
#         select(gene, HalfLife_cls_logPhase, crispr_ess) %>% 
#         group_by(HalfLife_cls_logPhase, crispr_ess) %>% 
#         summarise(freq = n())
# hl_ess_CDS_common %>% 
#         select(gene, HalfLife_cls_logPhase, crispr_ess) %>%
#         group_by(HalfLife_cls_logPhase, crispr_ess) %>% 
#         summarise(freq = n()) %>% 
#         ggplot(aes(x = HalfLife_cls_logPhase, y = freq, fill = HalfLife_cls_logPhase)) +
#         geom_bar(stat="identity", width = 0.7) +
#         scale_fill_manual(values = col_hl_cls) +
#         geom_text(aes(label = freq)) +
#         facet_wrap(~ crispr_ess, ncol = 1, scales = "free_y") +
#         theme_bw() +
#         theme(
#                 strip.background = element_blank(),
#                 panel.grid.major = element_blank(),
#                 panel.grid.minor = element_blank(),
#                 panel.border = element_blank(),
#                 axis.ticks = element_blank(),
#                 strip.text.x = element_blank(),
#                 axis.text = element_blank(),
#                 axis.title = element_blank()
#         )
# ggsave("./viz_essentialGeneInHalfLifeClass/CDS_halfLifeEss_logPhase_common.png", width = 3, height = 4, units = "cm", dpi = 600)

##### hypoxia
# hl_ess_CDS_common %>% 
#         select(gene, HalfLife_cls_hypoxia, crispr_ess) %>% 
#         group_by(HalfLife_cls_hypoxia, crispr_ess) %>% 
#         summarise(freq = n())
# hl_ess_CDS_common %>% 
#         select(HalfLife_cls_hypoxia, crispr_ess) %>% 
#         group_by(HalfLife_cls_hypoxia, crispr_ess) %>% 
#         summarise(freq = n()) %>% 
#         ggplot(aes(x = HalfLife_cls_hypoxia, y = freq, fill = HalfLife_cls_hypoxia)) +
#         geom_bar(stat="identity", width = 0.7) +
#         scale_fill_manual(values = col_hl_cls) +
#         geom_text(aes(label = freq)) +
#         facet_wrap(~ crispr_ess, ncol = 1, scales = "free_y") +
#         theme_bw() +
#         theme(
#                 strip.background = element_blank(),
#                 panel.grid.major = element_blank(),
#                 panel.grid.minor = element_blank(),
#                 panel.border = element_blank(),
#                 axis.ticks = element_blank(),
#                 strip.text.x = element_blank(),
#                 axis.text = element_blank(),
#                 axis.title = element_blank(),
#                 legend.position = "None"
#         )
# ggsave("./viz_essentialGeneInHalfLifeClass/CDS_halfLifeEss_hypoxia_common.png", width = 3, height = 4, units = "cm", dpi = 600)

##### half-life change from log phase to hypoxia 
# hl_ess_CDS_FC_common <- hl_ess_CDS %>% 
#         select(HalfLife_FCcls_hypoxiaToLogPhase, crispr_ess) %>% 
#         drop_na() %>% 
#         mutate(HalfLife_FCcls_hypoxiaToLogPhase = factor(HalfLife_FCcls_hypoxiaToLogPhase, 
#                                                          levels = c("Small", "Med-small", "Med-large", "Large")))
# hl_ess_CDS_FC_common %>% 
#         group_by(HalfLife_FCcls_hypoxiaToLogPhase, crispr_ess) %>% 
#         summarise(freq = n()) %>% 
#         ggplot(aes(x = HalfLife_FCcls_hypoxiaToLogPhase, y = freq, fill = HalfLife_FCcls_hypoxiaToLogPhase)) +
#         geom_bar(stat="identity", width = 0.7) +
#         scale_fill_manual(values = col_hl_cls) +
#         geom_text(aes(label = freq)) +
#         facet_wrap(~ crispr_ess, ncol = 1, scales = "free_y") +
#         theme_bw() +
#         theme(
#                 strip.background = element_blank(),
#                 panel.grid.major = element_blank(),
#                 panel.grid.minor = element_blank(),
#                 panel.border = element_blank(),
#                 axis.ticks = element_blank(),
#                 strip.text.x = element_blank(),
#                 axis.text = element_blank(),
#                 axis.title = element_blank(),
#                 legend.position = "None"
#         )
# ggsave("./viz_essentialGeneInHalfLifeClass/CDS_halfLifeFcEss_common.png", width = 3, height = 4, units = "cm", dpi = 600)

####### relative percentage of essential genes for half-life class of log phase, hypoxia, fold change
####### percentage is calculated by essential genes divided by total number of genes in each half-life class
####### with only available genes in both log phase and hypoxia class
hl_ess_CDS_common <- hl_ess_CDS %>% 
        drop_na()

###### log phase
hl_ess_CDS_common_logPhase_total <- hl_ess_CDS_common %>% 
        select(gene, HalfLife_cls_logPhase, crispr_ess) %>% 
        group_by(HalfLife_cls_logPhase) %>% 
        summarise(freq = n())

hl_ess_CDS_common_relative <- hl_ess_CDS_common %>% 
        select(gene, HalfLife_cls_logPhase, crispr_ess) %>%
        group_by(HalfLife_cls_logPhase, crispr_ess) %>% 
        summarise(freq = n()) %>% 
        filter(crispr_ess == "Essential") %>% 
        rename(freq_ess = freq) %>% 
        inner_join(hl_ess_CDS_common_logPhase_total, by = "HalfLife_cls_logPhase") %>% 
        rename(freq_total = freq) %>% 
        mutate(ess_percent = (freq_ess / freq_total) * 100) 

hl_ess_CDS_common_relative %>%      
        ggplot(aes(x = HalfLife_cls_logPhase, y = ess_percent, fill = HalfLife_cls_logPhase)) +
        geom_bar(stat = "identity", width = 0.7) +
        scale_fill_manual(values = col_hl_cls) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.ticks.x = element_blank(),
                axis.title = element_blank(),
                axis.text.y = element_text(size = 10, face = "bold"),
                axis.text.x = element_blank(),
                legend.text = element_text(size = 5),
                legend.title = element_blank(),
        ) + 
        ylim(0, 15)
ggsave("./viz_essentialGeneInHalfLifeClass/CDS_halfLifeEss_logPhase_common_relativeToEachClassTotal.png", width = 7, height = 4, units = "cm", dpi = 600)

##### Hypergeometric test for essential genes in half-life classes
#### enrichment 
p_fast <- phyper(89, 318, 3362, 1045, lower.tail = FALSE)
p_medFast <- phyper(96, 318, 3362, 934, lower.tail = FALSE)
p_medSlow <- phyper(55, 318, 3362, 876, lower.tail = FALSE)
p_slow <- phyper(74, 318, 3362, 825, lower.tail = FALSE)
pvals <- c(p_fast, p_medFast, p_medSlow, p_slow)
p.adjust(pvals, method = "fdr", n = length(pvals))

#### depletion
p_fast <- phyper(90, 318, 3362, 1045, lower.tail = TRUE)
p_medFast <- phyper(97, 318, 3362, 934, lower.tail = TRUE)
p_medSlow <- phyper(56, 318, 3362, 876, lower.tail = TRUE)
p_slow <- phyper(75, 318, 3362, 825, lower.tail = TRUE)
pvals <- c(p_fast, p_medFast, p_medSlow, p_slow)
p.adjust(pvals, method = "fdr", n = length(pvals))

###### hypoxia
hl_ess_CDS_common_hypoxia_total <- hl_ess_CDS_common %>% 
        select(gene, HalfLife_cls_hypoxia, crispr_ess) %>% 
        group_by(HalfLife_cls_hypoxia) %>% 
        summarise(freq = n())

hl_ess_CDS_common_relative <- hl_ess_CDS_common %>% 
        select(gene, HalfLife_cls_hypoxia, crispr_ess) %>%
        group_by(HalfLife_cls_hypoxia, crispr_ess) %>% 
        summarise(freq = n()) %>% 
        filter(crispr_ess == "Essential") %>% 
        rename(freq_ess = freq) %>% 
        inner_join(hl_ess_CDS_common_hypoxia_total, by = "HalfLife_cls_hypoxia") %>% 
        rename(freq_total = freq) %>% 
        mutate(ess_percent = (freq_ess / freq_total) * 100)

hl_ess_CDS_common_relative %>% 
        ggplot(aes(x = HalfLife_cls_hypoxia, y = ess_percent, fill = HalfLife_cls_hypoxia)) +
        geom_bar(stat="identity", width = 0.7) +
        scale_fill_manual(values = col_hl_cls) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.ticks.x = element_blank(),
                axis.title = element_blank(),
                axis.text.y = element_text(size = 10, face = "bold"),
                axis.text.x = element_blank(),
                legend.text = element_text(size = 5),
                legend.title = element_blank()
        ) + 
        ylim(0, 15)
ggsave("./viz_essentialGeneInHalfLifeClass/CDS_halfLifeEss_hypoxia_common_relativeToEachClassTotal.png", width = 7, height = 4, units = "cm", dpi = 600)

##### Hypergeometric test for essential genes in half-life classes
#### enrichment 
p_fast <- phyper(53, 318, 3362, 879, lower.tail = FALSE)
p_medFast <- phyper(83, 318, 3362, 949, lower.tail = FALSE)
p_medSlow <- phyper(78, 318, 3362, 964, lower.tail = FALSE)
p_slow <- phyper(100, 318, 3362, 888, lower.tail = FALSE)
pvals <- c(p_fast, p_medFast, p_medSlow, p_slow)
p.adjust(pvals, method = "fdr", n = length(pvals))

#### depletion
p_fast <- phyper(54, 318, 3362, 879, lower.tail = TRUE)
p_medFast <- phyper(84, 318, 3362, 949, lower.tail = TRUE)
p_medSlow <- phyper(79, 318, 3362, 964, lower.tail = TRUE)
p_slow <- phyper(101, 318, 3362, 888, lower.tail = TRUE)
pvals <- c(p_fast, p_medFast, p_medSlow, p_slow)
p.adjust(pvals, method = "fdr", n = length(pvals))

###### half-life change from log phase to hypoxia 
hl_ess_CDS_common_FC_total <- hl_ess_CDS_common %>% 
        select(HalfLife_FCcls_hypoxiaToLogPhase, crispr_ess) %>% 
        group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>% 
        summarise(freq = n())

hl_ess_CDS_common_relative <- hl_ess_CDS_common %>% 
        select(HalfLife_FCcls_hypoxiaToLogPhase, crispr_ess) %>%
        mutate(HalfLife_FCcls_hypoxiaToLogPhase = factor(HalfLife_FCcls_hypoxiaToLogPhase, 
                                                         levels = c("Small", "Med-small", "Med-large", "Large"))) %>%
        group_by(HalfLife_FCcls_hypoxiaToLogPhase, crispr_ess) %>% 
        summarise(freq = n()) %>% 
        filter(crispr_ess == "Essential") %>% 
        rename(freq_ess = freq) %>% 
        inner_join(hl_ess_CDS_common_FC_total, by = "HalfLife_FCcls_hypoxiaToLogPhase") %>% 
        rename(freq_total = freq) %>% 
        mutate(ess_percent = (freq_ess / freq_total) * 100)

hl_ess_CDS_common_relative %>% 
        ggplot(aes(x = HalfLife_FCcls_hypoxiaToLogPhase, y = ess_percent, fill = HalfLife_FCcls_hypoxiaToLogPhase)) +
        geom_bar(stat="identity", width = 0.7) +
        scale_fill_manual(values = col_hl_cls) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.ticks.x = element_blank(),
                axis.title = element_blank(),
                axis.text.y = element_text(size = 10, face = "bold"),
                axis.text.x = element_blank(),
                legend.text = element_text(size = 5),
                legend.title = element_blank()
        ) + 
        ylim(0, 15)
ggsave("./viz_essentialGeneInHalfLifeClass/CDS_halfLifeFcEss_common_relativeToEachClassTotal.png", width = 7, height = 4, units = "cm", dpi = 600)

##### Hypergeometric test for essential genes in half-life classes
#### enrichment 
p_small <- phyper(44, 318, 3362, 910, lower.tail = FALSE)
p_medSmall <- phyper(73, 318, 3362, 917, lower.tail = FALSE)
p_medLarge <- phyper(98, 318, 3362, 920, lower.tail = FALSE)
p_large <- phyper(99, 318, 3362, 933, lower.tail = FALSE)
pvals <- c(p_small, p_medSmall, p_medLarge, p_large)
p.adjust(pvals, method = "fdr", n = length(pvals))

#### depletion
p_small <- phyper(45, 318, 3362, 910, lower.tail = TRUE)
p_medSmall <- phyper(74, 318, 3362, 917, lower.tail = TRUE)
p_medLarge <- phyper(99, 318, 3362, 920, lower.tail = TRUE)
p_large <- phyper(100, 318, 3362, 933, lower.tail = TRUE)
pvals <- c(p_small, p_medSmall, p_medLarge, p_large)
p.adjust(pvals, method = "fdr", n = length(pvals))

####### relative percentage of essential genes for half-life class of log phase, hypoxia, fold change
###### percentage is calculated by essential genes divided by total number of essential or non-essential genes
###### with only available genes in all class
##### Vetoed
# hl_ess_CDS_common <- hl_ess_CDS %>% 
#         drop_na()
# 
# hl_ess_CDS_common_total <- hl_ess_CDS_common %>% 
#         select(gene, crispr_ess) %>% 
#         group_by(crispr_ess) %>% 
#         summarise(freq = n())

##### log phase
# hl_ess_CDS_common %>% 
#         select(gene, HalfLife_cls_logPhase, crispr_ess) %>%
#         group_by(HalfLife_cls_logPhase, crispr_ess) %>% 
#         summarise(freq = n()) %>% 
#         rename(freq_perHalfLifeCls = freq) %>% 
#         inner_join(hl_ess_CDS_common_total, by = "crispr_ess") %>% 
#         rename(freq_total = freq) %>% 
#         mutate(ess_percent = (freq_perHalfLifeCls / freq_total) * 100) %>% 
#         ggplot(aes(x = HalfLife_cls_logPhase, y = ess_percent, fill = HalfLife_cls_logPhase)) +
#         geom_bar(stat="identity", width = 0.7) +
#         scale_fill_manual(values = col_hl_cls) +
#         facet_wrap(~ crispr_ess, ncol = 1) +
#         theme_bw() +
#         theme(
#                 strip.background = element_blank(),
#                 panel.grid.major = element_blank(),
#                 panel.grid.minor = element_blank(),
#                 panel.border = element_blank(),
#                 axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
#                 axis.ticks.x = element_blank(),
#                 strip.text.x = element_blank(),
#                 axis.text.y = element_text(size = 7, face = "bold"),
#                 axis.text.x = element_blank(),
#                 axis.title = element_blank(),
#                 legend.position = "None"
#         ) +
#         ylim(0, 40)
# ggsave("./viz_essentialGeneInHalfLifeClass/CDS_halfLifeEss_logPhase_common_relativeToEssentialClassTotal.png", width = 7, height = 4, units = "cm", dpi = 600)

##### hypoxia
# hl_ess_CDS_common %>% 
#         select(gene, HalfLife_cls_hypoxia, crispr_ess) %>%
#         group_by(HalfLife_cls_hypoxia, crispr_ess) %>% 
#         summarise(freq = n()) %>% 
#         rename(freq_perHalfLifeCls = freq) %>% 
#         inner_join(hl_ess_CDS_common_total, by = "crispr_ess") %>% 
#         rename(freq_total = freq) %>% 
#         mutate(ess_percent = (freq_perHalfLifeCls / freq_total) * 100) %>% 
#         ggplot(aes(x = HalfLife_cls_hypoxia, y = ess_percent, fill = HalfLife_cls_hypoxia)) +
#         geom_bar(stat="identity", width = 0.7) +
#         scale_fill_manual(values = col_hl_cls) +
#         facet_wrap(~ crispr_ess, ncol = 1) +
#         theme_bw() +
#         theme(
#                 strip.background = element_blank(),
#                 panel.grid.major = element_blank(),
#                 panel.grid.minor = element_blank(),
#                 panel.border = element_blank(),
#                 axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
#                 axis.ticks.x = element_blank(),
#                 strip.text.x = element_blank(),
#                 axis.text.y = element_text(size = 7, face = "bold"),
#                 axis.text.x = element_blank(),
#                 axis.title = element_blank(),
#                 legend.position = "None"
#         ) +
#         ylim(0, 40)
# ggsave("./viz_essentialGeneInHalfLifeClass/CDS_halfLifeEss_hypoxia_common_relativeToEssentialClassTotal.png", width = 7, height = 4, units = "cm", dpi = 600)

##### half-life change from log phase to hypoxia 
# hl_ess_CDS_common %>% 
#         select(HalfLife_FCcls_hypoxiaToLogPhase, crispr_ess) %>%
#         mutate(HalfLife_FCcls_hypoxiaToLogPhase = factor(HalfLife_FCcls_hypoxiaToLogPhase, 
#                                                          levels = c("Small", "Med-small", "Med-large", "Large"))) %>%
#         group_by(HalfLife_FCcls_hypoxiaToLogPhase, crispr_ess) %>% 
#         summarise(freq = n()) %>% 
#         rename(freq_perHalfLifeCls = freq) %>% 
#         inner_join(hl_ess_CDS_common_total, by = "crispr_ess") %>% 
#         rename(freq_total = freq) %>% 
#         mutate(ess_percent = (freq_perHalfLifeCls / freq_total) * 100) %>% 
#         ggplot(aes(x = HalfLife_FCcls_hypoxiaToLogPhase, y = ess_percent, fill = HalfLife_FCcls_hypoxiaToLogPhase)) +
#         geom_bar(stat="identity", width = 0.7) +
#         scale_fill_manual(values = col_hl_cls) +
#         facet_wrap(~ crispr_ess, ncol = 1) +
#         theme_bw() +
#         theme(
#                 strip.background = element_blank(),
#                 panel.grid.major = element_blank(),
#                 panel.grid.minor = element_blank(),
#                 panel.border = element_blank(),
#                 axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
#                 axis.ticks.x = element_blank(),
#                 strip.text.x = element_blank(),
#                 axis.text.y = element_text(size = 7, face = "bold"),
#                 axis.text.x = element_blank(),
#                 axis.title = element_blank(),
#                 legend.position = "None"
#         ) +
#         ylim(0, 40)
# ggsave("./viz_essentialGeneInHalfLifeClass/CDS_halfLifeFcEss_common_relativeToEssentialClassTotal.png", width = 7, height = 4, units = "cm", dpi = 600)

####### relative percentage of essential genes for half-life class of log phase, hypoxia, fold change
####### percentage is calculated by essential genes divided by total number of genes in each half-life class
####### with all available leadered or leaderlss genes
hl_ess_CDS_common <- hl_ess_CDS %>% 
        drop_na()

###### get leadered and leaderless genes
gene_leadered <- read.table("../index/msmeg_CombinedAnnotation_CDS_leadered.bed", col.names = c("gene", "start", "end", "strand"))
gene_leaderless <- read.table("../index/msmeg_CombinedAnnotation_CDS_leaderless.bed", col.names = c("gene", "start", "end", "strand"))

hl_ess_CDS_common_leadered <- inner_join(hl_ess_CDS_common, gene_leadered, by = "gene")
hl_ess_CDS_common_leaderless <- inner_join(hl_ess_CDS_common, gene_leaderless, by = "gene")

###### leadered genes 
##### log phase
hl_ess_CDS_common_logPhase_total_leadered <- hl_ess_CDS_common_leadered %>% 
        select(gene, HalfLife_cls_logPhase, crispr_ess) %>% 
        group_by(HalfLife_cls_logPhase) %>% 
        summarise(freq = n())

hl_ess_CDS_common_relative_leadered <- hl_ess_CDS_common_leadered %>% 
        select(gene, HalfLife_cls_logPhase, crispr_ess) %>%
        group_by(HalfLife_cls_logPhase, crispr_ess) %>% 
        summarise(freq = n()) %>% 
        filter(crispr_ess == "Essential") %>% 
        rename(freq_ess = freq) %>% 
        inner_join(hl_ess_CDS_common_logPhase_total_leadered, by = "HalfLife_cls_logPhase") %>% 
        rename(freq_total = freq) %>% 
        mutate(ess_percent = (freq_ess / freq_total) * 100)

hl_ess_CDS_common_relative_leadered %>% 
        ggplot(aes(x = HalfLife_cls_logPhase, y = ess_percent, fill = HalfLife_cls_logPhase)) +
        geom_bar(stat="identity", width = 0.7) +
        scale_fill_manual(values = col_hl_cls) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.ticks.x = element_blank(),
                axis.title = element_blank(),
                axis.text.y = element_text(size = 10, face = "bold"),
                axis.text.x = element_blank(),
                legend.text = element_text(size = 5),
                legend.title = element_blank()
        ) + 
        ylim(0, 15)
ggsave("./viz_essentialGeneInHalfLifeClass/CDS_halfLifeEss_logPhase_common_relativeToEachClassTotal_leadered_complete.png", width = 7, height = 4, units = "cm", dpi = 600)

#### Hypergeometric test for essential genes in half-life classes
### enrichment 
p_fast <- phyper(46, 318, 3362, 417, lower.tail = FALSE)
p_medFast <- phyper(38, 318, 3362, 327, lower.tail = FALSE)
p_medSlow <- phyper(22, 318, 3362, 267, lower.tail = FALSE)
p_slow <- phyper(30, 318, 3362, 316, lower.tail = FALSE)
pvals <- c(p_fast, p_medFast, p_medSlow, p_slow)
p.adjust(pvals, method = "fdr", n = length(pvals))

### depletion
p_fast <- phyper(47, 318, 3362, 417, lower.tail = TRUE)
p_medFast <- phyper(39, 318, 3362, 327, lower.tail = TRUE)
p_medSlow <- phyper(23, 318, 3362, 267, lower.tail = TRUE)
p_slow <- phyper(31, 318, 3362, 316, lower.tail = TRUE)
pvals <- c(p_fast, p_medFast, p_medSlow, p_slow)
p.adjust(pvals, method = "fdr", n = length(pvals))

##### hypoxia
hl_ess_CDS_common_hypoxia_total_leadered <- hl_ess_CDS_common_leadered %>% 
        select(gene, HalfLife_cls_hypoxia, crispr_ess) %>% 
        group_by(HalfLife_cls_hypoxia) %>% 
        summarise(freq = n())

hl_ess_CDS_common_relative_leadered <- hl_ess_CDS_common_leadered %>% 
        select(gene, HalfLife_cls_hypoxia, crispr_ess) %>%
        group_by(HalfLife_cls_hypoxia, crispr_ess) %>% 
        summarise(freq = n()) %>% 
        filter(crispr_ess == "Essential") %>% 
        rename(freq_ess = freq) %>% 
        inner_join(hl_ess_CDS_common_hypoxia_total_leadered, by = "HalfLife_cls_hypoxia") %>% 
        rename(freq_total = freq) %>% 
        mutate(ess_percent = (freq_ess / freq_total) * 100)

hl_ess_CDS_common_relative_leadered %>% 
        ggplot(aes(x = HalfLife_cls_hypoxia, y = ess_percent, fill = HalfLife_cls_hypoxia)) +
        geom_bar(stat="identity", width = 0.7) +
        scale_fill_manual(values = col_hl_cls) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.ticks.x = element_blank(),
                axis.title = element_blank(),
                axis.text.y = element_text(size = 10, face = "bold"),
                axis.text.x = element_blank(),
                legend.text = element_text(size = 5),
                legend.title = element_blank()
        ) +
        ylim(0, 15)
ggsave("./viz_essentialGeneInHalfLifeClass/CDS_halfLifeEss_hypoxia_common_relativeToEachClassTotal_leadered_complete.png", width = 7, height = 4, units = "cm", dpi = 600)

#### Hypergeometric test for essential genes in half-life classes
### enrichment 
p_fast <- phyper(31, 318, 3362, 362, lower.tail = FALSE)
p_medFast <- phyper(35, 318, 3362, 362, lower.tail = FALSE)
p_medSlow <- phyper(32, 318, 3362, 338, lower.tail = FALSE)
p_slow <- phyper(38, 318, 3362, 265, lower.tail = FALSE)
pvals <- c(p_fast, p_medFast, p_medSlow, p_slow)
p.adjust(pvals, method = "fdr", n = length(pvals))

### depletion
p_fast <- phyper(32, 318, 3362, 362, lower.tail = TRUE)
p_medFast <- phyper(36, 318, 3362, 362, lower.tail = TRUE)
p_medSlow <- phyper(33, 318, 3362, 338, lower.tail = TRUE)
p_slow <- phyper(39, 318, 3362, 265, lower.tail = TRUE)
pvals <- c(p_fast, p_medFast, p_medSlow, p_slow)
p.adjust(pvals, method = "fdr", n = length(pvals))

##### half-life change from log phase to hypoxia 
hl_ess_CDS_common_FC_total_leadered <- hl_ess_CDS_common_leadered %>% 
        select(HalfLife_FCcls_hypoxiaToLogPhase, crispr_ess) %>% 
        group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>% 
        summarise(freq = n())

hl_ess_CDS_common_relative_leadered <- hl_ess_CDS_common_leadered %>% 
        select(HalfLife_FCcls_hypoxiaToLogPhase, crispr_ess) %>%
        mutate(HalfLife_FCcls_hypoxiaToLogPhase = factor(HalfLife_FCcls_hypoxiaToLogPhase, 
                                                         levels = c("Small", "Med-small", "Med-large", "Large"))) %>%
        group_by(HalfLife_FCcls_hypoxiaToLogPhase, crispr_ess) %>% 
        summarise(freq = n()) %>% 
        filter(crispr_ess == "Essential") %>% 
        rename(freq_ess = freq) %>% 
        inner_join(hl_ess_CDS_common_FC_total_leadered, by = "HalfLife_FCcls_hypoxiaToLogPhase") %>% 
        rename(freq_total = freq) %>% 
        mutate(ess_percent = (freq_ess / freq_total) * 100)

hl_ess_CDS_common_relative_leadered %>% 
        ggplot(aes(x = HalfLife_FCcls_hypoxiaToLogPhase, y = ess_percent, fill = HalfLife_FCcls_hypoxiaToLogPhase)) +
        geom_bar(stat="identity", width = 0.7) +
        scale_fill_manual(values = col_hl_cls) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.ticks.x = element_blank(),
                axis.title = element_blank(),
                axis.text.y = element_text(size = 10, face = "bold"),
                axis.text.x = element_blank(),
                legend.text = element_text(size = 5),
                legend.title = element_blank()
        ) + 
        ylim(0, 15)
ggsave("./viz_essentialGeneInHalfLifeClass/CDS_halfLifeFcEss_common_relativeToEachClassTotal_leadered_complete.png", width = 7, height = 4, units = "cm", dpi = 600)

#### Hypergeometric test for essential genes in half-life classes
### enrichment 
p_small <- phyper(22, 318, 3362, 366, lower.tail = FALSE)
p_medSmall <- phyper(35, 318, 3362, 347, lower.tail = FALSE)
p_medLarge <- phyper(40, 318, 3362, 333, lower.tail = FALSE)
p_large <- phyper(39, 318, 3362, 281, lower.tail = FALSE)
pvals <- c(p_small, p_medSmall, p_medLarge, p_large)
p.adjust(pvals, method = "fdr", n = length(pvals))

### depletion
p_small <- phyper(23, 318, 3362, 366, lower.tail = TRUE)
p_medSmall <- phyper(36, 318, 3362, 347, lower.tail = TRUE)
p_medLarge <- phyper(41, 318, 3362, 333, lower.tail = TRUE)
p_large <- phyper(40, 318, 3362, 281, lower.tail = TRUE)
pvals <- c(p_small, p_medSmall, p_medLarge, p_large)
p.adjust(pvals, method = "fdr", n = length(pvals))

###### leaderless genes 
##### log phase
hl_ess_CDS_common_logPhase_total_leaderless <- hl_ess_CDS_common_leaderless %>% 
        select(gene, HalfLife_cls_logPhase, crispr_ess) %>% 
        group_by(HalfLife_cls_logPhase) %>% 
        summarise(freq = n())

hl_ess_CDS_common_relative_leaderless <- hl_ess_CDS_common_leaderless %>% 
        select(gene, HalfLife_cls_logPhase, crispr_ess) %>%
        group_by(HalfLife_cls_logPhase, crispr_ess) %>% 
        summarise(freq = n()) %>% 
        filter(crispr_ess == "Essential") %>% 
        rename(freq_ess = freq) %>% 
        inner_join(hl_ess_CDS_common_logPhase_total_leaderless, by = "HalfLife_cls_logPhase") %>% 
        rename(freq_total = freq) %>% 
        mutate(ess_percent = (freq_ess / freq_total) * 100) 

hl_ess_CDS_common_relative_leaderless %>% 
        ggplot(aes(x = HalfLife_cls_logPhase, y = ess_percent, fill = HalfLife_cls_logPhase)) +
        geom_bar(stat="identity", width = 0.7) +
        scale_fill_manual(values = col_hl_cls) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.ticks.x = element_blank(),
                axis.title = element_blank(),
                axis.text.y = element_text(size = 10, face = "bold"),
                axis.text.x = element_blank(),
                legend.text = element_text(size = 5),
                legend.title = element_blank()
        ) + 
        ylim(0, 15)
ggsave("./viz_essentialGeneInHalfLifeClass/CDS_halfLifeEss_logPhase_common_relativeToEachClassTotal_leaderless_complete.png", width = 7, height = 4, units = "cm", dpi = 600)

#### Hypergeometric test for essential genes in half-life classes
### enrichment 
p_fast <- phyper(24, 318, 3362, 277, lower.tail = FALSE)
p_medFast <- phyper(30, 318, 3362, 212, lower.tail = FALSE)
p_medSlow <- phyper(11, 318, 3362, 188, lower.tail = FALSE)
p_slow <- phyper(10, 318, 3362, 116, lower.tail = FALSE)
pvals <- c(p_fast, p_medFast, p_medSlow, p_slow)
p.adjust(pvals, method = "fdr", n = length(pvals))

### depletion
p_fast <- phyper(25, 318, 3362, 277, lower.tail = TRUE)
p_medFast <- phyper(31, 318, 3362, 212, lower.tail = TRUE)
p_medSlow <- phyper(12, 318, 3362, 188, lower.tail = TRUE)
p_slow <- phyper(11, 318, 3362, 116, lower.tail = TRUE)
pvals <- c(p_fast, p_medFast, p_medSlow, p_slow)
p.adjust(pvals, method = "fdr", n = length(pvals))

##### hypoxia
hl_ess_CDS_common_hypoxia_total_leaderless <- hl_ess_CDS_common_leaderless %>% 
        select(gene, HalfLife_cls_hypoxia, crispr_ess) %>% 
        group_by(HalfLife_cls_hypoxia) %>% 
        summarise(freq = n())

hl_ess_CDS_common_relative_leaderless <- hl_ess_CDS_common_leaderless %>% 
        select(gene, HalfLife_cls_hypoxia, crispr_ess) %>%
        group_by(HalfLife_cls_hypoxia, crispr_ess) %>% 
        summarise(freq = n()) %>% 
        filter(crispr_ess == "Essential") %>% 
        rename(freq_ess = freq) %>% 
        inner_join(hl_ess_CDS_common_hypoxia_total_leaderless, by = "HalfLife_cls_hypoxia") %>% 
        rename(freq_total = freq) %>% 
        mutate(ess_percent = (freq_ess / freq_total) * 100)

hl_ess_CDS_common_relative_leaderless %>% 
        ggplot(aes(x = HalfLife_cls_hypoxia, y = ess_percent, fill = HalfLife_cls_hypoxia)) +
        geom_bar(stat="identity", width = 0.7) +
        scale_fill_manual(values = col_hl_cls) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.ticks.x = element_blank(),
                axis.title = element_blank(),
                axis.text.y = element_text(size = 10, face = "bold"),
                axis.text.x = element_blank(),
                legend.text = element_text(size = 5),
                legend.title = element_blank()
        ) +
        ylim(0, 15)
ggsave("./viz_essentialGeneInHalfLifeClass/CDS_halfLifeEss_hypoxia_common_relativeToEachClassTotal_leaderless_complete.png", width = 7, height = 4, units = "cm", dpi = 600)

#### Hypergeometric test for essential genes in half-life classes
### enrichment 
p_fast <- phyper(14, 318, 3362, 226, lower.tail = FALSE)
p_medFast <- phyper(28, 318, 3362, 261, lower.tail = FALSE)
p_medSlow <- phyper(23, 318, 3362, 192, lower.tail = FALSE)
p_slow <- phyper(10, 318, 3362, 114, lower.tail = FALSE)
pvals <- c(p_fast, p_medFast, p_medSlow, p_slow)
p.adjust(pvals, method = "fdr", n = length(pvals))

### depletion
p_fast <- phyper(15, 318, 3362, 226, lower.tail = TRUE)
p_medFast <- phyper(29, 318, 3362, 261, lower.tail = TRUE)
p_medSlow <- phyper(24, 318, 3362, 192, lower.tail = TRUE)
p_slow <- phyper(11, 318, 3362, 114, lower.tail = TRUE)
pvals <- c(p_fast, p_medFast, p_medSlow, p_slow)
p.adjust(pvals, method = "fdr", n = length(pvals))

##### half-life change from log phase to hypoxia 
hl_ess_CDS_common_FC_total_leaderless <- hl_ess_CDS_common_leaderless %>% 
        select(HalfLife_FCcls_hypoxiaToLogPhase, crispr_ess) %>% 
        group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>% 
        summarise(freq = n())

hl_ess_CDS_common_relative_leaderless <- hl_ess_CDS_common_leaderless %>% 
        select(HalfLife_FCcls_hypoxiaToLogPhase, crispr_ess) %>%
        mutate(HalfLife_FCcls_hypoxiaToLogPhase = factor(HalfLife_FCcls_hypoxiaToLogPhase, 
                                                         levels = c("Small", "Med-small", "Med-large", "Large"))) %>%
        group_by(HalfLife_FCcls_hypoxiaToLogPhase, crispr_ess) %>% 
        summarise(freq = n()) %>% 
        filter(crispr_ess == "Essential") %>% 
        rename(freq_ess = freq) %>% 
        inner_join(hl_ess_CDS_common_FC_total_leaderless, by = "HalfLife_FCcls_hypoxiaToLogPhase") %>% 
        rename(freq_total = freq) %>% 
        mutate(ess_percent = (freq_ess / freq_total) * 100)

hl_ess_CDS_common_relative_leaderless %>% 
        ggplot(aes(x = HalfLife_FCcls_hypoxiaToLogPhase, y = ess_percent, fill = HalfLife_FCcls_hypoxiaToLogPhase)) +
        geom_bar(stat="identity", width = 0.7) +
        scale_fill_manual(values = col_hl_cls) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.ticks.x = element_blank(),
                axis.title = element_blank(),
                axis.text.y = element_text(size = 10, face = "bold"),
                axis.text.x = element_blank(),
                legend.text = element_text(size = 5),
                legend.title = element_blank()
        ) + 
        ylim(0, 15)
ggsave("./viz_essentialGeneInHalfLifeClass/CDS_halfLifeFcEss_common_relativeToEachClassTotal_leaderless_complete.png", width = 7, height = 4, units = "cm", dpi = 600)

#### Hypergeometric test for essential genes in half-life classes
### enrichment 
p_small <- phyper(12, 318, 3362, 177, lower.tail = FALSE)
p_medSmall <- phyper(18, 318, 3362, 242, lower.tail = FALSE)
p_medLarge <- phyper(30, 318, 3362, 219, lower.tail = FALSE)
p_large <- phyper(15, 318, 3362, 155, lower.tail = FALSE)
pvals <- c(p_small, p_medSmall, p_medLarge, p_large)
p.adjust(pvals, method = "fdr", n = length(pvals))

### depletion
p_small <- phyper(13, 318, 3362, 177, lower.tail = TRUE)
p_medSmall <- phyper(19, 318, 3362, 242, lower.tail = TRUE)
p_medLarge <- phyper(31, 318, 3362, 219, lower.tail = TRUE)
p_large <- phyper(16, 318, 3362, 155, lower.tail = TRUE)
pvals <- c(p_small, p_medSmall, p_medLarge, p_large)
p.adjust(pvals, method = "fdr", n = length(pvals))
