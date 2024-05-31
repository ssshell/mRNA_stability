
####### start codon identities for leadered transcripts
####### check correlation with half-life and initial abundance

library(tidyverse)
library(ggbeeswarm)

####### get leadered transcripts
gene_leadered <- read.table("../index/msmeg_CombinedAnnotation_CDS_leadered.bed", col.names = c("gene", "start", "end", "strand"))

####### get start codon
CDS_startCodon <- read.table('../feature/FeatureSet_byType_Translation_leadered/msmeg_CDS_StartCodon.txt', header = TRUE)
CDS_startCodon$start_codon <- colnames(CDS_startCodon)[apply(CDS_startCodon, 1, which.max)]

####### start codon frequency of leadered transcript
CDS_startCodon_leadered <- CDS_startCodon %>%
        rownames_to_column('gene') %>% 
        inner_join(gene_leadered, by = 'gene')

CDS_startCodon_leadered_freq <- CDS_startCodon_leadered %>% 
        summarize(across(start_codon_AUG:start_codon_UUG, sum)) %>% 
        gather(startCodon, raw_freq) %>% 
        mutate(frequency_total = sum(raw_freq)) %>% 
        mutate(relative_freq = raw_freq / frequency_total)

CDS_startCodon_leadered_freq %>% 
        ggplot(aes(x = startCodon, y = relative_freq)) +
        geom_bar(stat="identity", width = 0.7, fill = "#424242") +
        scale_x_discrete(labels = c('AUG', 'GUG', 'UUG')) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 10, face = "bold"),
              axis.title =  element_text(size = 9, face = "bold"),
              legend.position = "None")
ggsave("./viz_startCodon_leadered/startCodon_leadered_freq.png", width = 6, height = 7, units = "cm", dpi = 600)

####### distributions of half-life and initial abundance for start codon (AUG & GUG)
###### half-life
hl_CDS_complete_class <- read.table('./half_life/CDS_halfLife_completeClass.txt', header = TRUE)
hl_CDS_logPhase <- hl_CDS_complete_class %>% 
        select(gene, HalfLife_vals_logPhase) %>% 
        drop_na()

CDS_startCodon_leadered_halfLife <- CDS_startCodon_leadered %>% 
        inner_join(hl_CDS_logPhase, by = 'gene')

CDS_startCodon_leadered_halfLife_AUG_GUG <- CDS_startCodon_leadered_halfLife %>% 
        filter(start_codon == 'start_codon_AUG' | start_codon == 'start_codon_GUG')

gini_median_halfLife_AUG_GUG <- CDS_startCodon_leadered_halfLife_AUG_GUG %>% 
        group_by(start_codon) %>% 
        summarise(median_val = median(HalfLife_vals_logPhase))

CDS_startCodon_leadered_halfLife_AUG_GUG %>%
        ggplot(aes(start_codon, HalfLife_vals_logPhase)) +
        geom_quasirandom(shape = 21, size = 2.5, color = "white",
                         fill = "#424242", width = 0.2) +
        geom_crossbar(data = gini_median_halfLife_AUG_GUG, aes(x = start_codon, y = median_val, ymin = median_val, ymax = median_val),
                      width = 0.5, colour = "#EB636A", size = 0.5) +
        ylim(0, 10) +
        scale_x_discrete(labels = c('AUG', 'GUG')) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 10, face = "bold"),
                axis.title =  element_text(size = 9, face = "bold"),
                legend.position = "None"
        )
ggsave("./viz_startCodon_leadered/startCodon_leadered_halfLife_AUG_GUG.png", width = 6, height = 7, units = "cm", dpi = 600)

CDS_startCodon_leadered_halfLife_AUG <- CDS_startCodon_leadered_halfLife %>% 
        filter(start_codon == 'start_codon_AUG')
CDS_startCodon_leadered_halfLife_GUG <- CDS_startCodon_leadered_halfLife %>% 
        filter(start_codon == 'start_codon_GUG')
wilcox.test(CDS_startCodon_leadered_halfLife_AUG$HalfLife_vals_logPhase, CDS_startCodon_leadered_halfLife_GUG$HalfLife_vals_logPhase)

###### initial abundance
initialAbundance_CDS_logPhase <- read.table('../feature/FeatureSet_byType_Others_logPhase/msmeg_initialAbundance_logPhase.txt', header = TRUE)
initialAbundance_CDS_logPhase_log2 <- initialAbundance_CDS_logPhase %>% 
        mutate(initialAbundance_logPhase_log2 = log2(initialAbundance_logPhase)) %>% 
        rownames_to_column('gene')

CDS_startCodon_leadered_initialAbundance_logPhaseLog2 <- CDS_startCodon_leadered %>% 
        inner_join(initialAbundance_CDS_logPhase_log2, by = 'gene')

CDS_startCodon_leadered_initialAbundance_logPhaseLog2_AUG_GUG <- CDS_startCodon_leadered_initialAbundance_logPhaseLog2 %>% 
        filter(start_codon == 'start_codon_AUG' | start_codon == 'start_codon_GUG')

gini_median_initialAbundance_AUG_GUG <- CDS_startCodon_leadered_initialAbundance_logPhaseLog2_AUG_GUG %>% 
        group_by(start_codon) %>% 
        summarise(median_val = median(initialAbundance_logPhase_log2))

CDS_startCodon_leadered_initialAbundance_logPhaseLog2_AUG_GUG %>%
        ggplot(aes(start_codon, initialAbundance_logPhase_log2)) +
        geom_quasirandom(shape = 21, size = 2.5, color = "white",
                         fill = "#424242", width = 0.2) +
        geom_crossbar(data = gini_median_initialAbundance_AUG_GUG, aes(x = start_codon, y = median_val, ymin = median_val, ymax = median_val),
                      width = 0.5, colour = "#EB636A", size = 0.5) +
        ylim(0, 25) +
        scale_x_discrete(labels = c('AUG', 'GUG')) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 10, face = "bold"),
                axis.title =  element_text(size = 9, face = "bold"),
                legend.position = "None"
        )
ggsave("./viz_startCodon_leadered/startCodon_leadered_initialAbundance_AUG_GUG.png", width = 6, height = 7, units = "cm", dpi = 600)

CDS_startCodon_leadered_initialAbundance_logPhaseLog2_AUG <- CDS_startCodon_leadered_initialAbundance_logPhaseLog2 %>% 
        filter(start_codon == 'start_codon_AUG')
CDS_startCodon_leadered_initialAbundance_logPhaseLog2_GUG <- CDS_startCodon_leadered_initialAbundance_logPhaseLog2 %>% 
        filter(start_codon == 'start_codon_GUG')
wilcox.test(CDS_startCodon_leadered_initialAbundance_logPhaseLog2_AUG$initialAbundance_logPhase_log2,
            CDS_startCodon_leadered_initialAbundance_logPhaseLog2_GUG$initialAbundance_logPhase_log2)
