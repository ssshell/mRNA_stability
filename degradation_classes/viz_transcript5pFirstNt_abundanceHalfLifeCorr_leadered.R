
####### 5'end first and second nucleotide of the transcript
####### check correlation with half-life and initial abundance
####### only for leadered defined TSS

library(tidyverse)
library(ggbeeswarm)

####### get leadered genes
gene_leadered <- read.table("../index/msmeg_CombinedAnnotation_CDS_leadered.bed", col.names = c("gene", "start", "end", "strand"))

####### 5'end first and second nucleotide
transcript_5p_nts <- read.table('../index/msmeg_transcript_5pFirstAndSecond_nt.txt', header = TRUE,
                                col.names = c('firstNT', 'secondNT')) %>% 
        rownames_to_column('gene')

####### nt frequency of 5'end of leadered transcript
transcript_5p_nts_leadered <- transcript_5p_nts %>%
        inner_join(gene_leadered, by = 'gene') %>% 
        mutate(firstNT = factor(firstNT, levels = c("A", "G", "C", "T")))

transcript_5p_nts_leadered_freq <- transcript_5p_nts_leadered %>%
        count(firstNT) %>%
        mutate(frequency_raw = n, frequency_percentage = n / sum(n))

transcript_5p_nts_leadered_freq %>% 
        ggplot(aes(x = firstNT, y = frequency_percentage)) +
        geom_bar(stat="identity", width = 0.7, fill = "#424242") +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 10, face = "bold"),
              axis.title =  element_text(size = 9, face = "bold"),
              legend.position = "None")
ggsave("./viz_transcript5pFirstNt_leadered/transcript5pFirstNt_leadered_freq.png", width = 6, height = 7, units = "cm", dpi = 600)

####### half-life and initial abundance of log phase
###### half-life
hl_CDS_complete_class <- read.table('./half_life/CDS_halfLife_completeClass.txt', header = TRUE)
hl_CDS_logPhase <- hl_CDS_complete_class %>% 
        select(gene, HalfLife_vals_logPhase) %>% 
        drop_na()

###### initial abundance
initialAbundance_CDS_logPhase <- read.table('../feature/FeatureSet_byType_Others_logPhase/msmeg_initialAbundance_logPhase.txt', header = TRUE)
initialAbundance_CDS_logPhase_log2 <- initialAbundance_CDS_logPhase %>% 
        mutate(initialAbundance_logPhase_log2 = log2(initialAbundance_logPhase)) %>% 
        rownames_to_column('gene')

####### distributions of half-life and initial abundance for 5'end first nt (A & G)
###### half-life
transcript_5p_nts_leadered_halfLife <- transcript_5p_nts_leadered %>% 
        inner_join(hl_CDS_logPhase, by = 'gene')

transcript_5p_nts_leadered_halfLife_AG <- transcript_5p_nts_leadered_halfLife %>% 
        filter(firstNT == 'A' | firstNT == 'G')

gini_median_halfLife_AG <- transcript_5p_nts_leadered_halfLife_AG %>% 
        group_by(firstNT) %>% 
        summarise(median_val = median(HalfLife_vals_logPhase))

transcript_5p_nts_leadered_halfLife_AG %>%
        ggplot(aes(firstNT, HalfLife_vals_logPhase)) +
        geom_quasirandom(shape = 21, size = 2.5, color = "white",
                         fill = "#424242", width = 0.2) +
        geom_crossbar(data = gini_median_halfLife_AG, aes(x = firstNT, y = median_val, ymin = median_val, ymax = median_val),
                      width = 0.5, colour = "#EB636A", size = 0.5) +
        ylim(0, 10) +
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
ggsave("./viz_transcript5pFirstNt_leadered/transcript5pFirstNt_leadered_halfLife_AG.png", width = 6, height = 7, units = "cm", dpi = 600)

transcript_5p_nts_leadered_halfLife_A <- transcript_5p_nts_leadered_halfLife %>% 
        filter(firstNT == 'A')
transcript_5p_nts_leadered_halfLife_G <- transcript_5p_nts_leadered_halfLife %>% 
        filter(firstNT == 'G')
wilcox.test(transcript_5p_nts_leadered_halfLife_A$HalfLife_vals_logPhase,
            transcript_5p_nts_leadered_halfLife_G$HalfLife_vals_logPhase)

###### initial abundance
transcript_5p_nts_leadered_initialAbundance_logPhaseLog2 <- transcript_5p_nts_leadered %>% 
        inner_join(initialAbundance_CDS_logPhase_log2, by = 'gene')

transcript_5p_nts_leadered_initialAbundance_logPhaseLog2_AG <- transcript_5p_nts_leadered_initialAbundance_logPhaseLog2 %>% 
        filter(firstNT == 'A' | firstNT == 'G')

gini_median_initialAbundance_AG <- transcript_5p_nts_leadered_initialAbundance_logPhaseLog2_AG %>% 
        group_by(firstNT) %>% 
        summarise(median_val = median(initialAbundance_logPhase_log2))

transcript_5p_nts_leadered_initialAbundance_logPhaseLog2_AG %>%
        ggplot(aes(firstNT, initialAbundance_logPhase_log2)) +
        geom_quasirandom(shape = 21, size = 2.5, color = "white",
                         fill = "#424242", width = 0.2) +
        geom_crossbar(data = gini_median_initialAbundance_AG, aes(x = firstNT, y = median_val, ymin = median_val, ymax = median_val),
                      width = 0.5, colour = "#EB636A", size = 0.5) +
        ylim(0, 25) +
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
ggsave("./viz_transcript5pFirstNt_leadered/transcript5pFirstNt_leadered_initialAbundance_AG.png", width = 6, height = 7, units = "cm", dpi = 600)

transcript_5p_nts_leadered_initialAbundance_logPhaseLog2_A <- transcript_5p_nts_leadered_initialAbundance_logPhaseLog2 %>% 
        filter(firstNT == 'A')
transcript_5p_nts_leadered_initialAbundance_logPhaseLog2_G <- transcript_5p_nts_leadered_initialAbundance_logPhaseLog2 %>% 
        filter(firstNT == 'G')
wilcox.test(transcript_5p_nts_leadered_initialAbundance_logPhaseLog2_A$initialAbundance_logPhase_log2,
            transcript_5p_nts_leadered_initialAbundance_logPhaseLog2_G$initialAbundance_logPhase_log2)
