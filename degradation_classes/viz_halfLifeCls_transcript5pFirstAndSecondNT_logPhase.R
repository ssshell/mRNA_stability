
####### 5'end first and second nucleotide of the transcript
####### explore potential correlation with half-life classes
####### only for leadered and leaderless with defined TSS in log phase

library(tidyverse)

####### 5'end first and second nucleotide
transcript_5p_nts <- read.table('../index/msmeg_transcript_5pFirstAndSecond_nt.txt', header = TRUE,
                                col.names = c('firstNT', 'secondNT')) %>% 
        rownames_to_column('gene')

####### complete half-life class
hl_CDS_complete_class <- read.table('./half_life/CDS_halfLife_completeClass.txt', header = TRUE)

####### overall distribution of nts over half-life classes in log phase
transcript_5p_nts_logPhase <- hl_CDS_complete_class %>% 
        select(gene, HalfLife_cls_logPhase) %>%
        inner_join(transcript_5p_nts, by = "gene") %>%
        mutate(HalfLife_cls_logPhase = factor(HalfLife_cls_logPhase)) %>% 
        mutate(firstNT = factor(firstNT)) %>% 
        mutate(secondNT = factor(secondNT)) %>% 
        drop_na()

###### first nt
transcript_5p_nts_logPhase_firstNT <- transcript_5p_nts_logPhase %>% 
        select(-secondNT) %>% 
        count(firstNT, HalfLife_cls_logPhase) %>% 
        group_by(HalfLife_cls_logPhase) %>% 
        mutate(relative_freq = n / sum(n))

transcript_5p_nts_logPhase_firstNT %>% 
        ggplot(aes(HalfLife_cls_logPhase, relative_freq)) + 
        geom_bar(stat = "identity", width = 0.7) +
        facet_wrap(~ firstNT) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text.y = element_text(size = 20, face = "bold"),
              axis.text.x = element_text(size = 20, face = "bold", angle = 30, hjust = 1),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None",
              strip.background = element_blank(),
              strip.text.x = element_text(size = 15, face = "bold")
              )

###### second nt
transcript_5p_nts_logPhase_secondNT <- transcript_5p_nts_logPhase %>% 
        select(-firstNT) %>% 
        count(secondNT, HalfLife_cls_logPhase) %>% 
        group_by(HalfLife_cls_logPhase) %>%
        mutate(relative_freq = n / sum(n))       

transcript_5p_nts_logPhase_secondNT %>% 
        ggplot(aes(HalfLife_cls_logPhase, relative_freq)) + 
        geom_bar(stat = "identity", width = 0.7) +
        facet_wrap(~ secondNT) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text.y = element_text(size = 20, face = "bold"),
              axis.text.x = element_text(size = 20, face = "bold", angle = 30, hjust = 1),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None",
              strip.background = element_blank(),
              strip.text.x = element_text(size = 15, face = "bold")
              )

####### transcript type specific distribution of nts over half-life classes in log phase
###### get leadered and leaderless genes
gene_leadered <- read.table("../index/msmeg_CombinedAnnotation_CDS_leadered.bed", col.names = c("gene", "start", "end", "strand"))
gene_leaderless <- read.table("../index/msmeg_CombinedAnnotation_CDS_leaderless.bed", col.names = c("gene", "start", "end", "strand"))

###### for leadered transcripts
transcript_5p_nts_logPhase_leadered <- transcript_5p_nts_logPhase %>% 
        inner_join(gene_leadered, by = "gene")

##### first nt
transcript_5p_nts_logPhase_leadered_firstNT <- transcript_5p_nts_logPhase_leadered %>% 
        select(-secondNT) %>% 
        count(firstNT, HalfLife_cls_logPhase) %>%
        group_by(HalfLife_cls_logPhase) %>%
        mutate(relative_freq = n / sum(n))

transcript_5p_nts_logPhase_leadered_firstNT %>% 
        ggplot(aes(HalfLife_cls_logPhase, relative_freq)) + 
        geom_bar(stat = "identity", width = 0.7) +
        facet_wrap(~ firstNT) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text.y = element_text(size = 20, face = "bold"),
              axis.text.x = element_text(size = 20, face = "bold", angle = 30, hjust = 1),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None",
              strip.background = element_blank(),
              strip.text.x = element_text(size = 15, face = "bold")
              )

##### second nt
transcript_5p_nts_logPhase_leadered_secondNT <- transcript_5p_nts_logPhase_leadered %>% 
        select(-firstNT) %>% 
        count(secondNT, HalfLife_cls_logPhase) %>%
        group_by(HalfLife_cls_logPhase) %>%
        mutate(relative_freq = n / sum(n))       

transcript_5p_nts_logPhase_leadered_secondNT %>% 
        ggplot(aes(HalfLife_cls_logPhase, relative_freq)) + 
        geom_bar(stat = "identity", width = 0.7) +
        facet_wrap(~ secondNT) +
        ylim(0, 0.4) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text.y = element_text(size = 20, face = "bold"),
              axis.text.x = element_text(size = 20, face = "bold", angle = 30, hjust = 1),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None",
              strip.background = element_blank(),
              strip.text.x = element_text(size = 15, face = "bold")
              )

###### for leaderless transcripts
transcript_5p_nts_logPhase_leaderless <- transcript_5p_nts_logPhase %>% 
        inner_join(gene_leaderless, by = "gene")

##### first nt
transcript_5p_nts_logPhase_leaderless_firstNT_p1 <- transcript_5p_nts_logPhase_leaderless %>% 
        select(-secondNT) %>% 
        count(firstNT, HalfLife_cls_logPhase) %>% 
        add_row(firstNT = 'T', HalfLife_cls_logPhase = 'Med-fast', n = 0) %>% 
        add_row(firstNT = 'T', HalfLife_cls_logPhase = 'Med-slow', n = 0) %>% 
        add_row(firstNT = 'T', HalfLife_cls_logPhase = 'Slow', n = 0) %>%
        add_row(firstNT = 'C', HalfLife_cls_logPhase = 'Fast', n = 0) %>%
        add_row(firstNT = 'C', HalfLife_cls_logPhase = 'Med-fast', n = 0) %>% 
        add_row(firstNT = 'C', HalfLife_cls_logPhase = 'Med-slow', n = 0) %>% 
        add_row(firstNT = 'C', HalfLife_cls_logPhase = 'Slow', n = 0)
        
transcript_5p_nts_logPhase_leaderless_firstNT <- transcript_5p_nts_logPhase_leaderless_firstNT_p1 %>% 
        group_by(HalfLife_cls_logPhase) %>%
        mutate(relative_freq = n / sum(n))

transcript_5p_nts_logPhase_leaderless_firstNT %>% 
        ggplot(aes(HalfLife_cls_logPhase, relative_freq)) + 
        geom_bar(stat = "identity", width = 0.7) +
        facet_wrap(~ firstNT) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text.y = element_text(size = 20, face = "bold"),
              axis.text.x = element_text(size = 20, face = "bold", angle = 30, hjust = 1),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None",
              strip.background = element_blank(),
              strip.text.x = element_text(size = 15, face = "bold")
              )

##### second nt
transcript_5p_nts_logPhase_leaderless_secondNT_p1 <- transcript_5p_nts_logPhase_leaderless %>% 
        select(-firstNT) %>% 
        count(secondNT, HalfLife_cls_logPhase) %>%
        add_row(secondNT = 'A', HalfLife_cls_logPhase = 'Fast', n = 0) %>%
        add_row(secondNT = 'A', HalfLife_cls_logPhase = 'Med-fast', n = 0) %>% 
        add_row(secondNT = 'A', HalfLife_cls_logPhase = 'Med-slow', n = 0) %>% 
        add_row(secondNT = 'A', HalfLife_cls_logPhase = 'Slow', n = 0) %>% 
        add_row(secondNT = 'C', HalfLife_cls_logPhase = 'Fast', n = 0) %>%
        add_row(secondNT = 'C', HalfLife_cls_logPhase = 'Med-fast', n = 0) %>% 
        add_row(secondNT = 'C', HalfLife_cls_logPhase = 'Med-slow', n = 0) %>% 
        add_row(secondNT = 'C', HalfLife_cls_logPhase = 'Slow', n = 0) %>% 
        add_row(secondNT = 'G', HalfLife_cls_logPhase = 'Fast', n = 0) %>%
        add_row(secondNT = 'G', HalfLife_cls_logPhase = 'Med-fast', n = 0) %>% 
        add_row(secondNT = 'G', HalfLife_cls_logPhase = 'Med-slow', n = 0) %>% 
        add_row(secondNT = 'G', HalfLife_cls_logPhase = 'Slow', n = 0)

transcript_5p_nts_logPhase_leaderless_secondNT <- transcript_5p_nts_logPhase_leaderless_secondNT_p1 %>% 
        group_by(HalfLife_cls_logPhase) %>%
        mutate(relative_freq = n / sum(n))       

transcript_5p_nts_logPhase_leaderless_secondNT %>% 
        ggplot(aes(HalfLife_cls_logPhase, relative_freq)) + 
        geom_bar(stat = "identity", width = 0.7) +
        facet_wrap(~ secondNT) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text.y = element_text(size = 20, face = "bold"),
              axis.text.x = element_text(size = 20, face = "bold", angle = 30, hjust = 1),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None",
              strip.background = element_blank(),
              strip.text.x = element_text(size = 15, face = "bold")
              )
