
####### compare half-life correlation of ribosome profile of last 20 percent CDS vs. CDS length

library(tidyverse)

####### get leadered and leaderless genes
gene_leadered <- read.table("../index/msmeg_CombinedAnnotation_CDS_leadered.bed", col.names = c("gene", "start", "end", "strand"))
gene_leaderless <- read.table("../index/msmeg_CombinedAnnotation_CDS_leaderless.bed", col.names = c("gene", "start", "end", "strand"))

####### CDS length 
CDS_length <- read.table('../index/msmeg_CDS_length.txt', header = TRUE) %>% 
        rownames_to_column('gene')

####### last 20 percent CDS ribosome profile
ribo_last20percent <- read.table('../feature/ribosome_profiling/ribo_normalizedBy_totalmRNA_3pLast20percent.txt',
                                 header = FALSE, col.names = c("gene", "riboProfile_last20percent"))

###### overall correlation
CDS_length_ribo <- CDS_length %>% 
        inner_join(ribo_last20percent, by = "gene") %>% 
        drop_na() %>% 
        mutate(riboProfile_last20percent_log2 = log2(riboProfile_last20percent)) %>% 
        mutate(CDS_length_log2 = log2(CDS_length))

ggplot(CDS_length_ribo, aes(x = CDS_length_log2, y = riboProfile_last20percent_log2)) + 
        geom_point(size = 5) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 30, face = "bold"),
              legend.position = "None")
cor.test(CDS_length_ribo$CDS_length_log2, CDS_length_ribo$riboProfile_last20percent_log2, method = "spearman")

###### correlation for leadered genes
CDS_length_ribo_leadered <- CDS_length_ribo %>% 
        inner_join(gene_leadered, by = "gene")

ggplot(CDS_length_ribo_leadered, aes(x = CDS_length_log2, y = riboProfile_last20percent_log2)) + 
        geom_point(size = 5) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 30, face = "bold"),
              legend.position = "None")
cor.test(CDS_length_ribo_leadered$CDS_length_log2, CDS_length_ribo_leadered$riboProfile_last20percent_log2, method = "spearman")

###### correlation for leaderless genes
CDS_length_ribo_leaderless <- CDS_length_ribo %>% 
        inner_join(gene_leaderless, by = "gene")

ggplot(CDS_length_ribo_leaderless, aes(x = CDS_length_log2, y = riboProfile_last20percent_log2)) + 
        geom_point(size = 5) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 30, face = "bold"),
              legend.position = "None")
cor.test(CDS_length_ribo_leaderless$CDS_length_log2, CDS_length_ribo_leaderless$riboProfile_last20percent_log2, method = "spearman")


####### CDS ribosome profile
ribo_full <- read.table('../feature/ribosome_profiling/msmeg_riboProfiling_CDS.txt') %>% 
        rownames_to_column('gene')

CDS_length_ribo_full <- CDS_length %>% 
        inner_join(ribo_full, by = "gene") %>% 
        drop_na() %>% 
        mutate(riboProfile_5p_excl_log2 = log2(ribo_5p_excl)) %>% 
        mutate(CDS_length_log2 = log2(CDS_length))

ggplot(CDS_length_ribo_full, aes(x = CDS_length_log2, y = riboProfile_5p_excl_log2)) + 
        geom_point(size = 5) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 30, face = "bold"),
              legend.position = "None")
cor.test(CDS_length_ribo_full$CDS_length_log2, CDS_length_ribo_full$riboProfile_5p_excl_log2, method = "spearman")

###### correlation for leadered genes
CDS_length_ribo_full_leadered <- CDS_length_ribo_full %>% 
        inner_join(gene_leadered, by = "gene")

ggplot(CDS_length_ribo_full_leadered, aes(x = CDS_length_log2, y = riboProfile_5p_excl_log2)) + 
        geom_point(size = 5) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 30, face = "bold"),
              legend.position = "None")
cor.test(CDS_length_ribo_full_leadered$CDS_length_log2, CDS_length_ribo_full_leadered$riboProfile_5p_excl_log2, method = "spearman")

###### correlation for leaderless genes
CDS_length_ribo_full_leaderless <- CDS_length_ribo_full %>% 
        inner_join(gene_leaderless, by = "gene")

ggplot(CDS_length_ribo_full_leaderless, aes(x = CDS_length_log2, y = riboProfile_5p_excl_log2)) + 
        geom_point(size = 5) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 30, face = "bold"),
              legend.position = "None")
cor.test(CDS_length_ribo_full_leaderless$CDS_length_log2, CDS_length_ribo_full_leaderless$riboProfile_5p_excl_log2, method = "spearman")

