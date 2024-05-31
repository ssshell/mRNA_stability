
####### check correlation between CDS half-life and length for different gene sets
####### gene sets are defined by the combinations of gene scope, gene type, condition and half-life restriction
####### gene scope: genes with entire gene half-life; genes with 5'end half-life
####### condition: log phase; hypoxia
####### gene type: all leadered; all leaderless (all leadered in ml models; all leaderless in ml models maybe later)
####### half-life restriction: try with only similar half-life coverage ratio (coverage subset G3 or G4; CDS length > 300 nt)

library(tidyverse)

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

##### log phase
hl_CDS_logPhase <- hl_CDS_complete_class %>% 
        select(gene, HalfLife_vals_logPhase) %>% 
        drop_na() %>% 
        rename(HalfLife_vals = HalfLife_vals_logPhase)

##### hypoxia 
hl_CDS_hypoxia <- hl_CDS_complete_class %>% 
        select(gene, HalfLife_vals_hypoxia) %>% 
        drop_na() %>% 
        rename(HalfLife_vals = HalfLife_vals_hypoxia)

####### all genes with similar half-life ratio (G3 or G4)
###### log phase
##### half-life of entire gene
hl_CDS_logPhase_CDSlength_G3 <- hl_CDS_logPhase %>%
        inner_join(halfLifeSubset_G3_logPhase, by = "gene") %>% 
        select(HalfLife_vals, CDS_length) %>% 
        mutate(CDSlength_logPhase_log2 = log2(CDS_length)) %>% 
        mutate(HalfLife_logPhase_log2 = log2(HalfLife_vals))
cor.test(hl_CDS_logPhase_CDSlength_G3$HalfLife_logPhase_log2, hl_CDS_logPhase_CDSlength_G3$CDSlength_logPhase_log2, method = "spearman")

hl_CDS_logPhase_CDSlength_G3 %>% 
        ggplot(aes(x = HalfLife_logPhase_log2, y = CDSlength_logPhase_log2)) + 
        geom_point(size = 4, color = "#036D9C", 
                   alpha = 0.3,
                   stroke = 0) +
        ylim(5, 15) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None")

##### half-life of 5'end
halfLifeSubset_G3_logPhase_CDSlength <- halfLifeSubset_G3_logPhase %>%
        select(HalfLife_5p300nt_logPhase, CDS_length) %>% 
        mutate(CDSlength_logPhase_log2 = log2(CDS_length)) %>% 
        mutate(HalfLife5p_logPhase_log2 = log2(HalfLife_5p300nt_logPhase))
cor.test(halfLifeSubset_G3_logPhase_CDSlength$HalfLife5p_logPhase_log2, halfLifeSubset_G3_logPhase_CDSlength$CDSlength_logPhase_log2, method = "spearman")

halfLifeSubset_G3_logPhase_CDSlength %>% 
        ggplot(aes(x = HalfLife5p_logPhase_log2, y = CDSlength_logPhase_log2)) + 
        geom_point(size = 4, color = "#036D9C", 
                   alpha = 0.3,
                   stroke = 0) +
        ylim(5, 15) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None")

#### adding regression fit
halfLifeSubset_G3_logPhase_CDSlength %>% 
        ggplot(aes(x = HalfLife5p_logPhase_log2, y = CDS_length)) + 
        geom_point(size = 4, fill = "#036D9C",
                   color = "white",
                   shape = 21) +
        geom_smooth(method = lm, se = FALSE, fullrange = TRUE, colour = "#EB636A") +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None") +
        xlim(-1.5, 2.5) +
        ylim(75, 2000)
ggsave("./viz_CDSlengthCorr_byHalfLifeRatio/HalfLifeRatio_logPhaseG3_allGenes_lengthCorrWith5pEndHalfLife.png", width = 15, height = 12, units = "cm", dpi = 600)
model <- lm(halfLifeSubset_G3_logPhase_CDSlength$CDS_length ~ halfLifeSubset_G3_logPhase_CDSlength$HalfLife5p_logPhase_log2)
summary(model)
summary(model)$r.squared

###### hypoxia
##### half-life of entire gene
hl_CDS_hypoxia_CDSlength_G4 <- hl_CDS_hypoxia %>%
        inner_join(halfLifeSubset_G4_hypoxia, by = "gene") %>% 
        select(HalfLife_vals, CDS_length) %>% 
        mutate(CDSlength_hypoxia_log2 = log2(CDS_length)) %>% 
        mutate(HalfLife_hypoxia_log2 = log2(HalfLife_vals))
cor.test(hl_CDS_hypoxia_CDSlength_G4$HalfLife_hypoxia_log2, hl_CDS_hypoxia_CDSlength_G4$CDSlength_hypoxia_log2, method = "spearman")

hl_CDS_hypoxia_CDSlength_G4 %>% 
        ggplot(aes(x = HalfLife_hypoxia_log2, y = CDSlength_hypoxia_log2)) + 
        geom_point(size = 4, color = "#83B7DF", 
                   alpha = 0.3,
                   stroke = 0) +
        ylim(5, 15) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None")

##### half-life of 5'end
halfLifeSubset_G4_hypoxia_CDSlength <- halfLifeSubset_G4_hypoxia %>%
        select(HalfLife_5p300nt_hypoxia, CDS_length) %>% 
        mutate(CDSlength_hypoxia_log2 = log2(CDS_length)) %>% 
        mutate(HalfLife5p_hypoxia_log2 = log2(HalfLife_5p300nt_hypoxia))
cor.test(halfLifeSubset_G4_hypoxia_CDSlength$HalfLife5p_hypoxia_log2, halfLifeSubset_G4_hypoxia_CDSlength$CDSlength_hypoxia_log2, method = "spearman")

halfLifeSubset_G4_hypoxia_CDSlength %>% 
        ggplot(aes(x = HalfLife5p_hypoxia_log2, y = CDSlength_hypoxia_log2)) + 
        geom_point(size = 4, color = "#83B7DF", 
                   alpha = 0.3,
                   stroke = 0) +
        ylim(5, 15) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None")

#### adding regression fit
halfLifeSubset_G4_hypoxia_CDSlength %>% 
        ggplot(aes(x = HalfLife5p_hypoxia_log2, y = CDS_length)) + 
        geom_point(size = 4, fill = "#83B7DF", 
                   color = "white",
                   shape = 21) +
        geom_smooth(method = lm, se = FALSE, fullrange = TRUE, colour = "#EB636A") +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None") +
        xlim(3, 6) +
        ylim(75, 2000)
ggsave("./viz_CDSlengthCorr_byHalfLifeRatio/HalfLifeRatio_hypoxiaG4_allGenes_lengthCorrWith5pEndHalfLife.png", width = 15, height = 12, units = "cm", dpi = 600)
model <- lm(halfLifeSubset_G4_hypoxia_CDSlength$CDS_length ~ halfLifeSubset_G4_hypoxia_CDSlength$HalfLife5p_hypoxia_log2)
summary(model)
summary(model)$r.squared

####### all leadered genes with similar half-life ratio (G3 or G4)
###### log phase
##### half-life of entire gene
hl_CDS_logPhase_CDSlength_G3_leadered <- hl_CDS_logPhase %>%
        inner_join(halfLifeSubset_G3_logPhase, by = "gene") %>%
        inner_join(gene_leadered, by = "gene") %>%
        select(HalfLife_vals, CDS_length) %>% 
        mutate(CDSlength_logPhase_log2 = log2(CDS_length)) %>% 
        mutate(HalfLife_logPhase_log2 = log2(HalfLife_vals)) %>% 
        drop_na()
cor.test(hl_CDS_logPhase_CDSlength_G3_leadered$HalfLife_logPhase_log2, hl_CDS_logPhase_CDSlength_G3_leadered$CDSlength_logPhase_log2, method = "spearman")

hl_CDS_logPhase_CDSlength_G3_leadered %>% 
        ggplot(aes(x = HalfLife_logPhase_log2, y = CDSlength_logPhase_log2)) + 
        geom_point(size = 4, color = "#036D9C", 
                   alpha = 0.3,
                   stroke = 0) +
        ylim(5, 15) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None")

##### half-life of 5'end
halfLifeSubset_G3_logPhase_CDSlength_leadered <- halfLifeSubset_G3_logPhase %>%
        select(gene, HalfLife_5p300nt_logPhase, CDS_length) %>%
        inner_join(gene_leadered, by = "gene") %>%
        mutate(CDSlength_logPhase_log2 = log2(CDS_length)) %>% 
        mutate(HalfLife5p_logPhase_log2 = log2(HalfLife_5p300nt_logPhase)) %>% 
        drop_na()
cor.test(halfLifeSubset_G3_logPhase_CDSlength_leadered$HalfLife5p_logPhase_log2, halfLifeSubset_G3_logPhase_CDSlength_leadered$CDSlength_logPhase_log2, method = "spearman")

halfLifeSubset_G3_logPhase_CDSlength_leadered %>% 
        ggplot(aes(x = HalfLife5p_logPhase_log2, y = CDSlength_logPhase_log2)) + 
        geom_point(size = 4, color = "#036D9C", 
                   alpha = 0.3,
                   stroke = 0) +
        ylim(5, 15) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None")

#### adding regression fit
halfLifeSubset_G3_logPhase_CDSlength_leadered %>% 
        ggplot(aes(x = HalfLife5p_logPhase_log2, y = CDS_length)) + 
        geom_point(size = 4, fill = "#036D9C", 
                   color = "white",
                   shape = 21) +
        geom_smooth(method = lm, se = FALSE, fullrange = TRUE, colour = "#EB636A") +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None") +
        xlim(-1.5, 2.5) +
        ylim(75, 2000)
ggsave("./viz_CDSlengthCorr_byHalfLifeRatio/HalfLifeRatio_logPhaseG3_leadered_lengthCorrWith5pEndHalfLife.png", width = 15, height = 12, units = "cm", dpi = 600)
model <- lm(halfLifeSubset_G3_logPhase_CDSlength_leadered$CDS_length ~ halfLifeSubset_G3_logPhase_CDSlength_leadered$HalfLife5p_logPhase_log2)
summary(model)

###### hypoxia
##### half-life of entire gene
hl_CDS_hypoxia_CDSlength_G4_leadered <- hl_CDS_hypoxia %>%
        inner_join(halfLifeSubset_G4_hypoxia, by = "gene") %>% 
        inner_join(gene_leadered, by = "gene") %>%
        select(HalfLife_vals, CDS_length) %>% 
        mutate(CDSlength_hypoxia_log2 = log2(CDS_length)) %>% 
        mutate(HalfLife_hypoxia_log2 = log2(HalfLife_vals)) %>% 
        drop_na()
cor.test(hl_CDS_hypoxia_CDSlength_G4_leadered$HalfLife_hypoxia_log2, hl_CDS_hypoxia_CDSlength_G4_leadered$CDSlength_hypoxia_log2, method = "spearman")

hl_CDS_hypoxia_CDSlength_G4_leadered %>% 
        ggplot(aes(x = HalfLife_hypoxia_log2, y = CDSlength_hypoxia_log2)) + 
        geom_point(size = 4, color = "#83B7DF", 
                   alpha = 0.3,
                   stroke = 0) +
        ylim(5, 15) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None")

##### half-life of 5'end
halfLifeSubset_G4_hypoxia_CDSlength_leadered <- halfLifeSubset_G4_hypoxia %>%
        select(gene, HalfLife_5p300nt_hypoxia, CDS_length) %>% 
        inner_join(gene_leadered, by = "gene") %>%
        mutate(CDSlength_hypoxia_log2 = log2(CDS_length)) %>% 
        mutate(HalfLife5p_hypoxia_log2 = log2(HalfLife_5p300nt_hypoxia)) %>% 
        drop_na()
cor.test(halfLifeSubset_G4_hypoxia_CDSlength_leadered$HalfLife5p_hypoxia_log2, halfLifeSubset_G4_hypoxia_CDSlength_leadered$CDSlength_hypoxia_log2, method = "spearman")

halfLifeSubset_G4_hypoxia_CDSlength_leadered %>% 
        ggplot(aes(x = HalfLife5p_hypoxia_log2, y = CDSlength_hypoxia_log2)) + 
        geom_point(size = 4, color = "#83B7DF", 
                   alpha = 0.3,
                   stroke = 0) +
        ylim(5, 15) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None")

#### adding regression fit
halfLifeSubset_G4_hypoxia_CDSlength_leadered %>% 
        ggplot(aes(x = HalfLife5p_hypoxia_log2, y = CDS_length)) + 
        geom_point(size = 4, fill = "#83B7DF", 
                   color = "white",
                   shape = 21) +
        geom_smooth(method = lm, se = FALSE, fullrange = TRUE, colour = "#EB636A") +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None") +
        xlim(3, 6) +
        ylim(75, 2000)
ggsave("./viz_CDSlengthCorr_byHalfLifeRatio/HalfLifeRatio_hypoxiaG4_leadered_lengthCorrWith5pEndHalfLife.png", width = 15, height = 12, units = "cm", dpi = 600)
model <- lm(halfLifeSubset_G4_hypoxia_CDSlength_leadered$CDS_length ~ halfLifeSubset_G4_hypoxia_CDSlength_leadered$HalfLife5p_hypoxia_log2)
summary(model)

####### all leaderless genes with similar half-life ratio (G3 or G4)
###### log phase
##### half-life of entire gene
hl_CDS_logPhase_CDSlength_G3_leaderless <- hl_CDS_logPhase %>%
        inner_join(halfLifeSubset_G3_logPhase, by = "gene") %>%
        inner_join(gene_leaderless, by = "gene") %>%
        select(HalfLife_vals, CDS_length) %>% 
        mutate(CDSlength_logPhase_log2 = log2(CDS_length)) %>% 
        mutate(HalfLife_logPhase_log2 = log2(HalfLife_vals)) %>% 
        drop_na()
cor.test(hl_CDS_logPhase_CDSlength_G3_leaderless$HalfLife_logPhase_log2, hl_CDS_logPhase_CDSlength_G3_leaderless$CDSlength_logPhase_log2, method = "spearman")

hl_CDS_logPhase_CDSlength_G3_leaderless %>% 
        ggplot(aes(x = HalfLife_logPhase_log2, y = CDSlength_logPhase_log2)) + 
        geom_point(size = 4, color = "#036D9C", 
                   alpha = 0.3,
                   stroke = 0) +
        ylim(5, 15) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None")

##### half-life of 5'end
halfLifeSubset_G3_logPhase_CDSlength_leaderless <- halfLifeSubset_G3_logPhase %>%
        select(gene, HalfLife_5p300nt_logPhase, CDS_length) %>%
        inner_join(gene_leaderless, by = "gene") %>%
        mutate(CDSlength_logPhase_log2 = log2(CDS_length)) %>% 
        mutate(HalfLife5p_logPhase_log2 = log2(HalfLife_5p300nt_logPhase)) %>% 
        drop_na()
cor.test(halfLifeSubset_G3_logPhase_CDSlength_leaderless$HalfLife5p_logPhase_log2, halfLifeSubset_G3_logPhase_CDSlength_leaderless$CDSlength_logPhase_log2, method = "spearman")

halfLifeSubset_G3_logPhase_CDSlength_leaderless %>% 
        ggplot(aes(x = HalfLife5p_logPhase_log2, y = CDSlength_logPhase_log2)) + 
        geom_point(size = 4, color = "#036D9C", 
                   alpha = 0.3,
                   stroke = 0) +
        ylim(5, 15) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None")

#### adding regression fit
halfLifeSubset_G3_logPhase_CDSlength_leaderless %>% 
        ggplot(aes(x = HalfLife5p_logPhase_log2, y = CDS_length)) + 
        geom_point(size = 4, fill = "#036D9C", 
                   color = "white",
                   shape = 21) +
        geom_smooth(method = lm, se = FALSE, fullrange = TRUE, colour = "#EB636A") +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None") +
        xlim(-1.5, 2.5) +
        ylim(75, 2000)
ggsave("./viz_CDSlengthCorr_byHalfLifeRatio/HalfLifeRatio_logPhaseG3_leaderless_lengthCorrWith5pEndHalfLife.png", width = 15, height = 12, units = "cm", dpi = 600)
model <- lm(halfLifeSubset_G3_logPhase_CDSlength_leaderless$CDS_length ~ halfLifeSubset_G3_logPhase_CDSlength_leaderless$HalfLife5p_logPhase_log2)
summary(model)

###### hypoxia
##### half-life of entire gene
hl_CDS_hypoxia_CDSlength_G4_leaderless <- hl_CDS_hypoxia %>%
        inner_join(halfLifeSubset_G4_hypoxia, by = "gene") %>% 
        inner_join(gene_leaderless, by = "gene") %>%
        select(HalfLife_vals, CDS_length) %>% 
        mutate(CDSlength_hypoxia_log2 = log2(CDS_length)) %>% 
        mutate(HalfLife_hypoxia_log2 = log2(HalfLife_vals)) %>% 
        drop_na()
cor.test(hl_CDS_hypoxia_CDSlength_G4_leaderless$HalfLife_hypoxia_log2, hl_CDS_hypoxia_CDSlength_G4_leaderless$CDSlength_hypoxia_log2, method = "spearman")

hl_CDS_hypoxia_CDSlength_G4_leaderless %>% 
        ggplot(aes(x = HalfLife_hypoxia_log2, y = CDSlength_hypoxia_log2)) + 
        geom_point(size = 4, color = "#83B7DF", 
                   alpha = 0.3,
                   stroke = 0) +
        ylim(5, 15) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None")

##### half-life of 5'end
halfLifeSubset_G4_hypoxia_CDSlength_leaderless <- halfLifeSubset_G4_hypoxia %>%
        select(gene, HalfLife_5p300nt_hypoxia, CDS_length) %>% 
        inner_join(gene_leaderless, by = "gene") %>%
        mutate(CDSlength_hypoxia_log2 = log2(CDS_length)) %>% 
        mutate(HalfLife5p_hypoxia_log2 = log2(HalfLife_5p300nt_hypoxia)) %>% 
        drop_na()
cor.test(halfLifeSubset_G4_hypoxia_CDSlength_leaderless$HalfLife5p_hypoxia_log2, halfLifeSubset_G4_hypoxia_CDSlength_leaderless$CDSlength_hypoxia_log2, method = "spearman")

halfLifeSubset_G4_hypoxia_CDSlength_leaderless %>% 
        ggplot(aes(x = HalfLife5p_hypoxia_log2, y = CDSlength_hypoxia_log2)) + 
        geom_point(size = 4, color = "#83B7DF", 
                   alpha = 0.3,
                   stroke = 0) +
        ylim(5, 15) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None")

#### adding regression fit
halfLifeSubset_G4_hypoxia_CDSlength_leaderless %>% 
        ggplot(aes(x = HalfLife5p_hypoxia_log2, y = CDS_length)) + 
        geom_point(size = 4, fill = "#83B7DF", 
                   color = "white",
                   shape = 21) +
        geom_smooth(method = lm, se = FALSE, fullrange = TRUE, colour = "#EB636A") +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None") +
        xlim(3, 6) +
        ylim(75, 2000)
ggsave("./viz_CDSlengthCorr_byHalfLifeRatio/HalfLifeRatio_hypoxiaG4_leaderless_lengthCorrWith5pEndHalfLife.png", width = 15, height = 12, units = "cm", dpi = 600)
model <- lm(halfLifeSubset_G4_hypoxia_CDSlength_leaderless$CDS_length ~ halfLifeSubset_G4_hypoxia_CDSlength_leaderless$HalfLife5p_hypoxia_log2)
summary(model)

####### genes with similar half-life ratio (G4) in hypoxia corresponding in log phase
###### all genes
hl_CDS_hypoxia_CDSlength_G4_forLogPhase <- hl_CDS_hypoxia %>%
        inner_join(halfLifeSubset_G4_hypoxia, by = "gene") %>% 
        select(gene, HalfLife_entireGene_hypoxia)

hl_CDS_logPhase_CDSlength_G4_byHypoxia <- hl_CDS5p300nt_logPhase %>%
        inner_join(hl_CDS_hypoxia_CDSlength_G4_forLogPhase, by = "gene") %>% 
        inner_join(CDS_length, by = "gene") %>% 
        mutate(CDSlength_logPhase_log2 = log2(CDS_length)) %>% 
        mutate(HalfLife5p_logPhase_log2 = log2(HalfLife_5p300nt_logPhase))

cor.test(hl_CDS_logPhase_CDSlength_G4_byHypoxia$HalfLife5p_logPhase_log2, hl_CDS_logPhase_CDSlength_G4_byHypoxia$CDSlength_logPhase_log2, method = "spearman")

hl_CDS_logPhase_CDSlength_G4_byHypoxia %>% 
        ggplot(aes(x = HalfLife5p_logPhase_log2, y = CDSlength_logPhase_log2)) + 
        geom_point(size = 4, color = "#036D9C", 
                   alpha = 0.3,
                   stroke = 0) +
        ylim(5, 15) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None")

#### adding regression fit
hl_CDS_logPhase_CDSlength_G4_byHypoxia %>% 
        ggplot(aes(x = HalfLife5p_logPhase_log2, y = CDS_length)) + 
        geom_point(size = 4, fill = "#036D9C", 
                   color = "white",
                   shape = 21) +
        geom_smooth(method = lm, se = FALSE, fullrange = TRUE, colour = "#EB636A") +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None") +
        xlim(-1.5, 2.5) +
        ylim(75, 2000)
ggsave("./viz_CDSlengthCorr_byHalfLifeRatio/HalfLifeRatio_logPhaseByHypoxiaG4_allGenes_lengthCorrWith5pEndHalfLife.png", width = 15, height = 12, units = "cm", dpi = 600)
model <- lm(hl_CDS_logPhase_CDSlength_G4_byHypoxia$CDS_length ~ hl_CDS_logPhase_CDSlength_G4_byHypoxia$HalfLife5p_logPhase_log2)
summary(model)

###### all leadered genes
hl_CDS_hypoxia_CDSlength_G4_forLogPhase_leadered <- hl_CDS_hypoxia %>%
        inner_join(halfLifeSubset_G4_hypoxia, by = "gene") %>% 
        inner_join(gene_leadered, by = "gene") %>% 
        select(gene, HalfLife_entireGene_hypoxia)

hl_CDS_logPhase_CDSlength_G4_byHypoxia_leadered <- hl_CDS5p300nt_logPhase %>%
        inner_join(hl_CDS_hypoxia_CDSlength_G4_forLogPhase_leadered, by = "gene") %>% 
        inner_join(CDS_length, by = "gene") %>% 
        mutate(CDSlength_logPhase_log2 = log2(CDS_length)) %>% 
        mutate(HalfLife5p_logPhase_log2 = log2(HalfLife_5p300nt_logPhase))

cor.test(hl_CDS_logPhase_CDSlength_G4_byHypoxia_leadered$HalfLife5p_logPhase_log2, hl_CDS_logPhase_CDSlength_G4_byHypoxia_leadered$CDSlength_logPhase_log2, method = "spearman")

hl_CDS_logPhase_CDSlength_G4_byHypoxia_leadered %>% 
        ggplot(aes(x = HalfLife5p_logPhase_log2, y = CDSlength_logPhase_log2)) + 
        geom_point(size = 4, color = "#036D9C", 
                   alpha = 0.3,
                   stroke = 0) +
        ylim(5, 15) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None")

#### adding regression fit
hl_CDS_logPhase_CDSlength_G4_byHypoxia_leadered %>% 
        ggplot(aes(x = HalfLife5p_logPhase_log2, y = CDS_length)) + 
        geom_point(size = 4, fill = "#036D9C", 
                   color = "white",
                   shape = 21) +
        geom_smooth(method = lm, se = FALSE, fullrange = TRUE, colour = "#EB636A") +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None") +
        xlim(-1.5, 2.5) +
        ylim(75, 2000)
ggsave("./viz_CDSlengthCorr_byHalfLifeRatio/HalfLifeRatio_logPhaseByHypoxiaG4_leadered_lengthCorrWith5pEndHalfLife.png", width = 15, height = 12, units = "cm", dpi = 600)
model <- lm(hl_CDS_logPhase_CDSlength_G4_byHypoxia_leadered$CDS_length ~ hl_CDS_logPhase_CDSlength_G4_byHypoxia_leadered$HalfLife5p_logPhase_log2)
summary(model)

###### all leaderless genes
hl_CDS_hypoxia_CDSlength_G4_forLogPhase_leaderless <- hl_CDS_hypoxia %>%
        inner_join(halfLifeSubset_G4_hypoxia, by = "gene") %>% 
        inner_join(gene_leaderless, by = "gene") %>% 
        select(gene, HalfLife_entireGene_hypoxia)

hl_CDS_logPhase_CDSlength_G4_byHypoxia_leaderless <- hl_CDS5p300nt_logPhase %>%
        inner_join(hl_CDS_hypoxia_CDSlength_G4_forLogPhase_leaderless, by = "gene") %>% 
        inner_join(CDS_length, by = "gene") %>% 
        mutate(CDSlength_logPhase_log2 = log2(CDS_length)) %>% 
        mutate(HalfLife5p_logPhase_log2 = log2(HalfLife_5p300nt_logPhase))

cor.test(hl_CDS_logPhase_CDSlength_G4_byHypoxia_leaderless$HalfLife5p_logPhase_log2, hl_CDS_logPhase_CDSlength_G4_byHypoxia_leaderless$CDSlength_logPhase_log2, method = "spearman")

hl_CDS_logPhase_CDSlength_G4_byHypoxia_leaderless %>% 
        ggplot(aes(x = HalfLife5p_logPhase_log2, y = CDSlength_logPhase_log2)) + 
        geom_point(size = 4, color = "#036D9C", 
                   alpha = 0.3,
                   stroke = 0) +
        ylim(5, 15) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None")

#### adding regression fit
hl_CDS_logPhase_CDSlength_G4_byHypoxia_leaderless %>% 
        ggplot(aes(x = HalfLife5p_logPhase_log2, y = CDS_length)) + 
        geom_point(size = 4, fill = "#036D9C", 
                   color = "white",
                   shape = 21) +
        geom_smooth(method = lm, se = FALSE, fullrange = TRUE, colour = "#EB636A") +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), 
              axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 20, face = "bold"),
              axis.title =  element_text(size = 25, face = "bold"),
              legend.position = "None") +
        xlim(-1.5, 2.5) +
        ylim(75, 2000)
ggsave("./viz_CDSlengthCorr_byHalfLifeRatio/HalfLifeRatio_logPhaseByHypoxiaG4_leaderless_lengthCorrWith5pEndHalfLife.png", width = 15, height = 12, units = "cm", dpi = 600)
model <- lm(hl_CDS_logPhase_CDSlength_G4_byHypoxia_leaderless$CDS_length ~ hl_CDS_logPhase_CDSlength_G4_byHypoxia_leaderless$HalfLife5p_logPhase_log2)
summary(model)
