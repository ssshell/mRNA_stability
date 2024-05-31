
####### get degradation pattern classes
####### class defined by hierarchical clustering
####### log phase

library(tidyverse)
library(RColorBrewer)
library(scales)
library(ggdendro)
library(dendextend)

####### run following script first to load function for tree branch color
####### viz_treeColorScript.R

set.seed(7)

####### get color scheme
paired_hl_logPhase_firstRound <- c("#7DC462", "#0D95D0", "#774FA0", "#E72F52")

####### get CDS annotations
gene_CDS <- read.table("../index/msmeg_CombinedAnnotation_CDS.bed", col.names = c("gene", "start", "end", "strand"))
gene_leadered <- read.table("../index/msmeg_CombinedAnnotation_CDS_leadered.bed", col.names = c("gene", "start", "end", "strand"))
gene_leaderless <- read.table("../index/msmeg_CombinedAnnotation_CDS_leaderless.bed", col.names = c("gene", "start", "end", "strand"))

####### get decay profile
raw_mean_log2 <- read.table("../normalization/NormCoverage_ReplicateAvg_log2.txt")
raw_mean_log2_logPhase <- raw_mean_log2[, c(7:13)]
colnames(raw_mean_log2_logPhase) <- c("0", "1", "2", "4", "8", "16", "32")

####### get CV of replicates
cv_all <- read.table('../normalization/NormCoverage_CV.txt')
cv_logPhase <- cv_all[, c(7:13)]
colnames(cv_logPhase) <- c("0", "1", "2", "4", "8", "16", "32")
cv_logPhase_cut <- cv_logPhase %>% 
        rownames_to_column("gene") %>%
        select("gene", "0", "1", "2", "4") %>% 
        filter_at(vars("0", "1", "2", "4"), all_vars(. < 0.75))

####### get normalized coverage after filtering
###### coverage relative to T0
raw_norm2zero <- read.table("../normalization/norm2first_NormCoverage.txt")
raw_norm2zero_logPhase <- raw_norm2zero[, c(7:13)]
colnames(raw_norm2zero_logPhase) <- c("0", "1", "2", "4", "8", "16", "32")
raw_norm2zero_logPhase_RAIcut <- raw_norm2zero_logPhase[apply(raw_norm2zero_logPhase, 1, max) <= 1.5, ]

norm2zero_logPhase <- raw_norm2zero_logPhase_RAIcut[(rownames(raw_norm2zero_logPhase_RAIcut) %in% cv_logPhase_cut$gene), ]
norm2zero_logPhase_CDS <- norm2zero_logPhase[(rownames(norm2zero_logPhase) %in% gene_CDS$gene), ]

###### coverage relative to T1
raw_norm2second <- read.table("../normalization/norm2second_NormCoverage.txt")
raw_norm2second_logPhase <- raw_norm2second[, c(7:13)]
colnames(raw_norm2second_logPhase) <- c("0", "1", "2", "4", "8", "16", "32")
raw_norm2second_logPhase_RAIcut <- raw_norm2second_logPhase[apply(raw_norm2second_logPhase, 1, max) <= 1.5, ]

norm2second_logPhase <- raw_norm2second_logPhase_RAIcut[(rownames(raw_norm2second_logPhase_RAIcut) %in% cv_logPhase_cut$gene), ]
norm2second_logPhase_CDS <- norm2second_logPhase[(rownames(norm2second_logPhase) %in% gene_CDS$gene), ]

###### combine kept genes by both T0 and T1
norm2zero_logPhase_CDS_kept <- norm2zero_logPhase_CDS %>% 
        rownames_to_column("gene")
norm2second_logPhase_CDS_kept <- norm2second_logPhase_CDS %>% 
        rownames_to_column("gene")
kept_by_both <- inner_join(norm2zero_logPhase_CDS_kept, 
                           norm2second_logPhase_CDS_kept, by = "gene")
norm2zero_logPhase_CDS_afterFilter <- norm2zero_logPhase_CDS[(rownames(norm2zero_logPhase_CDS) %in% kept_by_both$gene), ]
norm2second_logPhase_CDS_afterFilter <- norm2second_logPhase_CDS[(rownames(norm2second_logPhase_CDS) %in% kept_by_both$gene), ] %>% 
        select(c("1", "2", "4", "8", "16", "32"))
norm2zero_logPhase_CDS_log2 <- raw_mean_log2_logPhase[(rownames(raw_mean_log2_logPhase) %in% rownames(norm2zero_logPhase_CDS_afterFilter)), ]
norm2second_logPhase_CDS_log2 <- raw_mean_log2_logPhase[(rownames(raw_mean_log2_logPhase) %in% rownames(norm2second_logPhase_CDS_afterFilter)), ]

####### decay pattern first round clustering. To get rid of delay class 
###### relative to T0
hc_firstRound_norm2zero <- hclust(dist(norm2zero_logPhase_CDS_afterFilter, method = "euclidean"), method = "ward.D2")

hc_firstRound_norm2zero_order <- norm2zero_logPhase_CDS_afterFilter %>% 
        rownames_to_column("gene") %>%
        mutate(hc_cluster = factor(cutree(hc_firstRound_norm2zero, 4))) %>% 
        select(gene, hc_cluster)

###### relative to T1
hc_firstRound_norm2second <- hclust(dist(norm2second_logPhase_CDS_afterFilter, method = "euclidean"), method = "ward.D2")

hc_firstRound_norm2second_order <- norm2second_logPhase_CDS_afterFilter %>% 
        rownames_to_column("gene") %>%
        mutate(hc_cluster = factor(cutree(hc_firstRound_norm2second, 4))) %>% 
        select(gene, hc_cluster)

####### get degradation patterns to identify delay class
###### relative to T0
##### class 4 is the delay class
norm2zero_logPhase_CDS_log2 %>%
        # select(-c("16", "32")) %>%
        rownames_to_column("gene") %>%  
        inner_join(hc_firstRound_norm2zero_order, by = "gene") %>% 
        gather(time_points, log2_abundance_T0, -gene, -hc_cluster) %>%
        mutate_at("time_points", as.numeric) %>%
        group_by(time_points, hc_cluster) %>% 
        summarise(
                sd = sd(log2_abundance_T0), 
                log2_abundance_T0 = mean(log2_abundance_T0)
        ) %>% 
        ggplot(aes(time_points, log2_abundance_T0)) + 
        scale_colour_manual(values = paired_hl_logPhase_firstRound) +
        geom_line(aes(group = hc_cluster, color = hc_cluster),
                  position = position_dodge(0.3), size = 5) +
        geom_pointrange(aes(ymin = log2_abundance_T0 - sd, ymax = log2_abundance_T0 + sd, 
                            color = hc_cluster), position = position_dodge(0.3)) +
        ylim(0, 17) +
        xlab("Minutes after adding rifampicin") + 
        ylab("Log2 mRNA abundance") +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold"),
        )

###### relative to T1
##### class 4 is the delay class
norm2second_logPhase_CDS_log2 %>%
        # select(-c("16", "32")) %>%
        rownames_to_column("gene") %>%  
        inner_join(hc_firstRound_norm2second_order, by = "gene") %>% 
        gather(time_points, log2_abundance_T1, -gene, -hc_cluster) %>%
        mutate_at("time_points", as.numeric) %>%
        group_by(time_points, hc_cluster) %>% 
        summarise(
                sd = sd(log2_abundance_T1), 
                log2_abundance_T1 = mean(log2_abundance_T1)
        ) %>% 
        ggplot(aes(time_points, log2_abundance_T1)) + 
        scale_colour_manual(values = paired_hl_logPhase_firstRound) +
        geom_line(aes(group = hc_cluster, color = hc_cluster),
                  position = position_dodge(0.3), size = 5) +
        geom_pointrange(aes(ymin = log2_abundance_T1 - sd, ymax = log2_abundance_T1 + sd, 
                            color = hc_cluster), position = position_dodge(0.3)) +
        ylim(0, 17) +
        xlab("Minutes after adding rifampicin") + 
        ylab("Log2 mRNA abundance") +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold")
        )

####### get dendrogram of first round clustering
###### relative to T0
hcdata_firstRound_norm2zero <- dendro_data_k(hc_firstRound_norm2zero, 4)
paired_firstRound_norm2zero <- c("#999999", "#7DC462", "#774FA0", "#0D95D0", "#E72F52")

plot_ggdendro(hcdata_firstRound_norm2zero,
              direction   = "tb",
              scale.color = paired_firstRound_norm2zero,
              label.size  = 0,
              branch.size = 0.5,
              expand.y    = 0) +
        theme_void()

###### relative to T1
hcdata_firstRound_norm2second <- dendro_data_k(hc_firstRound_norm2second, 4)
paired_firstRound_norm2second <- c("#999999", "#7DC462", "#774FA0", "#0D95D0", "#E72F52")

plot_ggdendro(hcdata_firstRound_norm2second,
              direction   = "tb",
              scale.color = paired_firstRound_norm2second,
              label.size  = 0,
              branch.size = 0.5,
              expand.y    = 0) +
        theme_void()

####### get log2 abundance and normalized coverage after excluding delay class
###### relative to T0
norm2zero_logPhase_CDS_kept_exclDelay <- norm2zero_logPhase_CDS_afterFilter %>% 
        rownames_to_column("gene") %>%
        inner_join(hc_firstRound_norm2zero_order, by = "gene") %>% 
        filter(hc_cluster != 4)

###### relative to T1
norm2second_logPhase_CDS_kept_exclDelay <- norm2second_logPhase_CDS_afterFilter %>%
        rownames_to_column("gene") %>%
        inner_join(hc_firstRound_norm2second_order, by = "gene") %>% 
        filter(hc_cluster != 4)

###### combine kept genes by both T0 and T1
kept_by_both <- inner_join(norm2zero_logPhase_CDS_kept_exclDelay, 
                           norm2second_logPhase_CDS_kept_exclDelay, by = "gene")

####### decay pattern second round clustering
####### To separate pattern driven by degradation rate
####### use profile relative to T1
norm2second_logPhase_CDS_exclDelay <- norm2second_logPhase_CDS_afterFilter[(rownames(norm2second_logPhase_CDS_afterFilter) %in% kept_by_both$gene), ]

hc_secondRound_norm2second <- hclust(dist(norm2second_logPhase_CDS_exclDelay, method = "euclidean"), method = "ward.D2")

hc_secondRound_norm2second_order <- norm2second_logPhase_CDS_exclDelay %>% 
        rownames_to_column("gene") %>%
        mutate(hc_cluster = factor(cutree(hc_secondRound_norm2second, 4))) %>% 
        select(gene, hc_cluster)

decayPattern_logPhase_CDS <- hc_secondRound_norm2second_order %>% 
        mutate(Pattern_cls_logPhase = case_when(hc_cluster == 1 ~ "Fast",
                                                hc_cluster == 3 ~ "Med-fast",
                                                hc_cluster == 4 ~ "Med-slow",
                                                hc_cluster == 2 ~ "Slow")) %>% 
        select(gene, Pattern_cls_logPhase)
summary(factor(decayPattern_logPhase_CDS$Pattern_cls_logPhase))

write.table(decayPattern_logPhase_CDS, "./pattern/CDS_pattern_logPhaseClass.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

####### get degradation patterns after excluding delay class
####### use profile relative to T1
norm2second_logPhase_CDS_log2_exclDelay <- norm2second_logPhase_CDS_log2[(rownames(norm2second_logPhase_CDS_log2) %in% rownames(norm2second_logPhase_CDS_exclDelay)), ]

paired_hl_logPhase_secondRound <- c("#7DC462", "#0D95D0", "#774FA0", "#E72F52")

norm2second_logPhase_CDS_log2_exclDelay %>%
        # select(-c("16", "32")) %>%
        rownames_to_column("gene") %>%  
        inner_join(decayPattern_logPhase_CDS, by = "gene") %>% 
        gather(time_points, log2_abundance_T1, -gene, -Pattern_cls_logPhase) %>%
        mutate_at("time_points", as.numeric) %>%
        group_by(time_points, Pattern_cls_logPhase) %>% 
        summarise(
                sd = sd(log2_abundance_T1), 
                log2_abundance_T1 = mean(log2_abundance_T1)
        ) %>% 
        ggplot(aes(time_points, log2_abundance_T1)) + 
        scale_colour_manual(values = paired_hl_logPhase_secondRound) +
        geom_line(aes(group = Pattern_cls_logPhase, color = Pattern_cls_logPhase),
                  position = position_dodge(0.3), size = 5) +
        geom_pointrange(aes(ymin = log2_abundance_T1 - sd, ymax = log2_abundance_T1 + sd, 
                            color = Pattern_cls_logPhase), position = position_dodge(0.3)) +
        ylim(0, 17) +
        xlab("Minutes after adding rifampicin") + 
        ylab("Log2 mRNA abundance") +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold")
        )
ggsave("./viz_decayPattern/CDS_hierarchicalClassDecayPattern_complete_logPhase.png", width = 30, height = 20, units = "cm", dpi = 600)

##### with only complete leadered OR leaderless genes
#### leadered genes
norm2second_logPhase_CDS_log2_exclDelay_leadered <- norm2second_logPhase_CDS_log2_exclDelay %>% 
        rownames_to_column("gene") %>% 
        inner_join(gene_leadered, by = "gene") %>% 
        select(-start, -end, -strand) %>% 
        inner_join(decayPattern_logPhase_CDS, by = "gene")
summary(factor(norm2second_logPhase_CDS_log2_exclDelay_leadered$Pattern_cls_logPhase))

norm2second_logPhase_CDS_log2_exclDelay_leadered %>%
        # select(-c("16", "32")) %>%
        gather(time_points, log2_abundance_T1, -gene, -Pattern_cls_logPhase) %>%
        mutate_at("time_points", as.numeric) %>%
        group_by(time_points, Pattern_cls_logPhase) %>% 
        summarise(
                sd = sd(log2_abundance_T1), 
                log2_abundance_T1 = mean(log2_abundance_T1)
        ) %>% 
        ggplot(aes(time_points, log2_abundance_T1)) + 
        scale_colour_manual(values = paired_hl_logPhase_secondRound) +
        geom_line(aes(group = Pattern_cls_logPhase, color = Pattern_cls_logPhase),
                  position = position_dodge(0.3), size = 5) +
        geom_pointrange(aes(ymin = log2_abundance_T1 - sd, ymax = log2_abundance_T1 + sd, 
                            color = Pattern_cls_logPhase), position = position_dodge(0.3)) +
        ylim(0, 17) +
        xlab("Minutes after adding rifampicin") + 
        ylab("Log2 mRNA abundance") +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold")
        )
ggsave("./viz_decayPattern/CDS_hierarchicalClassDecayPattern_leadered_complete_logPhase.png", width = 30, height = 20, units = "cm", dpi = 600)

#### leaderless genes
norm2second_logPhase_CDS_log2_exclDelay_leaderless <- norm2second_logPhase_CDS_log2_exclDelay %>% 
        rownames_to_column("gene") %>% 
        inner_join(gene_leaderless, by = "gene") %>% 
        select(-start, -end, -strand) %>% 
        inner_join(decayPattern_logPhase_CDS, by = "gene")
summary(factor(norm2second_logPhase_CDS_log2_exclDelay_leaderless$Pattern_cls_logPhase))

norm2second_logPhase_CDS_log2_exclDelay_leaderless %>%
        # select(-c("16", "32")) %>%
        gather(time_points, log2_abundance_T1, -gene, -Pattern_cls_logPhase) %>%
        mutate_at("time_points", as.numeric) %>%
        group_by(time_points, Pattern_cls_logPhase) %>% 
        summarise(
                sd = sd(log2_abundance_T1), 
                log2_abundance_T1 = mean(log2_abundance_T1)
        ) %>% 
        ggplot(aes(time_points, log2_abundance_T1)) + 
        scale_colour_manual(values = paired_hl_logPhase_secondRound) +
        geom_line(aes(group = Pattern_cls_logPhase, color = Pattern_cls_logPhase),
                  position = position_dodge(0.3), size = 5) +
        geom_pointrange(aes(ymin = log2_abundance_T1 - sd, ymax = log2_abundance_T1 + sd, 
                            color = Pattern_cls_logPhase), position = position_dodge(0.3)) +
        ylim(0, 17) +
        xlab("Minutes after adding rifampicin") + 
        ylab("Log2 mRNA abundance") +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold")
        )
ggsave("./viz_decayPattern/CDS_hierarchicalClassDecayPattern_leaderless_complete_logPhase.png", width = 30, height = 20, units = "cm", dpi = 600)

####### get dendrogram of second round clustering
###### relative to T1
hcdata_secondRound_norm2second <- dendro_data_k(hc_secondRound_norm2second, 4)
paired_secondRound_norm2second <- c("#999999", "#7DC462", "#E72F52", "#0D95D0", "#774FA0")

plot_ggdendro(hcdata_secondRound_norm2second,
              direction   = "tb",
              scale.color = paired_secondRound_norm2second,
              label.size  = 0,
              branch.size = 0.5,
              expand.y    = 0) +
        theme_void()
ggsave("./viz_decayPattern/CDS_hierarchicalClassDendrogram_complete_logPhase.png", width = 30, height = 20, units = "cm", dpi = 600, bg = "white")
