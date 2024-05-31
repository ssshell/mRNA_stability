
####### get degradation pattern classes
####### class defined by hierarchical clustering
####### hypoxia

library(tidyverse)
library(RColorBrewer)
library(scales)
library(ggdendro)
library(dendextend)

####### run following script first to load function for tree branch color
####### viz_treeColorScript.R

set.seed(7)

####### get color scheme
paired_hl_hypoxia <- c("#7DC462", "#0D95D0", "#774FA0", "#E72F52")

####### get CDS annotations
gene_CDS <- read.table("../index/msmeg_CombinedAnnotation_CDS.bed", col.names = c("gene", "start", "end", "strand"))
gene_leadered <- read.table("../index/msmeg_CombinedAnnotation_CDS_leadered.bed", col.names = c("gene", "start", "end", "strand"))
gene_leaderless <- read.table("../index/msmeg_CombinedAnnotation_CDS_leaderless.bed", col.names = c("gene", "start", "end", "strand"))

####### get decay profile
raw_mean_log2 <- read.table("../normalization/NormCoverage_ReplicateAvg_log2.txt")
raw_mean_log2_hypoxia <- raw_mean_log2[, c(26:32)]
colnames(raw_mean_log2_hypoxia) <- c("0", "3", "6", "9", "15", "30", "60")

####### get CV of replicates
cv_all <- read.table('../normalization/NormCoverage_CV.txt')
cv_hypoxia <- cv_all[, c(26:32)]
colnames(cv_hypoxia) <- c("0", "3", "6", "9", "15", "30", "60")
cv_hypoxia_cut <- cv_hypoxia %>% 
        rownames_to_column("gene") %>%
        select("gene", "0", "3", "6", "9") %>% 
        filter_at(vars("0", "3", "6", "9"), all_vars(. < 0.75))

####### get normalized coverage after filtering
###### coverage relative to T0
raw_norm2zero <- read.table("../normalization/norm2first_NormCoverage.txt")
raw_norm2zero_hypoxia <- raw_norm2zero[, c(26:32)]
colnames(raw_norm2zero_hypoxia) <- c("0", "3", "6", "9", "15", "30", "60")
raw_norm2zero_hypoxia_RAIcut <- raw_norm2zero_hypoxia[apply(raw_norm2zero_hypoxia, 1, max) <= 1.5, ]

norm2zero_hypoxia <- raw_norm2zero_hypoxia_RAIcut[(rownames(raw_norm2zero_hypoxia_RAIcut) %in% cv_hypoxia_cut$gene), ]
norm2zero_hypoxia_CDS <- norm2zero_hypoxia[(rownames(norm2zero_hypoxia) %in% gene_CDS$gene), ]
norm2zero_hypoxia_CDS_log2 <- raw_mean_log2_hypoxia[(rownames(raw_mean_log2_hypoxia) %in% rownames(norm2zero_hypoxia_CDS)), ]

####### decay pattern clustering.
###### relative to T0
hc_norm2zero <- hclust(dist(norm2zero_hypoxia_CDS, method = "euclidean"), method = "ward.D2")

hc_norm2zero_order <- norm2zero_hypoxia_CDS %>% 
        rownames_to_column("gene") %>%
        mutate(hc_cluster = factor(cutree(hc_norm2zero, 4))) %>% 
        select(gene, hc_cluster)

decayPattern_hypoxia_CDS <- hc_norm2zero_order %>% 
        mutate(Pattern_cls_hypoxia = case_when(hc_cluster == 4 ~ "Fast",
                                               hc_cluster == 1 ~ "Med-fast",
                                               hc_cluster == 2 ~ "Med-slow",
                                               hc_cluster == 3 ~ "Slow"))%>% 
        select(gene, Pattern_cls_hypoxia)
summary(factor(decayPattern_hypoxia_CDS$Pattern_cls_hypoxia))

write.table(decayPattern_hypoxia_CDS, "./pattern/CDS_pattern_hypoxiaClass.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

####### get degradation patterns
###### relative to T0
norm2zero_hypoxia_CDS_log2 %>%
        rownames_to_column("gene") %>%  
        inner_join(decayPattern_hypoxia_CDS, by = "gene") %>% 
        gather(time_points, log2_abundance_T0, -gene, -Pattern_cls_hypoxia) %>%
        mutate_at("time_points", as.numeric) %>%
        group_by(time_points, Pattern_cls_hypoxia) %>% 
        summarise(
                sd = sd(log2_abundance_T0), 
                log2_abundance_T0 = mean(log2_abundance_T0)
        ) %>% 
        ggplot(aes(time_points, log2_abundance_T0)) + 
        scale_colour_manual(values = paired_hl_hypoxia) +
        geom_line(aes(group = Pattern_cls_hypoxia, color = Pattern_cls_hypoxia),
                  position = position_dodge(0.3), size = 5) +
        geom_pointrange(aes(ymin = log2_abundance_T0 - sd, ymax = log2_abundance_T0 + sd, 
                            color = Pattern_cls_hypoxia), position = position_dodge(0.3)) +
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
ggsave("./viz_decayPattern/CDS_hierarchicalClassDecayPattern_complete_hypoxia.png", width = 30, height = 20, units = "cm", dpi = 600)

##### with only complete leadered OR leaderless genes
#### leadered genes
norm2zero_hypoxia_CDS_log2_leadered <- norm2zero_hypoxia_CDS_log2 %>% 
        rownames_to_column("gene") %>% 
        inner_join(gene_leadered, by = "gene") %>% 
        select(-start, -end, -strand) %>% 
        inner_join(decayPattern_hypoxia_CDS, by = "gene")
summary(factor(norm2zero_hypoxia_CDS_log2_leadered$Pattern_cls_hypoxia))

norm2zero_hypoxia_CDS_log2_leadered %>%
        # select(-c("16", "32")) %>%
        gather(time_points, log2_abundance_T1, -gene, -Pattern_cls_hypoxia) %>%
        mutate_at("time_points", as.numeric) %>%
        group_by(time_points, Pattern_cls_hypoxia) %>% 
        summarise(
                sd = sd(log2_abundance_T1), 
                log2_abundance_T1 = mean(log2_abundance_T1)
        ) %>% 
        ggplot(aes(time_points, log2_abundance_T1)) + 
        scale_colour_manual(values = paired_hl_hypoxia) +
        geom_line(aes(group = Pattern_cls_hypoxia, color = Pattern_cls_hypoxia),
                  position = position_dodge(0.3), size = 5) +
        geom_pointrange(aes(ymin = log2_abundance_T1 - sd, ymax = log2_abundance_T1 + sd, 
                            color = Pattern_cls_hypoxia), position = position_dodge(0.3)) +
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
ggsave("./viz_decayPattern/CDS_hierarchicalClassDecayPattern_leadered_complete_hypoxia.png", width = 30, height = 20, units = "cm", dpi = 600)

#### leaderless genes
norm2zero_hypoxia_CDS_log2_leaderless <- norm2zero_hypoxia_CDS_log2 %>% 
        rownames_to_column("gene") %>% 
        inner_join(gene_leaderless, by = "gene") %>% 
        select(-start, -end, -strand) %>% 
        inner_join(decayPattern_hypoxia_CDS, by = "gene")
summary(factor(norm2zero_hypoxia_CDS_log2_leaderless$Pattern_cls_hypoxia))

norm2zero_hypoxia_CDS_log2_leaderless %>%
        # select(-c("16", "32")) %>%
        gather(time_points, log2_abundance_T1, -gene, -Pattern_cls_hypoxia) %>%
        mutate_at("time_points", as.numeric) %>%
        group_by(time_points, Pattern_cls_hypoxia) %>% 
        summarise(
                sd = sd(log2_abundance_T1), 
                log2_abundance_T1 = mean(log2_abundance_T1)
        ) %>% 
        ggplot(aes(time_points, log2_abundance_T1)) + 
        scale_colour_manual(values = paired_hl_hypoxia) +
        geom_line(aes(group = Pattern_cls_hypoxia, color = Pattern_cls_hypoxia),
                  position = position_dodge(0.3), size = 5) +
        geom_pointrange(aes(ymin = log2_abundance_T1 - sd, ymax = log2_abundance_T1 + sd, 
                            color = Pattern_cls_hypoxia), position = position_dodge(0.3)) +
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
ggsave("./viz_decayPattern/CDS_hierarchicalClassDecayPattern_leaderless_complete_hypoxia.png", width = 30, height = 20, units = "cm", dpi = 600)

####### get dendrogram of first round clustering
###### relative to T0
hcdata_norm2zero <- dendro_data_k(hc_norm2zero, 4)
paired_norm2zero <- c("#999999", "#0D95D0", "#774FA0", "#E72F52", "#7DC462")
plot_ggdendro(hcdata_norm2zero,
              direction   = "tb",
              scale.color = paired_norm2zero,
              label.size  = 0,
              branch.size = 0.5,
              expand.y    = 0) +
        theme_void()
ggsave("./viz_decayPattern/CDS_hierarchicalClassDendrogram_complete_hypoxia.png", width = 30, height = 20, units = "cm", dpi = 600, bg = "white")
