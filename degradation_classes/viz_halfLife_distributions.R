
####### plot distributions of half-life & half-life classes

library(tidyverse)
library(ggalluvial)
library(RColorBrewer)
library(scales)
library(ggbeeswarm)
library(ggbreak)

####### with all half-life available genes in each condition  
####### condition overview log phase vs. hypoxia

# display.brewer.all(type = 'qual')
# display.brewer.pal(n = 12, name = 'Paired')
# brewer.pal(n = 12, name = 'Paired')
# show_col("#1F78B4")

####### get color scheme for half-life classes
col_hl_cls <- c("#7DC462", "#0D95D0", "#774FA0", "#E72F52")
col_FC_hypoxiaToLogPhase <- c("#E72F52", "#774FA0", "#0D95D0", "#7DC462")

####### get gene list
gene_leadered <- read.table("../index/msmeg_CombinedAnnotation_CDS_leadered.bed", col.names = c("gene", "start", "end", "strand"))
gene_leaderless <- read.table("../index/msmeg_CombinedAnnotation_CDS_leaderless.bed", col.names = c("gene", "start", "end", "strand"))

####### get half-life values
hl_CDS_complete_class <- read.table('./half_life/CDS_halfLife_completeClass.txt', header = TRUE)

hl_CDS_logPhase <- hl_CDS_complete_class %>% 
        select(gene, HalfLife_vals_logPhase) %>% 
        drop_na() %>% 
        mutate(condition = "logPhase") %>% 
        rename(HalfLife_vals = HalfLife_vals_logPhase)
summary(hl_CDS_logPhase$HalfLife_vals)

hl_CDS_hypoxia <- hl_CDS_complete_class %>% 
        select(gene, HalfLife_vals_hypoxia) %>% 
        drop_na() %>% 
        mutate(condition = "hypoxia") %>% 
        rename(HalfLife_vals = HalfLife_vals_hypoxia)
summary(hl_CDS_hypoxia$HalfLife_vals)

####### distributions for log phase and hypoxia 
###### with complete available genes in each condition
hl_complete <- bind_rows(hl_CDS_logPhase, hl_CDS_hypoxia)

##### using log2 half-life values
hl_complete %>% 
        ggplot(aes(x = log2(HalfLife_vals), fill = condition, color = condition)) +
        scale_color_manual(values = c("#83B7DF", "#036D9C")) +
        scale_fill_manual(values = c("#83B7DF", "#036D9C")) +
        geom_histogram(binwidth = 0.1, position = 'identity') +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold")
                ) +
        xlab("Log2(mRNA Half-life, mins)") + 
        ylab("Gene count")
ggsave("./viz_halfLife/CDS_halfLife_complete.png", width = 30, height = 20, units = "cm", dpi = 600)

##### using raw half-life values
hl_complete %>% 
        ggplot(aes(x = HalfLife_vals, fill = condition, color = condition)) +
        scale_color_manual(values = c("#83B7DF", "#036D9C")) +
        scale_fill_manual(values = c("#83B7DF", "#036D9C")) +
        geom_histogram(binwidth = 0.1, position = 'identity') +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold")
        ) +
        xlim(0, 70) + 
        xlab("mRNA Half-life, mins") + 
        ylab("Gene count")
ggsave("./viz_halfLife/CDS_halfLife_complete_rawHalfLifeValues.png", width = 30, height = 20, units = "cm", dpi = 600)

##### using raw half-life values with y axis limit
hl_complete %>% 
        ggplot(aes(x = HalfLife_vals, fill = condition, color = condition)) +
        scale_color_manual(values = c("#83B7DF", "#036D9C")) +
        scale_fill_manual(values = c("#83B7DF", "#036D9C")) +
        geom_histogram(binwidth = 0.1, position = 'identity') +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold")
        ) +
        xlim(0, 70) + 
        ylim(0, 70) +
        xlab("mRNA Half-life, mins") + 
        ylab("Gene count")
ggsave("./viz_halfLife/CDS_halfLife_complete_rawHalfLifeValues_withYlimit.png", width = 30, height = 20, units = "cm", dpi = 600)

##### using raw half-life values with break on y axis
hl_complete %>% 
        ggplot(aes(x = HalfLife_vals, fill = condition, color = condition)) +
        scale_color_manual(values = c("#83B7DF", "#036D9C")) +
        scale_fill_manual(values = c("#83B7DF", "#036D9C")) +
        geom_histogram(binwidth = 0.1, position = 'identity') +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold")
        ) +
        scale_y_break(breaks = c(70, 200), scales = 0.5) +
        xlim(0, 70) + 
        xlab("mRNA Half-life, mins") + 
        ylab("Gene count")
ggsave("./viz_halfLife/CDS_halfLife_complete_rawHalfLifeValues_withYbreaks.png", width = 30, height = 20, units = "cm", dpi = 600)

###### with only complete leadered OR leaderless genes in each condition
##### leadered genes
hl_CDS_logPhase_leadered <- hl_CDS_logPhase %>% 
        inner_join(gene_leadered, by = "gene")
summary(hl_CDS_logPhase_leadered$HalfLife_vals)
hl_CDS_hypoxia_leadered <- hl_CDS_hypoxia %>% 
        inner_join(gene_leadered, by = "gene")
summary(hl_CDS_hypoxia_leadered$HalfLife_vals)

hl_leadered <- bind_rows(hl_CDS_logPhase_leadered, hl_CDS_hypoxia_leadered)

#### using log2 half-life values
hl_leadered %>% 
        ggplot(aes(x = log2(HalfLife_vals), fill = condition, color = condition)) +
        scale_color_manual(values = c("#83B7DF", "#036D9C")) +
        scale_fill_manual(values = c("#83B7DF", "#036D9C")) +
        geom_histogram(binwidth = 0.1, position = 'identity') +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold")
        ) +
        xlab("Log2(mRNA Half-life, mins)") + 
        ylab("Leadered gene count")
ggsave("./viz_halfLife/CDS_halfLife_leadered_complete.png", width = 30, height = 20, units = "cm", dpi = 600)

#### using raw half-life values
hl_leadered %>% 
        ggplot(aes(x = HalfLife_vals, fill = condition, color = condition)) +
        scale_color_manual(values = c("#83B7DF", "#036D9C")) +
        scale_fill_manual(values = c("#83B7DF", "#036D9C")) +
        geom_histogram(binwidth = 0.1, position = 'identity') +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold")
        ) +
        xlim(0, 70) +
        xlab("mRNA Half-life, mins") + 
        ylab("Leadered gene count")
ggsave("./viz_halfLife/CDS_halfLife_leadered_complete_rawHalfLifeValues.png", width = 30, height = 20, units = "cm", dpi = 600)

##### leaderless genes
hl_CDS_logPhase_leaderless <- hl_CDS_logPhase %>% 
        inner_join(gene_leaderless, by = "gene")
summary(hl_CDS_logPhase_leaderless$HalfLife_vals)
hl_CDS_hypoxia_leaderless <- hl_CDS_hypoxia %>% 
        inner_join(gene_leaderless, by = "gene")
summary(hl_CDS_hypoxia_leaderless$HalfLife_vals)

hl_leaderless <- bind_rows(hl_CDS_logPhase_leaderless, hl_CDS_hypoxia_leaderless)

#### using log2 half-life values
hl_leaderless %>% 
        ggplot(aes(x = log2(HalfLife_vals), fill = condition, color = condition)) +
        scale_color_manual(values = c("#83B7DF", "#036D9C")) +
        scale_fill_manual(values = c("#83B7DF", "#036D9C")) +
        geom_histogram(binwidth = 0.1, position = 'identity') +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold")
        ) +
        xlab("Log2(mRNA Half-life, mins)") + 
        ylab("Leaderless gene count")
ggsave("./viz_halfLife/CDS_halfLife_leaderless_complete.png", width = 30, height = 20, units = "cm", dpi = 600)

#### using raw half-life values
hl_leaderless %>% 
        ggplot(aes(x = HalfLife_vals, fill = condition, color = condition)) +
        scale_color_manual(values = c("#83B7DF", "#036D9C")) +
        scale_fill_manual(values = c("#83B7DF", "#036D9C")) +
        geom_histogram(binwidth = 0.1, position = 'identity') +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold")
        ) +
        xlim(0, 70) +
        ylim(0, 125) +
        xlab("mRNA Half-life, mins") + 
        ylab("Leaderless gene count")
ggsave("./viz_halfLife/CDS_halfLife_leaderless_complete_rawHalfLifeValues.png", width = 30, height = 20, units = "cm", dpi = 600)

###### directly comparing half-life between leadered and leaderless in two conditions  
set.seed(7)

##### log phase
hl_CDS_logPhase_leadered_status <- hl_CDS_logPhase_leadered %>% 
        mutate(leader_status = "Leadered")
hl_CDS_logPhase_leaderless_status <- hl_CDS_logPhase_leaderless %>% 
        mutate(leader_status = "Leaderless")

hl_CDS_logPhase_leadered_status %>% 
        bind_rows(hl_CDS_logPhase_leaderless_status) %>% 
        mutate(leader_status = factor(leader_status)) %>% 
        ggplot(aes(x = HalfLife_vals, y = leader_status, fill = leader_status)) +
        geom_violin(color = NA) +
        scale_fill_manual(values = c("#EB636A", "#F1C36B")) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold")
        ) +
        xlab("mRNA Half-life, mins") + 
        ylab("Leader status")
ggsave("./viz_halfLife/CDS_halfLife_distViolin_leaderedVSleaderless_logPhase.png", width = 35, height = 20, units = "cm", dpi = 600)

hl_CDS_logPhase_median <- hl_CDS_logPhase_leadered_status %>% 
        bind_rows(hl_CDS_logPhase_leaderless_status) %>% 
        mutate(leader_status = factor(leader_status)) %>% 
        group_by(leader_status) %>% 
        summarise(median_val = median(HalfLife_vals))

hl_CDS_logPhase_leadered_status %>% 
        bind_rows(hl_CDS_logPhase_leaderless_status) %>% 
        mutate(leader_status = factor(leader_status)) %>%
        ggplot(aes(leader_status, HalfLife_vals, fill = leader_status)) +
        geom_quasirandom(shape = 21, size = 2.5, color = "white",
                         width = 0.2) +
        geom_crossbar(data = hl_CDS_logPhase_median, aes(x = leader_status, y = median_val, ymin = median_val, ymax = median_val),
                      width = 0.5, colour = "#424242", size = 0.5) +
        scale_fill_manual(values = c("#EB636A", "#F1C36B")) +
        ylim(0, 5) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
                axis.text.y = element_text(size = 10),
                axis.title.x = element_blank(),
        )
ggsave("./viz_halfLife/CDS_halfLife_distPoint_leaderedVSleaderless_logPhase.png", width = 10, height = 7, units = "cm", dpi = 600)
wilcox.test(hl_CDS_logPhase_leadered$HalfLife_vals, hl_CDS_logPhase_leaderless$HalfLife_vals)

hl_CDS_logPhase_leadered_status %>% 
        bind_rows(hl_CDS_logPhase_leaderless_status) %>% 
        mutate(leader_status = factor(leader_status)) %>% 
        ggplot(aes(x = HalfLife_vals)) +
        stat_ecdf(aes(colour = leader_status), size = 2) +
        scale_color_manual(values = c("#EB636A", "#F1C36B")) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 27, face = "bold"),
                axis.title =  element_text(size = 33, face = "bold")
        ) +
        xlab("mRNA Half-life, mins") +
        ylab("Cumulative fraction of transcripts")
ggsave("./viz_halfLife/CDS_halfLife_distCDF_leaderedVSleaderless_logPhase.png", width = 30, height = 22, units = "cm", dpi = 600)
ks.test(hl_CDS_logPhase_leadered$HalfLife_vals, hl_CDS_logPhase_leaderless$HalfLife_vals)

##### hypoxia
hl_CDS_hypoxia_leadered_status <- hl_CDS_hypoxia_leadered %>% 
        mutate(leader_status = "Leadered")
hl_CDS_hypoxia_leaderless_status <- hl_CDS_hypoxia_leaderless %>% 
        mutate(leader_status = "Leaderless")

hl_CDS_hypoxia_leadered_status %>% 
        bind_rows(hl_CDS_hypoxia_leaderless_status) %>% 
        mutate(leader_status = factor(leader_status)) %>% 
        ggplot(aes(x = HalfLife_vals, y = leader_status, fill = leader_status)) +
        geom_violin(color = NA) +
        scale_fill_manual(values = c("#EB636A", "#F1C36B")) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold")
        ) +
        xlab("mRNA Half-life, mins") + 
        ylab("Leader status") +
        xlim(0, 150)
ggsave("./viz_halfLife/CDS_halfLife_distViolin_leaderedVSleaderless_hypoxia.png", width = 35, height = 20, units = "cm", dpi = 600)

hl_CDS_hypoxia_median <- hl_CDS_hypoxia_leadered_status %>% 
        bind_rows(hl_CDS_hypoxia_leaderless_status) %>% 
        mutate(leader_status = factor(leader_status)) %>% 
        group_by(leader_status) %>% 
        summarise(median_val = median(HalfLife_vals))

hl_CDS_hypoxia_leadered_status %>% 
        bind_rows(hl_CDS_hypoxia_leaderless_status) %>% 
        mutate(leader_status = factor(leader_status)) %>%
        ggplot(aes(leader_status, HalfLife_vals, fill = leader_status)) +
        geom_quasirandom(shape = 21, size = 2.5, color = "white",
                         width = 0.2) +
        geom_crossbar(data = hl_CDS_hypoxia_median, aes(x = leader_status, y = median_val, ymin = median_val, ymax = median_val),
                      width = 0.5, colour = "#424242", size = 0.5) +
        scale_fill_manual(values = c("#EB636A", "#F1C36B")) +
        ylim(0, 100) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
                axis.text.y = element_text(size = 10),
                axis.title.x = element_blank(),
        )
ggsave("./viz_halfLife/CDS_halfLife_distPoint_leaderedVSleaderless_hypoxia.png", width = 10, height = 7, units = "cm", dpi = 600)
wilcox.test(hl_CDS_hypoxia_leadered$HalfLife_vals, hl_CDS_hypoxia_leaderless$HalfLife_vals)

hl_CDS_hypoxia_leadered_status %>% 
        bind_rows(hl_CDS_hypoxia_leaderless_status) %>% 
        mutate(leader_status = factor(leader_status)) %>% 
        ggplot(aes(x = HalfLife_vals)) +
        stat_ecdf(aes(colour = leader_status), size = 2) +
        scale_color_manual(values = c("#EB636A", "#F1C36B")) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 27, face = "bold"),
                axis.title =  element_text(size = 33, face = "bold")
        ) +
        xlab("mRNA Half-life, mins") +
        ylab("Cumulative fraction of transcripts") +
        xlim(0, 150) 
ggsave("./viz_halfLife/CDS_halfLife_distCDF_leaderedVSleaderless_hypoxia.png", width = 30, height = 22, units = "cm", dpi = 600)
ks.test(hl_CDS_hypoxia_leadered$HalfLife_vals, hl_CDS_hypoxia_leaderless$HalfLife_vals)


####### half-life classes for log phase and hypoxia
###### log phase
table(hl_CDS_complete_class$HalfLife_cls_logPhase)
hl_CDS_complete_class %>% 
        select(gene, HalfLife_vals_logPhase, HalfLife_cls_logPhase) %>% 
        drop_na() %>% 
        ggplot(aes(x = HalfLife_vals_logPhase, fill = HalfLife_cls_logPhase)) +
        geom_histogram(binwidth = 0.02, position = 'identity') +
        scale_fill_manual(values = col_hl_cls) +
        xlim(0, 5) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold")
        ) +
        xlab("mRNA Half-life in log phase, mins") + 
        ylab("Gene count")
ggsave("./viz_halfLife/CDS_halfLifeClass_logPhase.png", width = 30, height = 20, units = "cm", dpi = 600)

###### hypoxia
table(hl_CDS_complete_class$HalfLife_cls_hypoxia)
hl_CDS_complete_class %>% 
        select(gene, HalfLife_vals_hypoxia, HalfLife_cls_hypoxia) %>% 
        drop_na() %>% 
        ggplot(aes(x = HalfLife_vals_hypoxia, fill = HalfLife_cls_hypoxia)) +
        geom_histogram(binwidth = 0.3, position = 'identity') +
        scale_fill_manual(values = col_hl_cls) +
        xlim(0, 70) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold")
        ) +
        xlab("mRNA Half-life in hypoxia, mins") + 
        ylab("Gene count")
ggsave("./viz_halfLife/CDS_halfLifeClass_hypoxia.png", width = 30, height = 20, units = "cm", dpi = 600)

####### half-life fold change classes
summary(hl_CDS_complete_class$HalfLife_FC_hypoxiaToLogPhase)
hlFC_CDS_hypoxiaToLogPhase <- hl_CDS_complete_class %>% 
        select(gene, HalfLife_FC_hypoxiaToLogPhase, HalfLife_FCcls_hypoxiaToLogPhase) %>% 
        drop_na()
summary(hlFC_CDS_hypoxiaToLogPhase$HalfLife_FC_hypoxiaToLogPhase)
table(hlFC_CDS_hypoxiaToLogPhase$HalfLife_FCcls_hypoxiaToLogPhase)
hl_CDS_complete_class %>% 
        drop_na() %>% 
        select(gene, HalfLife_FC_hypoxiaToLogPhase, HalfLife_FCcls_hypoxiaToLogPhase) %>% 
        ggplot(aes(x = HalfLife_FC_hypoxiaToLogPhase, fill = HalfLife_FCcls_hypoxiaToLogPhase)) +
        geom_histogram(binwidth = 0.2, position = 'identity') +
        scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
        xlim(0, 70) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold")
        ) +
        xlab("mRNA Half-life fold change") + 
        ylab("Gene count")
ggsave("./viz_halfLife/CDS_halfLifeFcClass_hypoxiaToLogPhase.png", width = 30, height = 20, units = "cm", dpi = 600)

####### half-life change from log phase to hypoxia 
###### with complete genes have half-life available in both conditions
hl_CDS_FC <- hl_CDS_complete_class %>% 
        drop_na() %>% 
        select(gene, HalfLife_cls_logPhase, HalfLife_cls_hypoxia)
table(hl_CDS_FC$HalfLife_cls_logPhase)
table(hl_CDS_FC$HalfLife_cls_hypoxia)

hl_CDS_FC_alluvial <- hl_CDS_FC %>% 
        mutate(across(c(HalfLife_cls_logPhase, HalfLife_cls_hypoxia), factor)) %>% 
        group_by(HalfLife_cls_logPhase, HalfLife_cls_hypoxia) %>% 
        summarise(freq = n())

############## test session
####### check dataframe format
is_alluvia_form(hl_CDS_FC_alluvial, axes = 1:3, silent = TRUE)

####### double check frequency count
sum(hl_CDS_FC_alluvial$freq)
test <- hl_CDS_FC %>% 
        # filter(HalfLife_cls_logPhase == 'Fast' & HalfLife_cls_hypoxia == 'Fast')
        # filter(HalfLife_cls_logPhase == 'Fast' & HalfLife_cls_hypoxia == 'Med-fast')
        # filter(HalfLife_cls_logPhase == 'Fast' & HalfLife_cls_hypoxia == 'Med-slow')
        # filter(HalfLife_cls_logPhase == 'Fast' & HalfLife_cls_hypoxia == 'Slow')
        # filter(HalfLife_cls_logPhase == 'Med-fast' & HalfLife_cls_hypoxia == 'Fast')
        # filter(HalfLife_cls_logPhase == 'Med-fast' & HalfLife_cls_hypoxia == 'Med-fast')
        # filter(HalfLife_cls_logPhase == 'Med-fast' & HalfLife_cls_hypoxia == 'Med-slow')
        # filter(HalfLife_cls_logPhase == 'Med-fast' & HalfLife_cls_hypoxia == 'Slow')
        # filter(HalfLife_cls_logPhase == 'Med-slow' & HalfLife_cls_hypoxia == 'Fast')
        # filter(HalfLife_cls_logPhase == 'Med-slow' & HalfLife_cls_hypoxia == 'Med-fast')
        # filter(HalfLife_cls_logPhase == 'Med-slow' & HalfLife_cls_hypoxia == 'Med-slow')
        # filter(HalfLife_cls_logPhase == 'Med-slow' & HalfLife_cls_hypoxia == 'Slow')
        # filter(HalfLife_cls_logPhase == 'Slow' & HalfLife_cls_hypoxia == 'Fast')
        # filter(HalfLife_cls_logPhase == 'Slow' & HalfLife_cls_hypoxia == 'Med-fast')
        # filter(HalfLife_cls_logPhase == 'Slow' & HalfLife_cls_hypoxia == 'Med-slow')
        filter(HalfLife_cls_logPhase == 'Slow' & HalfLife_cls_hypoxia == 'Slow')
############## end of test session

##### to eliminate warning message b/c log phase and hypoxia have same class value
hl_CDS_FC_alluvial$HalfLife_cls_logPhase <- paste(hl_CDS_FC_alluvial$HalfLife_cls_logPhase, "logPhase", sep = "_")
hl_CDS_FC_alluvial$HalfLife_cls_hypoxia <- paste(hl_CDS_FC_alluvial$HalfLife_cls_hypoxia, "hypoxia", sep = "_")

hl_CDS_FC_alluvial %>% 
        ggplot(aes(y = freq, axis1 = HalfLife_cls_logPhase, axis2 = HalfLife_cls_hypoxia)) +
        geom_alluvium(aes(fill = HalfLife_cls_logPhase), width = 1/12) +
        geom_stratum(width = 1/12) +
        # geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
        scale_x_discrete(expand = c(.07, .07)) +
        scale_fill_manual(values = col_hl_cls) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank()
              )
ggsave("./viz_halfLife/CDS_halfLifeFC_complete.png", width = 30, height = 20, units = "cm", dpi = 600)

###### with only complete leadered OR leaderless genes in each condition
##### leadered genes
hl_leadered_FC <- hl_CDS_complete_class %>% 
        drop_na() %>% 
        select(gene, HalfLife_cls_logPhase, HalfLife_cls_hypoxia) %>% 
        inner_join(gene_leadered, by = "gene")
table(hl_leadered_FC$HalfLife_cls_logPhase)
table(hl_leadered_FC$HalfLife_cls_hypoxia)

hl_leadered_FC_alluvial <- hl_leadered_FC %>% 
        mutate(across(c(HalfLife_cls_logPhase, HalfLife_cls_hypoxia), factor)) %>% 
        group_by(HalfLife_cls_logPhase, HalfLife_cls_hypoxia) %>% 
        summarise(freq = n())

##### check dataframe format
is_alluvia_form(hl_leadered_FC_alluvial, axes = 1:3, silent = TRUE)

#### to eliminate warning message b/c log phase and hypoxia have same class value
hl_leadered_FC_alluvial$HalfLife_cls_logPhase <- paste(hl_leadered_FC_alluvial$HalfLife_cls_logPhase, "logPhase", sep = "_")
hl_leadered_FC_alluvial$HalfLife_cls_hypoxia <- paste(hl_leadered_FC_alluvial$HalfLife_cls_hypoxia, "hypoxia", sep = "_")

hl_leadered_FC_alluvial %>% 
        ggplot(aes(y = freq, axis1 = HalfLife_cls_logPhase, axis2 = HalfLife_cls_hypoxia)) +
        geom_alluvium(aes(fill = HalfLife_cls_logPhase), width = 1/12) +
        geom_stratum(width = 1/12) +
        # geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
        scale_x_discrete(expand = c(.07, .07)) +
        scale_fill_manual(values = col_hl_cls) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank()
              )
ggsave("./viz_halfLife/CDS_halfLifeFC_leadered_complete.png", width = 30, height = 20, units = "cm", dpi = 600)

##### leaderless genes
hl_leaderless_FC <- hl_CDS_complete_class %>% 
        drop_na() %>% 
        select(gene, HalfLife_cls_logPhase, HalfLife_cls_hypoxia) %>% 
        inner_join(gene_leaderless, by = "gene")
table(hl_leaderless_FC$HalfLife_cls_logPhase)
table(hl_leaderless_FC$HalfLife_cls_hypoxia)

hl_leaderless_FC_alluvial <- hl_leaderless_FC %>% 
        mutate(across(c(HalfLife_cls_logPhase, HalfLife_cls_hypoxia), factor)) %>% 
        group_by(HalfLife_cls_logPhase, HalfLife_cls_hypoxia) %>% 
        summarise(freq = n())

##### check dataframe format
is_alluvia_form(hl_leaderless_FC_alluvial, axes = 1:3, silent = TRUE)

#### to eliminate warning message b/c log phase and hypoxia have same class value
hl_leaderless_FC_alluvial$HalfLife_cls_logPhase <- paste(hl_leaderless_FC_alluvial$HalfLife_cls_logPhase, "logPhase", sep = "_")
hl_leaderless_FC_alluvial$HalfLife_cls_hypoxia <- paste(hl_leaderless_FC_alluvial$HalfLife_cls_hypoxia, "hypoxia", sep = "_")

hl_leaderless_FC_alluvial %>% 
        ggplot(aes(y = freq, axis1 = HalfLife_cls_logPhase, axis2 = HalfLife_cls_hypoxia)) +
        geom_alluvium(aes(fill = HalfLife_cls_logPhase), width = 1/12) +
        geom_stratum(width = 1/12) +
        # geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
        scale_x_discrete(expand = c(.07, .07)) +
        scale_fill_manual(values = col_hl_cls) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank()
              )
ggsave("./viz_halfLife/CDS_halfLifeFC_leaderless_complete.png", width = 30, height = 20, units = "cm", dpi = 600)

####### degradation pattern according to half-life class definition
####### for log phase and hypoxia with complete available genes in each condition
####### and for complete available leadered and leaderless genes
###### get decay profile
decayProfile_mean_log2 <- read.table("../normalization/NormCoverage_ReplicateAvg_log2.txt")

###### get decay pattern by half-life class
##### log phase
decayProfile_logPhase <- decayProfile_mean_log2[, c(7:13)]
colnames(decayProfile_logPhase) <- c("0", "1", "2", "4", "8", "16", "32")

#### get CDS half-life class
hl_CDS_logPhase <- hl_CDS_complete_class %>% 
        select(gene, HalfLife_cls_logPhase, HalfLife_vals_logPhase) %>% 
        drop_na()

### with complete available genes
## combine decay profile and half-life class
decayProfile_logPhase_CDS <- decayProfile_logPhase %>% 
        rownames_to_column("gene") %>% 
        inner_join(hl_CDS_logPhase, by = "gene")
table(decayProfile_logPhase_CDS$HalfLife_cls_logPhase)
summary(decayProfile_logPhase_CDS$HalfLife_vals_logPhase)

decayProfile_logPhase_CDS %>%
        select(-HalfLife_vals_logPhase) %>% 
        gather(time_points, log2_abundance, -gene, -HalfLife_cls_logPhase) %>%
        mutate_at("time_points", as.numeric) %>%
        group_by(time_points, HalfLife_cls_logPhase) %>% 
        summarise(
                sd = sd(log2_abundance), 
                log2_abundance = mean(log2_abundance)
        ) %>% 
        ggplot(aes(time_points, log2_abundance)) + 
        scale_colour_manual(values = col_hl_cls) +
        geom_line(aes(group = HalfLife_cls_logPhase, color = HalfLife_cls_logPhase),
                  position = position_dodge(0.3), size = 5) +
        geom_pointrange(aes(ymin = log2_abundance - sd, ymax = log2_abundance + sd, 
                            color = HalfLife_cls_logPhase), position = position_dodge(0.3)) +
        ylim(0, 17) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold")
        ) +
        xlab("Minutes after adding rifampicin") + 
        ylab("Log2 mRNA abundance")
ggsave("./viz_halfLife/CDS_halfLifeClassDecayPattern_complete_logPhase.png", width = 30, height = 20, units = "cm", dpi = 600)

### with only complete leadered OR leaderless genes in each condition
## leadered genes
# combine decay profile and half-life class
decayProfile_logPhase_CDS_leadered <- decayProfile_logPhase %>% 
        rownames_to_column("gene") %>% 
        inner_join(hl_CDS_logPhase, by = "gene") %>% 
        inner_join(gene_leadered, by = "gene") %>% 
        select(-start, -end, -strand)
table(decayProfile_logPhase_CDS_leadered$HalfLife_cls_logPhase)
summary(decayProfile_logPhase_CDS_leadered$HalfLife_vals_logPhase)

decayProfile_logPhase_CDS_leadered %>%
        select(-HalfLife_vals_logPhase) %>% 
        gather(time_points, log2_abundance, -gene, -HalfLife_cls_logPhase) %>%
        mutate_at("time_points", as.numeric) %>%
        group_by(time_points, HalfLife_cls_logPhase) %>% 
        summarise(
                sd = sd(log2_abundance), 
                log2_abundance = mean(log2_abundance)
        ) %>% 
        ggplot(aes(time_points, log2_abundance)) + 
        scale_colour_manual(values = col_hl_cls) +
        geom_line(aes(group = HalfLife_cls_logPhase, color = HalfLife_cls_logPhase),
                  position = position_dodge(0.3), size = 5) +
        geom_pointrange(aes(ymin = log2_abundance - sd, ymax = log2_abundance + sd, 
                            color = HalfLife_cls_logPhase), position = position_dodge(0.3)) +
        ylim(0, 17) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold")
        ) +
        xlab("Minutes after adding rifampicin") + 
        ylab("Log2 mRNA abundance")
ggsave("./viz_halfLife/CDS_halfLifeClassDecayPattern_leadered_complete_logPhase.png", width = 30, height = 20, units = "cm", dpi = 600)

## leaderless genes
# combine decay profile and half-life class
decayProfile_logPhase_CDS_leaderless <- decayProfile_logPhase %>% 
        rownames_to_column("gene") %>% 
        inner_join(hl_CDS_logPhase, by = "gene") %>% 
        inner_join(gene_leaderless, by = "gene") %>% 
        select(-start, -end, -strand)
table(decayProfile_logPhase_CDS_leaderless$HalfLife_cls_logPhase)
summary(decayProfile_logPhase_CDS_leaderless$HalfLife_vals_logPhase)

decayProfile_logPhase_CDS_leaderless %>%
        select(-HalfLife_vals_logPhase) %>% 
        gather(time_points, log2_abundance, -gene, -HalfLife_cls_logPhase) %>%
        mutate_at("time_points", as.numeric) %>%
        group_by(time_points, HalfLife_cls_logPhase) %>% 
        summarise(
                sd = sd(log2_abundance), 
                log2_abundance = mean(log2_abundance)
        ) %>% 
        ggplot(aes(time_points, log2_abundance)) + 
        scale_colour_manual(values = col_hl_cls) +
        geom_line(aes(group = HalfLife_cls_logPhase, color = HalfLife_cls_logPhase),
                  position = position_dodge(0.3), size = 5) +
        geom_pointrange(aes(ymin = log2_abundance - sd, ymax = log2_abundance + sd, 
                            color = HalfLife_cls_logPhase), position = position_dodge(0.3)) +
        ylim(0, 17) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold")
        ) +
        xlab("Minutes after adding rifampicin") + 
        ylab("Log2 mRNA abundance")
ggsave("./viz_halfLife/CDS_halfLifeClassDecayPattern_leaderless_complete_logPhase.png", width = 30, height = 20, units = "cm", dpi = 600)

##### hypoxia
decayProfile_hypoxia <- decayProfile_mean_log2[, c(26:32)]
colnames(decayProfile_hypoxia) <- c("0", "3", "6", "9", "15", "30", "60")

#### get CDS half-life class
hl_CDS_hypoxia <- hl_CDS_complete_class %>% 
        select(gene, HalfLife_cls_hypoxia, HalfLife_vals_hypoxia) %>% 
        drop_na()

### with complete available genes
## combine decay profile and half-life class
decayProfile_hypoxia_CDS <- decayProfile_hypoxia %>% 
        rownames_to_column("gene") %>% 
        inner_join(hl_CDS_hypoxia, by = "gene")
table(decayProfile_hypoxia_CDS$HalfLife_cls_hypoxia)
summary(decayProfile_hypoxia_CDS$HalfLife_vals_hypoxia)

decayProfile_hypoxia_CDS %>%
        select(-HalfLife_vals_hypoxia) %>% 
        gather(time_points, log2_abundance, -gene, -HalfLife_cls_hypoxia) %>%
        mutate_at("time_points", as.numeric) %>%
        group_by(time_points, HalfLife_cls_hypoxia) %>% 
        summarise(
                sd = sd(log2_abundance), 
                log2_abundance = mean(log2_abundance)
        ) %>% 
        ggplot(aes(time_points, log2_abundance)) + 
        scale_colour_manual(values = col_hl_cls) +
        geom_line(aes(group = HalfLife_cls_hypoxia, color = HalfLife_cls_hypoxia),
                  position = position_dodge(0.3), size = 5) +
        geom_pointrange(aes(ymin = log2_abundance - sd, ymax = log2_abundance + sd, 
                            color = HalfLife_cls_hypoxia), position = position_dodge(0.3)) +
        ylim(0, 17) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold")
        ) +
        xlab("Minutes after adding rifampicin") + 
        ylab("Log2 mRNA abundance")
ggsave("./viz_halfLife/CDS_halfLifeClassDecayPattern_complete_hypoxia.png", width = 30, height = 20, units = "cm", dpi = 600)

### with only complete leadered OR leaderless genes in each condition
## leadered genes
# combine decay profile and half-life class
decayProfile_hypoxia_CDS_leadered <- decayProfile_hypoxia %>% 
        rownames_to_column("gene") %>% 
        inner_join(hl_CDS_hypoxia, by = "gene") %>% 
        inner_join(gene_leadered, by = "gene") %>% 
        select(-start, -end, -strand)
table(decayProfile_hypoxia_CDS_leadered$HalfLife_cls_hypoxia)
summary(decayProfile_hypoxia_CDS_leadered$HalfLife_vals_hypoxia)

decayProfile_hypoxia_CDS_leadered %>%
        select(-HalfLife_vals_hypoxia) %>% 
        gather(time_points, log2_abundance, -gene, -HalfLife_cls_hypoxia) %>%
        mutate_at("time_points", as.numeric) %>%
        group_by(time_points, HalfLife_cls_hypoxia) %>% 
        summarise(
                sd = sd(log2_abundance), 
                log2_abundance = mean(log2_abundance)
        ) %>% 
        ggplot(aes(time_points, log2_abundance)) + 
        scale_colour_manual(values = col_hl_cls) +
        geom_line(aes(group = HalfLife_cls_hypoxia, color = HalfLife_cls_hypoxia),
                  position = position_dodge(0.3), size = 5) +
        geom_pointrange(aes(ymin = log2_abundance - sd, ymax = log2_abundance + sd, 
                            color = HalfLife_cls_hypoxia), position = position_dodge(0.3)) +
        ylim(0, 17) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold")
        ) +
        xlab("Minutes after adding rifampicin") + 
        ylab("Log2 mRNA abundance")
ggsave("./viz_halfLife/CDS_halfLifeClassDecayPattern_leadered_complete_hypoxia.png", width = 30, height = 20, units = "cm", dpi = 600)

## leaderless genes
# combine decay profile and half-life class
decayProfile_hypoxia_CDS_leaderless <- decayProfile_hypoxia %>% 
        rownames_to_column("gene") %>% 
        inner_join(hl_CDS_hypoxia, by = "gene") %>% 
        inner_join(gene_leaderless, by = "gene") %>% 
        select(-start, -end, -strand)
table(decayProfile_hypoxia_CDS_leaderless$HalfLife_cls_hypoxia)
summary(decayProfile_hypoxia_CDS_leaderless$HalfLife_vals_hypoxia)

decayProfile_hypoxia_CDS_leaderless %>%
        select(-HalfLife_vals_hypoxia) %>% 
        gather(time_points, log2_abundance, -gene, -HalfLife_cls_hypoxia) %>%
        mutate_at("time_points", as.numeric) %>%
        group_by(time_points, HalfLife_cls_hypoxia) %>% 
        summarise(
                sd = sd(log2_abundance), 
                log2_abundance = mean(log2_abundance)
        ) %>% 
        ggplot(aes(time_points, log2_abundance)) + 
        scale_colour_manual(values = col_hl_cls) +
        geom_line(aes(group = HalfLife_cls_hypoxia, color = HalfLife_cls_hypoxia),
                  position = position_dodge(0.3), size = 5) +
        geom_pointrange(aes(ymin = log2_abundance - sd, ymax = log2_abundance + sd, 
                            color = HalfLife_cls_hypoxia), position = position_dodge(0.3)) +
        ylim(0, 17) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), 
                axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 30, face = "bold"),
                axis.title =  element_text(size = 35, face = "bold")
        ) +
        xlab("Minutes after adding rifampicin") + 
        ylab("Log2 mRNA abundance")
ggsave("./viz_halfLife/CDS_halfLifeClassDecayPattern_leaderless_complete_hypoxia.png", width = 30, height = 20, units = "cm", dpi = 600)
