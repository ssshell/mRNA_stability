
####### plot overview of degradation profiling in log phase and hypoxia
####### global expression of CDS in two dimension with all the time points

library(tidyverse)
library(ggsci)
library(umap)
library(RColorBrewer)
library(scales)

# display.brewer.all(type = 'qual')
# display.brewer.pal(n = 12, name = 'Paired')
# brewer.pal(n = 12, name = 'Paired')
# show_col("#1F78B4")

####### get CDS annotations
gene_CDS <- read.table("../index/msmeg_CombinedAnnotation_CDS.bed", col.names = c("gene", "start", "end", "strand"))

#######  overview viz
set.seed(7) 
coverage_norm_complete <- read.table('../normalization/NormCoverage_allReplicates_log2.txt')
coverage_norm_CDS <- coverage_norm_complete[which(rownames(coverage_norm_complete) %in% gene_CDS$gene), ]
coverage_norm_CDS_logPhaseAndHypoxia <- coverage_norm_CDS[, c(19:39, 76:96)]
mydat_CDS_logPhaseAndHypoxia <- as.data.frame(t(coverage_norm_CDS_logPhaseAndHypoxia))

custom.config = umap.defaults
custom.config$n_neighbors = 20
custom.config$min_dist = 0.25
custom.config$n_components = 2
custom.config$random_state = 77

############## test session
####### try scale per column/feature (normalize features) before using umap
####### umap split is basically the same
####### conclusion: scale is not necessary. Especially when features are not independent
###### or have meaningful relationships. Like here, the abundance change over time

# test <- scale(mydat_CDS_logPhaseAndHypoxia)
# umap_CDS_logPhaseAndHypoxia <- umap(test, custom.config)

umap_CDS_logPhaseAndHypoxia <- umap(mydat_CDS_logPhaseAndHypoxia, custom.config)
####### double check parameters
umap_CDS_logPhaseAndHypoxia$config

umap_CDS_logPhaseAndHypoxia_sample <- as.data.frame(umap_CDS_logPhaseAndHypoxia$layout)
umap_CDS_logPhaseAndHypoxia_sample <- tibble::rownames_to_column(umap_CDS_logPhaseAndHypoxia_sample, "sample")
umap_CDS_logPhaseAndHypoxia_sample$condition <- c(rep("logPhase", 21), rep("hypoxia", 21))
umap_CDS_logPhaseAndHypoxia_sample$condition <- factor(umap_CDS_logPhaseAndHypoxia_sample$condition, 
                                                       levels = c("hypoxia", "logPhase"))

####### without time points
umap_CDS_logPhaseAndHypoxia_sample %>% 
        rename(UMAP1 = "V1", UMAP2 = "V2") %>% 
        ggplot(aes(x = UMAP1, y = UMAP2, color = condition)) +  
        geom_point(size = 4, stroke = 0) +
        scale_color_manual(values = c("#83B7DF", "#036D9C")) +
        xlim(-4, 5) +
        ylim(-4, 5) +
        theme_bw() +
        theme(
              panel.border = element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 10, face = "bold"),
              axis.title = element_blank(),
              legend.position = "None",
              )
ggsave("./viz_decayProfile/CDS_globalDecayProfile_condition.png", width = 8, height = 8, units = "cm", dpi = 600)

sample_timePoints <- c()
sampleName_split <- str_split(umap_CDS_logPhaseAndHypoxia_sample$sample, "_")
for (i in 1:42){
        time_point_i <- rev(unlist(sampleName_split[[i]]))[2]
        sample_timePoints <- c(sample_timePoints, time_point_i)
}
time_points <- noquote(sample_timePoints)
umap_CDS_logPhaseAndHypoxia_sample$time_point <- factor(time_points, levels = c("60", "32", "30", "16", "15", "9",
                                                                                "8", "6", "4", "3", "2", "1", "0"))
####### with time points
umap_CDS_logPhaseAndHypoxia_sample %>% 
        rename(UMAP1 = "V1", UMAP2 = "V2") %>% 
        ggplot(aes(x = UMAP1, y = UMAP2)) +  
        geom_point(aes(alpha = time_point), size = 4, stroke = 0) +
        xlim(-4, 5) +
        ylim(-4, 5) +
        theme_bw() +
        theme(
              panel.border = element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
              axis.text = element_text(size = 10, face = "bold"),
              axis.title = element_blank(),
              legend.position = "None",
        )
ggsave("./viz_decayProfile/CDS_globalDecayProfile_timePoints.png", width = 8, height = 8, units = "cm", dpi = 600)
