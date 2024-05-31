
####### get half-life value and percentage classes

library(tidyverse)

####### get CDS annotations
gene_CDS <- read.table("../index/msmeg_CombinedAnnotation_CDS.bed", col.names = c("gene", "start", "end", "strand"))

####### log phase
half_life_logPhase <- read.table("./half_life/half_life_logPhase/msmeg_halfLife_logPhase.txt")
colnames(half_life_logPhase) <- c("gene", "HalfLife_vals_logPhase")

half_life_logPhase_CDS <- half_life_logPhase[which(half_life_logPhase$gene %in% gene_CDS$gene), ]
summary(half_life_logPhase_CDS)

half_life_logPhase_class <- half_life_logPhase_CDS %>%
        mutate(HalfLife_cls_logPhase = case_when(HalfLife_vals_logPhase <= 0.7182 ~ "Fast",
                                      HalfLife_vals_logPhase > 0.7182 & HalfLife_vals_logPhase <= 0.9511 ~ "Med-fast",
                                      HalfLife_vals_logPhase > 0.9511 & HalfLife_vals_logPhase <= 1.3737 ~ "Med-slow",
                                      HalfLife_vals_logPhase > 1.3737 ~ "Slow"))

####### hypoxia
half_life_hypoxia <- read.table("./half_life/half_life_hypoxia/msmeg_halfLife_hypoxia.txt")
colnames(half_life_hypoxia) <- c("gene", "HalfLife_vals_hypoxia")

half_life_hypoxia_CDS <- half_life_hypoxia[which(half_life_hypoxia$gene %in% gene_CDS$gene), ]
summary(half_life_hypoxia_CDS)

half_life_hypoxia_class <- half_life_hypoxia_CDS %>%
        mutate(HalfLife_cls_hypoxia = case_when(HalfLife_vals_hypoxia <= 17.25 ~ "Fast",
                                      HalfLife_vals_hypoxia > 17.25 & HalfLife_vals_hypoxia <= 24.27 ~ "Med-fast",
                                      HalfLife_vals_hypoxia > 24.27 & HalfLife_vals_hypoxia <= 38.78 ~ "Med-slow",
                                      HalfLife_vals_hypoxia > 38.78 ~ "Slow"))

####### half-life fold increase in hypoxia
half_life_FC <- read.table("./half_life/half_life_hyp2logPhaseFC/msmeg_halfLife_hyp2logPhaseFC.txt")
colnames(half_life_FC) <- c("gene", "HalfLife_FC_hypoxiaToLogPhase")

half_life_FC_CDS <- half_life_FC[which(half_life_FC$gene %in% gene_CDS$gene), ]
summary(half_life_FC_CDS)

half_life_FC_class <- half_life_FC_CDS %>%
        mutate(HalfLife_FCcls_hypoxiaToLogPhase = case_when(HalfLife_FC_hypoxiaToLogPhase <= 18.54 ~ "Small",
                                            HalfLife_FC_hypoxiaToLogPhase > 18.54 & HalfLife_FC_hypoxiaToLogPhase <= 26.45 ~ "Med-small",
                                            HalfLife_FC_hypoxiaToLogPhase > 26.45 & HalfLife_FC_hypoxiaToLogPhase <= 40.41 ~ "Med-large",
                                            HalfLife_FC_hypoxiaToLogPhase > 40.41 ~ "Large"))

####### combine classes into one table
###### round log phase half-life and calculate fold change by rounded half-life
##### checked new class assignment. majority remains the same as above
#### add mse of best interval to calculate half-life in log phase and hypoxia
half_life_CDS_full_class_p1 <- half_life_logPhase_class %>%
        mutate(HalfLife_vals_logPhase_round = round(HalfLife_vals_logPhase, 2)) %>% 
        mutate(HalfLife_cls_logPhase_round = case_when(HalfLife_vals_logPhase_round <= 0.72 ~ "Fast",
                                                  HalfLife_vals_logPhase_round > 0.72 & HalfLife_vals_logPhase_round <= 0.95 ~ "Med-fast",
                                                  HalfLife_vals_logPhase_round > 0.95 & HalfLife_vals_logPhase_round <= 1.37 ~ "Med-slow",
                                                  HalfLife_vals_logPhase_round > 1.37 ~ "Slow")) %>% 
        full_join(half_life_hypoxia_class, by = "gene") %>%
        mutate(HalfLife_FC_hypoxiaToLogPhase_round_p1 = HalfLife_vals_hypoxia / HalfLife_vals_logPhase_round) %>% 
        mutate(HalfLife_FC_hypoxiaToLogPhase_round = round(HalfLife_FC_hypoxiaToLogPhase_round_p1, 2)) %>% 
        mutate(HalfLife_FCcls_hypoxiaToLogPhase_round = case_when(HalfLife_FC_hypoxiaToLogPhase_round <= 18.55 ~ "Small",
                                                           HalfLife_FC_hypoxiaToLogPhase_round > 18.55 & HalfLife_FC_hypoxiaToLogPhase_round <= 26.49 ~ "Med-small",
                                                           HalfLife_FC_hypoxiaToLogPhase_round > 26.49 & HalfLife_FC_hypoxiaToLogPhase_round <= 40.37 ~ "Med-large",
                                                           HalfLife_FC_hypoxiaToLogPhase_round > 40.37 ~ "Large")) %>% 
        full_join(half_life_FC_class, by = "gene") 

half_life_CDS_full_class_p2 <- half_life_CDS_full_class_p1 %>% 
        select(gene, HalfLife_vals_logPhase_round, HalfLife_cls_logPhase_round,
               HalfLife_vals_hypoxia, HalfLife_cls_hypoxia, 
               HalfLife_FC_hypoxiaToLogPhase_round, HalfLife_FCcls_hypoxiaToLogPhase_round) %>% 
        rename(HalfLife_vals_logPhase = HalfLife_vals_logPhase_round, 
               HalfLife_cls_logPhase = HalfLife_cls_logPhase_round, 
               HalfLife_FC_hypoxiaToLogPhase = HalfLife_FC_hypoxiaToLogPhase_round, 
               HalfLife_FCcls_hypoxiaToLogPhase = HalfLife_FCcls_hypoxiaToLogPhase_round) %>% 
        arrange(gene)

#### loading MSE for log phase and hypoxia
mse_logPhase <- read.table("./half_life/half_life_logPhase/msmeg_halfLife_MSE_logPhase.txt",
                           col.names = c("gene", "HalfLife_BestIntervalMSE_logPhase_p1"))
mse_logPhase_CDS <- mse_logPhase[which(mse_logPhase$gene %in% gene_CDS$gene), ]
mse_logPhase_CDS$HalfLife_BestIntervalMSE_logPhase <- round(mse_logPhase_CDS$HalfLife_BestIntervalMSE_logPhase_p1, 4)

mse_hypoxia <- read.table("./half_life/half_life_hypoxia/msmeg_halfLife_MSE_hypoxia.txt",
                           col.names = c("gene", "HalfLife_BestIntervalMSE_hypoxia_p1"))
mse_hypoxia_CDS <- mse_hypoxia[which(mse_hypoxia$gene %in% gene_CDS$gene), ]
mse_hypoxia_CDS$HalfLife_BestIntervalMSE_hypoxia <- round(mse_hypoxia_CDS$HalfLife_BestIntervalMSE_hypoxia_p1, 4)

half_life_CDS_full_class <- half_life_CDS_full_class_p2 %>%
        full_join(mse_logPhase_CDS, by = "gene") %>% 
        full_join(mse_hypoxia_CDS, by = "gene") %>% 
        select(gene, HalfLife_vals_logPhase, HalfLife_cls_logPhase, HalfLife_BestIntervalMSE_logPhase,
               HalfLife_vals_hypoxia, HalfLife_cls_hypoxia, HalfLife_BestIntervalMSE_hypoxia,
               HalfLife_FC_hypoxiaToLogPhase, HalfLife_FCcls_hypoxiaToLogPhase)

write.table(half_life_CDS_full_class, "./half_life/CDS_halfLife_completeClass.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

############# test session
############ checking if the full table consistent with data from individual file
test <- half_life_CDS_full_class_p1 %>% 
        # drop_na()
        # filter(HalfLife_cls_logPhase != HalfLife_cls_logPhase_round)
        filter(HalfLife_FCcls_hypoxiaToLogPhase != HalfLife_FCcls_hypoxiaToLogPhase_round)

summary(half_life_CDS_full_class_p1$HalfLife_vals_logPhase_round)
summary(half_life_CDS_full_class_p1$HalfLife_vals_hypoxia)
summary(half_life_CDS_full_class_p1$HalfLife_FC_hypoxiaToLogPhase_round)

half_life_full_class_logPhase <- half_life_CDS_full_class_p1 %>% 
        select(gene, HalfLife_vals_logPhase_round, HalfLife_cls_logPhase_round) %>% 
        drop_na() %>% 
        arrange(gene)
summary(half_life_full_class_logPhase$HalfLife_vals_logPhase_round)
half_life_logPhase_CDS_sort <- half_life_logPhase_CDS %>% 
        arrange(gene)
table(half_life_full_class_logPhase$gene == half_life_logPhase_CDS_sort$gene)

half_life_full_class_hypoxia <- half_life_CDS_full_class_p1 %>% 
        select(gene, HalfLife_vals_hypoxia, HalfLife_cls_hypoxia) %>% 
        drop_na() %>% 
        arrange(gene)
summary(half_life_full_class_hypoxia$HalfLife_vals_hypoxia)
half_life_hypoxia_CDS_sort <- half_life_hypoxia_CDS %>% 
        arrange(gene)
table(half_life_full_class_hypoxia$gene == half_life_hypoxia_CDS_sort$gene)
table(half_life_full_class_hypoxia$HalfLife_vals_hypoxia == half_life_hypoxia_CDS_sort$HalfLife_vals_hypoxia)

half_life_full_class_FC <- half_life_CDS_full_class_p1 %>% 
        select(gene, HalfLife_FC_hypoxiaToLogPhase_round, HalfLife_FCcls_hypoxiaToLogPhase_round) %>% 
        drop_na() %>% 
        arrange(gene)
half_life_FC_CDS_sort <- half_life_FC_CDS %>% 
        arrange(gene)
summary(half_life_full_class_FC$HalfLife_FC_hypoxiaToLogPhase_round)
table(half_life_full_class_FC$gene == half_life_FC_CDS_sort$gene)

####### checking half-life values loading from different files.
####### summary table with half-life classes and fold stabilization in hypoxia
###### compare separate result tables with loading altogether from one single file
##### genes in this single file only have hypoxia half-life for those with fold change in hypoxia
#### for global half-life distribution in individual condition use separate files with the most complete available genes
### checked. half-life & fold change are the same for 3988 genes from these two ways of loading
summary_halfLife <- read.csv("./half_life/half_life_hyp2logPhaseFC/half_life_summary.csv", col.names = c("gene", "HalfLife_vals_logPhase", "HalfLife_vals_hypoxia", "HalfLife_FC_hypoxiaToLogPhase"))
summary_halfLife_CDS <- summary_halfLife[which(summary_halfLife$gene %in% gene_CDS$gene), ]

summary_halfLife_CDS_logPhase <- summary_halfLife_CDS %>% 
        select(gene, HalfLife_vals_logPhase) %>% 
        drop_na() %>% 
        arrange(gene)
summary(summary_halfLife_CDS_logPhase)
half_life_logPhase_CDS_sort <- half_life_logPhase_CDS %>% 
        arrange(gene)
table(half_life_logPhase_CDS_sort$HalfLife_vals_logPhase == summary_halfLife_CDS_logPhase$HalfLife_vals_logPhase)

summary_halfLife_CDS_hypoxia <- summary_halfLife_CDS %>% 
        select(gene, HalfLife_vals_hypoxia) %>% 
        drop_na() %>% 
        arrange(gene) %>% 
        mutate(HalfLife_vals_hypoxia_round = round(HalfLife_vals_hypoxia, 2))
summary(summary_halfLife_CDS_hypoxia)
half_life_hypoxia_CDS_comm <- half_life_hypoxia_CDS[which(half_life_hypoxia_CDS$gene %in% summary_halfLife_CDS_hypoxia$gene), ]
half_life_hypoxia_CDS_comm_sort <- half_life_hypoxia_CDS_comm %>% 
        arrange(gene)
table(half_life_hypoxia_CDS_comm_sort$HalfLife_vals_hypoxia == summary_halfLife_CDS_hypoxia$HalfLife_vals_hypoxia_round)

summary_halfLife_CDS_FC <- summary_halfLife_CDS %>% 
        select(gene, HalfLife_FC_hypoxiaToLogPhase) %>% 
        filter(HalfLife_FC_hypoxiaToLogPhase != 0) %>% 
        arrange(gene)
half_life_FC_CDS_sort <- half_life_FC_CDS %>% 
        arrange(gene)
table(half_life_FC_CDS_sort$HalfLife_FC_hypoxiaToLogPhase == summary_halfLife_CDS_FC$HalfLife_FC_hypoxiaToLogPhase)

####### checking half-life MSE values loading from different files.
summary_mse_logPhase <- read.csv("./half_life/half_life_hyp2logPhaseFC/half_life_summary_MSE_logPhase.csv",
                                 col.names = c("gene", "mse_logPhase"))
summary_mse_logPhase_CDS <- summary_mse_logPhase[which(summary_mse_logPhase$gene %in% gene_CDS$gene), ]
summary_mse_logPhase_CDS_sort <- arrange(summary_mse_logPhase_CDS, gene)
mse_logPhase_CDS_sort <- arrange(mse_logPhase_CDS, gene)
table(mse_logPhase_CDS_sort$gene == summary_mse_logPhase_CDS_sort$gene)
table(mse_logPhase_CDS_sort$HalfLife_BestIntervalMSE_logPhase_p1 == summary_mse_logPhase_CDS_sort$mse_logPhase)

summary_mse_hypoxia <- read.table("./half_life/half_life_hyp2logPhaseFC/half_life_summary_MSE_hypoxia.txt",
                                 col.names = c("gene", "mse_hypoxia"))
summary_mse_hypoxia_CDS <- summary_mse_hypoxia[which(summary_mse_hypoxia$gene %in% gene_CDS$gene), ]
summary_mse_hypoxia_CDS_sort <- arrange(summary_mse_hypoxia_CDS, gene)
mse_hypoxia_CDS_sort <- arrange(mse_hypoxia_CDS, gene)
table(mse_hypoxia_CDS_sort$gene == summary_mse_hypoxia_CDS_sort$gene)
table(mse_hypoxia_CDS_sort$HalfLife_BestIntervalMSE_hypoxia_p1 == summary_mse_hypoxia_CDS_sort$mse_hypoxia)

